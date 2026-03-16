#' Plot Gene Dynamics Trends Along Trajectory Tubes
#'
#' @description
#' Generates a high-clarity plot showing smoothed maturation trends for specific genes.
#' Automatically maps cells to their respective overlapping biological tubes (e.g., Tube_Cd4, Tube_Cd8a)
#' and plots the expression dynamics across Tocky Time. Cells with dual identities properly
#' contribute to the shared origins of multiple trajectories.
#'
#' @param object An mCanonicalTockyObj containing tube assignments.
#' @param seurat_obj The original Seurat object to fetch gene expression from.
#' @param features Character vector of genes to plot.
#' @param ncol Integer. Number of columns for the plot grid (default 2).
#' @param ... Additional arguments passed to methods. 
#' @export
#' @importFrom ggplot2 ggplot aes geom_smooth scale_color_manual scale_fill_manual theme_minimal labs theme element_text element_blank element_line
#' @importFrom patchwork wrap_plots
#' @importFrom Seurat FetchData
#' @importFrom stats setNames
#' @importFrom rlang .data
setGeneric("plotLineageGeneDynamics", function(object, seurat_obj, features, ...) standardGeneric("plotLineageGeneDynamics"))

#' @rdname plotLineageGeneDynamics
#' @export
setMethod("plotLineageGeneDynamics", signature(object = "mCanonicalTockyObj"),
          function(object, seurat_obj, features, ncol = 2) {
            
  df <- object@TockyData
  
  tube_cols <- grep("^Tube_", colnames(df), value = TRUE)
  if (length(tube_cols) == 0) {
    stop("No biological Tube definitions found. Please run mExtractTrajectoryTubes() first.")
  }
  
  if (!"angle" %in% colnames(df)) stop("Tocky angle missing.")
  
  valid_features <- intersect(features, rownames(seurat_obj))
  if (length(valid_features) == 0) {
    stop("None of the requested features were found in the Seurat object.")
  }
  
  valid_cells <- rownames(df)
  expr_data <- Seurat::FetchData(seurat_obj, vars = valid_features, cells = valid_cells)
  
  long_data_list <- list()
  
  for (t_col in tube_cols) {
    trajectory_name <- gsub("^Tube_", "", t_col)
    in_tube_cells <- rownames(df)[df[[t_col]] == TRUE]
    
    if (length(in_tube_cells) == 0) next
    
    tube_df <- data.frame(
      cell_barcode = in_tube_cells,
      angle = df[in_tube_cells, "angle"],
      Trajectory = trajectory_name
    )
    
    tube_df <- cbind(tube_df, expr_data[in_tube_cells, valid_features, drop = FALSE])
    long_data_list[[trajectory_name]] <- tube_df
  }
  
  plot_data <- do.call(rbind, long_data_list)
  if (is.null(plot_data) || nrow(plot_data) == 0) stop("No cells found within any trajectory tubes.")
  
  trajectories <- unique(plot_data$Trajectory)
  path_colors <- c("#000066", "#8B0000", "#006400", "#4B0082", "#FF8C00")
  
  if(length(trajectories) <= length(path_colors)) {
    color_mapping <- setNames(path_colors[1:length(trajectories)], trajectories)
  } else {
    color_mapping <- setNames(scales::hue_pal()(length(trajectories)), trajectories)
  }
  
  plot_list <- lapply(valid_features, function(gene) {
    ggplot(plot_data, aes(x = .data$angle, y = .data[[gene]],
                                            color = .data$Trajectory, fill = .data$Trajectory)) +
      geom_smooth(method = "gam",
                           formula = y ~ s(x, bs = "cs"),
                           linewidth = 1.8,
                           alpha = 0.2,
                           se = TRUE) +
      scale_color_manual(values = color_mapping) +
      scale_fill_manual(values = color_mapping) +
      theme_minimal() +
      labs(title = gene,
                    x = "Tocky Time",
                    y = "Expression",
                    color = "Trajectory", fill = "Trajectory") +
      theme(
        plot.title = element_text(face = "bold", size = 14),
        legend.position = "bottom",
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")
      )
  })
  
  return(patchwork::wrap_plots(plot_list, ncol = ncol))
})
