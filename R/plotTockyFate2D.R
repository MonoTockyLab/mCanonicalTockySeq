

#' Plot 2D Fate Space Maturation
#'
#' @description Generates a 2D scatter plot of lineage divergence,
#' colored continuously by Tocky Time (Blue -> Red).
#'
#' @param object An mCanonicalTockyObj with fate scores calculated.
#' @param ... Additional arguments passed to methods.
#' @export
#' @importFrom ggplot2 ggplot aes geom_point scale_color_identity geom_hline geom_vline labs theme_minimal theme element_text
#' @importFrom scales col_numeric
setGeneric("plotTockyFate2D", function(object, ...) standardGeneric("plotTockyFate2D"))

#' @rdname plotTockyFate2D
#' @export
setMethod("plotTockyFate2D", signature(object = "mCanonicalTockyObj"),
          function(object) {
            
  df <- object@TockyData
  score_cols <- grep("_Score$", colnames(df), value = TRUE)
  
  if (length(score_cols) != 2 || !"angle" %in% colnames(df)) {
    stop("Required columns missing. Run mGetFateScores() and mGradientTockySeq() first.")
  }
  
  x_col <- score_cols[1]
  y_col <- score_cols[2]
  destinations <- gsub("_Score$", "", score_cols)
  
  plot_df <- df[!is.na(df$angle), ]
  
  tocky_pal <- scales::col_numeric(palette = c("blue", "purple", "red"), domain = c(0, 90))
  plot_df$PointColor <- tocky_pal(plot_df$angle)

  
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]], color = .data$PointColor)) +
    ggplot2::geom_point(alpha = 0.3, size = 1.2) +
    ggplot2::scale_color_identity() +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::labs(
      title = "Lineage Divergence (2D Fate Space)",
      subtitle = "Color: Tocky Maturation (Blue=New, Red=Aged) | Alpha = 0.3",
      x = paste(destinations[1], "Commitment Score"),
      y = paste(destinations[2], "Commitment Score")
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
  
  return(p)
})



#' Plot Canonical Tocky Locus Features
#'
#' @description Generates a matrix plot showing Tocky progression and gene expression
#' across binned maturation stages in fate space.
#'
#' @param object An mCanonicalTockyObj.
#' @param seurat_obj The original Seurat object to fetch expression from.
#' @param features Character vector of genes to plot.
#' @param bins Number of angle bins (default 5).
#' @param ... Additional arguments passed to methods. 
#' @export
#' @importFrom ggplot2 ggplot aes geom_point scale_color_identity scale_color_gradient geom_hline geom_vline facet_wrap labs theme_minimal theme element_text element_rect element_blank unit
#' @importFrom patchwork wrap_plots plot_annotation
#' @importFrom dplyr mutate
#' @importFrom utils tail
#' @importFrom Seurat FetchData
setGeneric("plotTockyLocusFeatures", function(object, seurat_obj, features, ...) standardGeneric("plotTockyLocusFeatures"))

#' @rdname plotTockyLocusFeatures
#' @export
setMethod("plotTockyLocusFeatures", signature(object = "mCanonicalTockyObj"),
          function(object, seurat_obj, features, bins = 5) {
            
  df <- object@TockyData
  score_cols <- grep("_Score$", colnames(df), value = TRUE)
  
  if (length(score_cols) != 2 || !"angle" %in% colnames(df)) {
    stop("Required columns missing. Run mGetFateScores() first.")
  }
  
  x_axis <- score_cols[1]
  y_axis <- score_cols[2]
  destinations <- gsub("_Score$", "", score_cols)
  
  plot_df <- df[!is.na(df$angle), ]
  
  tocky_pal <- scales::col_numeric(palette = c("blue", "purple", "red"), domain = c(0, 90))
  plot_df$PointColor <- tocky_pal(plot_df$angle)
  
  valid_features <- intersect(features, rownames(seurat_obj))
  if(length(valid_features) == 0) stop("None of the features were found in the Seurat object.")
  
  expr_data <- Seurat::FetchData(seurat_obj, vars = valid_features, cells = rownames(plot_df))
  plot_df <- cbind(plot_df, expr_data[rownames(plot_df), , drop = FALSE])
  
  breaks <- seq(0, 90, length.out = bins + 1)
  labels <- paste0(round(head(breaks, -1)), "-", round(tail(breaks, -1)))
  
  plot_df$Angle_Category <- cut(plot_df$angle, breaks = breaks, include.lowest = TRUE, labels = labels)
  
  base_theme <- ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 11),
      strip.background = ggplot2::element_rect(fill = "grey95", color = NA),
      strip.text = ggplot2::element_text(face = "bold"),
      panel.spacing = ggplot2::unit(1, "lines"),
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank()
    )

  p_top <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data[[x_axis]], y = .data[[y_axis]], color = .data$PointColor)) +
    ggplot2::geom_point(alpha = 0.3, size = 0.8) +
    ggplot2::scale_color_identity() +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::facet_wrap(~Angle_Category, nrow = 1) +
    ggplot2::labs(title = "1. Tocky Time", y = paste(destinations[2], "Score")) +
    base_theme

  gene_plots <- lapply(seq_along(valid_features), function(i) {
    gene <- valid_features[i]
    is_last <- i == length(valid_features)
    
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data[[x_axis]], y = .data[[y_axis]], color = .data[[gene]])) +
      ggplot2::geom_point(alpha = 0.4, size = 0.8) +
      ggplot2::scale_color_gradient(low = "grey90", high = "blue", name = "Exp") +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
      ggplot2::facet_wrap(~Angle_Category, nrow = 1) +
      ggplot2::labs(title = paste0(i + 1, ". ", gene, " Expression"), y = paste(destinations[2], "Score")) +
      base_theme +
      ggplot2::theme(strip.background = ggplot2::element_blank(), strip.text = ggplot2::element_blank())
    
    if (is_last) {
      p <- p + ggplot2::theme(axis.title.x = ggplot2::element_text(), axis.text.x = ggplot2::element_text()) +
               ggplot2::labs(x = paste(destinations[1], "Commitment Score"))
    }
    return(p)
  })

  final_plot <- patchwork::wrap_plots(c(list(p_top), gene_plots), ncol = 1) +
    patchwork::plot_annotation(title = "mCanonicalTockySeq: Developmental Locus Map",
                               theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 14, hjust = 0.5)))
  
  return(final_plot)
})


#' Plot Smooth 3D Fate Trajectory (Expression-Anchored)
#'
#' @description
#' Generates an interactive 3D Plotly visualization of cells flowing through Tocky Time
#' alongside dynamically calculated LOESS-smoothed lineage paths. Automatically detects
#' biologically anchored tubes generated by mExtractTrajectoryTubes().
#'
#' @param object An mCanonicalTockyObj with fate scores and extracted tubes in @TockyData.
#' @param color_by Character. Either "time" (default) or "lineage".
#' @param window_size Numeric. Width of the sliding angle slice (default 10).
#' @param step_size Numeric. Step size for the sliding window (default 1).
#' @param span Numeric. LOESS smoothing span (default 0.4).
#' @param ... Additional arguments passed to methods. 
#' @export
#' @importFrom dplyr filter mutate group_by summarise arrange ungroup bind_rows n
#' @importFrom plotly plot_ly add_trace
#' @importFrom stats predict loess
#' @importFrom grDevices rainbow
setGeneric("plotTockyFate3D", function(object, ...) standardGeneric("plotTockyFate3D"))

#' @rdname plotTockyFate3D
#' @export
setMethod("plotTockyFate3D", signature(object = "mCanonicalTockyObj"),
          function(object, color_by = "time", window_size = 10, step_size = 1, span = 0.4) {
            
  color_by <- match.arg(color_by, c("time", "lineage"))
            
  df_3d <- object@TockyData
  
  score_cols <- grep("_Score$", colnames(df_3d), value = TRUE)
  if (length(score_cols) < 2) {
    stop("Could not auto-detect fate scores. Ensure mGetFateScores() was run first.")
  }
  
  dest1_col <- score_cols[1]
  dest2_col <- score_cols[2]
  destinations <- gsub("_Score$", "", score_cols)
  
  if (!"angle" %in% colnames(df_3d)) {
    stop("Tocky angle missing. Run mGradientTockySeq() first.")
  }
  
  if (!"Lineage_Identity" %in% colnames(df_3d)) {
    stop("Biological trajectories missing. Run mExtractTrajectoryTubes() first to anchor the manifold.")
  }
  
  df_3d <- df_3d[!is.na(df_3d$angle), ]
  
  tube_cols <- grep("^Tube_", colnames(df_3d), value = TRUE)
  if (length(tube_cols) == 0) stop("No biological Tube definitions found in the object.")
  
   slide_starts <- seq(0, 90 - window_size, by = step_size)
  path_data_list <- list()
  
  for (t_col in tube_cols) {
    gene_name <- gsub("^Tube_", "", t_col)
    
    tube_cells <- df_3d[df_3d[[t_col]] == TRUE, ]
    
    if (nrow(tube_cells) < 10) next
    
    centroids <- do.call(rbind, lapply(slide_starts, function(start) {
      slice <- tube_cells[tube_cells$angle >= start & tube_cells$angle < (start + window_size), ]
      if(nrow(slice) < 5) return(NULL)
      data.frame(
        Tocky_Time = start + (window_size / 2),
        s1 = mean(slice[[dest1_col]], na.rm = TRUE),
        s2 = mean(slice[[dest2_col]], na.rm = TRUE),
        n_cells = nrow(slice)
      )
    }))
    
    if (!is.null(centroids) && nrow(centroids) > 5) {
      centroids$s1 <- stats::predict(stats::loess(s1 ~ Tocky_Time, data = centroids, span = span))
      centroids$s2 <- stats::predict(stats::loess(s2 ~ Tocky_Time, data = centroids, span = span))
      path_data_list[[gene_name]] <- centroids
    }
  }

  p <- plotly::plot_ly()
  
  if (color_by == "time") {
     p <- plotly::add_trace(p, data = df_3d,
                           x = ~get(dest1_col), y = ~get(dest2_col), z = ~angle,
                           type = 'scatter3d', mode = 'markers',
                           marker = list(
                             size = 1.5, opacity = 0.15,
                             color = ~angle,
                             colorscale = list(c(0, "blue"), c(0.5, "purple"), c(1, "red")),
                             cmin = 0, cmax = max(df_3d$angle, na.rm = TRUE),
                             showscale = TRUE,
                             colorbar = list(
                               title = "Tocky Time",
                               len = 0.4,
                               thickness = 20
                             )
                           ),
                           showlegend = FALSE,
                           hoverinfo = "text",
                           text = ~paste("Tocky Time:", round(angle, 1),
                                         "<br>Identity:", Lineage_Identity))
  } else {
    identities <- unique(df_3d$Lineage_Identity)
    id_colors <- setNames(grDevices::rainbow(length(identities)), identities)
    if ("Unassigned" %in% names(id_colors)) id_colors["Unassigned"] <- "#E0E0E0"
    
    p <- plotly::add_trace(p, data = df_3d,
                           x = ~get(dest1_col), y = ~get(dest2_col), z = ~angle,
                           color = ~Lineage_Identity, colors = id_colors,
                           type = 'scatter3d', mode = 'markers',
                           marker = list(size = 1.5, opacity = 0.25),
                           showlegend = TRUE,
                           hoverinfo = "text",
                           text = ~paste("Tocky Time:", round(angle, 1),
                                         "<br>Identity:", Lineage_Identity))
  }

  path_colors <- c("#000066", "#8B0000", "#006400", "#4B0082", "#FF8C00")
  i <- 1
  
  for (gene_path in names(path_data_list)) {
    p_data <- path_data_list[[gene_path]]
    line_col <- path_colors[((i - 1) %% length(path_colors)) + 1]
    
    p <- plotly::add_trace(p, data = p_data,
                           x = ~s1, y = ~s2, z = ~Tocky_Time,
                           type = 'scatter3d', mode = 'lines',
                           line = list(color = line_col, width = 10),
                           name = paste(gene_path, 'Trajectory'))
    i <- i + 1
  }
  
  p <- plotly::layout(p,
    title = "Expression-Anchored 3D Thymic Trajectories",
    legend = list(itemsizing = 'constant', orientation = "v", x = 1.1),
    scene = list(
      xaxis = list(title = paste(destinations[1], 'Commitment')),
      yaxis = list(title = paste(destinations[2], 'Commitment')),
      zaxis = list(title = 'Tocky Time (0-90)')
    )
  )
  
  return(p)
})
