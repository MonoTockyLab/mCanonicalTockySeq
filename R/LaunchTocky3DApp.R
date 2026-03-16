#' Launch Interactive 3D Fate Explorer App
#'
#' @description
#' Launches a local Shiny web application to explore the 3D Tocky trajectory.
#' Utilizes dynamically anchored overlapping multi-identity tubes from mExtractTrajectoryTubes().
#'
#' @param object An mCanonicalTockyObj.
#' @param seurat_obj Original Seurat object to fetch gene expression from.
#' @param window_size Numeric. Width of the sliding angle slice (default 10).
#' @param step_size Numeric. Step size for the sliding window (default 1).
#' @param span Numeric. LOESS smoothing span (default 0.4).
#' @export
#' @importFrom shiny fluidPage sidebarLayout sidebarPanel mainPanel radioButtons checkboxInput textInput textOutput conditionalPanel shinyApp renderText reactive req hr helpText titlePanel
#' @importFrom plotly plot_ly add_trace plotlyOutput renderPlotly layout
#' @importFrom Seurat GetAssayData
#' @importFrom stats predict loess setNames
LaunchTocky3DApp <- function(object, seurat_obj, window_size = 10, step_size = 1, span = 0.4) {
            
  if (!requireNamespace("shiny", quietly = TRUE)) stop("Package 'shiny' is required.")
  

  df_3d <- object@TockyData
  score_cols <- grep("_Score$", colnames(df_3d), value = TRUE)
  if (length(score_cols) < 2) stop("Need at least two fate scores in TockyData.")
  
  dest1_col <- score_cols[1]
  dest2_col <- score_cols[2]
  destinations <- gsub("_Score$", "", score_cols)
  
  if (!all(c("angle", "Lineage_Identity") %in% colnames(df_3d))) {
    stop("Lineage_Identity missing. Run mExtractTrajectoryTubes() first.")
  }
  
  df_3d <- df_3d[!is.na(df_3d$angle), ]
  df_3d$cell_id <- rownames(df_3d)
  
  # Identify all dynamic tube columns
  tube_cols <- grep("^Tube_", colnames(df_3d), value = TRUE)
  if (length(tube_cols) == 0) stop("No Tube definitions found.")
  

  cat("Pre-computing 3D LOESS Manifolds for each independent biological tube...\n")
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
  
  tryCatch({ expr_mat <- Seurat::GetAssayData(seurat_obj, layer = "data") },
           error = function(e) { expr_mat <- Seurat::GetAssayData(seurat_obj, slot = "data") })
  available_genes <- rownames(expr_mat)

  identities <- unique(df_3d$Lineage_Identity)
  id_colors <- setNames(grDevices::rainbow(length(identities)), identities)
  if ("Unassigned" %in% names(id_colors)) id_colors["Unassigned"] <- "#E0E0E0"


  ui <- shiny::fluidPage(
    shiny::titlePanel(paste("mCanonicalTockySeq: 3D Fate Explorer -", destinations[1], "vs", destinations[2])),
    
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        shiny::radioButtons("color_mode", "Color Cells By:",
                            choices = c("Tocky Time" = "time",
                                        "Gene-Associated Identities" = "lineage",
                                        "Gene Expression" = "gene")),
        
        shiny::conditionalPanel(
          condition = "input.color_mode == 'gene'",
          shiny::textInput("selected_gene", "Enter Gene Symbol (e.g., Cd4):", value = ""),
          shiny::textOutput("gene_status")
        ),
        
        shiny::hr(),
        shiny::checkboxInput("tube_filter", "Hide 'Unassigned' Noise Cells", value = FALSE),
        shiny::hr(),
        shiny::helpText("LOESS-smoothed trajectories are calculated independently for each gene-associated tube defined by mExtractTrajectoryTubes().")
      ),
      
      shiny::mainPanel(
        plotly::plotlyOutput("tocky3d", height = "700px")
      )
    )
  )


  server <- function(input, output, session) {
    
    output$gene_status <- shiny::renderText({
      shiny::req(input$color_mode == "gene")
      gene <- trimws(input$selected_gene)
      if (gene == "") return("Waiting for input...")
      if (!gene %in% available_genes) return(paste("Gene '", gene, "' not found."))
      return(paste("Showing expression for:", gene))
    })
    
    plot_data <- shiny::reactive({
      dat <- df_3d
      if (input$tube_filter) {
        dat <- dat[dat$Lineage_Identity != "Unassigned", ]
      }
      
      mode <- input$color_mode
      if (mode == "time") {
        dat$plot_color <- dat$angle
      } else if (mode == "lineage") {
        dat$plot_color <- dat$Lineage_Identity
      } else if (mode == "gene") {
        gene <- trimws(input$selected_gene)
        if (gene != "" && gene %in% available_genes) {
          dat$plot_color <- as.numeric(expr_mat[gene, dat$cell_id])
        } else {
          dat$plot_color <- 0
        }
      }
      return(dat)
    })
    
    output$tocky3d <- plotly::renderPlotly({
      dat <- plot_data()
      mode <- input$color_mode
      
      p <- plotly::plot_ly()
      
      if (mode == "time") {
              p <- plotly::add_trace(p, data = dat,
                                     x = ~get(dest1_col), y = ~get(dest2_col), z = ~angle,
                                     type = 'scatter3d', mode = 'markers',
                                     marker = list(
                                       color = ~plot_color,

                                       colorscale = list(c(0, "blue"), c(0.5, "purple"), c(1, "red")),
                                       cmin = 0,
                                       cmax = 90,
                                       size = 2,
                                       opacity = 0.2,
                                       showscale = TRUE,
                                       colorbar = list(
                                         title = "Tocky Time",
                                         len = 0.3,
                                         thickness = 15
                                       )
                                     ),
                                     name = "Cells")
              
        } else if (mode == "lineage") {

        identity_counts <- table(dat$Lineage_Identity)
        
        dat$legend_labels <- paste0(dat$Lineage_Identity, " (n=", identity_counts[dat$Lineage_Identity], ")")
        
        dynamic_colors <- id_colors[names(identity_counts)]
        names(dynamic_colors) <- paste0(names(identity_counts), " (n=", identity_counts, ")")
        
        p <- plotly::add_trace(p, data = dat, x = ~get(dest1_col), y = ~get(dest2_col), z = ~angle,
                               color = ~legend_labels, colors = dynamic_colors,
                               type = 'scatter3d', mode = 'markers',
                               marker = list(size = 2, opacity = 0.6))
        
      } else if (mode == "gene") {
              gene_expr <- dat$plot_color
              cmax_val <- max(gene_expr, na.rm = TRUE)
              if (cmax_val == 0) cmax_val <- 1
              
              p <- plotly::add_trace(p,
                                     data = dat,
                                     x = ~get(dest1_col), y = ~get(dest2_col), z = ~angle,
                                     type = 'scatter3d',
                                     mode = 'markers',
                                     marker = list(
                                       color = ~plot_color,
                                       colorscale = list(c(0, "#D3D3D3"), c(1, "#0000FF")),
                                       cmin = 0,
                                       cmax = cmax_val,
                                       size = 2.5,
                                       opacity = 0.3,
                                       showscale = TRUE,
                                       colorbar = list(
                                         title = "Log Exp",
                                         len = 0.2,
                                         thickness = 15,
                                         y = 0.75
                                       )
                                     ),
                                     name = "Expression",
                                     text = ~paste("Expression:", round(plot_color, 2)),
                                     hoverinfo = "text")
            }
      
      path_colors <- c("#000066", "#8B0000", "#006400", "#4B0082", "#FF8C00")
      i <- 1
      for (gene_path in names(path_data_list)) {
        p_data <- path_data_list[[gene_path]]
        line_col <- path_colors[((i - 1) %% length(path_colors)) + 1]
        
        p <- plotly::add_trace(p, data = p_data, x = ~s1, y = ~s2, z = ~Tocky_Time,
                               type = 'scatter3d', mode = 'lines',
                               line = list(color = line_col, width = 10),
                               name = paste(gene_path, 'Trajectory'))
        i <- i + 1
      }
      
      plotly::layout(p,
                     scene = list(
                       xaxis = list(title = paste(destinations[1], 'differentiation score')),
                       yaxis = list(title = paste(destinations[2], 'differentiation score')),
                       zaxis = list(title = 'Tocky Time (0-90)')
                     ),
                     legend = list(
                       itemsizing = 'constant',
                       font = list(size = 14)
                     ))
    })
  }
  
  shiny::shinyApp(ui = ui, server = server)
}
