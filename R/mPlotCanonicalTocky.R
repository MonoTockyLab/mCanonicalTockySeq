#' Plot Canonical Tocky Ordination with Continuous Gradient Legend
#'
#' @description Visually evaluates the constrained ordination space. Automatically detects
#' the number of available RDA axes and generates a multi-panel layout showing all
#' sequential axis pairs (1v2, 2v3, 3v4, etc.) followed by a continuous Tocky Time legend.
#'
#' @param object An mCanonicalTockyObj.
#' @param alpha_level Numeric transparency level for cells (0-1). Default is 0.2.
#' @param ... Additional arguments passed to methods.
#'
#' @export
#' @importFrom grDevices colorRampPalette rgb as.raster col2rgb
#' @importFrom graphics plot abline arrows text par rasterImage rect title layout segments
setGeneric("mPlotCanonicalTocky", function(object, ...) standardGeneric("mPlotCanonicalTocky"))

#' @rdname mPlotCanonicalTocky
#' @export
setMethod("mPlotCanonicalTocky", signature(object = "mCanonicalTockyObj"),
          function(object, alpha_level = 0.2) {

  if (!"angle" %in% colnames(object@TockyData)) {
    warning("Tocky Time (angle) not found. Cells will be colored grey. Run mGradientTockySeq() for coloring.")
    tocky_time <- rep(NA, ncol(object@X))
    names(tocky_time) <- colnames(object@X)
  } else {
    tocky_time <- object@TockyData$angle
    names(tocky_time) <- rownames(object@TockyData)
  }
  
  cell_scores <- object@cell_scores
  biplot_vecs <- object@biplot

  num_dims <- ncol(cell_scores)
  if (num_dims < 2) stop("At least 2 dimensions are required to plot.")
  num_pairs <- num_dims - 1

  graphics::layout(matrix(1:(num_pairs + 1), nrow = 1), widths = c(rep(1, num_pairs), 0.6))

  n_colors <- 100
  rbPal <- grDevices::colorRampPalette(c("blue", "purple", "red"))
  color_map <- rbPal(n_colors)
  rgb_vals <- grDevices::col2rgb(color_map)
  transparent_map <- grDevices::rgb(rgb_vals[1,], rgb_vals[2,], rgb_vals[3,],
                                    maxColorValue = 255, alpha = alpha_level * 255)
  grey_transparent <- grDevices::rgb(0.8, 0.8, 0.8, alpha = alpha_level)

  cell_colors <- rep(grey_transparent, nrow(cell_scores))
  names(cell_colors) <- rownames(cell_scores)
  
  valid_cells <- intersect(names(tocky_time)[!is.na(tocky_time)], names(cell_colors))
  if (length(valid_cells) > 0) {
    vals <- pmax(0, pmin(90, tocky_time[valid_cells]))
    idx <- floor((vals / 90) * (n_colors - 1)) + 1
    cell_colors[valid_cells] <- transparent_map[idx]
  }


  for (i in 1:num_pairs) {
    ax_x <- i
    ax_y <- i + 1
    
    graphics::par(mar = c(5, 4, 4, 1))
    coords <- cell_scores[, c(ax_x, ax_y)]
    vecs   <- biplot_vecs[, c(ax_x, ax_y)]
    
    scale_fac <- (max(abs(coords)) / max(abs(vecs))) * 0.5
    vecs_sc <- vecs * scale_fac
    
    graphics::plot(coords, col = cell_colors, pch = 16, cex = 1.0,
                   xlab = paste("Axis", ax_x), ylab = paste("Axis", ax_y),
                   main = paste("Projection:", ax_x, "vs", ax_y),
                   las = 1, bty = "l")
    graphics::abline(h = 0, v = 0, lty = 2, col = "grey80")
    graphics::arrows(0, 0, vecs_sc[,1], vecs_sc[,2], length = 0.1, lwd = 2, col = "black")
    graphics::text(vecs_sc[,1]*1.2, vecs_sc[,2]*1.2, rownames(vecs), col="black", font=2, cex=1.0)
  }

  graphics::par(mar = c(5, 0, 4, 2))
  graphics::plot(c(0, 1), c(0, 1), type = "n", axes = FALSE, xlab = "", ylab = "", main = "")
  graphics::title("Tocky Time", line = 1)
  
  legend_raster <- grDevices::as.raster(matrix(rev(transparent_map), ncol=1))
  graphics::rasterImage(legend_raster, 0.2, 0.2, 0.3, 0.8)
  
  tick_y <- c(0.2, 0.5, 0.8)
  tick_lab <- c("0 (New)", "45 (Persistent)", "90 (Arrested)")
  tick_col <- c("blue", "purple", "red")
  
  graphics::segments(0.3, tick_y, 0.35, tick_y, lwd = 1)
  graphics::text(0.38, tick_y, labels = tick_lab, pos = 4, col = tick_col, cex = 1.1, font=2)
  
  graphics::rect(0.2, 0.05, 0.3, 0.1, col = grey_transparent, border = NA)
  graphics::text(0.38, 0.075, labels = "Timer Neg", pos = 4, col = "grey50", cex = 1.1)
})



#' Plot Gene Expression Dynamics by Group
#'
#' @description Visualizes single-cell gene expression trends along the Tocky Time axis.
#' Accepts multiple genes and automatically generates a grid layout.
#' Automatically detects Lineage_Bias from the mCanonicalTockyObj.
#'
#' @param object An mCanonicalTockyObj.
#' @param seurat_obj Original Seurat object.
#' @param features Character vector of genes to plot.
#' @param group_by Optional metadata column in Seurat or TockyData to color by. Defaults to "Lineage_Bias".
#' @param span Numeric. Smoothing span for LOESS (default 0.8).
#' @param ncol Integer. Number of columns for the plot grid (default 2).
#' @param jitter_amount Numeric. Amount of vertical jitter applied to expression values to reduce overplotting (default 0.2).
#' @param pt_alpha Numeric. Transparency level for the plotted cells, ranging from 0 (transparent) to 1 (opaque) (default 0.2).
#' @param m Numeric. Multiplier for the maximum y-axis limit to ensure curves and data points fit comfortably within the plot (default 1.15).
#' @param ... Additional arguments passed to methods.
#' @export
#' @importFrom graphics plot rect lines legend grid par
#' @importFrom stats loess predict setNames
#' @importFrom grDevices rgb adjustcolor rainbow
#' @importFrom Seurat GetAssayData
setGeneric("mPlotGeneDynamics", function(object, seurat_obj, features, ...) standardGeneric("mPlotGeneDynamics"))

#' @rdname mPlotGeneDynamics
#' @export
setMethod("mPlotGeneDynamics", signature(object = "mCanonicalTockyObj"),
          function(object, seurat_obj, features, group_by = "Lineage_Bias",
                   span = 0.8, jitter_amount = 0.2, pt_alpha = 0.2, m = 1.15, ncol = 2) {
            
  if (!"angle" %in% colnames(object@TockyData)) stop("Run mGradientTockySeq() first.")
  
  tocky_time <- object@TockyData$angle
  names(tocky_time) <- rownames(object@TockyData)
  
  tryCatch({
    expr_mat <- Seurat::GetAssayData(seurat_obj, layer = "data")
  }, error = function(e) {
    expr_mat <- Seurat::GetAssayData(seurat_obj, slot = "data")
  })
  
  valid_features <- intersect(features, rownames(expr_mat))
  if (length(valid_features) == 0) stop("None of the requested features were found.")
  
  common_cells <- intersect(colnames(expr_mat), names(tocky_time))

  if (!is.null(group_by)) {
    if (group_by %in% colnames(object@TockyData)) {
      groups <- object@TockyData[[group_by]]
      names(groups) <- rownames(object@TockyData)
    } else if (group_by %in% colnames(seurat_obj@meta.data)) {
      groups <- seurat_obj@meta.data[[group_by]]
      names(groups) <- colnames(seurat_obj)
    } else {
      warning("Group not found. Plotting all cells in black.")
      groups <- NULL
    }
    if (!is.null(groups)) common_cells <- intersect(common_cells, names(groups))
  } else {
    groups <- NULL
  }

  valid_cells <- common_cells[!is.na(tocky_time[common_cells])]
  if (!is.null(groups)) valid_cells <- valid_cells[!is.na(groups[valid_cells])]
  
  plot_x <- as.numeric(tocky_time[valid_cells])

  if (!is.null(groups)) {
    plot_group <- as.factor(groups[valid_cells])
    unique_groups <- levels(plot_group)
    pt_colors_base <- setNames(grDevices::rainbow(length(unique_groups)), unique_groups)

    if (all(c("CD4_Path", "CD8_Path") %in% unique_groups)) pt_colors_base[c("CD4_Path", "CD8_Path")] <- c("#228B22", "#FF8C00")
    pt_colors <- grDevices::adjustcolor(pt_colors_base[as.character(plot_group)], alpha.f = pt_alpha)
  } else {
    plot_group <- NULL
    pt_colors <- grDevices::rgb(0.2, 0.2, 0.2, pt_alpha)
    unique_groups <- "Cells"
  }
  
  n_plots <- length(valid_features)
  if (n_plots > 1) {
    plot_cols <- min(ncol, n_plots)
    plot_rows <- ceiling(n_plots / plot_cols)
    old_par <- graphics::par(mfrow = c(plot_rows, plot_cols))
    on.exit(graphics::par(old_par))
  }

  for (gene in valid_features) {
    plot_y <- as.numeric(expr_mat[gene, valid_cells])
    plot_y_jittered <- pmax(0, jitter(plot_y, amount = jitter_amount))
    ymax <- max(plot_y_jittered, na.rm=TRUE) * m
    
    graphics::plot(plot_x, plot_y_jittered, pch = 16, cex = 0.6, col = pt_colors,
                   xlab = "Gradient Tocky Time (Angle)", ylab = paste(gene, "Expression (Log)"),
                   main = gene, las = 1, xlim = c(0, 90), ylim = c(0, max(ymax, 1)))
    
    graphics::rect(0, -10, 30, ymax*2,  col = grDevices::rgb(0, 0, 1, 0.05), border = NA)
    graphics::rect(30, -10, 60, ymax*2, col = grDevices::rgb(0.5, 0, 0.5, 0.05), border = NA)
    graphics::rect(60, -10, 90, ymax*2, col = grDevices::rgb(1, 0, 0, 0.05), border = NA)
    
    draw_curve <- function(x_vec, y_vec, col_hex, label_txt) {
      if (length(x_vec) < 10) return(NULL)
      fit <- tryCatch(stats::loess(y_vec ~ x_vec, span = span), error = function(e) NULL)
      if (is.null(fit)) return(NULL)
      x_seq <- seq(min(x_vec), max(x_vec), length.out = 100)
      y_pred <- pmax(0, stats::predict(fit, newdata = data.frame(x_vec = x_seq)))
      graphics::lines(x_seq, y_pred, col = col_hex, lwd = 3)
      return(paste0(label_txt, " (", round(sum(y_vec > 0)/length(y_vec)*100), "%)"))
    }
    
    if (is.null(plot_group)) {
      lab <- draw_curve(plot_x, plot_y, "black", "Cells")
      if(!is.null(lab)) graphics::legend("topleft", legend=lab, col="black", lwd=3, bty="n", cex=0.8)
    } else {
      legend_text <- c(); legend_cols_vec <- c()
      for (grp in unique_groups) {
        idx <- plot_group == grp
        lab <- draw_curve(plot_x[idx], plot_y[idx], pt_colors_base[grp], grp)
        if (!is.null(lab)) { legend_text <- c(legend_text, lab); legend_cols_vec <- c(legend_cols_vec, pt_colors_base[grp]) }
      }
      if(length(legend_text) > 0) graphics::legend("topleft", legend = legend_text, col = legend_cols_vec, lwd = 3, bty = "n", cex = 0.8)
    }
    graphics::grid()
  }
  
  return(invisible(NULL))
})


#' Plot Tocky-Time Heatmap
#'
#' @description Generates a binned and smoothed heatmap ordered by Tocky Time.
#'
#' @param object An mCanonicalTockyObj.
#' @param seurat_obj Original Seurat object.
#' @param genes Character vector of genes to plot.
#' @param n_bins Integer. Number of bins (default 100).
#' @param ordering_method Character. "peak" (sort by max time) or "cluster".
#' @param span Numeric. The smoothing parameter (alpha) for LOESS smoothing across Tocky Time bins. Default is 0.5.
#' @param scale_rows Logical. Whether to perform row-wise Z-score scaling on the smoothed expression matrix. Default is TRUE.
#' @param ... Additional arguments passed to methods.
#' @export
#' @importFrom stats loess predict hclust dist sd complete.cases
#' @importFrom pheatmap pheatmap
#' @importFrom grDevices colorRampPalette
#' @importFrom Seurat GetAssayData
setGeneric("mPlotTockyHeatmap", function(object, seurat_obj, genes, ...) standardGeneric("mPlotTockyHeatmap"))

#' @rdname mPlotTockyHeatmap
#' @export
setMethod("mPlotTockyHeatmap", signature(object = "mCanonicalTockyObj"),
          function(object, seurat_obj, genes, n_bins = 100, ordering_method = c("peak", "cluster"), span = 0.5, scale_rows = TRUE) {
            
  if (!requireNamespace("pheatmap", quietly = TRUE)) stop("Please install 'pheatmap'.")
  ordering_method <- match.arg(ordering_method)
  
  if (!"angle" %in% colnames(object@TockyData)) stop("Run mGradientTockySeq() first.")
  tocky_time <- object@TockyData$angle
  names(tocky_time) <- rownames(object@TockyData)
  
  tryCatch({ counts <- Seurat::GetAssayData(seurat_obj, layer = "data") },
           error = function(e) { counts <- Seurat::GetAssayData(seurat_obj, slot = "data") })
  
  common_cells <- intersect(colnames(counts), names(tocky_time))
  cells_use <- common_cells[!is.na(tocky_time[common_cells])]
  
  if (length(cells_use) < 10) stop("Too few valid cells.")
  
  valid_genes <- intersect(genes, rownames(counts))
  expr_mat <- as.matrix(counts[valid_genes, cells_use])
  time_vec <- tocky_time[cells_use]
  
  breaks <- seq(0, 90, length.out = n_bins + 1)
  bin_centers <- (breaks[1:n_bins] + breaks[2:(n_bins+1)]) / 2
  cell_bins <- cut(time_vec, breaks = breaks, labels = FALSE, include.lowest = TRUE)
  
  binned_mat <- matrix(NA, nrow = nrow(expr_mat), ncol = n_bins)
  rownames(binned_mat) <- rownames(expr_mat)
  colnames(binned_mat) <- round(bin_centers, 1)
  
  for (i in 1:n_bins) {
    cells_in_bin <- which(cell_bins == i)
    if (length(cells_in_bin) == 1) binned_mat[, i] <- expr_mat[, cells_in_bin]
    else if (length(cells_in_bin) > 1) binned_mat[, i] <- rowMeans(expr_mat[, cells_in_bin])
  }
  
  smoothed_mat <- binned_mat
  if (span > 0) {
    x_grid <- 1:n_bins
    for (i in 1:nrow(binned_mat)) {
      y_vals <- binned_mat[i, ]
      valid_bins <- !is.na(y_vals)
      if (sum(valid_bins) >= 5) {
        fit <- tryCatch(stats::loess(y_vals[valid_bins] ~ x_grid[valid_bins], span = span), error = function(e) NULL)
        if (!is.null(fit)) smoothed_mat[i, ] <- stats::predict(fit, newdata = x_grid)
      }
    }
  }
  
  smoothed_mat <- smoothed_mat[stats::complete.cases(smoothed_mat), , drop = FALSE]
  if (nrow(smoothed_mat) < 2) return(NULL)

  plot_mat <- smoothed_mat
  if (scale_rows) {
    row_means <- rowMeans(plot_mat)
    row_sds   <- apply(plot_mat, 1, stats::sd)
    plot_mat  <- (plot_mat - row_means) / replace(row_sds, row_sds==0, 1)
    plot_mat[plot_mat > 3]  <- 3
    plot_mat[plot_mat < -3] <- -3
  }
  
  if (ordering_method == "peak") {
    peak_col <- apply(plot_mat, 1, which.max)
    plot_mat <- plot_mat[order(peak_col, decreasing = FALSE), ]
    cluster_rows_flag <- FALSE
    main_title <- "Tocky Cascade (Ordered by Peak Time)"
  } else {
    cluster_rows_flag <- TRUE
    main_title <- "Tocky Heatmap (Hierarchical Clustering)"
  }

  heat_colors <- grDevices::colorRampPalette(c("#7b3294", "#f7f7f7", "#d95f0e"))(100)
  pheatmap::pheatmap(plot_mat, cluster_rows = cluster_rows_flag, cluster_cols = FALSE,
                     show_colnames = FALSE, main = main_title, color = heat_colors,
                     border_color = NA, fontsize_row = 8)
})
