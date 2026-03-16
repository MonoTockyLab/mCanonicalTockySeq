#' Run TradeSeq Differential Dynamics Analysis (Optimized for 2 Lineages)
#'
#' @description
#' Evaluates differential gene expression across exactly two developmental tubes.
#' Employs aggressive proportional thresholding to drop uninformative genes,
#' bypassing severe RAM and processing bottlenecks in tradeSeq.
#'
#' @param object An mCanonicalTockyObj.
#' @param seurat_obj The Seurat object to fetch counts from.
#' @param tubes Character vector of exactly two genes defining the tubes to compare (e.g., c("Cd4", "Cd8a")). If NULL, auto-detects if exactly two exist.
#' @param min_pct Numeric. Minimum fraction of cells a gene must be expressed in (default 0.05, i.e., 5%).
#' @param n_knots Integer. Number of knots for the GAM (default 5).
#' @param n_cores Integer. Number of cores for parallel processing (default 1).
#' @param p_adjust_method Character. Correction method: "fdr" (default), "bonferroni", "holm", etc.
#' @param ... Additional arguments passed to methods. 
#'
#' @export
#' @importFrom stats p.adjust
#' @importFrom Seurat GetAssayData
#' @importFrom tradeSeq fitGAM
#' @importFrom Matrix rowSums
#' @importFrom BiocParallel SnowParam
setGeneric("mRunTradeSeqTocky", function(object, seurat_obj, ...) standardGeneric("mRunTradeSeqTocky"))

#' @rdname mRunTradeSeqTocky
#' @export
setMethod("mRunTradeSeqTocky", signature(object = "mCanonicalTockyObj"),
          function(object, seurat_obj, tubes = NULL, min_pct = 0.05, n_knots = 5, n_cores = 1, p_adjust_method = "fdr") {
            
  if (!requireNamespace("tradeSeq", quietly = TRUE)) stop("Package 'tradeSeq' is required.")
  if (!requireNamespace("BiocParallel", quietly = TRUE)) stop("Package 'BiocParallel' is required.")
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Package 'Matrix' is required.")
  
  df <- object@TockyData

  if (!is.null(tubes)) {
    tube_cols <- paste0("Tube_", tubes)
    if (!all(tube_cols %in% colnames(df))) stop("One or both specified tubes not found in object.")
  } else {
    tube_cols <- grep("^Tube_", colnames(df), value = TRUE)
  }
  
  if (length(tube_cols) != 2) {
    stop(sprintf("tradeSeq optimization requires exactly 2 tubes. Found %d. Please specify using the 'tubes' argument (e.g., tubes = c('Runx3', 'Zbtb7b')).", length(tube_cols)))
  }
  
  if (!"angle" %in% colnames(df)) stop("Tocky angle missing.")

  in_any_tube <- rowSums(df[, tube_cols, drop = FALSE], na.rm = TRUE) > 0
  valid_df <- df[!is.na(df$angle) & in_any_tube, ]
  valid_cells <- rownames(valid_df)
  
  valid_cells <- intersect(valid_cells, colnames(seurat_obj))
  if (length(valid_cells) < 10) stop("Fewer than 10 valid cells found. Cannot run GAM.")
  
  cat(sprintf("Found %d cells actively participating in %s vs %s...\n",
              length(valid_cells), tube_cols[1], tube_cols[2]))

  counts_mat <- tryCatch({
    Seurat::GetAssayData(seurat_obj, layer = "counts")[, valid_cells, drop = FALSE]
  }, error = function(e) {
    Seurat::GetAssayData(seurat_obj, slot = "counts")[, valid_cells, drop = FALSE]
  })

  min_cells <- ceiling(length(valid_cells) * min_pct)
  cat(sprintf("Filtering out genes expressed in fewer than %d cells (%.1f%% threshold)...\n", min_cells, min_pct * 100))
  
  gene_cell_counts <- Matrix::rowSums(counts_mat > 0)
  valid_genes <- names(gene_cell_counts[gene_cell_counts >= min_cells])
  
  if (length(valid_genes) == 0) stop("No genes passed the min_pct expression threshold.")
  y_counts <- as.matrix(counts_mat[valid_genes, ])
  t_time   <- valid_df[valid_cells, "angle"]
  
  cat(sprintf("Running tradeSeq GAM on %d robustly expressed genes across %d cells...\n",
              nrow(y_counts), ncol(y_counts)))

  pseudotime_mat <- matrix(rep(t_time, 2), ncol = 2)
  rownames(pseudotime_mat) <- valid_cells
  colnames(pseudotime_mat) <- gsub("^Tube_", "", tube_cols)
  
  weights_mat <- as.matrix(valid_df[valid_cells, tube_cols])
  mode(weights_mat) <- "numeric"
  rownames(weights_mat) <- valid_cells
  colnames(weights_mat) <- gsub("^Tube_", "", tube_cols)
  
  bpparam <- if (n_cores > 1) {
    if (.Platform$OS.type == "windows") BiocParallel::SnowParam(workers = n_cores)
    else BiocParallel::MulticoreParam(workers = n_cores)
  } else {
    BiocParallel::SerialParam()
  }
  
  gam_fit <- tradeSeq::fitGAM(counts = y_counts,
                              pseudotime = pseudotime_mat,
                              cellWeights = weights_mat,
                              nknots = n_knots,
                              parallel = (n_cores > 1),
                              BPPARAM = bpparam,
                              verbose = TRUE)
  
  cat("Testing for differential expression patterns across the 2 tubes (patternTest)...\n")
  cond_res <- tradeSeq::patternTest(gam_fit, global = TRUE)
  
  cat(sprintf("Adjusting p-values using the '%s' method...\n", p_adjust_method))
  cond_res$padj <- stats::p.adjust(cond_res$pvalue, method = p_adjust_method)
  
  cond_res <- cond_res[order(cond_res$padj), ]
  
    tube1_name <- colnames(weights_mat)[1]
    tube2_name <- colnames(weights_mat)[2]
    cat(sprintf("Running Seurat overall comparison: %s vs %s...\n", tube1_name, tube2_name))
    
    cell_tube_idx <- max.col(weights_mat, ties.method = "first")
    dominant_tubes <- colnames(weights_mat)[cell_tube_idx]
    names(dominant_tubes) <- valid_cells
    
    seurat_subset <- subset(seurat_obj, cells = valid_cells)
    Seurat::Idents(seurat_subset) <- dominant_tubes[colnames(seurat_subset)]
    
    seurat_markers <- Seurat::FindMarkers(seurat_subset,
                                          ident.1 = tube1_name,
                                          ident.2 = tube2_name,
                                          features = rownames(cond_res),
                                          logfc.threshold = 0,
                                          min.pct = 0,
                                          verbose = FALSE)
    
    cat("Classifying trajectory patterns...\n")
    
    cond_res$Overall_Log2FC <- seurat_markers[rownames(cond_res), "avg_log2FC"]
    
    cond_res$Pattern_Class <- "Not Dynamic"
    sig_genes <- which(cond_res$padj < 0.05)
    
    cond_res$Pattern_Class[sig_genes] <- ifelse(
      cond_res$Overall_Log2FC[sig_genes] > 0.25, paste("Overall High:", tube1_name),
      ifelse(cond_res$Overall_Log2FC[sig_genes] < -0.25, paste("Overall High:", tube2_name),
             "Complex Dynamics")
    )
    
    cond_res <- cond_res[order(cond_res$padj), ]
  
  return(list(results = cond_res, gam_fit = gam_fit))
})
