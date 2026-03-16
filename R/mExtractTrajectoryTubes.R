#' Extract Independent Biological Trajectory Tubes
#'
#' @description
#' Models independent developmental trajectories for an arbitrary number of genes.
#' For each gene, it isolates the "Core" expressing cells to calculate an expression-weighted
#' barycenter and covariance matrix, bypassing the gravitational pull of ambient RNA noise.
#' It then draws a perfect geometric ellipse using the Chi-Square distribution. Cells falling
#' within the confidence interval of a gene's trajectory are marked as TRUE for that specific tube,
#' allowing for overlapping dual identities.
#'
#' @param object An mCanonicalTockyObj.
#' @param seurat_obj Original Seurat object.
#' @param features Character vector of genes to model trajectories for (e.g., c("Cd4", "Cd8a")).
#' @param ci Numeric. Confidence interval for the ellipses (default 0.90).
#' @param window_size Numeric. Width of the sliding angle slice (default 10).
#' @param span Numeric. LOESS smoothing span for the barycenters (default 0.4).
#' @param reset Logical. If TRUE (default), deletes previously created lineage tubes. If FALSE, appends new tubes and updates the composite Lineage_Identity.
#' @param ... Additional arguments passed to methods. 
#'
#' @export
#' @importFrom stats cov.wt mahalanobis quantile loess predict qchisq cov
#' @importFrom Seurat GetAssayData
#' @importFrom MASS kde2d
setGeneric("mExtractTrajectoryTubes", function(object, ...) standardGeneric("mExtractTrajectoryTubes"))

#' @rdname mExtractTrajectoryTubes
#' @export
setMethod("mExtractTrajectoryTubes", signature(object = "mCanonicalTockyObj"),
          function(object, seurat_obj, features, ci = 0.90, window_size = 10, span = 0.4, reset = TRUE) {
            
  if (reset) {
    cat("Resetting previously computed lineage data...\n")
    old_cols <- c(grep("^Tube_", colnames(object@TockyData), value = TRUE),
                  "Lineage_Identity", "Lineage_Bias", "In_Tube", "Distance_to_Path")
    for (col in old_cols) {
      if (col %in% colnames(object@TockyData)) object@TockyData[[col]] <- NULL
    }
  }
            
  df <- object@TockyData
  df$cell_barcode <- rownames(df)
  
  score_cols <- grep("_Score$", colnames(df), value = TRUE)
  if (length(score_cols) < 2) stop("Need at least two fate scores in TockyData.")
  
  dest1_col <- score_cols[1]
  dest2_col <- score_cols[2]
  
  if (!"angle" %in% colnames(df)) stop("Run mGradientTockySeq() first.")
  df <- df[!is.na(df$angle), ]
  
  tryCatch({ expr_mat <- Seurat::GetAssayData(seurat_obj, layer = "data") },
           error = function(e) { expr_mat <- Seurat::GetAssayData(seurat_obj, slot = "data") })
  
  valid_features <- intersect(features, rownames(expr_mat))
  if (length(valid_features) == 0) stop("None of the specified features were found.")
  
  slide_starts <- seq(0, 90 - window_size, by = 1)
  
  coords <- as.matrix(df[, c(dest1_col, dest2_col)])
  
  tube_results <- list()
  
  for (gene in valid_features) {
    cat(sprintf("Modeling independent trajectory tube for: %s...\n", gene))
    
    expr_vec <- as.numeric(expr_mat[gene, df$cell_barcode])
    gene_assignments <- matrix(NA, nrow = nrow(df), ncol = length(slide_starts))
    rownames(gene_assignments) <- df$cell_barcode
    
    for (i in seq_along(slide_starts)) {
      start <- slide_starts[i]
      window_idx <- which(df$angle >= start & df$angle < (start + window_size))
      
      if (length(window_idx) < 10) next
      
      sub_expr <- expr_vec[window_idx]
      sub_coords <- coords[window_idx, , drop = FALSE]
      
      positive_idx <- sub_expr > 0
      positive_expr <- sub_expr[positive_idx]
            
      if (length(positive_expr) < 10) next
            
      core_cutoff <- stats::quantile(positive_expr, probs = 0.50, na.rm = TRUE)
      core_idx <- sub_expr >= core_cutoff
            
      core_coords <- sub_coords[core_idx, , drop = FALSE]
            
      if (nrow(unique(core_coords)) < 3) next
            
      dens <- MASS::kde2d(core_coords[, 1], core_coords[, 2], n = 50)
      max_idx <- which(dens$z == max(dens$z), arr.ind = TRUE)
            
      peak_x <- dens$x[max_idx[1, 1]]
      peak_y <- dens$y[max_idx[1, 2]]
      barycenter <- c(peak_x, peak_y)
            
      cov_mat <- stats::cov(core_coords)
      diag(cov_mat) <- diag(cov_mat) + 1e-6
            
      mah_dist <- stats::mahalanobis(sub_coords, center = barycenter, cov = cov_mat)
      threshold <- stats::qchisq(ci, df = 2)
            
      gene_assignments[window_idx, i] <- mah_dist <= threshold
    }
    
    eval_counts <- rowSums(!is.na(gene_assignments))
    pass_counts <- rowSums(gene_assignments == TRUE, na.rm = TRUE)
    
    final_in_tube <- ifelse(eval_counts > 0, (pass_counts / eval_counts) > 0.5, FALSE)
    tube_results[[paste0("Tube_", gene)]] <- final_in_tube
  }
  
  tube_df <- as.data.frame(tube_results)
  rownames(tube_df) <- df$cell_barcode
  
  if (!reset) {
    existing_tubes <- grep("^Tube_", colnames(object@TockyData), value = TRUE)
    for (et in existing_tubes) {
      if (!et %in% colnames(tube_df)) {
        tube_df[[et]] <- object@TockyData[rownames(tube_df), et]
      }
    }
  }
  
  all_tube_cols <- grep("^Tube_", colnames(tube_df), value = TRUE)
  all_tube_genes <- gsub("^Tube_", "", all_tube_cols)
  
  identity_list <- lapply(rownames(tube_df), function(cell) {
    cell_vals <- tube_df[cell, all_tube_cols]
    cell_vals[is.na(cell_vals)] <- FALSE
    
    active_genes <- all_tube_genes[as.logical(cell_vals)]
    if (length(active_genes) == 0) return("Unassigned")
    return(paste0(active_genes, "+", collapse = ""))
  })
  
  tube_df$Lineage_Identity <- unlist(identity_list)
  
  matching_cells <- intersect(rownames(object@TockyData), rownames(tube_df))
  
  for (col in colnames(tube_df)) {
    object@TockyData[[col]] <- NA
    object@TockyData[matching_cells, col] <- tube_df[matching_cells, col]
  }
  
  cat("Trajectory extraction complete. Generated mathematically rigorous, overlapping multi-identity tubes.\n")
  print(table(object@TockyData$Lineage_Identity))
  
  return(object)
})



#' Calculate Fate Commitment Scores
#'
#' @description
#' Calculates the dot product projection of cells towards the lineage endpoints.
#' Automatically detects the lineage destinations established during mCanonicalTockySeq().
#'
#' @param object An mCanonicalTockyObj.
#' @param ... Additional arguments passed to methods. 
#' @export
setGeneric("mGetFateScores", function(object, ...) standardGeneric("mGetFateScores"))

#' @rdname mGetFateScores
#' @export
setMethod("mGetFateScores", signature(object = "mCanonicalTockyObj"),
          function(object) {
            
  destinations <- object@marker_info$lineage_names
  
  if (is.null(destinations)) {
    destinations <- object@marker_info$lineage_idents
  }
  
  if (is.null(destinations) || length(destinations) != 2) {
    stop("Could not auto-detect exactly two lineage destinations. Ensure mCanonicalTockySeq() was run with two lineage endpoints.")
  }
  
  cell_coords <- as.matrix(object@cell_scores)
  available_landmarks <- rownames(object@biplot)
  missing <- setdiff(destinations, available_landmarks)
  
  if (length(missing) > 0) {
    stop(paste("The following destinations were not found in the RDA biplot:",
               paste(missing, collapse = ", ")))
  }
  
  fate_matrix <- t(object@biplot[destinations, , drop = FALSE])
  fate_scores <- cell_coords %*% fate_matrix
  colnames(fate_scores) <- paste0(destinations, "_Score")
  
  if (nrow(object@TockyData) == 0) stop("Please run mGradientTockySeq() first to generate Tocky Time.")
  
  object@TockyData[colnames(fate_scores)] <- fate_scores
  
  cat(sprintf("Fate scores for %s & %s automatically added to @TockyData.\n",
              destinations[1], destinations[2]))
  
  return(object)
})

