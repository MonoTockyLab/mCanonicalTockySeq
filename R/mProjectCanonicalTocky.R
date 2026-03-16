#' Project New ScRNA-seq Data into an Established Canonical Tocky Space
#'
#' @description
#' Performs Redundancy Analysis (RDA) on a new expression matrix (e.g., human orthologs)
#' using the exact explanatory variable matrix (Z) and internal mappings from a reference
#' Tocky model (e.g., mouse).
#'
#' @param X_mat A Matrix (dense or sparse) of expression data for the new dataset.
#' @param ref_obj The reference mCanonicalTockyObj.
#' @param seurat_obj The Seurat object corresponding to the cells in X_mat.
#' @param scale_Z Logical. Whether to scale the Z matrix (default TRUE).
#' @param ... Additional arguments passed to methods.

#'
#' @return An mCanonicalTockyObj containing the projected cells.
#' @importFrom irlba irlba
#' @importFrom methods new
#' @export
setGeneric("mProjectCanonicalTocky", function(X_mat, ref_obj, seurat_obj, ...) standardGeneric("mProjectCanonicalTocky"))

#' @rdname mProjectCanonicalTocky
#' @export
setMethod("mProjectCanonicalTocky", signature(X_mat = "Matrix"),
          function(X_mat, ref_obj, seurat_obj, scale_Z = TRUE) {
            
  if (!requireNamespace("irlba", quietly = TRUE)) stop("Package 'irlba' is required.")

  Z_mat <- ref_obj@Z
  if (is.null(Z_mat)) stop("Reference Z matrix is missing from the provided object.")

  common_genes <- intersect(rownames(X_mat), rownames(Z_mat))
  if (length(common_genes) < 10) {
    stop("Too few overlapping genes between the new dataset and the reference Z matrix.")
  }
  
  cat(sprintf("Aligning matrices: Projecting %d cells across %d shared orthologous genes...\n",
              ncol(X_mat), length(common_genes)))
  
  X <- as.matrix(X_mat[common_genes, , drop = FALSE])
  Z <- Z_mat[common_genes, , drop = FALSE]
  
  Z_scaled <- if(scale_Z) scale(Z, center = FALSE, scale = TRUE) else Z
  S <- scale(X)

  inv_ZtZ <- solve(t(Z_scaled) %*% Z_scaled)
  projection_coeff <- Z_scaled %*% inv_ZtZ

  ZtS <- t(Z_scaled) %*% S
  S_star <- projection_coeff %*% ZtS
  
  cat("Calculating Singular Value Decomposition (SVD) on constrained projection space...\n")
  k_comps <- min(ncol(Z), nrow(S_star), ncol(S_star))
  svd_decomp <- irlba::irlba(S_star, nv = k_comps, nu = k_comps)
  
  U <- svd_decomp$u
  D_alpha <- diag(svd_decomp$d, nrow = length(svd_decomp$d))
  V <- svd_decomp$v

  gene_expression_scores <- (U %*% D_alpha)
  fitted_cell_scores <- scale(V %*% D_alpha, center = TRUE, scale = FALSE)
  cell_scores <- scale(crossprod(S, U), center = TRUE, scale = FALSE)
  biplot_values <- solve(t(Z) %*% Z) %*% t(Z) %*% gene_expression_scores
  
  fitted_cell_scores <- as.data.frame(fitted_cell_scores)
  cell_scores        <- as.data.frame(cell_scores)
  gene_expression_scores <- as.data.frame(gene_expression_scores)
  
  rownames(gene_expression_scores) <- rownames(X)
  rownames(fitted_cell_scores) <- colnames(X)
  rownames(cell_scores)        <- colnames(X)
  
  axis_names <- paste0("Axis", 1:k_comps)
  colnames(gene_expression_scores) <- axis_names
  colnames(fitted_cell_scores) <- axis_names
  colnames(cell_scores)        <- axis_names
  
  biplot_limit <- min(ncol(Z), ncol(biplot_values))
  biplot_values <- as.data.frame(biplot_values[, 1:biplot_limit, drop = FALSE])
  colnames(biplot_values) <- paste0("Axis", 1:ncol(biplot_values))
  
  rownames(biplot_values) <- unname(ref_obj@marker_info$mapping[colnames(Z)])
  
  common_cells <- intersect(colnames(X), colnames(seurat_obj))
  meta_df <- if (length(common_cells) > 0) seurat_obj@meta.data[common_cells, , drop=FALSE] else data.frame(row.names = colnames(X))
  
  human_tocky_obj <- new("mCanonicalTockyObj",
                         X = X,
                         Z = Z,
                         metadata = meta_df,
                         genes = rownames(X),
                         expression_scores = gene_expression_scores,
                         fitted_cell_scores = fitted_cell_scores,
                         cell_scores = cell_scores,
                         biplot = biplot_values,
                         marker_info = ref_obj@marker_info,
                         TockyData = data.frame()
  )
  
  cat("Successfully generated Projected mCanonicalTockyObj.\n")
  return(human_tocky_obj)
})
