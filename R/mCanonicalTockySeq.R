#' Canonical Redundancy Analysis for Tocky Differentiation
#'
#' @param object A Seurat object.
#' @param temporal_col Character. Metadata column for temporal anchors.
#' @param b_ident Character. Identity for the Blue (New) landmark (default "B").
#' @param br_ident Character. Identity for the Blue-Red (Persistent) landmark (default "BR").
#' @param r_ident Character. Identity for the Red (Arrested) landmark (default "R").
#' @param lineage_col Character. Metadata column for lineage endpoints.
#' @param lineage_idents Character vector of identities in Seurat object for lineage endpoints.
#' @param lineage_names Optional character vector. Biological names for the lineages (e.g., c("CD4", "CD8")).
#' @param top_n Integer. Number of marker genes to extract per group (default 100).
#' @param custom_genes Optional character vector to skip marker calculation.
#' @param layer Character. The Seurat layer to use (default "data").
#' @param scale_Z Logical. Whether to scale the Z matrix (default TRUE).
#'
#' @return An object of class \code{mCanonicalTockyObj}.
#' @importFrom irlba irlba
#' @importFrom utils head
#' @importFrom methods new
#' @importFrom Seurat GetAssayData FindMarkers Idents
#' @importFrom Matrix rowMeans
#' @export
mCanonicalTockySeq <- function(object,
                               temporal_col, b_ident = "B", br_ident = "BR", r_ident = "R",
                               lineage_col, lineage_idents, lineage_names = NULL, # <-- NEW PARAMETER
                               top_n = 100, custom_genes = NULL,
                               layer = "data", scale_Z = TRUE){
    
    if (!requireNamespace("irlba", quietly = TRUE)) {
      stop("Package 'irlba' is required for fast SVD. Please install it.")
    }

    if (is.null(lineage_names)) {
      lineage_names <- lineage_idents
    } else if (length(lineage_names) != length(lineage_idents)) {
      stop("Error: 'lineage_names' must be the exact same length as 'lineage_idents'.")
    }

    temporal_idents <- c(b_ident, br_ident, r_ident)

    expr_mat <- tryCatch({
      Seurat::GetAssayData(object, layer = layer)
    }, error = function(e) {
      Seurat::GetAssayData(object, slot = layer)
    })
    if (is.null(custom_genes)) {
      cat("Calculating markers for", length(temporal_idents) + length(lineage_idents), "groups...\n")
      all_markers <- list()
        
      Seurat::Idents(object) <- temporal_col
      for (id in temporal_idents) {
        cat("  Finding markers for temporal group:", id, "\n")
        m <- Seurat::FindMarkers(object, ident.1 = id, only.pos = TRUE, verbose = FALSE)
        all_markers[[id]] <- rownames(head(m, top_n))
      }
        
      Seurat::Idents(object) <- lineage_col
      for (id in lineage_idents) {
        cat("  Finding markers for lineage group:", id, "\n")
        m <- Seurat::FindMarkers(object, ident.1 = id, only.pos = TRUE, verbose = FALSE)
        all_markers[[id]] <- rownames(head(m, top_n))
      }
        
      sig_genes <- unique(unlist(all_markers))
      cat("Identified", length(sig_genes), "unique marker genes.\n")
    } else {
      sig_genes <- intersect(custom_genes, rownames(expr_mat))
      cat("Using", length(sig_genes), "user-provided genes.\n")
    }
      
    all_groups_fetch <- c(temporal_idents, lineage_idents)
    all_groups_display <- c(temporal_idents, lineage_names)
    
    Z_mat <- matrix(0, nrow = length(sig_genes), ncol = length(all_groups_fetch))
    rownames(Z_mat) <- sig_genes
    colnames(Z_mat) <- all_groups_display
      
    get_avg <- function(col_name, val) {
      cells <- colnames(object)[object[[col_name]] == val]
      if (length(cells) == 0) stop(paste("No cells found for", val, "in", col_name))
      Matrix::rowMeans(expr_mat[sig_genes, cells, drop = FALSE])
    }
      
    cat("Populating Z Matrix averages...\n")
    for (i in seq_along(temporal_idents)) { Z_mat[, temporal_idents[i]] <- get_avg(temporal_col, temporal_idents[i]) }
    for (i in seq_along(lineage_idents))  { Z_mat[, lineage_names[i]]   <- get_avg(lineage_col, lineage_idents[i]) }
      
    X_mat <- expr_mat[sig_genes, ]
      
    standard_names <- c(paste0("Temp_", 1:length(temporal_idents)),
                        paste0("Lin_", 1:length(lineage_names)))
    
    marker_mapping <- setNames(all_groups_display, standard_names)
    colnames(Z_mat) <- standard_names

    X <- X_mat
    Z <- Z_mat
    
    Z_scaled <- if(scale_Z) scale(Z, center = FALSE, scale = TRUE) else Z
    S <- scale(X)
    
    inv_ZtZ <- solve(t(Z_scaled) %*% Z_scaled)
    projection_coeff <- Z_scaled %*% inv_ZtZ
    ZtS <- t(Z_scaled) %*% S
    S_star <- projection_coeff %*% ZtS
    
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
    
    rownames(biplot_values) <- unname(marker_mapping[rownames(biplot_values)])
    
    tocky_obj <- new("mCanonicalTockyObj",
                     X = X_mat,
                     Z = Z_mat,
                     metadata = object@meta.data[colnames(X_mat), ],
                     genes = sig_genes,
                     expression_scores = gene_expression_scores,
                     fitted_cell_scores = fitted_cell_scores,
                     cell_scores = cell_scores,
                     biplot = biplot_values,
                     marker_info = list(mapping = marker_mapping,
                                        temporal_col = temporal_col,
                                        lineage_col = lineage_col,
                                        lineage_idents = lineage_idents,
                                        lineage_names = lineage_names,
                                        tocky_landmarks = list(B = b_ident, BR = br_ident, R = r_ident)),
                     TockyData = data.frame()
    )
    
    cat("mCanonicalTockyObj created successfully.\n")
    return(tocky_obj)
}
