# Copyright 2026 Masahiro Ono
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#' S4 Class to hold mCanonicalTockySeq Results
#' @slot X Matrix. The single-cell expression matrix (Genes x Cells), accepts either dense or sparse matrices.
#' @slot Z matrix. The reference landmark matrix (Genes x Vectors).
#' @slot metadata data.frame. Cell-level metadata.
#' @slot genes character. The features used for the RDA.
#' @slot expression_scores data.frame. Gene loadings.
#' @slot fitted_cell_scores data.frame. Fitted cell scores.
#' @slot cell_scores data.frame. Projected cell scores.
#' @slot biplot data.frame. Constraint vectors.
#' @slot marker_info list. Provenance of landmarks.
#' @slot TockyData data.frame. Gradient mapping and lineage trajectory results.
#' @export
setClass("mCanonicalTockyObj",
         slots = list(
           X = "Matrix",
           Z = "matrix",
           metadata = "data.frame",
           genes = "character",
           expression_scores = "data.frame",
           fitted_cell_scores = "data.frame",
           cell_scores = "data.frame",
           biplot = "data.frame",
           marker_info = "list",
           TockyData = "data.frame" 
         ))

#' Show method for mCanonicalTockyObj
#'
#' @param object An mCanonicalTockyObj
#' @importFrom methods show
#' @export
setMethod("show", "mCanonicalTockyObj", function(object) {
  
  cat("An object of class \"mCanonicalTockyObj\"\n")
  cat("=========================================\n")
  
  num_genes <- nrow(object@X)
  num_cells <- ncol(object@X)
  cat(sprintf("Single-cell data: %d genes across %d cells\n", num_genes, num_cells))

  num_vectors <- ncol(object@Z)
  cat(sprintf("Reference Landmarks (Z): %d vectors\n", num_vectors))
  if (length(object@marker_info$mapping) > 0) {
    cat("  - Groups:", paste(names(object@marker_info$mapping), collapse = ", "), "\n")
  }
  
  cat("-----------------------------------------\n")
  
  cat("Tocky Manifold Anchors:\n")
  landmarks <- object@marker_info$tocky_landmarks
  if (!is.null(landmarks) && length(landmarks) > 0) {
    cat(sprintf("  - Blue (New): %s\n", landmarks$B))
    cat(sprintf("  - Blue-Red (Persistent): %s\n", landmarks$BR))
    cat(sprintf("  - Red (Arrested): %s\n", landmarks$R))
  } else {
    cat("  - (Not defined in marker_info)\n")
  }
  
  lin_idents <- object@marker_info$lineage_idents
  lin_names <- object@marker_info$lineage_names
  
  if (!is.null(lin_idents) && !is.null(lin_names)) {
    cat("Lineage Endpoints:\n")
    for (i in seq_along(lin_idents)) {
      cat(sprintf("  - %s (Mapped from cluster: %s)\n", lin_names[i], lin_idents[i]))
    }
  }
  
  cat("-----------------------------------------\n")
  
  if (nrow(object@TockyData) > 0) {
    cat(sprintf("Gradient Mapping (Slerp): Completed (%d valid cells mapped)\n", nrow(object@TockyData)))
  } else {
    cat("Gradient Mapping (Slerp): Not yet performed. Run mGradientTockySeq().\n")
  }
  
  cat("=========================================\n")
})
