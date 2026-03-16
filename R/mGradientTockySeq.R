#' Gradient Tocky Sequence (Slerp Model)
#'
#' @description
#' Projects cells onto a piecewise Slerp manifold defined by the landmarks B -> BR -> R.
#' Automatically extracts the landmark vectors from the provided `mCanonicalTockyObj`.
#' Cells with negative correlations to all three landmarks are classified as "Timer Negative" (NA),
#' as they lie in the ordination space opposite to the Tocky trajectory.
#'
#' @param object An object of class \code{mCanonicalTockyObj}.
#' @param filter_negative Logical. If TRUE, assigns NA to cells with negative projection to all three landmarks.
#' @param intensity_q Numeric. Quantile threshold (0 to 1) for normalized intensity.
#'        Cells falling below this threshold are treated as non-participating and assigned NA for angle. Default is 0.75.
#' @param ... Additional arguments.
#'
#' @return An updated \code{mCanonicalTockyObj} with the Slerp mapping results stored in the \code{@TockyData} slot.
#' @export
#' @importFrom stats optimize quantile

setGeneric("mGradientTockySeq", function(object, ...) standardGeneric("mGradientTockySeq"))

#' @rdname mGradientTockySeq
#' @export
setMethod("mGradientTockySeq", signature(object = "mCanonicalTockyObj"),
          function(object, filter_negative = TRUE, intensity_q = 0.75) {
            
  landmarks <- object@marker_info$tocky_landmarks
  if (is.null(landmarks$B) || is.null(landmarks$BR) || is.null(landmarks$R)) {
    stop("Tocky landmarks (B, BR, R) not found in the object. Ensure it was created with the updated mCanonicalTockySeq.")
  }
  
  B  <- as.numeric(object@biplot[landmarks$B, ])
  BR <- as.numeric(object@biplot[landmarks$BR, ])
  R  <- as.numeric(object@biplot[landmarks$R, ])
  
  n <- length(B)
  cell_mat <- as.matrix(object@cell_scores[, 1:n])
  n_cells <- nrow(cell_mat)
  
  is_timer_neg <- rep(FALSE, n_cells)
  if (filter_negative) {
    is_timer_neg <- (cell_mat %*% B < 0) & (cell_mat %*% BR < 0) & (cell_mat %*% R < 0)
  }

  norm_B <- .get_norm(B); unit_B <- B / norm_B
  norm_BR <- .get_norm(BR); unit_BR <- BR / norm_BR
  norm_R <- .get_norm(R); unit_R <- R / norm_R
  
  omega1 <- acos(pmax(pmin(sum(unit_B * unit_BR), 1), -1))
  omega2 <- acos(pmax(pmin(sum(unit_BR * unit_R), 1), -1))

  t_values <- rep(NA, n_cells)
  sim_values <- rep(NA, n_cells)
  valid_idx <- which(!is_timer_neg)
  
  if (length(valid_idx) > 0) {
    subset_res <- apply(cell_mat[valid_idx, , drop=FALSE], 1, .get_best_t_for_cell,
                        B = B, BR = BR, R = R, omega1 = omega1, omega2 = omega2)
    
    t_values[valid_idx] <- subset_res["t", ]
    sim_values[valid_idx] <- subset_res["sim", ]
  }

  cell_norms <- sqrt(rowSums(cell_mat^2))
  expected_mags <- sapply(t_values, function(t) {
    if (is.na(t)) return(NA)
    ref <- if(t <= 0.5) .fast_slerp(B, BR, t*2, omega1) else .fast_slerp(BR, R, (t-0.5)*2, omega2)
    .get_norm(ref)
  })
  
  norm_intensities <- cell_norms / expected_mags

  num_low_intensity <- 0
  if (!is.null(intensity_q) && intensity_q > 0) {
    thresh <- stats::quantile(norm_intensities, probs = intensity_q, na.rm = TRUE)
    low_int_mask <- !is.na(norm_intensities) & (norm_intensities < thresh)

    t_values[low_int_mask] <- NA
    sim_values[low_int_mask] <- NA
    num_low_intensity <- sum(low_int_mask)
  }

  valid_mapped <- sum(!is.na(t_values))
  cat(sprintf("Mapped %d cells to the Tocky manifold.\n", valid_mapped))
  cat(sprintf("Filtered %d Timer Negative cells and %d Low-Intensity cells (q=%.2f).\n",
              sum(is_timer_neg), num_low_intensity, intensity_q))

  out_df <- data.frame(
    angle = t_values * 90,
    intensity = cell_norms,
    norm_intensity = norm_intensities,
    similarity = sim_values,
    row.names = rownames(cell_mat)
  )

  object@TockyData <- out_df
    
  return(object)
})
