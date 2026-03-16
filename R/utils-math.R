#' @importFrom stats optimize
NULL

#' Internal Slerp Interpolation
#' @keywords internal
#' @noRd
.fast_slerp <- function(v1, v2, h, omega) {
  if (abs(omega) < 1e-6) return(v1)
  o <- 1/sin(omega)
  (sin((1 - h) * omega) * o) * v1 + (sin(h * omega) * o) * v2
}

#' Internal Euclidean Norm
#' @keywords internal
#' @noRd
.get_norm <- function(x) sqrt(sum(x^2))

#' Find best T for a single cell
#' @keywords internal
#' @noRd
.get_best_t_for_cell <- function(cell_vec, B, BR, R, omega1, omega2) {
  
  obj_fun <- function(h, v_start, v_end, omega_val) {
    ref_vec <- .fast_slerp(v_start, v_end, h, omega_val)
    sim <- sum(cell_vec * ref_vec) / .get_norm(ref_vec)
    return(-sim)
  }
  
  opt1 <- stats::optimize(obj_fun, interval = c(0, 1),
                          v_start = B, v_end = BR, omega_val = omega1)
  
  opt2 <- stats::optimize(obj_fun, interval = c(0, 1),
                          v_start = BR, v_end = R, omega_val = omega2)
  
  if (opt1$objective <= opt2$objective) {
    return(c(t = opt1$minimum * 0.5, sim = -opt1$objective))
  } else {
    return(c(t = 0.5 + opt2$minimum * 0.5, sim = -opt2$objective))
  }
}
