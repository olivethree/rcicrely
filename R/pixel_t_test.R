#' Vectorised pixel-wise Welch t-test
#'
#' @description
#' Computes Welch's t (unequal variances, Student's is **not**
#' appropriate, conditions can differ in both N and variance) per
#' pixel between two condition signal matrices. Fully vectorised: no
#' per-row `apply()`. Returns a length-`n_pixels` numeric vector.
#'
#' Use this when you want the raw t-map (e.g. for custom plotting or
#' threshold sensitivity analysis) without the cluster-level
#' inference machinery of [rel_cluster_test()].
#'
#' @section Reading the result:
#' Numeric vector, one t-value per pixel. Pixels with zero variance
#' in both conditions get `0` rather than `NaN` so cluster utilities
#' don't have to special-case them.
#'
#' @section Reliability metrics expect raw masks:
#' Welch t is variance-based and sensitive to scaling. See
#' `vignette("tutorial", package = "rcicrely")` chapter 3.
#'
#' @param signal_matrix_a,signal_matrix_b Pixels x participants,
#'   base-subtracted. Row counts must match.
#' @return Numeric vector of length `nrow(signal_matrix_a)`.
#' @seealso [rel_cluster_test()]
#' @export
pixel_t_test <- function(signal_matrix_a, signal_matrix_b) {
  validate_two_signal_matrices(signal_matrix_a, signal_matrix_b)
  n_a <- ncol(signal_matrix_a)
  n_b <- ncol(signal_matrix_b)

  mean_a <- rowMeans(signal_matrix_a)
  mean_b <- rowMeans(signal_matrix_b)

  # row variances, vectorised
  var_a <- rowSums((signal_matrix_a - mean_a)^2) / (n_a - 1L)
  var_b <- rowSums((signal_matrix_b - mean_b)^2) / (n_b - 1L)

  se <- sqrt(var_a / n_a + var_b / n_b)
  t_vec <- (mean_a - mean_b) / se

  # pixels with zero variance in both conditions -> t is 0/0; report NaN
  # but don't let a few pixels contaminate downstream (cluster_utils
  # treats NaN as below-threshold).
  t_vec[!is.finite(t_vec)] <- 0
  t_vec
}
