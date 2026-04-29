#' Run every between-condition discriminability metric
#'
#' @description
#' Convenience orchestrator that runs [rel_cluster_test()] and
#' [rel_dissimilarity()] on two condition signal matrices and wraps
#' both results in an `rcicrely_report` for joint printing / plotting.
#'
#' Use this when you want both the spatial-pattern test and the
#' overall magnitude test in one call.
#'
#' @section Reading the result:
#' `$results$cluster_test` and `$results$dissimilarity`, one result
#' object each, fields as in the standalone functions.
#' `$method = "between"`.
#'
#' @section Reliability metrics expect raw masks:
#' Both downstream metrics are scale-sensitive: the cluster test
#' uses variance-based Welch t, and Euclidean distance in
#' `rel_dissimilarity()` is sensitive to any scaling. See
#' `vignette("tutorial", package = "rcicrely")` chapter 3.
#'
#' @param signal_matrix_a,signal_matrix_b Pixels x participants,
#'   base-subtracted. Row counts must match.
#' @param img_dims Integer `c(nrow, ncol)`. If `NULL`, inferred from
#'   the `img_dims` attribute on `signal_matrix_a`, or from
#'   `sqrt(n_pixels)` if that's a whole number.
#' @param n_permutations Passed to [rel_cluster_test()]. Default 2000.
#' @param n_boot Passed to [rel_dissimilarity()]. Default 2000.
#' @param cluster_threshold Passed to [rel_cluster_test()]. Default 2.0.
#' @param alpha Passed to [rel_cluster_test()]. Default 0.05.
#' @param ci_level Passed to [rel_dissimilarity()]. Default 0.95.
#' @param mask Optional logical vector of length
#'   `nrow(signal_matrix_a)`. Threaded through to both
#'   `rel_cluster_test()` (zero-out pattern, preserving 2D structure)
#'   and `rel_dissimilarity()` (drop pattern, both matrices subsetted
#'   identically). See `vignette("tutorial")` §6.6.
#' @param seed Optional integer.
#' @param progress Show `cli` progress bars.
#' @param acknowledge_scaling Logical. Forwarded to both
#'   [rel_cluster_test()] and [rel_dissimilarity()]; when `FALSE`
#'   (default), known-rendered inputs error.
#' @return Object of class `rcicrely_report` with `$results` =
#'   named list of two `rcicrely_*` objects (`cluster_test`,
#'   `dissimilarity`) and `$method = "between"`.
#' @seealso [rel_cluster_test()], [rel_dissimilarity()], [run_within()]
#' @export
run_between <- function(signal_matrix_a,
                        signal_matrix_b,
                        img_dims            = NULL,
                        n_permutations      = 2000L,
                        n_boot              = 2000L,
                        cluster_threshold   = 2.0,
                        alpha               = 0.05,
                        ci_level            = 0.95,
                        mask                = NULL,
                        seed                = NULL,
                        progress            = TRUE,
                        acknowledge_scaling = FALSE) {
  validate_two_signal_matrices(signal_matrix_a, signal_matrix_b)
  if (is.null(img_dims)) {
    img_dims <- attr(signal_matrix_a, "img_dims")
  }

  results <- list(
    cluster_test = rel_cluster_test(
      signal_matrix_a, signal_matrix_b,
      img_dims            = img_dims,
      n_permutations      = n_permutations,
      cluster_threshold   = cluster_threshold,
      alpha               = alpha,
      mask                = mask,
      seed                = seed,
      progress            = progress,
      acknowledge_scaling = acknowledge_scaling
    ),
    dissimilarity = rel_dissimilarity(
      signal_matrix_a, signal_matrix_b,
      n_boot              = n_boot,
      ci_level            = ci_level,
      mask                = mask,
      seed                = seed,
      progress            = progress,
      acknowledge_scaling = acknowledge_scaling
    )
  )
  new_rcicrely_report(
    results, method = "between",
    img_dims = results$cluster_test$img_dims
  )
}
