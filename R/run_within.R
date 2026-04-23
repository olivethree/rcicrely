#' Run every within-condition reliability metric
#'
#' @description
#' Convenience orchestrator that runs [rel_split_half()], [rel_loo()],
#' and [rel_icc()] on a single condition's signal matrix and wraps
#' the three results in an `rcicrely_report` for joint printing /
#' plotting.
#'
#' Use this when you want the full within-condition reliability
#' report in one call. Pass each metric individually if you need to
#' tune arguments per metric.
#'
#' @section Reading the result:
#' `$results$split_half`, `$results$loo`, `$results$icc`, one
#' result object each, with the same fields as the standalone
#' functions. `$method = "within"`.
#'
#' @section Reliability metrics expect raw masks:
#' All three downstream metrics expect the raw mask. See
#' `vignette("tutorial", package = "rcicrely")` chapter 3.
#'
#' @param signal_matrix Pixels x participants, base-subtracted.
#' @param n_permutations Passed to [rel_split_half()]. Default 2000.
#' @param flag_threshold Passed to [rel_loo()]. Default 2.5.
#' @param flag_method Passed to [rel_loo()]. Default `"sd"`.
#' @param flag_threshold_sd Deprecated alias for `flag_threshold`.
#' @param icc_variants Passed to [rel_icc()].
#' @param seed Optional integer; used for the split-half permutations.
#' @param progress Show `cli` progress bars.
#' @return Object of class `rcicrely_report` with `$results` =
#'   named list of three `rcicrely_*` objects (`split_half`, `loo`,
#'   `icc`) and `$method = "within"`.
#' @seealso [rel_split_half()], [rel_loo()], [rel_icc()], [run_between()]
#' @export
run_within <- function(signal_matrix,
                       n_permutations    = 2000L,
                       flag_threshold    = 2.5,
                       flag_method       = c("sd", "mad"),
                       flag_threshold_sd = NULL,
                       icc_variants      = c("3_1", "3_k"),
                       seed              = NULL,
                       progress          = TRUE) {
  validate_signal_matrix(signal_matrix)
  flag_method <- match.arg(flag_method)
  if (!is.null(flag_threshold_sd)) flag_threshold <- flag_threshold_sd
  img_dims <- attr(signal_matrix, "img_dims")

  results <- list(
    split_half = rel_split_half(
      signal_matrix,
      n_permutations = n_permutations,
      seed           = seed,
      progress       = progress
    ),
    loo = rel_loo(
      signal_matrix,
      flag_threshold = flag_threshold,
      flag_method    = flag_method
    ),
    icc = rel_icc(
      signal_matrix,
      variants = icc_variants
    )
  )
  new_rcicrely_report(results, method = "within", img_dims = img_dims)
}
