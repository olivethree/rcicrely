#' Run every within-condition reliability metric
#'
#' Convenience orchestrator that runs [rel_split_half()], [rel_loo()],
#' and [rel_icc()] on a single condition's signal matrix and wraps
#' the three results in an `rcicrely_report`.
#'
#' @param signal_matrix Pixels x participants, base-subtracted.
#' @param n_permutations Passed to [rel_split_half()]. Default 2000.
#' @param flag_threshold_sd Passed to [rel_loo()]. Default 2.5.
#' @param icc_variants Passed to [rel_icc()].
#' @param seed Optional integer; used for the split-half permutations.
#' @param progress Show `cli` progress bars.
#' @return Object of class `rcicrely_report` with `$results` =
#'   named list of three `rcicrely_*` objects (`split_half`, `loo`,
#'   `icc`) and `$method = "within"`.
#' @seealso [rel_split_half()], [rel_loo()], [rel_icc()], [run_between()]
#' @export
run_within <- function(signal_matrix,
                       n_permutations     = 2000L,
                       flag_threshold_sd  = 2.5,
                       icc_variants       = c("3_1", "3_k"),
                       seed               = NULL,
                       progress           = TRUE) {
  validate_signal_matrix(signal_matrix)
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
      flag_threshold_sd = flag_threshold_sd
    ),
    icc = rel_icc(
      signal_matrix,
      variants = icc_variants
    )
  )
  new_rcicrely_report(results, method = "within", img_dims = img_dims)
}
