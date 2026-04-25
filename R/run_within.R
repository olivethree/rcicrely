#' Run every within-condition reliability metric
#'
#' @description
#' Convenience orchestrator that runs the reliability metrics proper
#' (split-half with Spearman-Brown correction and ICC) on a single
#' condition's signal matrix and wraps the two results in an
#' `rcicrely_report` for joint printing and plotting.
#'
#' Use this when you want the full within-condition reliability
#' report in one call. Pass each metric individually if you need to
#' tune arguments per metric.
#'
#' @section What is included (and what is not):
#' `run_within()` returns the two metrics that quantify the
#' reliability of the group-level classification image proper:
#' **split-half** (a permutation-based estimate of group-CI stability
#' with Spearman-Brown projection to the full sample) and **ICC(3,*)**
#' (the psychometric variance decomposition). These are
#' non-redundant — split-half is resolution-stable and
#' assumption-light; ICC decomposes variance and reports both
#' single-producer and average-measures reliability.
#'
#' Leave-one-out influence screening lives in [rel_loo()] and is
#' **not** bundled here. Its output is an influence diagnostic, not a
#' reliability statistic, and mixing it into the reliability report
#' invites mis-reading `r_loo` values (which are near 1 by
#' construction) as reliability.
#'
#' @section Reading the result:
#' `$results$split_half`, `$results$icc`, one result object each,
#' with the same fields as the standalone functions.
#' `$method = "within"`.
#'
#' @section Reliability metrics expect raw masks:
#' Both downstream metrics expect the raw mask. See
#' `vignette("tutorial", package = "rcicrely")` chapter 3.
#'
#' @param signal_matrix Pixels x participants, base-subtracted.
#' @param n_permutations Passed to [rel_split_half()]. Default 2000.
#' @param icc_variants Passed to [rel_icc()].
#' @param mask Optional logical vector of length `nrow(signal_matrix)`
#'   restricting all metrics to a region (e.g., from [face_mask()]).
#'   Threaded through to both `rel_split_half()` and `rel_icc()`.
#' @param seed Optional integer; used for the split-half permutations.
#' @param progress Show `cli` progress bars.
#' @return Object of class `rcicrely_report` with `$results` =
#'   named list of two `rcicrely_*` objects (`split_half`, `icc`)
#'   and `$method = "within"`.
#' @seealso [rel_split_half()], [rel_icc()], [rel_loo()] for the
#'   influence diagnostic, [run_between()].
#' @export
run_within <- function(signal_matrix,
                       n_permutations    = 2000L,
                       icc_variants      = c("3_1", "3_k"),
                       mask              = NULL,
                       seed              = NULL,
                       progress          = TRUE) {
  validate_signal_matrix(signal_matrix)
  img_dims <- attr(signal_matrix, "img_dims")

  results <- list(
    split_half = rel_split_half(
      signal_matrix,
      n_permutations = n_permutations,
      mask           = mask,
      seed           = seed,
      progress       = progress
    ),
    icc = rel_icc(
      signal_matrix,
      variants = icc_variants,
      mask     = mask
    )
  )
  new_rcicrely_report(results, method = "within", img_dims = img_dims)
}
