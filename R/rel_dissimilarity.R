#' Between-condition dissimilarity with bootstrap confidence intervals
#'
#' @description
#' Quantifies overall dissimilarity between two conditions' group-level
#' classification images. The primary statistic is **Euclidean
#' distance** between the two group-mean CIs, reported both raw and
#' normalised by `sqrt(n_pixels)` for cross-resolution comparability.
#' Percentile bootstrap 95% confidence intervals are computed by
#' resampling participants with replacement within each condition.
#'
#' Pair with [rel_cluster_test()] when you also want to know **where**
#' the two conditions differ.
#'
#' @details
#' For the observed statistics and each bootstrap replicate `i`:
#'
#' ```
#' mean_a      = rowMeans(signal_matrix_a)
#' mean_b      = rowMeans(signal_matrix_b)
#' observed_dist            = sqrt(sum((mean_a - mean_b)^2))
#' observed_dist_normalised = observed_dist / sqrt(n_pixels)
#' ```
#'
#' Percentile CI via base R `quantile()`; no `boot` dependency. BCa
#' intervals are on the roadmap for a future release.
#'
#' @section Why Euclidean and not Pearson correlation as the primary:
#' An earlier release emphasised Pearson correlation between the two
#' group-mean CIs as a co-primary metric. The current release retains
#' those fields **for backwards compatibility only**; they are
#' scheduled for removal in v0.3 and **should not be reported as
#' primary** in new work. The reason is a positive baseline: two
#' base-subtracted CIs share systematic image-domain spatial
#' structure (face shape, oval signal support, low-frequency
#' Gaussian-noise smoothness) that pushes their correlation above
#' zero even when the underlying mental representations are
#' unrelated. Absolute correlation values therefore do not cleanly
#' mean "these conditions are similar"; they conflate real
#' similarity with shared image-domain structure.
#'
#' Euclidean distance does not share this failure mode. It is an
#' honest magnitude summary of how far the two group-mean CIs sit
#' from each other in pixel-signal units. For cross-study comparisons,
#' use `$euclidean_normalised` (divided by `sqrt(n_pixels)` so the
#' metric is resolution-invariant).
#'
#' @param signal_matrix_a,signal_matrix_b Pixels x participants,
#'   base-subtracted. Row counts must match.
#' @param n_boot Bootstrap replicates. Default 2000.
#' @param ci_level Confidence level. Default 0.95.
#' @param seed Optional integer; RNG state restored on exit.
#' @param progress Show a `cli` progress bar.
#' @section Reading the result:
#' * `$euclidean`, observed Euclidean distance between group means
#'   on the full sample (primary statistic).
#' * `$euclidean_normalised`, `$euclidean / sqrt(n_pixels)`. Use this
#'   for cross-resolution comparisons.
#' * `$boot_dist`, `$ci_dist`, `$boot_se_dist`, bootstrap distribution,
#'   percentile CI, and SE of the Euclidean distance.
#' * `$correlation`, `$boot_cor`, `$ci_cor`, `$boot_se_cor`: Pearson
#'   correlation of the group means and its bootstrap summary.
#'   **Deprecated.** Retained for backwards compatibility with v0.1.x
#'   code and scheduled for removal in v0.3. Treat these as
#'   within-study descriptive statistics only; see **Why Euclidean
#'   and not Pearson correlation**.
#' * `$n_boot`, `$ci_level`, metadata.
#'
#' @section Common mistakes:
#' * Reporting `$correlation` near 1 as "the conditions look the
#'   same". Two CIs with identical spatial pattern but different
#'   magnitudes produce a correlation near 1; magnitude differences
#'   only show up in `$euclidean`.
#' * Comparing `$euclidean` across studies with different image
#'   resolutions. Use `$euclidean_normalised` for that.
#' * Reading `$correlation` as a standalone "similarity score".
#'   The baseline is above zero for independent CIs (see above).
#'
#' @section Reliability metrics expect raw masks:
#' Euclidean distance is sensitive to **any** scaling; Pearson r
#' survives a single uniform scaling but breaks under per-CI
#' `"matched"`-style scaling. If `signal_matrix_a` / `_b` came from
#' [read_cis()] / [extract_signal()] on rendered PNGs, treat results
#' as approximate. See
#' `vignette("tutorial", package = "rcicrely")` chapter 3.
#'
#' @return Object of class `rcicrely_dissim`. Fields described in
#'   **Reading the result** above.
#' @seealso [rel_cluster_test()], [run_between()]
#' @references
#' Efron, B., & Tibshirani, R. J. (1994). *An introduction to the
#' bootstrap*. Chapman & Hall / CRC.
#' @export
rel_dissimilarity <- function(signal_matrix_a,
                              signal_matrix_b,
                              n_boot   = 2000L,
                              ci_level = 0.95,
                              seed     = NULL,
                              progress = TRUE) {
  validate_two_signal_matrices(signal_matrix_a, signal_matrix_b)
  if (looks_scaled(signal_matrix_a) || looks_scaled(signal_matrix_b)) {
    warn_looks_scaled("signal_matrix_a / _b")
  }
  n_a      <- ncol(signal_matrix_a)
  n_b      <- ncol(signal_matrix_b)
  n_pixels <- nrow(signal_matrix_a)
  n_boot   <- as.integer(n_boot)

  mean_a <- rowMeans(signal_matrix_a)
  mean_b <- rowMeans(signal_matrix_b)

  observed_cor       <- stats::cor(mean_a, mean_b)
  observed_dist      <- sqrt(sum((mean_a - mean_b)^2))
  observed_dist_norm <- observed_dist / sqrt(n_pixels)

  boot_cor  <- numeric(n_boot)
  boot_dist <- numeric(n_boot)

  pid <- progress_start(n_boot, "dissimilarity bootstrap",
                        show = progress)
  on.exit(progress_done(pid), add = TRUE)

  with_seed(seed, {
    for (i in seq_len(n_boot)) {
      idx_a <- sample.int(n_a, n_a, replace = TRUE)
      idx_b <- sample.int(n_b, n_b, replace = TRUE)
      m_a <- rowMeans(signal_matrix_a[, idx_a, drop = FALSE])
      m_b <- rowMeans(signal_matrix_b[, idx_b, drop = FALSE])
      boot_cor[i]  <- stats::cor(m_a, m_b)
      boot_dist[i] <- sqrt(sum((m_a - m_b)^2))
      progress_tick(pid)
    }
  })

  tail <- (1 - ci_level) / 2
  ci_cor  <- stats::quantile(boot_cor,  c(tail, 1 - tail),
                             na.rm = TRUE, names = FALSE)
  ci_dist <- stats::quantile(boot_dist, c(tail, 1 - tail),
                             na.rm = TRUE, names = FALSE)

  new_rcicrely_dissim(
    correlation          = observed_cor,
    euclidean            = observed_dist,
    euclidean_normalised = observed_dist_norm,
    boot_cor             = boot_cor,
    boot_dist            = boot_dist,
    ci_cor               = ci_cor,
    ci_dist              = ci_dist,
    boot_se_cor          = stats::sd(boot_cor,  na.rm = TRUE),
    boot_se_dist         = stats::sd(boot_dist, na.rm = TRUE),
    n_boot               = n_boot,
    ci_level             = ci_level,
    n_pixels             = n_pixels
  )
}
