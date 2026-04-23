#' Representational dissimilarity with bootstrap CIs
#'
#' @description
#' Compute two scalar dissimilarity metrics between two conditions'
#' group-level signals (Pearson correlation and Euclidean distance),
#' with percentile bootstrap confidence intervals.
#'
#' Use this when you want a single number (or two) summarising "how
#' different are these two conditions' CIs?", with uncertainty.
#' Pair with [rel_cluster_test()] when you also want to know
#' **where** they differ.
#'
#' Bootstrap: within each condition, resample participants with
#' replacement; recompute the group mean and the two metrics.
#' `n_boot` replicates give percentile CIs via `quantile()`.
#' No `boot` dependency. BCa is future work.
#'
#' @param signal_matrix_a,signal_matrix_b Pixels x participants,
#'   base-subtracted. Row counts must match.
#' @param n_boot Bootstrap replicates. Default 2000.
#' @param ci_level Confidence level. Default 0.95.
#' @param seed Optional integer; RNG state restored on exit.
#' @param progress Show a `cli` progress bar.
#' @section Reading the result:
#' * `$correlation`, `$euclidean`, observed values on the full
#'   sample. Pearson r is scale-tolerant for a single uniform scaling;
#'   Euclidean distance is sensitive to any scaling.
#' * `$boot_cor`, `$boot_dist`, full bootstrap distributions.
#' * `$ci_cor`, `$ci_dist`, percentile CIs at `ci_level`.
#' * `$boot_se_cor`, `$boot_se_dist`, bootstrap SEs.
#'
#' @section Common mistakes:
#' * Reading `$correlation` near 1 as "the conditions look the same".
#'   Two CIs with identical spatial pattern but different magnitudes
#'   produce r near 1; the magnitude difference shows up in
#'   `$euclidean`.
#' * Treating `$euclidean` as comparable across pixel resolutions
#'   (it scales with sqrt(n_pixels)). Use it as a within-study metric.
#'
#' @section Reliability metrics expect raw masks:
#' Pearson r survives a single uniform linear scaling but breaks
#' under per-CI `"matched"`-style scaling; Euclidean distance is
#' sensitive to **any** scaling. If `signal_matrix_a` / `_b` came
#' from [read_cis()] / [extract_signal()] on rendered PNGs, treat
#' results as approximate. See
#' `vignette("tutorial", package = "rcicrely")` chapter 3.
#'
#' @return Object of class `rcicrely_dissim`. Fields described in
#'   **Reading the result** above.
#' @seealso [rel_cluster_test()], [run_between()]
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
  n_a <- ncol(signal_matrix_a)
  n_b <- ncol(signal_matrix_b)
  n_boot <- as.integer(n_boot)

  mean_a <- rowMeans(signal_matrix_a)
  mean_b <- rowMeans(signal_matrix_b)

  observed_cor  <- stats::cor(mean_a, mean_b)
  observed_dist <- sqrt(sum((mean_a - mean_b)^2))

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
    correlation   = observed_cor,
    euclidean     = observed_dist,
    boot_cor      = boot_cor,
    boot_dist     = boot_dist,
    ci_cor        = ci_cor,
    ci_dist       = ci_dist,
    boot_se_cor   = stats::sd(boot_cor,  na.rm = TRUE),
    boot_se_dist  = stats::sd(boot_dist, na.rm = TRUE),
    n_boot        = n_boot,
    ci_level      = ci_level
  )
}
