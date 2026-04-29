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
#' @param paired Logical. `FALSE` (default) for a between-subjects
#'   design: participants are resampled within A and B independently.
#'   `TRUE` for a within-subjects design: A and B share a single
#'   resample index per replicate so the paired covariance structure
#'   is preserved. Requires `ncol(a) == ncol(b)` and, if named,
#'   identical column names.
#' @param n_boot Bootstrap replicates. Default 2000.
#' @param ci_level Confidence level. Default 0.95.
#' @param mask Optional logical vector of length
#'   `nrow(signal_matrix_a)` restricting the Euclidean / correlation
#'   computation to a region. Both matrices are subsetted with the
#'   same mask. The reported `n_pixels` and
#'   `euclidean_normalised` (= `euclidean / sqrt(n_pixels)`) reflect
#'   the masked count.
#' @param null One of `"none"` (default) or `"permutation"`.
#'   `"permutation"` builds an empirical chance baseline for the
#'   Euclidean distance:
#'   * Between-subjects (`paired = FALSE`): stratified condition-
#'     label permutation across producers, preserving `(N_a, N_b)`.
#'   * Paired (`paired = TRUE`): random sign-flip on per-producer
#'     `A - B` differences (exact null under exchangeability of pair
#'     sign).
#' @param n_permutations Integer. Number of null iterations when
#'   `null = "permutation"`. Default 2000.
#' @param seed Optional integer; RNG state restored on exit.
#' @param progress Show a `cli` progress bar.
#' @section Reading the result:
#' * `$euclidean`, observed Euclidean distance between group means
#'   on the full sample (primary statistic).
#' * `$euclidean_normalised`, `$euclidean / sqrt(n_pixels)`. Use this
#'   for cross-resolution comparisons.
#' * `$boot_dist`, `$ci_dist`, `$boot_se_dist`, bootstrap distribution,
#'   percentile CI, and SE of the Euclidean distance.
#' * `$null` (character), the null mode used.
#' * `$null_distribution`, when `null != "none"`: numeric vector of
#'   per-iteration Euclidean distances under the chosen null. `NULL`
#'   otherwise.
#' * `$d_null_p95`, 95th percentile of the null distribution
#'   (cutoff for "above-chance" magnitude). `NA_real_` when
#'   `null = "none"`.
#' * `$d_z`, z-equivalent effect size: `(observed_d - mean(null)) /
#'   sd(null)`. `NA_real_` when `null = "none"`.
#' * `$d_ratio`, observed Euclidean over the null median.
#'   `NA_real_` when `null = "none"`.
#' * `$correlation`, `$boot_cor`, `$ci_cor`, `$boot_se_cor`: Pearson
#'   correlation of the group means and its bootstrap summary.
#'   **Deprecated.** Retained for backwards compatibility with v0.1.x
#'   code and scheduled for removal in v0.3. Treat these as
#'   within-study descriptive statistics only; see **Why Euclidean
#'   and not Pearson correlation**.
#' * `$n_boot`, `$ci_level`, `$paired`, metadata.
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
#' `"matched"`-style scaling. Inputs with
#' `attr(., "source") == "rendered"` (set automatically by Mode 1
#' readers like [extract_signal()]) error unless
#' `acknowledge_scaling = TRUE`. See
#' `vignette("tutorial", package = "rcicrely")` chapter 3.
#'
#' @param acknowledge_scaling Logical. When `FALSE` (default), the
#'   shared `assert_raw_signal()` helper errors on a known-rendered
#'   matrix on either side. Set `TRUE` to override.
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
                              paired              = FALSE,
                              n_boot              = 2000L,
                              ci_level            = 0.95,
                              null                = c("none",
                                                       "permutation"),
                              n_permutations      = 2000L,
                              mask                = NULL,
                              seed                = NULL,
                              progress            = TRUE,
                              acknowledge_scaling = FALSE) {
  validate_two_signal_matrices(signal_matrix_a, signal_matrix_b)
  if (isTRUE(paired)) {
    validate_paired_matrices(signal_matrix_a, signal_matrix_b)
  }
  assert_raw_signal(signal_matrix_a, acknowledge_scaling,
                    name = "signal_matrix_a")
  assert_raw_signal(signal_matrix_b, acknowledge_scaling,
                    name = "signal_matrix_b")
  signal_matrix_a <- apply_mask_to_signal(signal_matrix_a, mask,
                                           name = "signal_matrix_a")
  signal_matrix_b <- apply_mask_to_signal(signal_matrix_b, mask,
                                           name = "signal_matrix_b")
  null <- match.arg(null)
  n_a      <- ncol(signal_matrix_a)
  n_b      <- ncol(signal_matrix_b)
  n_pixels <- nrow(signal_matrix_a)
  n_boot   <- as.integer(n_boot)
  n_permutations <- as.integer(n_permutations)

  mean_a <- rowMeans(signal_matrix_a)
  mean_b <- rowMeans(signal_matrix_b)

  observed_cor       <- stats::cor(mean_a, mean_b)
  observed_dist      <- sqrt(sum((mean_a - mean_b)^2))
  observed_dist_norm <- observed_dist / sqrt(n_pixels)

  boot_cor  <- numeric(n_boot)
  boot_dist <- numeric(n_boot)

  total_iter <- n_boot + (if (null == "none") 0L else n_permutations)
  pid <- progress_start(total_iter, "dissimilarity bootstrap",
                        show = progress)
  on.exit(progress_done(pid), add = TRUE)

  with_seed(seed, {
    for (i in seq_len(n_boot)) {
      if (isTRUE(paired)) {
        # Shared resample index preserves producer-level covariance
        # between A and B.
        idx   <- sample.int(n_a, n_a, replace = TRUE)
        idx_a <- idx
        idx_b <- idx
      } else {
        idx_a <- sample.int(n_a, n_a, replace = TRUE)
        idx_b <- sample.int(n_b, n_b, replace = TRUE)
      }
      m_a <- rowMeans(signal_matrix_a[, idx_a, drop = FALSE])
      m_b <- rowMeans(signal_matrix_b[, idx_b, drop = FALSE])
      boot_cor[i]  <- stats::cor(m_a, m_b)
      boot_dist[i] <- sqrt(sum((m_a - m_b)^2))
      progress_tick(pid)
    }
  })

  d_null <- if (null == "none") {
    NULL
  } else {
    if (isTRUE(paired)) {
      # Sign-flip on per-producer (A - B) — exact null under
      # exchangeability of the pair sign for paired designs.
      diff_mat <- signal_matrix_a - signal_matrix_b
      with_seed(if (is.null(seed)) NULL else seed + 1L, {
        out <- numeric(n_permutations)
        for (i in seq_len(n_permutations)) {
          signs <- sample(c(-1, 1), n_a, replace = TRUE)
          d_flipped <- sweep(diff_mat, 2L, signs, `*`)
          # Under H0 (no condition difference), mean of A - B is 0.
          # We measure the null Euclidean of (mean of sign-flipped
          # diff), since the bootstrapped estimator is the
          # Euclidean of (mean_a - mean_b).
          m_d <- rowMeans(d_flipped)
          out[i] <- sqrt(sum(m_d * m_d))
          progress_tick(pid)
        }
        out
      })
    } else {
      # Stratified condition-label permutation across producers.
      combined <- cbind(signal_matrix_a, signal_matrix_b)
      n_total <- n_a + n_b
      with_seed(if (is.null(seed)) NULL else seed + 1L, {
        out <- numeric(n_permutations)
        for (i in seq_len(n_permutations)) {
          perm <- sample.int(n_total)
          idx_a <- perm[seq_len(n_a)]
          idx_b <- perm[(n_a + 1L):n_total]
          m_a   <- rowMeans(combined[, idx_a, drop = FALSE])
          m_b   <- rowMeans(combined[, idx_b, drop = FALSE])
          out[i] <- sqrt(sum((m_a - m_b)^2))
          progress_tick(pid)
        }
        out
      })
    }
  }

  tail <- (1 - ci_level) / 2
  ci_cor  <- stats::quantile(boot_cor,  c(tail, 1 - tail),
                             na.rm = TRUE, names = FALSE)
  ci_dist <- stats::quantile(boot_dist, c(tail, 1 - tail),
                             na.rm = TRUE, names = FALSE)

  if (is.null(d_null)) {
    d_null_p95 <- NA_real_
    d_z        <- NA_real_
    d_ratio    <- NA_real_
  } else {
    sd_null    <- stats::sd(d_null, na.rm = TRUE)
    mean_null  <- mean(d_null, na.rm = TRUE)
    d_null_p95 <- unname(stats::quantile(d_null, 0.95,
                                         na.rm = TRUE))
    d_z        <- if (sd_null > 0) {
      (observed_dist - mean_null) / sd_null
    } else NA_real_
    med_null   <- stats::median(d_null, na.rm = TRUE)
    d_ratio    <- if (med_null > 0) observed_dist / med_null else NA_real_
  }

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
    n_pixels             = n_pixels,
    null                 = null,
    null_distribution    = d_null,
    d_null_p95           = d_null_p95,
    d_z                  = d_z,
    d_ratio              = d_ratio,
    paired               = isTRUE(paired)
  )
}
