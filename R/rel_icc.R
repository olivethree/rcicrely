#' Intraclass correlation coefficients via direct mean squares
#'
#' @description
#' Reports the reliability of producer-level CI signals as a single
#' scalar (or two scalars) per condition. Use **ICC(3,1)** to ask
#' "how informative is one producer's CI as a noisy estimate of the
#' group pattern?" and **ICC(3,k)** to ask "how stable is the
#' group-mean CI this experiment produced?". Most papers report
#' ICC(3,k) as the headline.
#'
#' @details
#' Reports **ICC(3,1)** and **ICC(3,k)**, two-way mixed model with
#' pixels fixed and participants random. Pixels are a fixed
#' `img_size x img_size` grid (not a random sample), so ICC(2,\*) is
#' mis-specified even when numerically similar. ICC(2,1) and
#' ICC(2,k) are available via `variants` for users whose reviewers
#' explicitly ask.
#'
#' Computed directly from ANOVA mean squares (never via
#' `psych::ICC()`, which allocates intermediates that blow memory on
#' a 262,144 x 30 matrix).
#'
#' @section What this ICC is, and is not:
#' Most ICCs reported in the reverse-correlation literature are
#' **trait-rating** reliability, phase-2 naive raters scoring CIs on
#' trait dimensions (trustworthy, competent, ...). rcicrely's ICC is
#' structurally different: it operates on the pixel-level signal
#' produced by the original producers. No phase-2 rating study is
#' involved. This sidesteps the two-phase design Cone, Brown-Iannuzzi,
#' Lei, & Dotsch (2021) showed inflates Type I error.
#'
#' @param signal_matrix Pixels x participants (targets x raters),
#'   base-subtracted.
#' @param variants Character vector of which ICC variants to return.
#'   Subset of `c("3_1", "3_k", "2_1", "2_k")`. Defaults to
#'   `c("3_1", "3_k")`.
#' @section Reading the result:
#' * `$icc_3_1`, `$icc_3_k`, two-way-mixed ICCs for single rater
#'   and average rater respectively. `$icc_2_1` / `$icc_2_k` only
#'   present if requested via `variants`.
#' * `$ms_rows`, `$ms_cols`, `$ms_error`, the underlying ANOVA mean
#'   squares for transparency / reproducibility.
#' * `$n_raters` (= participants), `$n_targets` (= pixels),
#'   `$model` (description string), `$variants` (which were returned).
#'
#' @section Common mistakes:
#' * Asking for ICC(2,*) because a reviewer expects it. ICC(2,*)
#'   treats pixels as a random sample from a pixel population, which
#'   they aren't, the image grid is fixed. Numbers are usually
#'   close to ICC(3,*) at high pixel counts but the model is
#'   mis-specified. Use `variants = c("3_1", "3_k", "2_1", "2_k")` to
#'   report both side-by-side.
#' * Comparing this ICC to a phase-2 trait-rating ICC from a different
#'   paper. Different statistical objects (see **What this ICC is**
#'   section below).
#'
#' @section Reliability metrics expect raw masks:
#' ICC is a variance-based statistic. **Strongly** sensitive to any
#' scaling step. If `signal_matrix` was extracted from rendered
#' (scaled) PNGs via [read_cis()] / [extract_signal()], the variance
#' decomposition is computed on `scaling(mask)` rather than `mask` and
#' the result will not match what raw responses would produce. See
#' `vignette("tutorial", package = "rcicrely")` chapter 3.
#'
#' @return Object of class `rcicrely_icc`. Fields described in
#'   **Reading the result** above.
#' @seealso [rel_split_half()], [rel_loo()], [run_within()]
#' @references
#' Shrout, P. E., & Fleiss, J. L. (1979). Intraclass correlations:
#' uses in assessing rater reliability. *Psychological Bulletin*,
#' 86(2), 420-428. \doi{10.1037/0033-2909.86.2.420}
#'
#' McGraw, K. O., & Wong, S. P. (1996). Forming inferences about some
#' intraclass correlation coefficients. *Psychological Methods*,
#' 1(1), 30-46. \doi{10.1037/1082-989X.1.1.30}
#'
#' Cone, J., Brown-Iannuzzi, J. L., Lei, R., & Dotsch, R. (2021). Type
#' I error is inflated in the two-phase reverse correlation procedure.
#' *Social Psychological and Personality Science*, 12(5), 760-768.
#' \doi{10.1177/1948550620938616}
#' @export
rel_icc <- function(signal_matrix,
                    variants = c("3_1", "3_k")) {
  validate_signal_matrix(signal_matrix)
  if (looks_scaled(signal_matrix)) warn_looks_scaled("signal_matrix")
  variants <- match.arg(variants,
                        choices = c("3_1", "3_k", "2_1", "2_k"),
                        several.ok = TRUE)

  n <- nrow(signal_matrix)   # targets (pixels, fixed)
  k <- ncol(signal_matrix)   # raters (participants, random)

  # Large-image advisory for ICC(3,k) asymptote. Once per session.
  if ("3_k" %in% variants && n > 50000L) {
    warn_icc_large_image(n)
  }

  row_means  <- rowMeans(signal_matrix)
  col_means  <- colMeans(signal_matrix)
  grand_mean <- mean(signal_matrix)

  ss_rows  <- k * sum((row_means - grand_mean)^2)
  ss_cols  <- n * sum((col_means - grand_mean)^2)
  ss_total <- sum((signal_matrix - grand_mean)^2)
  ss_error <- ss_total - ss_rows - ss_cols

  ms_rows  <- ss_rows  / (n - 1)
  ms_cols  <- ss_cols  / (k - 1)
  ms_error <- ss_error / ((n - 1) * (k - 1))

  icc_3_1 <- if ("3_1" %in% variants) {
    (ms_rows - ms_error) / (ms_rows + (k - 1) * ms_error)
  } else NA_real_
  icc_3_k <- if ("3_k" %in% variants) {
    (ms_rows - ms_error) / ms_rows
  } else NA_real_
  icc_2_1 <- if ("2_1" %in% variants) {
    (ms_rows - ms_error) /
      (ms_rows + (k - 1) * ms_error + (k / n) * (ms_cols - ms_error))
  } else NA_real_
  icc_2_k <- if ("2_k" %in% variants) {
    (ms_rows - ms_error) / (ms_rows + (ms_cols - ms_error) / n)
  } else NA_real_

  model <-
    "ICC(3,*) / two-way mixed; pixels fixed, participants random"

  new_rcicrely_icc(
    icc_3_1   = icc_3_1,
    icc_3_k   = icc_3_k,
    icc_2_1   = icc_2_1,
    icc_2_k   = icc_2_k,
    ms_rows   = ms_rows,
    ms_cols   = ms_cols,
    ms_error  = ms_error,
    n_raters  = k,
    n_targets = n,
    model     = model,
    variants  = variants
  )
}
