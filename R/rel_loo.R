#' Leave-one-out influence screening
#'
#' @description
#' Screens for influential producers by asking: "if this producer is
#' removed from the sample, how much does the group classification
#' image change?" Producers whose removal moves the group pattern
#' disproportionately are flagged for inspection. This is an
#' **influence / outlier screening** tool, not a reliability
#' statistic — see below.
#'
#' @details
#' For each producer `i`, compute the Pearson correlation between the
#' full-sample group CI and the group CI recomputed without that
#' producer:
#' ```
#' full        <- rowMeans(signal_matrix)
#' r_loo[i]    <- cor(full, rowMeans(signal_matrix[, -i]))
#' ```
#'
#' Because the full-sample mean and the leave-one-out mean share
#' `(N - 1) / N` of their data, `r_loo` values are near 1 by
#' construction even on noisy data, typically `[0.95, 0.999]` at
#' `N = 30`. **The absolute level of `r_loo` is not informative.**
#' What is informative is the *relative ordering*: producers whose
#' `r_loo` sits clearly below the pack are candidates for inspection.
#' For this reason the function also returns a **z-scored** version of
#' `r_loo` in `$z_scores`, using the same centre / spread estimators
#' as the flagging rule. `$z_scores` is the recommended quantity to
#' plot or report.
#'
#' Two flagging rules:
#'
#' * `"sd"` (default): flag producers with
#'   `r_loo < mean(r) - flag_threshold * sd(r)`, equivalently
#'   `z_scores < -flag_threshold`. Standard convention; sensitive to
#'   the very outliers it is trying to detect, so the flagging
#'   threshold can be pulled into the outliers' region.
#' * `"mad"`: flag producers with
#'   `r_loo < median(r) - flag_threshold * mad(r)`, equivalently
#'   `z_scores < -flag_threshold` with a robust centre/spread. Robust
#'   to the few atypical producers that typically show up in RC
#'   datasets and to skewed correlation distributions.
#'
#' The default `flag_threshold = 2.5` is calibrated so that a
#' 30-producer dataset flags roughly 0.3 producers by chance under
#' `"sd"`, rather than the ~1.5 a 2-SD rule would produce. Under
#' `"mad"` it is roughly comparable thanks to MAD's 1.4826 consistency
#' factor.
#'
#' @section What this function is, and is not:
#' `rel_loo()` is an **influence-screening diagnostic**. It answers
#' "which producers disproportionately shape the group CI?", not "how
#' reliable is the group CI?". For reliability, use
#' [rel_split_half()] or [rel_icc()]. A flag does not mean the
#' producer is "bad"; it means the producer's individual CI sits far
#' enough from the group pattern that the data deserve a second look
#' (response coding, fatigue, task misunderstanding, or a genuinely
#' atypical mental representation).
#'
#' @section Reading the result:
#' * `$z_scores`, named numeric vector, per-producer standardised
#'   influence. This is the recommended quantity to plot, report, or
#'   threshold. Values near 0 = typical producer; values below
#'   `-flag_threshold` = flagged.
#' * `$correlations`, named numeric vector, raw per-producer
#'   `r_loo` values. Included for transparency; see the note above
#'   about why the raw level is not informative.
#' * `$mean_r`, `$sd_r`, `$median_r`, `$mad_r`, centre / spread of
#'   the correlation distribution under each rule.
#' * `$threshold`, the raw cutoff value on `r_loo` under the chosen
#'   `flag_method`.
#' * `$flagged`, character vector of producer ids below threshold.
#' * `$summary_df`, one row per producer with `correlation`,
#'   `z_score`, and `flag`, sorted by `z_score`.
#' * `$flag_method`, `$flag_threshold`, what was used.
#'
#' @section Common mistakes:
#' * Reading `r_loo` as a reliability. An `r_loo` of .98 does not
#'   mean the CI is 98% reliable; it means a single producer's
#'   removal changed the group mean by 2%, which is the expected
#'   scale at `N = 30`. Report reliability via [rel_split_half()] or
#'   [rel_icc()].
#' * Treating `$flagged` as "drop these producers". Investigate first
#'   (response coding, fatigue, atypical strategy). Cross-check with
#'   the `rcicrdiagnostics` companion package to rule out
#'   response-coding errors.
#' * Lowering `flag_threshold` below 2 to flag more producers. That
#'   trades real signal for noise; use `flag_method = "mad"` instead
#'   if the SD rule is dominated by the outliers it would otherwise
#'   catch.
#'
#' @section Reliability metrics expect raw masks:
#' Operates on the raw mask; results may be distorted if
#' `signal_matrix` was extracted from rendered (scaled) PNGs. See
#' `vignette("tutorial", package = "rcicrely")` chapter 3.
#'
#' @param signal_matrix Pixels x participants, base-subtracted.
#' @param flag_threshold Numeric multiplier on `sd` (or `mad`) below
#'   the centre, defining the outlier cutoff. Default 2.5. Was
#'   `flag_threshold_sd` before v0.1.1; the old name still works as
#'   an alias.
#' @param flag_method One of `"sd"` (default) or `"mad"`. See
#'   **Details** for which to choose.
#' @param flag_threshold_sd Deprecated alias for `flag_threshold`.
#'   Kept for backwards compatibility with v0.1.0.
#' @param mask Optional logical vector of length
#'   `nrow(signal_matrix)` restricting computation to a region (e.g.,
#'   from [face_mask()] or [load_face_mask()]).
#' @return Object of class `rcicrely_loo`. Fields described in
#'   **Reading the result** above.
#' @seealso [rel_loo_z()] for a tidy z-score accessor;
#'   [rel_split_half()], [rel_icc()] for reliability metrics proper;
#'   [run_within()].
#' @export
rel_loo <- function(signal_matrix,
                    flag_threshold    = 2.5,
                    flag_method       = c("sd", "mad"),
                    flag_threshold_sd = NULL,
                    mask              = NULL) {
  validate_signal_matrix(signal_matrix)
  signal_matrix <- apply_mask_to_signal(signal_matrix, mask)
  if (looks_scaled(signal_matrix)) warn_looks_scaled("signal_matrix")
  flag_method <- match.arg(flag_method)
  if (!is.null(flag_threshold_sd)) {
    flag_threshold <- flag_threshold_sd
  }
  if (!is.numeric(flag_threshold) || length(flag_threshold) != 1L ||
        !is.finite(flag_threshold) || flag_threshold <= 0) {
    cli::cli_abort(
      "{.arg flag_threshold} must be a positive finite numeric scalar."
    )
  }

  n <- ncol(signal_matrix)
  if (is.null(colnames(signal_matrix))) {
    colnames(signal_matrix) <- sprintf("p%03d", seq_len(n))
  }
  full <- rowMeans(signal_matrix)
  cors <- numeric(n)
  names(cors) <- colnames(signal_matrix)
  for (i in seq_len(n)) {
    held_out <- rowMeans(signal_matrix[, -i, drop = FALSE])
    cors[i]  <- stats::cor(full, held_out)
  }

  mean_r   <- mean(cors)
  sd_r     <- stats::sd(cors)
  median_r <- stats::median(cors)
  mad_r    <- stats::mad(cors)

  # z-scores use the same centre/spread as the chosen flagging rule,
  # so the cutoff is always `z < -flag_threshold`.
  z_scores <- if (flag_method == "sd") {
    if (sd_r > 0) (cors - mean_r) / sd_r else rep(0, n)
  } else {
    if (mad_r > 0) (cors - median_r) / mad_r else rep(0, n)
  }
  names(z_scores) <- names(cors)

  threshold <- if (flag_method == "sd") {
    mean_r - flag_threshold * sd_r
  } else {
    median_r - flag_threshold * mad_r
  }
  flagged <- names(cors)[cors < threshold]

  summary_df <- data.frame(
    participant_id = names(cors),
    correlation    = unname(cors),
    z_score        = unname(z_scores),
    flag           = cors < threshold,
    stringsAsFactors = FALSE
  )
  summary_df <- summary_df[order(summary_df$z_score), , drop = FALSE]
  rownames(summary_df) <- NULL

  new_rcicrely_loo(
    correlations   = cors,
    z_scores       = z_scores,
    mean_r         = mean_r,
    sd_r           = sd_r,
    median_r       = median_r,
    mad_r          = mad_r,
    threshold      = threshold,
    flagged        = flagged,
    summary_df     = summary_df,
    flag_method    = flag_method,
    flag_threshold = flag_threshold
  )
}


#' Z-scored leave-one-out influence (accessor)
#'
#' @description
#' Convenience accessor that returns a data frame of producer ids and
#' their z-scored LOO influence, ordered from most-influential
#' (lowest, most negative `z_score`) to least. Accepts either a
#' signal matrix (runs [rel_loo()] under the hood) or an existing
#' `rcicrely_loo` result object (cheap, no recomputation).
#'
#' Use this when you want the ordered ranking of producer influence
#' without the full `rcicrely_loo` list — e.g., for joining against
#' producer-level metadata, passing to `dplyr`/`ggplot2`, or a
#' tidy Supplementary table.
#'
#' @param x Either a `pixels x participants` signal matrix or an
#'   object of class `rcicrely_loo` (as returned by [rel_loo()]).
#' @param ... Passed to [rel_loo()] when `x` is a signal matrix
#'   (e.g. `flag_threshold`, `flag_method`). Ignored when `x` is
#'   already a result object.
#' @return A data frame with columns `participant_id`, `correlation`,
#'   `z_score`, `flag`, sorted by `z_score` ascending.
#' @seealso [rel_loo()]
#' @export
rel_loo_z <- function(x, ...) {
  if (inherits(x, "rcicrely_loo")) {
    return(x$summary_df)
  }
  if (is.matrix(x) && is.numeric(x)) {
    return(rel_loo(x, ...)$summary_df)
  }
  cli::cli_abort(c(
    "{.arg x} must be a numeric signal matrix or an \\
     {.cls rcicrely_loo} object.",
    "i" = "Got {.cls {class(x)}}."
  ))
}
