#' Leave-one-out sensitivity
#'
#' @description
#' For each producer, correlates the full-sample group CI with the
#' group CI computed without that producer. A producer whose LOO
#' correlation falls more than `flag_threshold` units below the
#' centre (mean or median, depending on `flag_method`) is flagged as
#' an influential outlier, their removal changes the group CI
#' noticeably more than the others.
#'
#' Use this to spot producers whose individual CI sits far from the
#' group pattern. A flag does not imply a bad producer; it could
#' equally indicate a genuinely atypical mental representation. Cross-
#' check against `rcicrdiagnostics` first to rule out response-coding
#' errors.
#'
#' @details
#' Two outlier rules are available:
#'
#' * `"sd"` (default): flag producers with `r_loo < mean(r) -
#'   flag_threshold * sd(r)`. Standard convention; sensitive to the
#'   very outliers it is trying to detect.
#' * `"mad"`: flag producers with `r_loo < median(r) -
#'   flag_threshold * mad(r)`. Robust to the few atypical producers
#'   that often show up in RC datasets and to skewed correlation
#'   distributions.
#'
#' The default `flag_threshold = 2.5` is calibrated so that a
#' 30-producer dataset flags roughly 0.3 producers by chance under
#' `"sd"`, rather than the ~1.5 a 2-SD rule would produce. Under
#' `"mad"` it is roughly comparable thanks to MAD's 1.4826 consistency
#' factor.
#'
#' @section Reading the result:
#' * `$correlations`, named numeric vector, per-producer correlation
#'   between the full-sample mean and the leave-one-out mean.
#'   Higher = the producer's CI looks like the rest.
#' * `$mean_r`, `$sd_r`, `$median_r`, `$mad_r`, centre / spread of
#'   the correlation distribution.
#' * `$threshold`, the cutoff value computed under the chosen
#'   `flag_method`.
#' * `$flagged`, character vector of producer ids below threshold.
#' * `$summary_df`, one row per producer, with `correlation` and
#'   `flag`.
#' * `$flag_method`, `$flag_threshold`, what was used.
#'
#' @section Common mistakes:
#' * Treating `$flagged` as "drop these producers". Investigate first
#'   (response coding, fatigue, atypical strategy).
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
#' @return Object of class `rcicrely_loo`. Fields described in
#'   **Reading the result** above.
#' @seealso [rel_split_half()], [rel_icc()], [run_within()]
#' @export
rel_loo <- function(signal_matrix,
                    flag_threshold    = 2.5,
                    flag_method       = c("sd", "mad"),
                    flag_threshold_sd = NULL) {
  validate_signal_matrix(signal_matrix)
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

  threshold <- if (flag_method == "sd") {
    mean_r - flag_threshold * sd_r
  } else {
    median_r - flag_threshold * mad_r
  }
  flagged <- names(cors)[cors < threshold]

  summary_df <- data.frame(
    participant_id = names(cors),
    correlation    = unname(cors),
    flag           = cors < threshold,
    stringsAsFactors = FALSE
  )

  new_rcicrely_loo(
    correlations   = cors,
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
