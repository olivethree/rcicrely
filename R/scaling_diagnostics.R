## Heuristic detection of "this signal matrix looks like a rendered
## (scaled) CI rather than a raw mask".
##
## Intentionally not exported. The check is a cheap sanity net wired
## into every `rel_*()` and `run_*()` entry, not a precise classifier.
## False positives are silenced with
## `options(rcicrely.silence_scaling_warning = TRUE)`.

#' Cheap heuristic: does this signal matrix look rendered (scaled)?
#'
#' Raw masks (from Brief-RC `genMask` or rcicr's `$ci`) are sums of
#' sign-weighted noise patterns divided by trial count, so values
#' typically have a small overall standard deviation (rough order
#' 0.02 - 0.05). Rendered CIs read from PNG via [read_cis()] /
#' [extract_signal()] carry per-pixel `scaling(mask)` magnitudes,
#' which are typically an order of magnitude larger (SD ~ 0.15+),
#' since the rendered CI fills the 0-to-1 PNG range with a base
#' that itself varies across the image.
#'
#' Heuristic: flag when the overall standard deviation exceeds a
#' generous raw-mask bound (default 0.15). Intentionally permissive:
#' legitimate raw masks with very high SNR may also exceed it, but
#' the warning is once-per-session and easy to silence. The
#' alternative (false negatives on rendered PNGs) is the failure mode
#' we care about.
#'
#' @param x A pixels x participants numeric matrix.
#' @param raw_sd_bound Numeric. Overall SD above which the matrix is
#'   flagged as suspected-rendered. Default 0.15.
#' @return Logical scalar.
#' @keywords internal
#' @noRd
looks_scaled <- function(x, raw_sd_bound = 0.15) {
  if (!is.matrix(x) || !is.numeric(x) || length(x) == 0L) return(FALSE)
  vals <- x[is.finite(x)]
  if (length(vals) < 2L) return(FALSE)
  isTRUE(stats::sd(vals) > raw_sd_bound)
}
