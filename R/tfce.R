## Threshold-free cluster enhancement (Smith & Nichols, 2009).
##
## TFCE replaces the single arbitrary cluster-forming threshold of
## classical cluster inference with an integral over thresholds. Each
## pixel's enhanced value weights the local cluster size at many
## heights, so the statistic is sensitive to both tall-narrow and
## short-broad effects without forcing the analyst to commit to a
## cutoff.
##
## The core formula (Smith & Nichols, 2009, eq. 3):
##
##   TFCE(p) = integral_{h=h0}^{t(p)} e(h, p)^E * h^H dh
##
## where e(h, p) is the extent (in pixels) of the connected component
## containing p at threshold h, H = 2.0 and E = 0.5 are the default
## height and extent exponents, and the integral is approximated on a
## dense threshold grid.
##
## Implementation is pure R. Four-connectivity (via label_components_4conn
## in R/cluster_utils.R) matches the threshold-based cluster path; the
## two methods can be compared head-to-head without confounding
## connectivity choice.

#' TFCE enhancement of a 2D statistic map
#'
#' Positive and negative tails are enhanced separately and combined
#' into a signed map. For a pure-positive or pure-negative use case,
#' pass a non-negative input.
#'
#' @param t_map Numeric vector of length `prod(img_dims)`.
#' @param img_dims Integer `c(nrow, ncol)`.
#' @param H Height exponent. Default 2.0 (Smith & Nichols default).
#' @param E Extent exponent. Default 0.5 (Smith & Nichols default).
#' @param n_steps Number of thresholds in the integration grid.
#'   Default 100. Finer grids are more accurate but scale the cost
#'   linearly.
#' @return A numeric vector the same length as `t_map`, containing
#'   TFCE-enhanced values with sign preserved.
#' @keywords internal
#' @noRd
tfce_enhance <- function(t_map, img_dims, H = 2.0, E = 0.5,
                         n_steps = 100L) {
  stopifnot(length(t_map) == prod(img_dims))
  tmat <- matrix(t_map, nrow = img_dims[1L], ncol = img_dims[2L])
  pos <- tfce_enhance_half( tmat, H, E, n_steps)
  neg <- tfce_enhance_half(-tmat, H, E, n_steps)
  as.vector(pos - neg)
}

#' One-sided TFCE: non-negative input, non-negative output
#'
#' @keywords internal
#' @noRd
tfce_enhance_half <- function(tmat, H, E, n_steps) {
  # Zero-out negatives so thresholds run from 0 up.
  tmat[tmat < 0] <- 0
  max_t <- max(tmat)
  if (!is.finite(max_t) || max_t <= 0) {
    return(matrix(0, nrow(tmat), ncol(tmat)))
  }
  dh         <- max_t / n_steps
  thresholds <- seq(dh, max_t, length.out = n_steps)

  out <- matrix(0, nrow(tmat), ncol(tmat))
  for (h in thresholds) {
    mask <- tmat >= h
    if (!any(mask)) next
    labels <- label_components_4conn(mask)
    n_lab  <- max(labels)
    if (n_lab == 0L) next
    sizes <- tabulate(labels[labels != 0L], nbins = n_lab)
    # Contribution at this height: size^E * h^H * dh, added to every
    # pixel inside each cluster.
    contrib_per_id <- sizes^E * h^H * dh
    # Vectorised update: for each labelled pixel, look up its
    # cluster contribution.
    lab_vec <- as.vector(labels)
    mask    <- lab_vec != 0L
    out_vec <- as.vector(out)
    out_vec[mask] <- out_vec[mask] + contrib_per_id[lab_vec[mask]]
    out <- matrix(out_vec, nrow(tmat), ncol(tmat))
  }
  out
}
