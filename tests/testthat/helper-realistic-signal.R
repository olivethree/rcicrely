## Realistic-scale fixture generators for rcicrely tests. Bigger and
## slower than the toy fixtures in helper-fixtures.R (128 x 128 = 16,384
## pixels, N = 30) so use `skip_on_cran()` on every test that calls
## these, or downsize `n_side` for local runs.

#' Spatially-contiguous signal at known SNR
#'
#' Generates a pixels x participants matrix with a circular
#' signal region in the centre, anti-aliased, at a known
#' signal-to-noise ratio. `snr` is defined as the ratio of
#' target-signal amplitude to noise amplitude, matching the archetype
#' definitions used elsewhere in the package (low ~ 0.1,
#' medium ~ 0.3, high ~ 0.5). Per-producer signal amplitudes
#' jitter around 1.0 so producers are not identical.
#'
#' @param n_side Image side (pixels). Default 128.
#' @param n_p Number of producers. Default 30.
#' @param snr Signal-to-noise ratio. Default 0.3 (medium).
#' @param signal_centre Named list with `row`, `col` in (0, 1) units,
#'   default centred.
#' @param signal_radius Radius as fraction of image side. Default 0.2.
#' @param seed RNG seed.
#' @return pixels x participants matrix with `img_dims` attr and a
#'   `signal_mask` attr (length-`n_side^2` logical for assertions).
make_realistic_sig <- function(n_side        = 128L,
                               n_p           = 30L,
                               snr           = 0.3,
                               signal_centre = list(row = 0.5, col = 0.5),
                               signal_radius = 0.2,
                               seed          = 1L) {
  n_pix <- n_side * n_side
  # build a smooth circular signal mask centred at signal_centre
  rr <- row(matrix(0, n_side, n_side)) / n_side - signal_centre$row
  cc <- col(matrix(0, n_side, n_side)) / n_side - signal_centre$col
  d  <- sqrt(rr^2 + cc^2)
  mask_mat <- pmax(0, 1 - d / signal_radius)         # smooth 1 -> 0 edge
  mask_vec <- as.vector(mask_mat)
  signal_mask <- mask_mat > 0

  set.seed(seed)
  noise <- matrix(stats::rnorm(n_pix * n_p, sd = 1.0), n_pix, n_p)
  # jitter per-producer signal amplitude around 1.0
  amps  <- stats::runif(n_p, 0.7, 1.3)
  sig   <- snr * outer(mask_vec, amps)
  out   <- sig + noise

  attr(out, "img_dims")    <- c(as.integer(n_side), as.integer(n_side))
  attr(out, "signal_mask") <- as.vector(signal_mask)
  out
}

#' Two conditions with disjoint spatial signals at matched SNR
#'
#' Condition A has a signal region in the top-left quadrant;
#' condition B has a signal region in the bottom-right quadrant.
#' Producers in A and B are independent (between-subjects design,
#' which is what the package's between-condition tests assume).
#'
#' @param n_side Image side. Default 128.
#' @param n_p Producers per condition. Default 30.
#' @param snr Signal-to-noise ratio. Default 0.3.
#' @param seed RNG seed.
#' @return List with `a`, `b` (pixels x N matrices with `img_dims`
#'   attr), `img_dims`, `signal_mask_a`, `signal_mask_b`.
make_realistic_pair <- function(n_side = 128L,
                                n_p    = 30L,
                                snr    = 0.3,
                                seed   = 1L) {
  a <- make_realistic_sig(
    n_side = n_side, n_p = n_p, snr = snr,
    signal_centre = list(row = 0.30, col = 0.30),
    signal_radius = 0.18,
    seed = seed
  )
  b <- make_realistic_sig(
    n_side = n_side, n_p = n_p, snr = snr,
    signal_centre = list(row = 0.70, col = 0.70),
    signal_radius = 0.18,
    seed = seed + 1000L
  )
  list(
    a             = a,
    b             = b,
    img_dims      = c(as.integer(n_side), as.integer(n_side)),
    signal_mask_a = attr(a, "signal_mask"),
    signal_mask_b = attr(b, "signal_mask")
  )
}
