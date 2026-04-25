#' Per-pixel agreement heatmap for a producer signal matrix
#'
#' @description
#' Visualises where producers in a single condition agree on the
#' direction of signal. For each pixel, computes a one-sample
#' t-statistic against zero across producers
#' (`mean / (sd / sqrt(N))`), then displays the resulting map with
#' a diverging colour palette (positive = agreement on positive
#' signal, negative = agreement on negative signal, zero = no
#' agreement). Saturation of the colour is the magnitude of the
#' agreement, *not* the value of the group-mean CI.
#'
#' Use this to answer "where do producers consistently *agree* the
#' target trait lives in the face?". Pair with the group-mean CI
#' image (raw mask or rendered) to see *direction* and *agreement*
#' side by side.
#'
#' @details
#' This is structurally a one-sample t-map (vs zero); pixels where
#' producers' contributions are large *and consistent in sign* get
#' high |t|, pixels where contributions are random get t near zero.
#' Cluster-permutation inference would normally accompany this for
#' formal pixel-level FWER control between conditions
#' ([rel_cluster_test()]); the agreement map is the **descriptive**
#' counterpart for a single condition.
#'
#' @param signal_matrix Pixels x participants raw mask (as returned
#'   by `ci_from_responses_*()` or `read_cis()` + `extract_signal()`).
#' @param img_dims Integer `c(nrow, ncol)`. If `NULL`, inferred from
#'   `attr(signal_matrix, "img_dims")` or from `sqrt(n_pixels)` if
#'   the latter is a whole number.
#' @param mask Optional logical vector of length `nrow(signal_matrix)`
#'   restricting display to a region (e.g.,
#'   `face_mask(img_dims, region = "eyes")`).
#' @param threshold Optional positive numeric. When supplied,
#'   pixels with `|t| < threshold` are rendered in the neutral
#'   (white) colour, making clusters of agreement stand out.
#'   Default `NULL` (full continuous map).
#' @param zlim Numeric `c(low, high)` for the colour scale. Default
#'   is symmetric around zero at `c(-max(|t|), max(|t|))` so the
#'   neutral colour aligns with t = 0.
#' @param palette Character. `"diverging"` (default; positive =
#'   blue, negative = red, neutral = white) or `"viridis"` (no
#'   neutral; for absolute-magnitude views).
#' @param main Title.
#' @param ... Passed to `graphics::image()`.
#' @return Invisibly, a list with `t_map` (numeric vector of t values
#'   per pixel), `n` (producer count), `img_dims`, and `mask` (if
#'   supplied) — useful for further analysis or replotting.
#' @seealso [face_mask()], [rel_cluster_test()] for inferential
#'   between-condition tests.
#' @export
#' @examples
#' \dontrun{
#' # Synthetic example
#' set.seed(1)
#' n_side <- 32L
#' n_pix  <- n_side * n_side
#' rr <- row(matrix(0, n_side, n_side)) / n_side - 0.5
#' cc <- col(matrix(0, n_side, n_side)) / n_side - 0.5
#' d  <- sqrt(rr^2 + cc^2)
#' mask_vec <- as.vector(pmax(0, 1 - d / 0.2))
#' signal <- 0.5 * outer(mask_vec, runif(20L, 0.7, 1.3)) +
#'           matrix(rnorm(n_pix * 20L), n_pix, 20L)
#' plot_agreement_map(signal, img_dims = c(n_side, n_side))
#' }
plot_agreement_map <- function(signal_matrix,
                               img_dims  = NULL,
                               mask      = NULL,
                               threshold = NULL,
                               zlim      = NULL,
                               palette   = c("diverging", "viridis"),
                               main      = "Per-pixel producer agreement (t-map)",
                               ...) {
  if (!is.matrix(signal_matrix) || !is.numeric(signal_matrix)) {
    cli::cli_abort("{.arg signal_matrix} must be a numeric matrix.")
  }
  palette <- match.arg(palette)

  n_pix <- nrow(signal_matrix)
  if (is.null(img_dims)) {
    a <- attr(signal_matrix, "img_dims")
    if (!is.null(a)) {
      img_dims <- a
    } else {
      side <- sqrt(n_pix)
      if (side != as.integer(side)) {
        cli::cli_abort(
          "Cannot infer {.arg img_dims}; pass it explicitly."
        )
      }
      img_dims <- c(as.integer(side), as.integer(side))
    }
  }

  n   <- ncol(signal_matrix)
  m   <- rowMeans(signal_matrix)
  v   <- rowSums((signal_matrix - m)^2) / (n - 1L)
  se  <- sqrt(v / n)
  t_map <- m / se
  t_map[!is.finite(t_map)] <- 0

  display <- t_map
  if (!is.null(mask)) {
    if (!is.logical(mask) || length(mask) != n_pix) {
      cli::cli_abort(
        "{.arg mask} must be a logical vector of length {n_pix}."
      )
    }
    display[!mask] <- NA
  }
  if (!is.null(threshold)) {
    display[abs(display) < threshold] <- 0
  }

  if (is.null(zlim)) {
    rng  <- max(abs(display), na.rm = TRUE)
    if (!is.finite(rng) || rng == 0) rng <- 1
    zlim <- c(-rng, rng)
  }

  col_vec <- if (palette == "diverging") {
    grDevices::hcl.colors(256L, "RdBu", rev = TRUE)
  } else {
    grDevices::hcl.colors(256L, "viridis")
  }

  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  graphics::par(mar = c(1, 1, 3, 6) + 0.1)
  mat <- matrix(display, img_dims[1L], img_dims[2L])

  graphics::image(
    seq_len(img_dims[2L]), seq_len(img_dims[1L]),
    t(mat[nrow(mat):1L, ]),
    col       = col_vec,
    zlim      = zlim,
    main      = main,
    axes      = FALSE,
    xlab      = "", ylab = "",
    asp       = img_dims[1L] / img_dims[2L],
    useRaster = TRUE,
    ...
  )
  graphics::box(col = "grey80", lwd = 0.5)

  add_colour_bar(zlim, col_vec, label = "t-value (one-sample vs 0)")

  graphics::mtext(
    sprintf("N = %d producers,  %d x %d pixels%s",
            n, img_dims[1L], img_dims[2L],
            if (!is.null(threshold))
              sprintf(",  thresholded at |t| > %.2f", threshold)
            else ""),
    side = 3, line = 0.3, cex = 0.85, col = "grey30"
  )

  invisible(list(t_map = t_map, n = n,
                 img_dims = img_dims, mask = mask))
}

#' Add a vertical colour bar in the right margin of the active plot
#'
#' @keywords internal
#' @noRd
add_colour_bar <- function(zlim, col, label = NULL,
                           bar_width_frac = 0.05) {
  usr <- graphics::par("usr")
  pin <- graphics::par("pin")
  # Place bar in the right margin
  x_left  <- usr[2] + (usr[2] - usr[1]) * 0.03
  x_right <- usr[2] + (usr[2] - usr[1]) * 0.07
  y_bot   <- usr[3]
  y_top   <- usr[4]
  n_col   <- length(col)
  ys      <- seq(y_bot, y_top, length.out = n_col + 1L)
  graphics::par(xpd = TRUE)
  on.exit(graphics::par(xpd = FALSE), add = TRUE)
  for (i in seq_len(n_col)) {
    graphics::rect(x_left, ys[i], x_right, ys[i + 1L],
                   col = col[i], border = NA)
  }
  graphics::rect(x_left, y_bot, x_right, y_top,
                 border = "grey60", lwd = 0.5)
  ticks <- pretty(zlim, n = 5L)
  ticks <- ticks[ticks >= zlim[1] & ticks <= zlim[2]]
  tick_y <- y_bot + (ticks - zlim[1]) / diff(zlim) * (y_top - y_bot)
  graphics::segments(x_right, tick_y,
                     x_right + (usr[2] - usr[1]) * 0.01,
                     tick_y, col = "grey60", lwd = 0.5)
  graphics::text(x_right + (usr[2] - usr[1]) * 0.015, tick_y,
                 labels = format(ticks), pos = 4, cex = 0.7,
                 col = "grey20")
  if (!is.null(label)) {
    graphics::text(x_right + (usr[2] - usr[1]) * 0.06,
                   (y_bot + y_top) / 2,
                   labels = label, srt = -90, cex = 0.75,
                   col = "grey30")
  }
}
