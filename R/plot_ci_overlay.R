#' Render a group CI as a translucent agreement-map overlay on a base image
#'
#' @description
#' Plots the producer-mean classification image (or the per-pixel
#' *t*-map from [agreement_map_test()]) as a diverging-palette
#' translucent layer composited over the base face image. Useful as
#' a publication-grade single figure that combines anatomical
#' context (the base) with the empirical signal (the CI / agreement
#' map). Optionally overlays significance contours from an
#' [agreement_map_test()] result.
#'
#' Pure base graphics — no `grid` / `ggplot2` dependency. Returns
#' invisibly the composed `nrow x ncol x 3` numeric raster so it
#' can be saved to disk via [png::writePNG()] / [jpeg::writeJPEG()]
#' if needed.
#'
#' @param signal_matrix Pixels x participants; the producer mean is
#'   computed and rendered as the heatmap layer. Alternatively, a
#'   numeric vector of length `prod(img_dims)` (e.g.
#'   `agreement_map_test(...)$observed_t`) which is rendered
#'   directly.
#' @param base_image Either a numeric matrix (`nrow x ncol`,
#'   grayscale, values in 0-1) or a path to a PNG/JPEG file. When
#'   it's a path, the image is read as grayscale (BT.709 luminance
#'   for RGB inputs) the same way as [read_cis()].
#' @param img_dims Integer `c(nrow, ncol)`. If `NULL`, inferred
#'   from `attr(signal_matrix, "img_dims")` or, when `signal_matrix`
#'   is a vector, the dimensions of `base_image`.
#' @param test Optional [agreement_map_test()] result. When
#'   supplied, significance contours (around `significant_mask`)
#'   are drawn on top of the heatmap.
#' @param mask Optional logical vector restricting the visible
#'   overlay region; pixels outside the mask render as the bare
#'   base image.
#' @param threshold Optional numeric. Pixels with `|signal| <
#'   threshold` are rendered as fully transparent (only the base
#'   shows through). Default `NULL` (no threshold).
#' @param alpha_max Numeric in `[0, 1]`. Maximum opacity at the
#'   peak of the signal magnitude. Default 0.7.
#' @param palette Character. Currently only `"diverging"` is
#'   implemented (red-white-blue). Future-proofed by argument.
#' @param contour_col,contour_lwd Significance-contour colour and
#'   line width when `test` is supplied. Defaults: `"black"`, `1.0`.
#' @param main Optional plot title.
#' @return Invisibly the composed `nrow x ncol x 3` raster. The
#'   plot is drawn on the active device as a side effect.
#' @seealso [agreement_map_test()], [face_mask()]
#' @export
plot_ci_overlay <- function(signal_matrix,
                            base_image,
                            img_dims    = NULL,
                            test        = NULL,
                            mask        = NULL,
                            threshold   = NULL,
                            alpha_max   = 0.7,
                            palette     = c("diverging"),
                            contour_col = "black",
                            contour_lwd = 1.0,
                            main        = NULL) {
  palette <- match.arg(palette)

  # ---- resolve base image ---------------------------------------------
  if (is.character(base_image) && length(base_image) == 1L) {
    base_mat <- read_image_as_gray(base_image)
  } else if (is.matrix(base_image) && is.numeric(base_image)) {
    base_mat <- base_image
  } else {
    cli::cli_abort(
      "{.arg base_image} must be a numeric matrix or a path to PNG/JPEG."
    )
  }
  base_dims <- as.integer(dim(base_mat))

  # ---- resolve signal vector + dims -----------------------------------
  if (is.matrix(signal_matrix)) {
    sig_vec <- rowMeans(signal_matrix)
    if (is.null(img_dims)) {
      img_dims <- attr(signal_matrix, "img_dims")
    }
  } else if (is.numeric(signal_matrix)) {
    sig_vec <- signal_matrix
  } else {
    cli::cli_abort(
      "{.arg signal_matrix} must be a numeric matrix or vector."
    )
  }
  if (is.null(img_dims)) img_dims <- base_dims
  img_dims <- as.integer(img_dims)
  if (length(sig_vec) != prod(img_dims)) {
    cli::cli_abort(c(
      "Pixel count of {.arg signal_matrix} does not match {.arg img_dims}.",
      "*" = "signal length: {length(sig_vec)}",
      "*" = "img_dims: {img_dims[1]} x {img_dims[2]} = {prod(img_dims)}"
    ))
  }
  if (!identical(base_dims, img_dims)) {
    cli::cli_abort(c(
      "Base image and {.arg img_dims} disagree.",
      "*" = "base: {base_dims[1]} x {base_dims[2]}",
      "*" = "img_dims: {img_dims[1]} x {img_dims[2]}"
    ))
  }

  # ---- normalise signal to [-1, 1] ------------------------------------
  rng <- max(abs(sig_vec), na.rm = TRUE)
  if (rng == 0 || !is.finite(rng)) rng <- 1
  norm_vec <- pmax(-1, pmin(1, sig_vec / rng))
  if (!is.null(threshold)) {
    norm_vec[abs(sig_vec) < threshold] <- 0
  }
  if (!is.null(mask)) {
    if (!is.logical(mask) || length(mask) != length(sig_vec)) {
      cli::cli_abort(
        "{.arg mask} must be a logical vector of length {length(sig_vec)}."
      )
    }
    norm_vec[!mask] <- 0
  }

  # ---- diverging palette: red (positive), blue (negative) -------------
  # Output: nrow x ncol x 3 numeric raster in [0, 1]
  base_layer <- array(rep(as.vector(base_mat), 3L),
                      dim = c(img_dims[1], img_dims[2], 3L))
  norm_mat <- matrix(norm_vec, nrow = img_dims[1], ncol = img_dims[2])
  alpha    <- abs(norm_mat) * alpha_max

  # color: r, g, b at full opacity
  pos <- norm_mat > 0
  neg <- norm_mat < 0
  fg_r <- ifelse(pos, 0.85, ifelse(neg, 0.10, 0))
  fg_g <- ifelse(pos, 0.10, ifelse(neg, 0.20, 0))
  fg_b <- ifelse(pos, 0.10, ifelse(neg, 0.85, 0))

  composed <- array(0, dim = c(img_dims[1], img_dims[2], 3L))
  composed[, , 1] <- (1 - alpha) * base_layer[, , 1] + alpha * fg_r
  composed[, , 2] <- (1 - alpha) * base_layer[, , 2] + alpha * fg_g
  composed[, , 3] <- (1 - alpha) * base_layer[, , 3] + alpha * fg_b
  # pmax/pmin drop the array dim; clamp without flattening.
  composed[composed < 0] <- 0
  composed[composed > 1] <- 1

  # ---- draw -----------------------------------------------------------
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  graphics::par(mar = c(0.5, 0.5, if (is.null(main)) 0.5 else 2.5,
                        0.5))
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, img_dims[2]),
                        ylim = c(0, img_dims[1]),
                        asp = 1, xaxs = "i", yaxs = "i")
  graphics::rasterImage(composed, 0, 0, img_dims[2], img_dims[1],
                        interpolate = FALSE)

  # ---- significance contours ------------------------------------------
  if (!is.null(test)) {
    if (!inherits(test, "rcicrely_agreement_map_test")) {
      cli::cli_abort(
        "{.arg test} must be a result from \\
         {.fn agreement_map_test}."
      )
    }
    sig_mat <- matrix(as.numeric(test$significant_mask),
                      nrow = img_dims[1], ncol = img_dims[2])
    # contour() requires monotonic x and y; rasterImage draws the
    # array top-down, so flip the matrix vertically and pass
    # increasing coordinates.
    graphics::contour(
      x = seq_len(img_dims[2]) - 0.5,
      y = seq_len(img_dims[1]) - 0.5,
      z = t(sig_mat[nrow(sig_mat):1L, , drop = FALSE]),
      levels = 0.5, drawlabels = FALSE, add = TRUE,
      col = contour_col, lwd = contour_lwd
    )
  }
  if (!is.null(main)) {
    graphics::title(main = main, line = 1, cex.main = 1.0,
                    font.main = 1)
  }
  invisible(composed)
}
