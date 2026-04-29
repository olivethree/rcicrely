#' Subtract the base image from each column of a CI matrix
#'
#' @description
#' The pixel values returned by [read_cis()] contain the shared base
#' image plus a (typically scaled) noise contribution. Reliability
#' computations must operate on the **noise** alone, otherwise the
#' shared base inflates inter-participant correlations. This function
#' reads the base image, converts to grayscale if needed, validates
#' dimensions, and returns `cis - base` column-wise.
#'
#' @section Reading the result:
#' Numeric matrix the same shape as `cis`. Column names propagate.
#' The `img_dims` attribute is preserved (or filled in from the base
#' image if missing).
#'
#' @section Common mistakes:
#' * Passing a base image with different dimensions, aborts with a
#'   clear pixel-count mismatch.
#' * Treating the result as the raw mask. The output is
#'   `scaling(mask)` if the input PNGs were rendered with any scaling.
#'   See **Raw vs. rendered CIs** below.
#'
#' @inheritSection read_cis Raw vs. rendered CIs
#'
#' @param cis A numeric matrix from [read_cis()] (pixels x participants).
#' @param base_image_path Path to the base face image. PNG or JPEG.
#' @param acknowledge_scaling Set to `TRUE` to silence the
#'   once-per-session warning. See [read_cis()].
#' @return A numeric matrix the same shape as `cis`, base-subtracted.
#' @seealso [read_cis()], [load_signal_matrix()]
#' @export
extract_signal <- function(cis,
                           base_image_path,
                           acknowledge_scaling = FALSE) {
  if (!is.matrix(cis) || !is.numeric(cis)) {
    cli::cli_abort(
      "{.arg cis} must be a numeric matrix from {.fn read_cis}."
    )
  }
  base <- read_image_as_gray(base_image_path)
  base_vec <- as.vector(base)
  if (length(base_vec) != nrow(cis)) {
    base_dims <- dim(base)
    cli::cli_abort(c(
      "Base image does not match CI pixel count.",
      "*" = "{.arg cis}: {nrow(cis)} pixels per column",
      "*" = "{.arg base_image_path}: \\
             {base_dims[1]} x {base_dims[2]} = {length(base_vec)} pixels"
    ))
  }
  out <- cis - base_vec
  # preserve img_dims attribute if read_cis set it
  if (!is.null(attr(cis, "img_dims"))) {
    attr(out, "img_dims") <- attr(cis, "img_dims")
  } else {
    attr(out, "img_dims") <- as.integer(dim(base))
  }
  # Mode 1 input boundary: PNG-derived signal is rendered, not raw.
  attr(out, "source") <- "rendered"
  warn_mode1_scaling(acknowledge_scaling = acknowledge_scaling)
  out
}
