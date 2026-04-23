#' Subtract the base image from each column of a CI matrix
#'
#' The raw CI pixels returned by [read_cis()] contain the shared base
#' image plus the participant's signal contribution. Reliability
#' computations must operate on the **signal** alone  --  otherwise the
#' shared base inflates inter-participant correlations. This function
#' reads the base image, converts to grayscale if needed, validates
#' dimensions, and returns `cis - base` column-wise.
#'
#' @param cis A numeric matrix from [read_cis()] (pixels x participants).
#' @param base_image_path Path to the base face image. PNG or JPEG.
#' @return A numeric matrix the same shape as `cis`, base-subtracted.
#' @seealso [read_cis()], [load_signal_matrix()]
#' @export
extract_signal <- function(cis, base_image_path) {
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
  out
}
