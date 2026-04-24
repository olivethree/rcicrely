#' Oval face-region mask for a square image
#'
#' @description
#' Returns a logical vector of length `prod(img_dims)` marking an
#' elliptical face region centred on the image. This is the mask
#' Schmitz, Rougier, & Yzerbyt (2024) apply before computing
#' informational value on Brief-RC classification images; passing it
#' to [infoval()] restricts both observed and reference Frobenius
#' norms to the same face region.
#'
#' The default ellipse roughly matches the Schmitz geometry:
#' half-width 35% of image side, half-height 45% of image side,
#' centred at the image midpoint. Use the arguments to tune for
#' your base image.
#'
#' @param img_dims Integer `c(nrow, ncol)`.
#' @param centre Numeric `c(row, col)` in (0, 1) coordinates. Default
#'   `c(0.5, 0.5)`.
#' @param half_width Ellipse horizontal half-axis, in (0, 1) fraction
#'   of image width. Default 0.35.
#' @param half_height Ellipse vertical half-axis, in (0, 1) fraction
#'   of image height. Default 0.45.
#' @return Logical vector of length `prod(img_dims)`, column-major
#'   to match the package's image vectorisation convention.
#' @seealso [infoval()]
#' @references
#' Schmitz, M., Rougier, M., & Yzerbyt, V. (2024). Introducing the
#' brief reverse correlation: an improved tool to assess visual
#' representations. *European Journal of Social Psychology*.
#' \doi{10.1002/ejsp.3100}
#' @export
#' @examples
#' fm <- face_mask(c(128L, 128L))
#' sum(fm) / length(fm)  # roughly pi * 0.35 * 0.45 ~ 0.49
face_mask <- function(img_dims,
                      centre      = c(0.5, 0.5),
                      half_width  = 0.35,
                      half_height = 0.45) {
  img_dims <- as.integer(img_dims)
  if (length(img_dims) == 1L) img_dims <- c(img_dims, img_dims)
  if (length(img_dims) != 2L || any(img_dims < 1L)) {
    cli::cli_abort("{.arg img_dims} must be a positive length-1 or 2 integer.")
  }
  nr <- img_dims[1L]; nc <- img_dims[2L]
  rr <- (row(matrix(0, nr, nc)) - 1L) / (nr - 1L) - centre[1L]
  cc <- (col(matrix(0, nr, nc)) - 1L) / (nc - 1L) - centre[2L]
  inside <- (rr / half_height)^2 + (cc / half_width)^2 <= 1
  as.vector(inside)
}
