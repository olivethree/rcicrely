#' Oval face-region mask for a square image
#'
#' @description
#' Returns a logical vector of length `prod(img_dims)` marking an
#' elliptical face region (or a face sub-region) centred on the
#' image. Pass to [infoval()] to restrict both observed and reference
#' Frobenius norms to the masked region; pass to any `rel_*()`
#' function via row-subsetting (`signal_matrix[mask, ]`) to compute
#' reliability or discriminability on a single anatomical region.
#'
#' Five regions are supported:
#' * `"full"` (default): the full face oval (Schmitz, Rougier, &
#'   Yzerbyt 2024 geometry).
#' * `"eyes"`: two small ellipses at the typical eye positions.
#' * `"nose"`: a narrow vertical ellipse along the midline.
#' * `"mouth"`: a wide-and-short ellipse below centre.
#' * `"upper_face"`, `"lower_face"`: the top and bottom halves of the
#'   full face oval, useful for crude top/bottom comparisons.
#'
#' All region geometries are heuristic approximations matched to a
#' typical centred face on a square base image (e.g., 256x256 male
#' face from the KDEF). For non-default base images, use `centre`,
#' `half_width`, `half_height` to tune the full-face oval; the
#' sub-region geometries scale relative to that ellipse.
#'
#' @param img_dims Integer `c(nrow, ncol)`.
#' @param region Character. One of `"full"` (default), `"eyes"`,
#'   `"nose"`, `"mouth"`, `"upper_face"`, `"lower_face"`.
#' @param centre Numeric `c(row, col)` in (0, 1) coordinates. Default
#'   `c(0.5, 0.5)`. Defines the centre of the full-face oval; sub-
#'   region positions are calibrated relative to this centre.
#' @param half_width Full-face ellipse horizontal half-axis, in (0, 1)
#'   fraction of image width. Default 0.35.
#' @param half_height Full-face ellipse vertical half-axis, in (0, 1)
#'   fraction of image height. Default 0.45.
#' @return Logical vector of length `prod(img_dims)`, column-major
#'   to match the package's image vectorisation convention.
#' @seealso [infoval()], [plot_agreement_map()]
#' @references
#' Schmitz, M., Rougier, M., & Yzerbyt, V. (2024). Introducing the
#' brief reverse correlation: an improved tool to assess visual
#' representations. *European Journal of Social Psychology*.
#' \doi{10.1002/ejsp.3100}
#' @export
#' @examples
#' full  <- face_mask(c(128L, 128L))                   # ~0.49 of image
#' eyes  <- face_mask(c(128L, 128L), region = "eyes")  # ~0.02 of image
#' mouth <- face_mask(c(128L, 128L), region = "mouth") # ~0.04 of image
#' c(full = mean(full), eyes = mean(eyes), mouth = mean(mouth))
face_mask <- function(img_dims,
                      region      = c("full", "eyes", "nose",
                                      "mouth", "upper_face",
                                      "lower_face"),
                      centre      = c(0.5, 0.5),
                      half_width  = 0.35,
                      half_height = 0.45) {
  region   <- match.arg(region)
  img_dims <- as.integer(img_dims)
  if (length(img_dims) == 1L) img_dims <- c(img_dims, img_dims)
  if (length(img_dims) != 2L || any(img_dims < 1L)) {
    cli::cli_abort("{.arg img_dims} must be a positive length-1 or 2 integer.")
  }

  nr <- img_dims[1L]; nc <- img_dims[2L]
  full_oval <- ellipse_mask(img_dims, centre[1L], centre[2L],
                            half_height, half_width)

  # Sub-region geometries, expressed as offsets / scales relative to the
  # full-face ellipse. Calibrated for a typical centred face (KDEF male,
  # 256x256). Tune via custom mask construction if your base image
  # geometry differs strongly.
  if (region == "full") {
    return(as.vector(full_oval))
  }

  result <- matrix(FALSE, nr, nc)

  if (region == "eyes") {
    eye_row    <- centre[1L] - 0.10 * (2 * half_height)
    eye_offset <- 0.40 * half_width
    eye_hh     <- 0.10 * half_height
    eye_hw     <- 0.16 * half_width
    result <- result |
      ellipse_mask(img_dims, eye_row, centre[2L] - eye_offset,
                   eye_hh, eye_hw) |
      ellipse_mask(img_dims, eye_row, centre[2L] + eye_offset,
                   eye_hh, eye_hw)
  } else if (region == "nose") {
    nose_row <- centre[1L] + 0.05 * (2 * half_height)
    nose_hh  <- 0.20 * half_height
    nose_hw  <- 0.12 * half_width
    result <- ellipse_mask(img_dims, nose_row, centre[2L],
                           nose_hh, nose_hw)
  } else if (region == "mouth") {
    mouth_row <- centre[1L] + 0.32 * (2 * half_height)
    mouth_hh  <- 0.10 * half_height
    mouth_hw  <- 0.30 * half_width
    result <- ellipse_mask(img_dims, mouth_row, centre[2L],
                           mouth_hh, mouth_hw)
  } else if (region == "upper_face") {
    rr <- (row(matrix(0, nr, nc)) - 1L) / (nr - 1L)
    result <- full_oval & (rr < centre[1L])
  } else if (region == "lower_face") {
    rr <- (row(matrix(0, nr, nc)) - 1L) / (nr - 1L)
    result <- full_oval & (rr >= centre[1L])
  }

  # Sub-regions are always intersected with the full face oval so they
  # never extend into the hair / background.
  result <- result & full_oval
  as.vector(result)
}

#' Internal: ellipse-mask helper
#'
#' Builds a logical matrix of pixels inside an ellipse with the
#' given centre and half-axes (all in (0, 1) image fractions).
#' Centre = (row, col) in (0, 1); half-axes are in (0, 1) of the
#' corresponding image side.
#'
#' @keywords internal
#' @noRd
ellipse_mask <- function(img_dims, c_row, c_col,
                         half_height, half_width) {
  nr <- img_dims[1L]; nc <- img_dims[2L]
  rr <- (row(matrix(0, nr, nc)) - 1L) / (nr - 1L) - c_row
  cc <- (col(matrix(0, nr, nc)) - 1L) / (nc - 1L) - c_col
  (rr / half_height)^2 + (cc / half_width)^2 <= 1
}
