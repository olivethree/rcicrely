#' Read a directory of CI images into a pixels x participants matrix
#'
#' Reads every PNG / JPEG in `dir` (non-recursive), converts each to
#' grayscale via ITU-R BT.709 luminance weights
#' (`0.2126 R + 0.7152 G + 0.0722 B`), and returns a numeric matrix
#' with one column per file, sorted alphabetically by filename.
#' Column names are the filename without the extension.
#'
#' All images must share the same dimensions; mixed sizes abort.
#' The returned matrix is the **raw CI** pixel data  --  to obtain the
#' base-subtracted signal matrix that every `rel_*()` function takes,
#' pass this through [extract_signal()] or use [load_signal_matrix()]
#' which composes the two steps.
#'
#' @param dir Directory containing the CI images.
#' @param pattern Regular expression used to select files. Defaults
#'   to any PNG / JPG / JPEG (case-insensitive).
#' @return A numeric matrix with `nrow = n_pixels` and
#'   `ncol = n_participants`. Row-major pixel ordering follows
#'   `as.vector(img)` of the image matrix.
#' @seealso [extract_signal()], [load_signal_matrix()]
#' @export
#' @examples
#' \dontrun{
#' cis <- read_cis("data/cis_condition_A/")
#' signal <- extract_signal(cis, base_image_path = "data/base.jpg")
#' }
read_cis <- function(dir,
                     pattern = "\\.(png|jpe?g)$") {
  if (!dir.exists(dir)) {
    cli::cli_abort("Directory not found: {.path {dir}}")
  }
  files <- list.files(
    dir,
    pattern     = pattern,
    ignore.case = TRUE,
    full.names  = TRUE
  )
  files <- sort(files)
  if (length(files) == 0L) {
    cli::cli_abort(c(
      "No image files matching {.val {pattern}} in {.path {dir}}.",
      "i" = "Pattern is case-insensitive by default."
    ))
  }

  first <- read_image_as_gray(files[1L])
  dims  <- dim(first)
  n_pix <- prod(dims)

  out <- matrix(NA_real_, nrow = n_pix, ncol = length(files))
  colnames(out) <- tools::file_path_sans_ext(basename(files))
  out[, 1L] <- as.vector(first)

  for (i in seq_along(files)[-1L]) {
    img <- read_image_as_gray(files[i])
    if (!identical(dim(img), dims)) {
      cli::cli_abort(c(
        "Image dimensions differ inside {.path {dir}}.",
        "*" = "{.path {basename(files[1L])}}: \\
               {dims[1]} x {dims[2]}",
        "*" = "{.path {basename(files[i])}}: \\
               {dim(img)[1]} x {dim(img)[2]}"
      ))
    }
    out[, i] <- as.vector(img)
  }

  attr(out, "img_dims") <- as.integer(dims)
  out
}


#' Read a directory of CI images and extract the base-subtracted signal
#'
#' Convenience wrapper composing [read_cis()] and [extract_signal()].
#' Returns the **signal matrix** (pixels x participants, base-subtracted)
#' that every `rel_*()` function takes.
#'
#' @param dir Directory of CI images (as for [read_cis()]).
#' @param base_image_path Path to the base face image used at
#'   stimulus-generation time.
#' @param pattern Regex for file selection (see [read_cis()]).
#' @return A numeric pixels x participants matrix, base-subtracted.
#' @seealso [read_cis()], [extract_signal()]
#' @export
#' @examples
#' \dontrun{
#' signal <- load_signal_matrix(
#'   dir             = "data/cis_condition_A/",
#'   base_image_path = "data/base.jpg"
#' )
#' }
load_signal_matrix <- function(dir,
                               base_image_path,
                               pattern = "\\.(png|jpe?g)$") {
  cis <- read_cis(dir, pattern = pattern)
  extract_signal(cis, base_image_path = base_image_path)
}
