#' Read a directory of CI images into a pixels x participants matrix
#'
#' @description
#' Reads every PNG / JPEG in `dir` (non-recursive), converts each to
#' grayscale via ITU-R BT.709 luminance weights
#' (`0.2126 R + 0.7152 G + 0.0722 B`), and returns a numeric matrix
#' with one column per file, sorted alphabetically by filename.
#'
#' Use this when you have one CI image per producer already saved on
#' disk and you want a numeric matrix you can pass to
#' [extract_signal()] or [load_signal_matrix()] for the base
#' subtraction step.
#'
#' @details
#' All images must share the same dimensions; mixed sizes abort.
#' Column names are the filename without the extension.
#'
#' @section Reading the result:
#' The returned matrix has `nrow = n_pixels` and
#' `ncol = n_participants`. Row-major pixel ordering follows
#' `as.vector(img)` of the underlying image matrix. The integer
#' `c(nrow, ncol)` of the original image is stamped on the result as
#' an `img_dims` attribute, which downstream functions read
#' automatically.
#'
#' @section Common mistakes:
#' * Forgetting to subtract the base. The matrix returned here
#'   contains the **rendered CI** (base + scaled noise). Pass it
#'   through [extract_signal()] before any reliability metric.
#' * Mixing CIs from different image sizes in one directory, this
#'   aborts with a clear error rather than silently misaligning.
#' * Trusting reliability numbers without first reading the **Note**
#'   below on raw vs. rendered CIs.
#'
#' @section Raw vs. rendered CIs:
#' PNG pixel values are necessarily what was rendered to disk:
#' `base + scaling(mask)`. Subtracting the base in [extract_signal()]
#' therefore yields `scaling(mask)`, not the raw mask. Variance-based
#' reliability metrics ([rel_icc()], Euclidean half of
#' [rel_dissimilarity()], [pixel_t_test()], [rel_cluster_test()]) are
#' sensitive to the scaling transform; only correlation-based metrics
#' ([rel_split_half()], [rel_loo()], correlation half of
#' [rel_dissimilarity()]) survive a single uniform linear scaling
#' unmodified, and even those break under per-CI `"matched"`-style
#' scaling. The standard 2IFC `infoVal` path
#' (`rcicr::computeInfoVal2IFC()`) extracts the raw `$ci` internally
#' from the rcicr CI-list and is **not** affected; hand-rolled
#' `infoVal` implementations (Brief-RC, custom code) **are**, and
#' should be fed the raw mask. Prefer [ci_from_responses_2ifc()] /
#' [ci_from_responses_briefrc()] when raw responses are available;
#' both return the raw mask. See
#' `vignette("tutorial", package = "rcicrely")` chapter 3 for the
#' full discussion.
#'
#' @param dir Directory containing the CI images.
#' @param pattern Regular expression used to select files. Defaults
#'   to any PNG / JPG / JPEG (case-insensitive).
#' @param acknowledge_scaling Set to `TRUE` to silence the
#'   once-per-session warning about Mode 1 returning rendered (scaled)
#'   pixel data. Recommended only when you generated the PNGs yourself
#'   with `rcicr::generateCI2IFC(scaling = "none")` (or equivalent)
#'   and the PNG truly carries unscaled noise.
#' @return A numeric matrix with `nrow = n_pixels` and
#'   `ncol = n_participants`. Carries an `img_dims` attribute.
#' @seealso [extract_signal()], [load_signal_matrix()],
#'   [ci_from_responses_2ifc()], [ci_from_responses_briefrc()]
#' @export
#' @examples
#' \dontrun{
#' cis <- read_cis("data/cis_condition_A/")
#' signal <- extract_signal(cis, base_image_path = "data/base.jpg")
#' }
read_cis <- function(dir,
                     pattern             = "\\.(png|jpe?g)$",
                     acknowledge_scaling = FALSE) {
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
  # PNG-derived raw values are pre-base-subtraction but still rendered.
  attr(out, "source") <- "rendered"
  warn_mode1_scaling(acknowledge_scaling = acknowledge_scaling)
  out
}


#' Read a directory of CI images and extract the base-subtracted signal
#'
#' @description
#' Convenience wrapper composing [read_cis()] and [extract_signal()].
#' Returns a pixels x participants signal matrix in one call. Use this
#' when you don't need to intervene between reading and extracting
#' (e.g. masking pixels, swapping the base).
#'
#' @section Reading the result:
#' Pixels x participants numeric matrix, base-subtracted. Carries an
#' `img_dims` attribute. Column names are the filename without the
#' extension.
#'
#' @section Common mistakes:
#' * Pointing `base_image_path` at the wrong file (e.g. a generated CI
#'   PNG instead of the base face used at stimulus-generation time).
#'   The dimensions check will catch a wholly wrong image but not a
#'   wrong-base-of-the-right-size.
#' * Treating the result as raw mask. See **Note** below.
#'
#' @inheritSection read_cis Raw vs. rendered CIs
#'
#' @param dir Directory of CI images (as for [read_cis()]).
#' @param base_image_path Path to the base face image used at
#'   stimulus-generation time.
#' @param pattern Regex for file selection (see [read_cis()]).
#' @param acknowledge_scaling Set to `TRUE` to silence the
#'   once-per-session warning. See [read_cis()].
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
                               pattern             = "\\.(png|jpe?g)$",
                               acknowledge_scaling = FALSE) {
  cis <- read_cis(dir, pattern = pattern,
                  acknowledge_scaling = acknowledge_scaling)
  extract_signal(cis, base_image_path = base_image_path,
                 acknowledge_scaling = TRUE)  # warning already fired
}
