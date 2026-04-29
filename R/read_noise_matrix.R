#' Load the noise matrix from an rcicr `.Rdata`, an `.rds`, or a text file
#'
#' @description
#' Produces the pixels x n_trials numeric matrix of individual noise
#' patterns that [ci_from_responses_briefrc()] needs. Use this when
#' you want to load (or reconstruct) the noise pool once and pass it
#' to multiple Brief-RC CI computations without re-loading.
#'
#' Three input formats are supported:
#'
#' 1. **rcicr `.Rdata`** (the output of
#'    [rcicr::generateStimuli2IFC()]). The rdata does not store the
#'    noise matrix directly; it stores `stimuli_params` (per-trial
#'    parameter rows) and `p` (the sinusoid / gabor basis). This loader
#'    reconstructs each trial's noise pattern via
#'    [rcicr::generateNoiseImage()] and stacks them. Requires `rcicr`
#'    to be installed.
#' 2. **R `.rds`** previously saved from
#'    `rcicr::generateStimuli2IFC(return_as_dataframe = TRUE)`. No
#'    reconstruction needed.
#' 3. **Plain text** (whitespace- or comma-delimited) matrix, one
#'    column per trial.
#'
#' The loader dispatches on file magic bytes (not extension; rcicr
#' writes lowercase `.Rdata` on some filesystems, which trips
#' case-sensitive matching).
#'
#' @section Reading the result:
#' Numeric pixels x n_trials matrix. Each column is one trial's noise
#' pattern (a sinusoidal or gabor combination of the basis `p`).
#' Pool ids referenced by Brief-RC `responses$stimulus` are column
#' indices into this matrix.
#'
#' @section Common mistakes:
#' * Pointing at the rcicr `.Rdata` and assuming it contains a
#'   `stimuli` object, it doesn't. The loader reconstructs by
#'   calling [rcicr::generateNoiseImage()] per row, which is correct
#'   but slow (~4s per 120 trials at 64x64).
#' * Caching the captured return of
#'   `rcicr::generateStimuli2IFC(return_as_dataframe = TRUE)` to
#'   `.rds` once and reusing it is much faster than reconstructing
#'   from the rdata each run.
#' * Multiple stale `.rdata` accumulating in a directory because
#'   rcicr timestamps its filenames; pick the most recent by mtime if
#'   you have several.
#' * On macOS the file is saved with a lowercase `.Rdata` extension;
#'   `list.files(pattern = "\\.RData$")` will miss it (use
#'   `ignore.case = TRUE`).
#'
#' @section What is in an rcicr `.Rdata`:
#' `load("rcic_stimuli.Rdata")` adds the following objects to your
#' environment. Names are short and not self-explanatory:
#'
#' \describe{
#'   \item{`base_face_files`}{Named list of file paths to the original
#'     face images. The list names (e.g. `"base"`) are the labels you
#'     pass downstream as `baseimage = "..."`.}
#'   \item{`base_faces`}{Base images themselves, loaded as numeric
#'     matrices of grayscale pixels in `[0, 1]`. One matrix per label;
#'     shape `img_size x img_size`.}
#'   \item{`img_size`}{Side length of the (square) images in pixels.}
#'   \item{`n_trials`}{Number of stimulus pairs (and unique noise
#'     patterns) generated.}
#'   \item{`noise_type`}{Type of noise basis. The rcicr default and only
#'     well-tested option is `"sinusoid"`.}
#'   \item{`p`}{The noise basis. A list with `$patches` (a stack of
#'     standard sinusoidal patterns at multiple scales, orientations,
#'     phases) and `$patchIdx` (an index that maps positions in a
#'     parameter vector to entries in `$patches`). Think of `$patches`
#'     as a dictionary of sinusoidal "ingredients".}
#'   \item{`stimuli_params`}{The per-trial recipe for combining `p`'s
#'     ingredients. A named list of matrices, one per base face. Each
#'     row is one trial; each entry is a contrast weight.
#'     `rcicr::generateNoiseImage(stimuli_params[[base]][i, ], p)`
#'     reconstructs trial `i`.}
#'   \item{`seed`}{Random seed used at generation, for reproducibility.}
#'   \item{`label`, `stimulus_path`, `trial`, `generator_version`,
#'     `use_same_parameters`}{Bookkeeping; not consumed by analysis.}
#'   \item{`reference_norms`}{Random-responder Frobenius norms. Not
#'     present at first; created and inserted in place the first time
#'     `rcicr::computeInfoVal2IFC()` is called on the file. Copy the
#'     rdata first if you want it untouched.}
#' }
#'
#' The load-bearing objects for analysis are `base_faces`,
#' `stimuli_params`, `p`, and `img_size`. Trial-level noise images are
#' not stored; this loader recomputes them on demand via
#' `generateNoiseImage()`.
#'
#' @param path Path to the file.
#' @param baseimage For rcicr `.Rdata` inputs: which base label to
#'   reconstruct noise for. Defaults to the first key in `base_faces`
#'   if there is exactly one.
#' @param stimuli_object For `.Rdata` files that contain a pre-saved
#'   `stimuli` object (rare, rcicr does not save one by default),
#'   the name to look up. Defaults to `"stimuli"`.
#' @return A numeric pixels x n_trials matrix.
#' @seealso [ci_from_responses_briefrc()]
#' @export
#' @examples
#' \dontrun{
#' # Plain text (Schmitz et al. 2024 OSF format):
#' nm <- read_noise_matrix("data/noise_matrix.txt")
#'
#' # Cached .rds (from rcicr's return-as-dataframe output):
#' nm <- read_noise_matrix("data/noise_matrix.rds")
#'
#' # rcicr .Rdata (reconstructs each trial; slow):
#' nm <- read_noise_matrix("data/rcicr_stimuli.Rdata")
#' }
read_noise_matrix <- function(path, baseimage = NULL,
                              stimuli_object = "stimuli") {
  if (!file.exists(path)) {
    cli::cli_abort("File not found: {.path {path}}")
  }

  # sniff the first few bytes; .Rdata / .rds start with gzip magic
  # (0x1f 0x8b) for compressed, or RDX/RDA for uncompressed.
  con <- file(path, open = "rb")
  magic <- readBin(con, "raw", n = 4L)
  close(con)

  is_gzip <- length(magic) >= 2L &&
    magic[1] == as.raw(0x1f) && magic[2] == as.raw(0x8b)
  is_rdx <- length(magic) >= 4L &&
    rawToChar(magic[1:4]) %in% c("RDX2", "RDX3", "RDA2", "RDA3")
  is_binary <- is_gzip || is_rdx

  if (is_binary) {
    # try rds first
    rds_obj <- tryCatch(readRDS(path), error = function(e) NULL)
    if (!is.null(rds_obj)) {
      nm <- as.matrix(rds_obj)
      dimnames(nm) <- NULL
      storage.mode(nm) <- "double"
      return(nm)
    }
    # fall through: RData
    env <- new.env(parent = emptyenv())
    load_ok <- tryCatch({ load(path, envir = env); TRUE },
                       error = function(e) FALSE)
    if (!load_ok) {
      # gzip magic without an R header may be a gzipped text matrix
      if (is_gzip && !is_rdx) {
        text_nm <- tryCatch({
          tbl <- data.table::fread(path, header = FALSE)
          if (ncol(tbl) > 0L &&
                all(vapply(tbl, is.numeric, logical(1L)))) {
            m <- as.matrix(tbl); dimnames(m) <- NULL
            storage.mode(m) <- "double"
            m
          } else NULL
        }, error = function(e) NULL)
        if (!is.null(text_nm)) return(text_nm)
      }
      cli::cli_abort(c(
        "Could not read {.path {path}} as .rds, .Rdata, or gzipped text.",
        "i" = "Magic bytes suggest a binary R file, but none of \\
               {.fn readRDS}, {.fn load}, or {.fn data.table::fread} \\
               could parse it."
      ))
    }
    objs <- ls(env)

    # Path 1: user pre-saved a `stimuli` object
    if (stimuli_object %in% objs) {
      nm <- as.matrix(env[[stimuli_object]])
      dimnames(nm) <- NULL
      storage.mode(nm) <- "double"
      return(nm)
    }

    # Path 2: rcicr rdata - reconstruct from stimuli_params + p
    if (all(c("stimuli_params", "p", "img_size") %in% objs)) {
      if (!requireNamespace("rcicr", quietly = TRUE)) {
        cli::cli_abort(c(
          "Reconstructing the noise matrix from an rcicr rdata \\
           requires the {.pkg rcicr} package.",
          "i" = "Install: \\
                 {.run remotes::install_github(\"rdotsch/rcicr\")}"
        ))
      }
      labels <- names(env$stimuli_params)
      if (is.null(baseimage)) {
        if (length(labels) != 1L) {
          cli::cli_abort(c(
            "Multiple base labels in rdata; pick one via \\
             {.arg baseimage}.",
            "i" = "Available: {.val {labels}}"
          ))
        }
        baseimage <- labels[1L]
      }
      if (!baseimage %in% labels) {
        cli::cli_abort(c(
          "{.arg baseimage} = {.val {baseimage}} not in rdata.",
          "i" = "Available: {.val {labels}}"
        ))
      }
      params <- env$stimuli_params[[baseimage]]
      p_basis <- env$p
      img_size <- env$img_size
      n_trials <- nrow(params)
      n_pix <- as.integer(img_size) * as.integer(img_size)
      nm <- matrix(NA_real_, nrow = n_pix, ncol = n_trials)
      for (i in seq_len(n_trials)) {
        nm[, i] <- as.vector(rcicr::generateNoiseImage(
          params[i, ], p_basis
        ))
      }
      return(nm)
    }

    cli::cli_abort(c(
      "Unrecognised objects in {.path {path}}.",
      "i" = "Expected either {.var stimuli} (pre-saved matrix) or \\
             {.var stimuli_params} + {.var p} + {.var img_size} \\
             (rcicr rdata).",
      "*" = "Found: {.val {objs}}"
    ))
  }

  # text file fallback
  tbl <- data.table::fread(path, header = FALSE)
  nm <- as.matrix(tbl)
  dimnames(nm) <- NULL
  storage.mode(nm) <- "double"
  nm
}
