#' Load the noise matrix from an rcicr `.Rdata` file or a text file
#'
#' rcicr's [rcicr::generateStimuli2IFC()] writes a timestamped `.Rdata`
#' containing a `stimuli` object (pixels x n_trials data frame / matrix
#' of noise patterns). For Brief-RC workflows this same pool of noise
#' patterns is the input to the native mask computation. This loader
#' accepts either:
#'
#' - an `.Rdata` file as produced by rcicr  --  `stimuli` is extracted and
#'   coerced to a numeric matrix;
#' - a whitespace- or comma-delimited text file  --  read with
#'   [data.table::fread()], then coerced.
#'
#' rcicr writes its file with a lowercase `.Rdata` extension on some
#' filesystems and `.RData` on others. The loader matches both by
#' examining the first magic bytes (gzip / RData markers) rather than
#' trusting the extension (CLAUDE.md sec.8.8).
#'
#' @param path Path to the file.
#' @param stimuli_object Name of the object inside the `.Rdata` that
#'   holds the noise patterns. Defaults to `"stimuli"`, the rcicr
#'   v1.0.1 convention.
#' @return A numeric pixels x n_trials matrix.
#' @export
#' @examples
#' \dontrun{
#' nm <- read_noise_matrix("data/rcicr_bogus_20250101_1200.Rdata")
#' }
read_noise_matrix <- function(path, stimuli_object = "stimuli") {
  if (!file.exists(path)) {
    cli::cli_abort("File not found: {.path {path}}")
  }
  # sniff the first few bytes; .Rdata / .RData start with gzip magic
  # (0x1f 0x8b) for compressed, or "RDX2" / "RDX3" for uncompressed.
  con <- file(path, open = "rb")
  on.exit(close(con), add = TRUE)
  magic <- readBin(con, "raw", n = 4L)
  close(con)
  on.exit()

  is_gzip <- length(magic) >= 2L &&
    magic[1] == as.raw(0x1f) && magic[2] == as.raw(0x8b)
  is_rdx <- length(magic) >= 4L &&
    rawToChar(magic[1:4]) %in% c("RDX2", "RDX3", "RDA2", "RDA3")
  is_rdata <- is_gzip || is_rdx

  if (is_rdata) {
    env <- new.env(parent = emptyenv())
    load(path, envir = env)
    if (!stimuli_object %in% ls(env)) {
      cli::cli_abort(c(
        "Object {.var {stimuli_object}} not found in {.path {path}}.",
        "i" = "Objects present: {.var {ls(env)}}."
      ))
    }
    nm <- as.matrix(env[[stimuli_object]])
    if (!is.numeric(nm)) {
      storage.mode(nm) <- "double"
    }
    return(nm)
  }

  # fall through: text file
  tbl <- data.table::fread(path, header = FALSE)
  nm <- as.matrix(tbl)
  dimnames(nm) <- NULL   # fread assigns V1, V2, ... ; drop for cleanliness
  if (!is.numeric(nm)) {
    storage.mode(nm) <- "double"
  }
  nm
}
