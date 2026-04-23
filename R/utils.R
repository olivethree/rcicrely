## Internal helpers: validators, the rcicr soft-dep attacher, etc.
## None of these are exported.

#' Attach a package without using `library()`
#'
#' `library()` inside package code triggers R CMD check's
#' "Dependence on R services not declared". We need `foreach`,
#' `tibble`, and `dplyr` **attached** (not just loaded) when
#' `rcicr::generateStimuli2IFC()` or `computeInfoVal2IFC()` run
#' because they use `%dopar%`, `tribble()`, `%>%` and `filter()` at
#' evaluation time without namespace prefixes. CLAUDE.md sec.8.3, sec.8.7.
#'
#' @param pkgs Character vector of package names to attach.
#' @return Invisibly `TRUE` if all attached; errors if any missing.
#' @keywords internal
#' @noRd
ensure_attached <- function(pkgs) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      cli::cli_abort(
        "Package {.pkg {pkg}} is required here but not installed."
      )
    }
    attached <- paste0("package:", pkg) %in% search()
    if (!attached) {
      suppressPackageStartupMessages(attachNamespace(pkg))
    }
  }
  invisible(TRUE)
}

#' Assert a pixels x participants signal matrix is well-formed
#'
#' Every `rel_*` function takes a signal matrix. This centralises the
#' common validation: numeric matrix, at least 4 participants (minimum
#' meaningful for split-half), warn at < 30 per Cone et al. (2020).
#'
#' @param x Candidate signal matrix.
#' @param name Argument name, for diagnostic messages.
#' @param min_participants Hard floor below which we abort.
#' @return Invisibly the input.
#' @keywords internal
#' @noRd
validate_signal_matrix <- function(x,
                                   name = "signal_matrix",
                                   min_participants = 4L) {
  if (!is.matrix(x) || !is.numeric(x)) {
    cli::cli_abort(c(
      "{.arg {name}} must be a numeric matrix.",
      "i" = "Got {.cls {class(x)}}."
    ))
  }
  if (ncol(x) < min_participants) {
    cli::cli_abort(c(
      "{.arg {name}} must have at least {min_participants} participants \\
       (columns).",
      "i" = "Got {ncol(x)}. Reliability metrics are undefined below this."
    ))
  }
  if (ncol(x) < 30L) {
    cli::cli_warn(c(
      "{.arg {name}} has fewer than 30 participants ({ncol(x)}).",
      "i" = "Reliability estimates will be noisy; Cone et al. (2020) \\
             recommend N >= 60 for stable assessment."
    ))
  }
  invisible(x)
}

#' Assert two signal matrices have compatible shape for between-condition
#'
#' @keywords internal
#' @noRd
validate_two_signal_matrices <- function(a, b,
                                         name_a = "signal_matrix_a",
                                         name_b = "signal_matrix_b") {
  validate_signal_matrix(a, name_a)
  validate_signal_matrix(b, name_b)
  if (nrow(a) != nrow(b)) {
    cli::cli_abort(c(
      "Pixel counts differ between conditions.",
      "*" = "{.arg {name_a}}: {nrow(a)} pixels",
      "*" = "{.arg {name_b}}: {nrow(b)} pixels"
    ))
  }
  invisible(TRUE)
}

#' Validate and normalise an `img_dims` argument
#'
#' Accepts `c(nrow, ncol)` or a single integer (interpreted as square
#' `c(n, n)`). Returns integer length-2.
#'
#' @keywords internal
#' @noRd
validate_img_dims <- function(img_dims, n_pixels,
                              name = "img_dims") {
  if (length(img_dims) == 1L) {
    img_dims <- c(img_dims, img_dims)
  }
  if (length(img_dims) != 2L || any(!is.finite(img_dims)) ||
        any(img_dims < 1L)) {
    cli::cli_abort(
      "{.arg {name}} must be a positive integer of length 1 or 2 \\
       (`c(nrow, ncol)`)."
    )
  }
  img_dims <- as.integer(img_dims)
  if (prod(img_dims) != n_pixels) {
    cli::cli_abort(c(
      "{.arg {name}} is inconsistent with the signal matrix.",
      "i" = "{img_dims[1]} x {img_dims[2]} = {prod(img_dims)} != \\
             {n_pixels} pixels."
    ))
  }
  img_dims
}

#' Optionally call `set.seed()` inside a function without leaking state
#'
#' Saves the caller's RNG state, sets the package seed for the
#' duration of `expr`, and restores the state on exit. If `seed` is
#' `NULL`, `expr` runs with the caller's current RNG.
#'
#' @keywords internal
#' @noRd
with_seed <- function(seed, expr) {
  if (is.null(seed)) return(force(expr))
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    old <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit(assign(".Random.seed", old, envir = .GlobalEnv), add = TRUE)
  } else {
    on.exit(
      rm(".Random.seed", envir = .GlobalEnv, inherits = FALSE),
      add = TRUE
    )
  }
  set.seed(seed)
  force(expr)
}

#' Optional progress ticker using `cli::cli_progress_*`
#'
#' @keywords internal
#' @noRd
progress_start <- function(total, label, show = TRUE) {
  if (!isTRUE(show)) return(NULL)
  cli::cli_progress_bar(label, total = total, clear = TRUE)
}

#' @keywords internal
#' @noRd
progress_tick <- function(id) {
  if (is.null(id)) return(invisible())
  cli::cli_progress_update(id = id)
}

#' @keywords internal
#' @noRd
progress_done <- function(id) {
  if (is.null(id)) return(invisible())
  cli::cli_progress_done(id = id)
}

#' Read and convert an image file to a grayscale numeric matrix
#'
#' Handles PNG and JPEG via the optional `png` / `jpeg` packages,
#' collapses RGB channels to luminance via ITU-R BT.709 weights if
#' needed, returns an `nrow x ncol` numeric matrix with values in
#' the 0-1 range. Internal helper used by `read_cis()` and
#' `extract_signal()`.
#'
#' @keywords internal
#' @noRd
read_image_as_gray <- function(path) {
  if (!file.exists(path)) {
    cli::cli_abort("Image file not found: {.path {path}}")
  }
  ext <- tolower(tools::file_ext(path))
  img <- switch(
    ext,
    png = {
      if (!requireNamespace("png", quietly = TRUE)) {
        cli::cli_abort(
          "Reading PNG files requires the {.pkg png} package."
        )
      }
      png::readPNG(path)
    },
    jpg = ,
    jpeg = {
      if (!requireNamespace("jpeg", quietly = TRUE)) {
        cli::cli_abort(
          "Reading JPEG files requires the {.pkg jpeg} package."
        )
      }
      jpeg::readJPEG(path)
    },
    cli::cli_abort(
      "Unsupported image extension {.val {ext}} for {.path {path}}."
    )
  )
  if (length(dim(img)) == 2L) {
    return(img)
  }
  nch <- dim(img)[3]
  if (nch >= 3L) {
    0.2126 * img[, , 1] + 0.7152 * img[, , 2] + 0.0722 * img[, , 3]
  } else {
    img[, , 1]
  }
}
