## Internal helpers: validators, the rcicr soft-dep attacher, etc.
## None of these are exported.

#' Attach a package without using `library()`
#'
#' `library()` inside package code triggers R CMD check's
#' "Dependence on R services not declared". We need `foreach`,
#' `tibble`, and `dplyr` **attached** (not just loaded) when
#' `rcicr::generateStimuli2IFC()` or `computeInfoVal2IFC()` run
#' because they use `%dopar%`, `tribble()`, `%>%` and `filter()` at
#' evaluation time without namespace prefixes.
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
#' meaningful for split-half), warn at < 30.
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
      "i" = "Reliability estimates will be noisy; a sample size of \\
             N >= 60 is recommended for stable assessment."
    ))
  }
  invisible(x)
}

#' Validate and apply a region mask to a signal matrix
#'
#' Centralised handling for the `mask` argument now accepted by every
#' `rel_*()` function. Validates that `mask` is a logical vector of
#' the right length (matching `nrow(signal_matrix)`), then row-
#' subsets the signal matrix to the masked pixels.
#'
#' Returns the (possibly-subsetted) matrix, with the
#' `img_dims` attribute preserved when no subsetting happens (when
#' `mask = NULL`) and dropped when the mask reshapes the rows
#' (subsetted matrices are no longer 2D-image-shaped). The
#' `source` attribute (`"raw"` / `"rendered"`) is preserved across
#' subsetting because masking pixels does not change whether the
#' values are raw or scaled.
#'
#' @param signal_matrix Pixels x participants numeric matrix.
#' @param mask Optional logical vector of length `nrow(signal_matrix)`.
#' @param name Argument name of the matrix at the call site, for
#'   error messages.
#' @return The signal matrix unchanged (`mask = NULL`) or
#'   `signal_matrix[mask, , drop = FALSE]`.
#' @keywords internal
#' @noRd
apply_mask_to_signal <- function(signal_matrix, mask = NULL,
                                 name = "signal_matrix") {
  if (is.null(mask)) return(signal_matrix)
  if (!is.logical(mask)) {
    cli::cli_abort(
      "{.arg mask} must be a logical vector (got {.cls {class(mask)}})."
    )
  }
  if (length(mask) != nrow(signal_matrix)) {
    cli::cli_abort(c(
      "{.arg mask} length must match {.code nrow({name})}.",
      "*" = "{.arg mask}: {length(mask)}",
      "*" = "{.code nrow({name})}: {nrow(signal_matrix)}"
    ))
  }
  if (sum(mask) < 4L) {
    cli::cli_abort(c(
      "{.arg mask} selects too few pixels ({sum(mask)}).",
      "i" = "Reliability statistics on a few-pixel subset are \\
             dominated by noise. Use a more permissive mask."
    ))
  }
  src <- attr(signal_matrix, "source")
  out <- signal_matrix[mask, , drop = FALSE]
  # img_dims attribute is no longer meaningful after row-subset, drop it
  attr(out, "img_dims") <- NULL
  if (!is.null(src)) attr(out, "source") <- src
  out
}

#' Assert that a signal matrix is raw (not rendered/scaled)
#'
#' The shared enforcement point for variance-based metrics
#' (`pixel_t_test()`, `rel_icc()`, `rel_dissimilarity()`).
#' `rel_cluster_test()` inherits enforcement via its internal
#' `pixel_t_test()` call.
#'
#' Decision tree:
#' * `attr(x, "source") == "raw"`: pass.
#' * `attr(x, "source") == "rendered"`: error unless
#'   `acknowledge_scaling = TRUE` (variance-based metrics on scaled
#'   data give wrong numbers).
#' * No `source` attribute: fall back to `looks_scaled()` heuristic.
#'   If flagged, emit the once-per-session warning (do not error,
#'   we cannot be sure). If clean, pass silently.
#'
#' Why not a single enforcement point in `pixel_t_test()`:
#' `rel_icc()` and `rel_dissimilarity()` do not call
#' `pixel_t_test()`, so a single-point design would silently fail
#' to enforce on those code paths. Each variance-based metric calls
#' this helper directly.
#'
#' @param x Pixels x participants signal matrix.
#' @param acknowledge_scaling Logical. When `TRUE`, suppress the
#'   error/warning for known-rendered or heuristic-flagged input.
#' @param name Argument name at the call site, for error messages.
#' @return Invisibly `TRUE` if pass; `cli::cli_abort()`s if the
#'   matrix is `"rendered"` and `acknowledge_scaling` is `FALSE`.
#' @keywords internal
#' @noRd
assert_raw_signal <- function(x,
                              acknowledge_scaling = FALSE,
                              name = "signal_matrix") {
  src <- attr(x, "source")
  if (identical(src, "raw")) return(invisible(TRUE))
  if (identical(src, "rendered") && !isTRUE(acknowledge_scaling)) {
    cli::cli_abort(c(
      "{.arg {name}} is a {.strong rendered} CI (PNG-derived); \\
       variance-based metrics give wrong numbers on rendered data.",
      "i" = "Mode 2 ({.fn ci_from_responses_2ifc} / \\
             {.fn ci_from_responses_briefrc}) returns the raw mask.",
      "i" = "If you understand the trade-off and want to proceed \\
             anyway, pass {.code acknowledge_scaling = TRUE}.",
      "i" = "See {.code vignette(\"tutorial\", package = \"rcicrely\")} \\
             chapter 3 for the full discussion."
    ))
  }
  if (identical(src, "rendered") && isTRUE(acknowledge_scaling)) {
    return(invisible(TRUE))
  }
  # No source attribute: fall back to heuristic.
  if (!isTRUE(acknowledge_scaling) && looks_scaled(x)) {
    warn_looks_scaled(name)
  }
  invisible(TRUE)
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
#' `cli` ties a progress bar's lifetime to the environment in which
#' it was created and looks bars up by id within that environment. If
#' `progress_start()` and `progress_tick()` resolve to different
#' frames — e.g. when the loop runs inside a `with_seed({...})`
#' captured expression, or under knitr — the lookup misses with
#' `Cannot find progress bar 'cli-XXXX-YY'`. Fix: pin every bar to a
#' dedicated package-private env, which lives for the whole R
#' session and is the same object for `start`, `tick`, and `done`.
#' The env's parent is `baseenv()` so `cli`'s `glue`-based template
#' formatting can resolve `::` and other base R names.
#'
#' @keywords internal
#' @noRd
.rcicrely_progress_env <- new.env(parent = baseenv())

#' @keywords internal
#' @noRd
progress_start <- function(total, label, show = TRUE) {
  if (!isTRUE(show)) return(NULL)
  cli::cli_progress_bar(label, total = total, clear = TRUE,
                        .envir = .rcicrely_progress_env)
}

#' @keywords internal
#' @noRd
progress_tick <- function(id) {
  if (is.null(id)) return(invisible())
  cli::cli_progress_update(id = id, .envir = .rcicrely_progress_env)
}

#' @keywords internal
#' @noRd
progress_done <- function(id) {
  if (is.null(id)) return(invisible())
  cli::cli_progress_done(id = id, .envir = .rcicrely_progress_env)
}

## ---- session-state warning helpers ----------------------------------
##
## The "raw vs rendered" caveat is the load-bearing distinction in this
## package: PNG-derived signals contain `scaling(noise)`, not raw `noise`.
## We emit the warning once per session so the user sees it without being
## hammered every call. State lives in a package-private env that the
## test suite can reset via `reset_session_warnings()`.

.rcicrely_session_state <- new.env(parent = emptyenv())
.rcicrely_session_state$mode1_warning_emitted        <- FALSE
.rcicrely_session_state$looks_scaled_warning_emitted <- FALSE
.rcicrely_session_state$icc_resolution_warning_emitted <- FALSE
.rcicrely_session_state$loo_sd_deprecated_emitted    <- FALSE

#' @keywords internal
#' @noRd
reset_session_warnings <- function() {
  .rcicrely_session_state$mode1_warning_emitted          <- FALSE
  .rcicrely_session_state$looks_scaled_warning_emitted   <- FALSE
  .rcicrely_session_state$icc_resolution_warning_emitted <- FALSE
  .rcicrely_session_state$loo_sd_deprecated_emitted      <- FALSE
  invisible()
}

#' Emit the once-per-session rel_loo(flag_method = "sd") deprecation
#'
#' Default switched to `"mad"` in v0.3; `"sd"` retained for
#' backwards compatibility and slated for removal in v0.4. The MAD
#' rule is robust to the influential producers LOO is designed to
#' detect; the SD rule is mask-prone (one outlier inflates `sd_r`
#' and pulls `mean_r`, hiding itself behind a `z` near `-1`).
#'
#' @keywords internal
#' @noRd
warn_loo_sd_deprecated <- function() {
  if (isTRUE(getOption("rcicrely.silence_loo_deprecation", FALSE))) {
    return(invisible())
  }
  if (isTRUE(.rcicrely_session_state$loo_sd_deprecated_emitted)) {
    return(invisible())
  }
  cli::cli_warn(c(
    "{.code rel_loo(flag_method = \"sd\")} is deprecated since v0.3 \\
     and will be removed in v0.4.",
    "i" = "The MAD/median rule (the new default) is robust to the \\
           influential producers LOO is designed to flag.",
    "*" = "Drop the {.arg flag_method} argument or pass \\
           {.code flag_method = \"mad\"} explicitly.",
    "i" = "Silence: {.code options(rcicrely.silence_loo_deprecation = TRUE)}."
  ))
  .rcicrely_session_state$loo_sd_deprecated_emitted <- TRUE
  invisible()
}

#' Emit the once-per-session ICC(3,k) resolution-asymptote warning
#'
#' Called by `rel_icc()` when the image has more than 50,000 pixels.
#' Flags that ICC(3,k) approaches 1 at large image sizes and is not
#' resolution-comparable.
#'
#' @keywords internal
#' @noRd
warn_icc_large_image <- function(n_targets) {
  if (isTRUE(getOption("rcicrely.silence_icc_warning", FALSE))) {
    return(invisible())
  }
  if (isTRUE(.rcicrely_session_state$icc_resolution_warning_emitted)) {
    return(invisible())
  }
  cli::cli_warn(c(
    "ICC(3,k) tends toward 1 at large image sizes ({n_targets} pixels).",
    "i" = "This is a property of ICC(*,k), not a reliability statement: \\
           the error mean square scales as 1/((n-1)(k-1)) while the \\
           row mean square scales as 1/(n-1), so ICC(3,k) asymptotes \\
           to 1 as n grows.",
    "*" = "Report {.strong ICC(3,1)} as the primary, resolution-comparable \\
           statistic when comparing across image sizes.",
    "i" = "Silence: {.code options(rcicrely.silence_icc_warning = TRUE)}."
  ))
  .rcicrely_session_state$icc_resolution_warning_emitted <- TRUE
  invisible()
}

#' Emit the once-per-session Mode-1 scaling warning
#'
#' Called by `read_cis()`, `extract_signal()`, `load_signal_matrix()`.
#' Silent if `acknowledge_scaling = TRUE`, if the option
#' `rcicrely.silence_scaling_warning` is `TRUE`, or if the warning has
#' already fired in this session.
#'
#' @keywords internal
#' @noRd
warn_mode1_scaling <- function(acknowledge_scaling = FALSE) {
  if (isTRUE(acknowledge_scaling)) return(invisible())
  if (isTRUE(getOption("rcicrely.silence_scaling_warning", FALSE))) {
    return(invisible())
  }
  if (isTRUE(.rcicrely_session_state$mode1_warning_emitted)) {
    return(invisible())
  }
  cli::cli_warn(c(
    "PNG-derived signal matrix contains the {.strong rendered} CI, \\
     not the raw mask.",
    "i" = "PNGs encode {.code base + scaling(mask)}, so \\
           {.code cis - base} recovers {.code scaling(mask)}, not the \\
           raw mask itself.",
    "*" = "Pearson-based metrics ({.fn rel_split_half}, \\
           {.fn rel_loo}, correlation half of {.fn rel_dissimilarity}) \\
           survive a uniform linear scaling but are still distorted by \\
           per-CI {.val matched} scaling.",
    "*" = "{.fn rel_icc}, Euclidean half of {.fn rel_dissimilarity}, \\
           {.fn pixel_t_test}, {.fn rel_cluster_test} and any \\
           hand-rolled {.code infoVal} computation (e.g. for Brief-RC, \\
           which has no upstream implementation) are sensitive to \\
           {.strong any} scaling.",
    "i" = "{.code rcicr::computeInfoVal2IFC()} extracts the raw \\
           {.field $ci} component internally, so the standard 2IFC \\
           infoVal path is unaffected.",
    "i" = "Prefer Mode 2 ({.fn ci_from_responses_2ifc} / \\
           {.fn ci_from_responses_briefrc}) when raw responses are \\
           available, those return the raw mask.",
    "i" = "Silence: pass {.code acknowledge_scaling = TRUE}, set \\
           {.code options(rcicrely.silence_scaling_warning = TRUE)}, \\
           or generate PNGs with {.code scaling = \"none\"} so the \\
           output is effectively raw.",
    "i" = "See {.code vignette(\"tutorial\", package = \"rcicrely\")} \\
           chapter 3 for the full discussion."
  ))
  .rcicrely_session_state$mode1_warning_emitted <- TRUE
  invisible()
}

#' Emit the once-per-session "looks scaled" warning
#'
#' Called by every `rel_*()` and `run_*()` on entry when
#' `looks_scaled()` flags the input. Quieter than the Mode-1 warning
#' (which targets the input boundary directly).
#'
#' @keywords internal
#' @noRd
warn_looks_scaled <- function(name = "signal_matrix") {
  if (isTRUE(getOption("rcicrely.silence_scaling_warning", FALSE))) {
    return(invisible())
  }
  if (isTRUE(.rcicrely_session_state$looks_scaled_warning_emitted)) {
    return(invisible())
  }
  cli::cli_warn(c(
    "{.arg {name}} looks like it may be a {.strong rendered} CI \\
     (per-column dynamic range is highly heterogeneous).",
    "i" = "If the matrix came from {.fn read_cis} / \\
           {.fn extract_signal} on rendered PNGs, downstream metrics \\
           may be distorted. Mode 2 ({.fn ci_from_responses_2ifc} / \\
           {.fn ci_from_responses_briefrc}) returns the raw mask.",
    "i" = "Heuristic only, silence with \\
           {.code options(rcicrely.silence_scaling_warning = TRUE)} \\
           if the matrix is genuinely raw."
  ))
  .rcicrely_session_state$looks_scaled_warning_emitted <- TRUE
  invisible()
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
