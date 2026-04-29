## Per-result-class regression list. When `print.rcicrely_*()` is
## called on a result object whose `attr(., "rcicrely_version")`
## matches a known-bad version for that class, emit a one-time
## per-session warning telling the user to recompute.
##
## Add an entry by calling `register_known_regression()` at package
## load time (see .onLoad in zzz.R), or by listing it directly in
## `.rcicrely_regressions` below.
##
## Each entry: list(class = <result class>, version = <package_version>,
## message = <user-facing string>).

.rcicrely_regressions <- list(
  list(
    class   = "rcicrely_infoval",
    version = package_version("0.2.0"),
    message = "Result was computed by rcicrely 0.2.0, which sampled \\
               reference stim ids with replacement. This biased \\
               z-scores downward (typically by ~1 unit at \\
               n_trials == n_pool). Recompute with rcicrely \\
               >= 0.2.1 for correct calibration."
  )
)

#' Lookup table of known regressions for a given result class
#'
#' Returns the list of `(version, message)` pairs registered for the
#' given S3 class, or `NULL` if none.
#'
#' @keywords internal
#' @noRd
known_regressions_for <- function(cls) {
  hits <- vapply(.rcicrely_regressions,
                 function(e) identical(e$class, cls),
                 logical(1L))
  if (!any(hits)) return(NULL)
  .rcicrely_regressions[hits]
}

#' Emit the once-per-session "result computed by buggy version" warning
#'
#' Called from `print.rcicrely_*()` after looking up the result's
#' `rcicrely_version` attribute against `.rcicrely_regressions`.
#'
#' @keywords internal
#' @noRd
warn_known_regression <- function(x) {
  if (isTRUE(getOption("rcicrely.silence_regression_warning",
                       FALSE))) {
    return(invisible())
  }
  cls <- class(x)[[1L]]
  v   <- attr(x, "rcicrely_version")
  if (is.null(v)) return(invisible())
  hits <- known_regressions_for(cls)
  if (is.null(hits)) return(invisible())
  for (h in hits) {
    if (identical(v, h$version)) {
      key <- paste0("regression_", cls, "_", as.character(v))
      emitted <- isTRUE(.rcicrely_session_state[[key]])
      if (emitted) next
      cli::cli_warn(c(
        "{.cls {cls}} result was computed by rcicrely \\
         {as.character(v)} (known regression).",
        "i" = h$message,
        "i" = "Silence: \\
               {.code options(rcicrely.silence_regression_warning = TRUE)}."
      ))
      .rcicrely_session_state[[key]] <- TRUE
    }
  }
  invisible()
}
