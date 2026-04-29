#' Compute individual 2IFC CIs from trial-level responses
#'
#' @description
#' Thin wrapper around [rcicr::batchGenerateCI2IFC()] that handles
#' the known rcicr-integration gotchas and exposes a uniform return
#' shape across the 2IFC and Brief-RC paths.
#'
#' Use this when you have 2IFC trial-level responses and want the
#' package to compute individual CIs ready for the reliability metrics.
#'
#' @details
#' What the wrapper does for you:
#' - Attaches `foreach`, `tibble`, `dplyr` at runtime
#'   (`%dopar%` / `tribble` / `%>%` are not namespace-prefixed inside
#'   rcicr).
#' - Matches the `.Rdata` extension case-insensitively (rcicr writes
#'   lowercase on some filesystems).
#' - Extracts the per-participant CI **noise component** (`$ci`) from
#'   the rcicr result list and stacks it into a pixels x participants
#'   signal matrix, already base-subtracted, ready for `rel_*()`.
#' - When `keep_rendered = TRUE`, also extracts the rendered
#'   `$combined` image and stacks it as `$rendered_ci` for
#'   visualisation (display only, not for downstream stats).
#'
#' `rcicr` must be installed (it's a Suggests dep; install with
#' `remotes::install_github("rdotsch/rcicr")` if needed).
#'
#' @section Reading the result:
#' * `$signal_matrix` is the **raw mask** (rcicr's `$ci` per producer).
#'   This is what every `rel_*` function expects.
#' * `$rendered_ci`, present when `keep_rendered = TRUE`, is the
#'   `base + scaling(mask)` image rcicr would have written to PNG.
#'   **Visualisation only.**
#' * `$rcicr_result` is the raw return value of
#'   [rcicr::batchGenerateCI2IFC()] in case you need anything else
#'   (e.g. the per-producer base, scaled-only ci).
#'
#' @section Common mistakes:
#' * Passing `$rendered_ci` to `rel_*` functions, it carries the
#'   scaling step, which distorts variance-based metrics.
#' * Forgetting that the warning about Mode 1 also applies if you
#'   later read `$rendered_ci` back from PNG instead of using
#'   `$signal_matrix` directly.
#'
#' @param responses A data frame with one row per trial. Must contain
#'   `participant_col` (producer id), the column named by `stimulus_col`
#'   (stimulus id, integer, index into the rcicr noise pool), and
#'   `response_col` with values in `{-1, +1}`.
#' @param rdata_path Path to the `.Rdata` file produced by
#'   [rcicr::generateStimuli2IFC()].
#' @param baseimage Base image label used at stimulus-generation time
#'   (`baseimage = "base"` is the convention in the test harness;
#'   whatever you passed to `generateStimuli2IFC(base_face_files =
#'   list(<label> = path))`). If `NULL`, tries to read the single base
#'   stored in `rdata_path`.
#' @param participant_col,stimulus_col,response_col Column names in
#'   `responses`.
#' @param scaling rcicr scaling option; one of `"autoscale"`,
#'   `"independent"`, `"constant"`, `"none"`. Passed through to
#'   `rcicr::batchGenerateCI2IFC()`. Scaling only affects the
#'   `$combined` (rendered) image, not `$ci`, so the returned
#'   `$signal_matrix` is the raw mask regardless of this argument.
#' @param keep_rendered If `TRUE`, also extract rcicr's `$combined`
#'   image (base + scaled noise) per producer and stack into
#'   `$rendered_ci`. Default `FALSE`. Use only for visualisation; do
#'   not feed `$rendered_ci` to `rel_*` functions.
#' @param targetpath Where rcicr writes PNGs (defaults to an
#'   auto-cleaned tempdir so the working directory isn't polluted).
#' @param save_as_png Whether rcicr writes individual CI PNGs.
#'   Defaults to `FALSE` for speed in a pure reliability pipeline.
#' @return A list with
#'   * `signal_matrix`, pixels x participants raw mask (always raw,
#'     regardless of `scaling`).
#'   * `rendered_ci`, pixels x participants, present only when
#'     `keep_rendered = TRUE`. Visualisation only.
#'   * `participants`, character vector of participant ids in the
#'     column order of `signal_matrix`.
#'   * `img_dims`, integer `c(nrow, ncol)`.
#'   * `scaling`, the `scaling` value passed to rcicr.
#'   * `rcicr_result`, the raw rcicr list (one element per producer).
#' @seealso [ci_from_responses_briefrc()], [rcicr::batchGenerateCI2IFC()]
#' @export
#' @examples
#' \dontrun{
#' res <- ci_from_responses_2ifc(
#'   responses  = my_responses,
#'   rdata_path = "data/rcicr_stimuli.Rdata",
#'   baseimage  = "base"
#' )
#' dim(res$signal_matrix)  # n_pixels x n_participants
#' }
ci_from_responses_2ifc <- function(responses,
                                   rdata_path,
                                   baseimage       = NULL,
                                   participant_col = "participant_id",
                                   stimulus_col    = "stimulus",
                                   response_col    = "response",
                                   scaling         = "autoscale",
                                   keep_rendered   = FALSE,
                                   targetpath      = tempfile("rcicrely_2ifc_"),
                                   save_as_png     = FALSE) {
  if (!requireNamespace("rcicr", quietly = TRUE)) {
    cli::cli_abort(c(
      "The 2IFC CI path requires the {.pkg rcicr} package.",
      "i" = "Install via: \\
             {.run remotes::install_github(\"rdotsch/rcicr\")}"
    ))
  }

  responses <- as.data.frame(responses)
  required <- c(participant_col, stimulus_col, response_col)
  missing_cols <- setdiff(required, colnames(responses))
  if (length(missing_cols) > 0L) {
    n_missing <- length(missing_cols)
    cli::cli_abort(c(
      "Missing {n_missing} column{?s} in {.arg responses}:",
      "*" = "{.val {missing_cols}}",
      "i" = "Have: {.val {colnames(responses)}}"
    ))
  }

  unique_resp <- sort(unique(as.numeric(responses[[response_col]])))
  if (!identical(unique_resp, c(-1, 1))) {
    msg <- c(
      "Column {.var {response_col}} must contain only values in \\
       {.val {c(-1, 1)}}.",
      "*" = "Got: {.val {unique_resp}}"
    )
    if (identical(unique_resp, c(0, 1)) || identical(unique_resp, c(0, 1L))) {
      msg <- c(
        msg,
        "i" = "Did you mean {.code -1 / +1}? The {.code {{0, 1}}} \\
               coding is the most common silent failure in 2IFC \\
               pipelines (often produced by experiment software \\
               that records 'left' / 'right' as 0 / 1).",
        "i" = "Recode in one line: \\
               {.code responses${response_col} <- 2 * responses${response_col} - 1}"
      )
    }
    cli::cli_abort(msg)
  }

  if (!file.exists(rdata_path)) {
    cli::cli_abort("rdata not found: {.path {rdata_path}}")
  }

  # determine baseimage label
  env <- new.env(parent = emptyenv())
  load(rdata_path, envir = env)
  if (is.null(baseimage)) {
    if (!"base_faces" %in% ls(env)) {
      cli::cli_abort(c(
        "Cannot infer {.arg baseimage}, no {.var base_faces} in \\
         {.path {rdata_path}}.",
        "i" = "Pass {.arg baseimage} explicitly."
      ))
    }
    labels <- names(env$base_faces)
    if (length(labels) != 1L) {
      cli::cli_abort(c(
        "Multiple base images in rdata; pick one via {.arg baseimage}.",
        "i" = "Available: {.val {labels}}"
      ))
    }
    baseimage <- labels[1L]
  }
  if (!baseimage %in% names(env$base_faces)) {
    cli::cli_abort(c(
      "{.arg baseimage} = {.val {baseimage}} not found in rdata.",
      "i" = "Available: {.val {names(env$base_faces)}}"
    ))
  }
  base_mat <- env$base_faces[[baseimage]]
  img_dims <- as.integer(dim(base_mat))

  # rcicr uses %dopar%, tribble, %>% without namespace prefixes
  ensure_attached(c("foreach", "tibble", "dplyr"))

  dir.create(targetpath, recursive = TRUE, showWarnings = FALSE)

  cis <- rcicr::batchGenerateCI2IFC(
    data        = responses,
    by          = participant_col,
    stimuli     = stimulus_col,
    responses   = response_col,
    baseimage   = baseimage,
    rdata       = rdata_path,
    save_as_png = save_as_png,
    targetpath  = targetpath,
    scaling     = scaling
  )

  # stack $ci (the base-less noise component) per producer
  participants <- unique(as.character(responses[[participant_col]]))
  n_pix <- prod(img_dims)
  signal_matrix <- matrix(
    NA_real_,
    nrow = n_pix,
    ncol = length(participants),
    dimnames = list(NULL, participants)
  )
  rendered <- if (isTRUE(keep_rendered)) {
    matrix(NA_real_, nrow = n_pix, ncol = length(participants),
           dimnames = list(NULL, participants))
  } else NULL

  for (pid in participants) {
    key <- grep(paste0(participant_col, "_", pid, "$"), names(cis),
                value = TRUE)
    if (length(key) == 0L) {
      cli::cli_abort(c(
        "rcicr produced no CI for producer {.val {pid}}.",
        "i" = "batchGenerateCI2IFC keys: {.val {names(cis)}}"
      ))
    }
    ci_mat <- cis[[key[1L]]]$ci
    if (is.null(ci_mat)) {
      cli::cli_abort("rcicr CI element missing {.var $ci} component.")
    }
    signal_matrix[, pid] <- as.vector(ci_mat)
    if (!is.null(rendered)) {
      combined <- cis[[key[1L]]]$combined
      if (is.null(combined)) {
        cli::cli_abort(
          "rcicr CI element missing {.var $combined} component, \\
           cannot fill {.var $rendered_ci}."
        )
      }
      rendered[, pid] <- as.vector(combined)
    }
  }
  attr(signal_matrix, "img_dims") <- img_dims
  attr(signal_matrix, "source") <- "raw"
  if (!is.null(rendered)) {
    attr(rendered, "img_dims") <- img_dims
    attr(rendered, "source") <- "rendered"
  }

  out <- list(
    signal_matrix = signal_matrix,
    participants  = participants,
    img_dims      = img_dims,
    scaling       = scaling,
    rcicr_result  = cis
  )
  if (!is.null(rendered)) out$rendered_ci <- rendered
  out
}
