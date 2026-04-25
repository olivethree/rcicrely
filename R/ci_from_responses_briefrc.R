#' Compute individual Brief-RC CIs from trial-level responses
#'
#' @description
#' Native implementation of Schmitz, Rougier & Yzerbyt (2024)'s Brief-RC
#' mask. **Does not** call any `rcicr::*_brief` function, those do
#' not exist in the upstream rcicr v1.0.1. Only rcicr's noise-pattern
#' pool (the `stimuli` object inside an `.Rdata` from
#' `generateStimuli2IFC()`) is reused; the mask is computed in pure R.
#'
#' Use this when you have Brief-RC 12 trial-level responses and want
#' the package to produce per-producer noise masks ready for the
#' reliability metrics.
#'
#' @details
#' Formula (Schmitz's `genMask()` exactly):
#' ```
#' X    <- data.table(response, stim)
#' X    <- X[, .(response = mean(response)), stim]  # collapse duplicates
#' mask <- (noiseMatrix[, X$stim] %*% X$response) / length(X$response)
#' ```
#' The `length(X$response)` denominator is the number of **unique** pool
#' ids chosen by that participant, not the raw trial count. If a
#' participant chooses the same stimulus on two trials with opposite
#' responses, those two cancel.
#'
#' @section Reading the result:
#' * `$signal_matrix` is the **raw mask** per producer, the object
#'   every `rel_*` function expects. Always pass this (and only this)
#'   to reliability metrics or to any external infoVal computation.
#' * `$rendered_ci`, when present, is `base + scaling(mask)` per
#'   producer. **Visualization only.** Save to PNG, plot for
#'   inspection, but do not feed into reliability metrics, the
#'   scaling step distorts variance-based statistics. Also do not
#'   feed it to a hand-rolled Brief-RC `infoVal` (no upstream
#'   implementation exists in either rcicr or rcicrely): the
#'   reference distribution would compare against a different
#'   scale.
#' * `$participants` and `$img_dims` are convenience metadata.
#'
#' @section Common mistakes:
#' * Passing the "expanded" 12-rows-per-trial response format. Brief-RC
#'   12 is one row per trial; `stimulus` is the chosen pool id and
#'   `response` is `+1` (oriented) or `-1` (inverted). Schmitz 2024
#'   sec.3.1.2.
#' * Using `$rendered_ci` for downstream stats. It exists only because
#'   the user often wants to save a PNG.
#' * Asking for `method = "briefrc20"`, that's reserved for a future
#'   release and aborts.
#'
#' @param responses Data frame with one row per trial. Must contain the
#'   columns named by `participant_col`, `stimulus_col`, `response_col`.
#'   `response` values must be in `{-1, +1}`. Brief-RC 12 trial
#'   structure: `stimulus` = chosen pool id; `response = +1` if the
#'   oriented version was chosen, `-1` if the inverted version was
#'   chosen. Unchosen faces are not recorded.
#' @param rdata_path,noise_matrix Exactly one must be supplied. Provide
#'   `rdata_path` to read the noise matrix from an rcicr `.Rdata`
#'   (`stimuli` object), or pass a pre-loaded `noise_matrix` directly.
#' @param base_image_path Path to the base face image (PNG / JPEG).
#'   Used to validate the noise matrix shape and (when `scaling` is
#'   not `"none"`) to render the visualisation-only `$rendered_ci`
#'   field.
#' @param participant_col,stimulus_col,response_col Column names.
#' @param method One of `"briefrc12"`. `"briefrc20"` is reserved for a
#'   future release and currently aborts with a clear message.
#' @param scaling Visualisation-only scaling for the optional
#'   `$rendered_ci` field. One of `"none"` (default; no rendered_ci
#'   returned), `"matched"` (stretch mask to base range, then add to
#'   base, Schmitz's default), or `"constant"` (multiply mask by
#'   `scaling_constant`, then add to base). The mathematical
#'   `$signal_matrix` returned by this function is always the raw
#'   unscaled mask, regardless of this argument.
#' @param scaling_constant Numeric multiplier used when
#'   `scaling = "constant"`. Ignored otherwise.
#' @return A list with
#'   * `signal_matrix`, pixels x participants raw mask (always raw,
#'     even when `scaling != "none"`).
#'   * `rendered_ci`, pixels x participants, present only when
#'     `scaling != "none"`. Visualisation only.
#'   * `participants`, character vector of participant ids.
#'   * `img_dims`, integer `c(nrow, ncol)`.
#'   * `scaling`, the scaling option that was used.
#' @seealso [ci_from_responses_2ifc()], [run_within()], [run_between()]
#' @references
#' Schmitz, M., Rougier, M., & Yzerbyt, V. (2024). Introducing the
#' brief reverse correlation: an improved tool to assess visual
#' representations. *European Journal of Social Psychology*.
#' \doi{10.1002/ejsp.3100}
#' @export
#' @examples
#' \dontrun{
#' res <- ci_from_responses_briefrc(
#'   responses       = my_responses,
#'   rdata_path      = "data/rcicr_stimuli.Rdata",
#'   base_image_path = "data/base.jpg"
#' )
#' # signal_matrix goes into rel_*; rendered_ci is for plotting only.
#' rel_split_half(res$signal_matrix, seed = 1L)
#' }
ci_from_responses_briefrc <- function(responses,
                                      rdata_path       = NULL,
                                      noise_matrix     = NULL,
                                      base_image_path,
                                      participant_col  = "participant_id",
                                      stimulus_col     = "stimulus",
                                      response_col     = "response",
                                      method           = c("briefrc12",
                                                           "briefrc20"),
                                      scaling          = c("none",
                                                           "matched",
                                                           "constant"),
                                      scaling_constant = NULL) {
  method  <- match.arg(method)
  scaling <- match.arg(scaling)
  if (method == "briefrc20") {
    cli::cli_abort(c(
      "Brief-RC 20 is not supported in this release.",
      "i" = "rcicrely v0.1 supports {.val briefrc12} only. \\
             v0.2 is planned to add it."
    ))
  }
  if (scaling == "constant") {
    if (is.null(scaling_constant) ||
          !is.numeric(scaling_constant) ||
          length(scaling_constant) != 1L ||
          !is.finite(scaling_constant)) {
      cli::cli_abort(
        "{.arg scaling_constant} must be a finite numeric scalar when \\
         {.code scaling = \"constant\"}."
      )
    }
  }

  responses <- as.data.frame(responses)
  required <- c(participant_col, stimulus_col, response_col)
  missing_cols <- setdiff(required, colnames(responses))
  if (length(missing_cols) > 0L) {
    n_missing <- length(missing_cols)
    cli::cli_abort(c(
      "Missing {n_missing} column{?s} in {.arg responses}:",
      "*" = "{.val {missing_cols}}"
    ))
  }

  # force numeric: data.table::fread reads {-1, 1} as integer, which
  # then fails identical() checks against c(-1, 1) (double).
  resp_values <- as.numeric(responses[[response_col]])
  if (!all(is.finite(resp_values))) {
    cli::cli_abort(
      "Column {.var {response_col}} contains non-finite values."
    )
  }
  uniq <- sort(unique(resp_values))
  if (!identical(uniq, c(-1, 1))) {
    msg <- c(
      "Column {.var {response_col}} must contain only values in \\
       {.val {c(-1, 1)}}.",
      "*" = "Got: {.val {uniq}}"
    )
    if (identical(uniq, c(0, 1)) || identical(uniq, c(0, 1L))) {
      msg <- c(
        msg,
        "i" = "Did you mean {.code -1 / +1}? The {.code {{0, 1}}} \\
               coding is the most common silent failure in RC \\
               pipelines (often produced by experiment software \\
               that records 'left' / 'right' as 0 / 1).",
        "i" = "Recode in one line: \\
               {.code responses${response_col} <- 2 * responses${response_col} - 1}"
      )
    }
    cli::cli_abort(msg)
  }

  # load noise matrix
  if (is.null(noise_matrix) && is.null(rdata_path)) {
    cli::cli_abort(
      "Pass either {.arg rdata_path} or {.arg noise_matrix}."
    )
  }
  if (!is.null(rdata_path) && !is.null(noise_matrix)) {
    cli::cli_warn(
      "Both {.arg rdata_path} and {.arg noise_matrix} supplied, \\
       using {.arg noise_matrix}."
    )
  }
  if (is.null(noise_matrix)) {
    noise_matrix <- read_noise_matrix(rdata_path)
  } else {
    if (!is.matrix(noise_matrix)) {
      noise_matrix <- as.matrix(noise_matrix)
    }
    storage.mode(noise_matrix) <- "double"
  }

  # read base image for img_dims validation (and rendering, if requested)
  base_img <- read_image_as_gray(base_image_path)
  img_dims <- as.integer(dim(base_img))
  base_vec <- as.vector(base_img)
  if (prod(img_dims) != nrow(noise_matrix)) {
    cli::cli_abort(c(
      "Base image and noise matrix pixel counts disagree.",
      "*" = "Base: {img_dims[1]} x {img_dims[2]} = \\
             {prod(img_dims)}",
      "*" = "Noise matrix: {nrow(noise_matrix)} rows"
    ))
  }
  n_pool <- ncol(noise_matrix)

  participants <- unique(as.character(responses[[participant_col]]))
  signal_matrix <- matrix(
    NA_real_,
    nrow = nrow(noise_matrix),
    ncol = length(participants),
    dimnames = list(NULL, participants)
  )

  response_char <- as.character(responses[[participant_col]])
  stim_all      <- as.integer(responses[[stimulus_col]])
  resp_all      <- resp_values

  if (any(stim_all < 1L | stim_all > n_pool)) {
    bad <- range(stim_all)
    cli::cli_abort(c(
      "Column {.var {stimulus_col}} has ids outside the pool range.",
      "*" = "Range in data: [{bad[1]}, {bad[2]}]",
      "*" = "Noise matrix pool size: {n_pool}"
    ))
  }

  for (pid in participants) {
    idx <- which(response_char == pid)
    x <- data.table::data.table(
      response = resp_all[idx],
      stim     = stim_all[idx]
    )
    x <- x[, list(response = mean(response)), by = "stim"]
    mask <- (noise_matrix[, x$stim, drop = FALSE] %*% x$response) /
      nrow(x)
    signal_matrix[, pid] <- as.numeric(mask)
  }
  attr(signal_matrix, "img_dims") <- img_dims

  out <- list(
    signal_matrix = signal_matrix,
    participants  = participants,
    img_dims      = img_dims,
    scaling       = scaling
  )

  if (scaling != "none") {
    out$rendered_ci <- render_brief_ci(
      signal_matrix    = signal_matrix,
      base_vec         = base_vec,
      scaling          = scaling,
      scaling_constant = scaling_constant
    )
    attr(out$rendered_ci, "img_dims") <- img_dims
  }

  out
}

#' Render a raw mask matrix into base+scaling(mask) for visualisation
#'
#' @keywords internal
#' @noRd
render_brief_ci <- function(signal_matrix, base_vec,
                            scaling, scaling_constant = NULL) {
  base_rng <- diff(range(base_vec, finite = TRUE))
  out <- signal_matrix
  for (j in seq_len(ncol(signal_matrix))) {
    mask <- signal_matrix[, j]
    scaled <- switch(
      scaling,
      matched  = {
        mr <- diff(range(mask, finite = TRUE))
        if (mr == 0) rep(0, length(mask)) else mask * (base_rng / mr)
      },
      constant = mask * scaling_constant,
      mask
    )
    out[, j] <- base_vec + scaled
  }
  out
}
