#' Compute individual Brief-RC CIs from trial-level responses
#'
#' Native implementation of Schmitz, Rougier & Yzerbyt (2024)'s Brief-RC
#' mask, per CLAUDE.md sec.5.2. **Does not** call any `rcicr::*_brief`
#' function  --  those do not exist in canonical rcicr v1.0.1. Only rcicr's
#' noise-pattern pool (the `stimuli` object inside an `.Rdata` from
#' `generateStimuli2IFC()`) is reused; the mask is computed in pure R.
#'
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
#'   Used only to stamp `$img_dims` on the result and to validate the
#'   noise matrix shape.
#' @param participant_col,stimulus_col,response_col Column names.
#' @param method One of `"briefrc12"`. `"briefrc20"` is reserved for a
#'   future release and currently aborts with a clear message.
#' @return A list with
#'   * `signal_matrix`  --  pixels x participants numeric matrix (mask per
#'     producer).
#'   * `participants`  --  character vector of participant ids.
#'   * `img_dims`  --  integer `c(nrow, ncol)`.
#' @seealso [ci_from_responses_2ifc()]
#' @references
#' Schmitz, M., Rougier, M., & Yzerbyt, V. (2024). Introducing the
#' brief reverse correlation: an improved tool to assess visual
#' representations. *Behavior Research Methods*.
#' @export
#' @examples
#' \dontrun{
#' res <- ci_from_responses_briefrc(
#'   responses       = my_responses,
#'   rdata_path      = "data/rcicr_stimuli.Rdata",
#'   base_image_path = "data/base.jpg"
#' )
#' }
ci_from_responses_briefrc <- function(responses,
                                      rdata_path      = NULL,
                                      noise_matrix    = NULL,
                                      base_image_path,
                                      participant_col = "participant_id",
                                      stimulus_col    = "stimulus",
                                      response_col    = "response",
                                      method          = c("briefrc12",
                                                          "briefrc20")) {
  method <- match.arg(method)
  if (method == "briefrc20") {
    cli::cli_abort(c(
      "Brief-RC 20 is not supported in this release.",
      "i" = "rcicrely v0.1 supports {.val briefrc12} only \\
             (CLAUDE.md sec.5.1). v0.2 is planned to add it."
    ))
  }

  responses <- as.data.frame(responses)
  required <- c(participant_col, stimulus_col, response_col)
  missing_cols <- setdiff(required, colnames(responses))
  if (length(missing_cols) > 0L) {
    cli::cli_abort(c(
      "Missing column{?s} in {.arg responses}:",
      "*" = "{.val {missing_cols}}"
    ))
  }

  # force numeric to catch the fread int-vs-double trap (CLAUDE.md sec.9.2)
  resp_values <- as.numeric(responses[[response_col]])
  if (!all(is.finite(resp_values))) {
    cli::cli_abort(
      "Column {.var {response_col}} contains non-finite values."
    )
  }
  uniq <- sort(unique(resp_values))
  if (!identical(uniq, c(-1, 1))) {
    cli::cli_abort(c(
      "Column {.var {response_col}} must contain only values in \\
       {.val {c(-1, 1)}}.",
      "*" = "Got: {.val {uniq}}"
    ))
  }

  # load noise matrix
  if (is.null(noise_matrix) && is.null(rdata_path)) {
    cli::cli_abort(
      "Pass either {.arg rdata_path} or {.arg noise_matrix}."
    )
  }
  if (!is.null(rdata_path) && !is.null(noise_matrix)) {
    cli::cli_warn(
      "Both {.arg rdata_path} and {.arg noise_matrix} supplied  --  \\
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

  # read base image for img_dims validation
  base_img <- read_image_as_gray(base_image_path)
  img_dims <- as.integer(dim(base_img))
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

  list(
    signal_matrix = signal_matrix,
    participants  = participants,
    img_dims      = img_dims
  )
}
