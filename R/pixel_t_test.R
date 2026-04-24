#' Vectorised pixel-wise t-test (independent or paired)
#'
#' @description
#' At every pixel, tests whether condition A's mean signal differs
#' from condition B's. Two modes:
#'
#' * `paired = FALSE` (default): independent-samples Welch t per
#'   pixel. Correct when producers in A and B are different people
#'   (between-subjects design).
#' * `paired = TRUE`: paired t per pixel on the per-producer
#'   difference `A - B`. Correct when the same producers contributed
#'   to both conditions (within-subjects design). Paired t is
#'   strictly more powerful than Welch t whenever producer-level
#'   variance is present.
#'
#' Returns a numeric vector of t-values, length `n_pixels`. Used as
#' an intermediate by [rel_cluster_test()]; not intended as a
#' standalone inferential test.
#'
#' @param signal_matrix_a,signal_matrix_b Pixels x participants,
#'   base-subtracted. Row counts must match. When `paired = TRUE`
#'   the column counts must also match, and column names must
#'   correspond to the same producer across matrices.
#' @param paired Logical. `FALSE` (default) uses independent
#'   Welch t; `TRUE` uses paired t.
#' @return Numeric vector of length `nrow(signal_matrix_a)`.
#' @section Reading the result:
#' Numeric vector, one t-value per pixel. Pixels with zero variance
#' get `0` rather than `NaN` so cluster utilities don't have to
#' special-case them.
#'
#' @section Reliability metrics expect raw masks:
#' Welch t and paired t are variance-based and sensitive to scaling.
#' See `vignette("tutorial", package = "rcicrely")` chapter 3.
#'
#' @seealso [rel_cluster_test()]
#' @export
pixel_t_test <- function(signal_matrix_a, signal_matrix_b,
                         paired = FALSE) {
  validate_two_signal_matrices(signal_matrix_a, signal_matrix_b)

  if (isTRUE(paired)) {
    validate_paired_matrices(signal_matrix_a, signal_matrix_b)
    diff <- signal_matrix_a - signal_matrix_b
    n    <- ncol(diff)
    mean_d <- rowMeans(diff)
    var_d  <- rowSums((diff - mean_d)^2) / (n - 1L)
    se     <- sqrt(var_d / n)
    t_vec  <- mean_d / se
    t_vec[!is.finite(t_vec)] <- 0
    return(t_vec)
  }

  n_a <- ncol(signal_matrix_a)
  n_b <- ncol(signal_matrix_b)

  mean_a <- rowMeans(signal_matrix_a)
  mean_b <- rowMeans(signal_matrix_b)

  var_a <- rowSums((signal_matrix_a - mean_a)^2) / (n_a - 1L)
  var_b <- rowSums((signal_matrix_b - mean_b)^2) / (n_b - 1L)

  se <- sqrt(var_a / n_a + var_b / n_b)
  t_vec <- (mean_a - mean_b) / se

  t_vec[!is.finite(t_vec)] <- 0
  t_vec
}

#' Validate that two signal matrices are paired (same producers)
#'
#' @keywords internal
#' @noRd
validate_paired_matrices <- function(a, b,
                                     name_a = "signal_matrix_a",
                                     name_b = "signal_matrix_b") {
  if (ncol(a) != ncol(b)) {
    cli::cli_abort(c(
      "{.arg {name_a}} and {.arg {name_b}} have different \\
       column counts.",
      "*" = "{.arg {name_a}}: {ncol(a)} producers",
      "*" = "{.arg {name_b}}: {ncol(b)} producers",
      "i" = "For a paired test, the two matrices must carry the \\
             same producers."
    ))
  }
  if (!is.null(colnames(a)) && !is.null(colnames(b))) {
    if (!identical(colnames(a), colnames(b))) {
      cli::cli_abort(c(
        "Column names of {.arg {name_a}} and {.arg {name_b}} do not \\
         match.",
        "i" = "For a paired test, column j of {.arg {name_a}} must \\
               correspond to the same producer as column j of \\
               {.arg {name_b}}."
      ))
    }
  }
  invisible(TRUE)
}
