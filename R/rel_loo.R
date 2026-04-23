#' Leave-one-out sensitivity
#'
#' For each participant `i`, correlates the full-sample group CI
#' (`rowMeans(signal_matrix)`) with the group CI computed after
#' removing participant `i` (`rowMeans(signal_matrix[, -i])`). A
#' participant whose LOO correlation falls more than
#' `flag_threshold_sd` standard deviations below the mean is flagged
#' as an influential outlier  --  their removal changes the group CI
#' noticeably more than the others.
#'
#' The default `flag_threshold_sd = 2.5` is calibrated so that a
#' 30-participant dataset flags roughly 0.3 participants by chance
#' rather than the ~1.5 that a 2-SD rule would produce.
#'
#' @param signal_matrix Pixels x participants, base-subtracted.
#' @param flag_threshold_sd Numeric threshold in standard deviations
#'   below the mean LOO correlation.
#' @return Object of class `rcicrely_loo`. Fields: `$correlations`
#'   (named vector), `$mean_r`, `$sd_r`, `$threshold`, `$flagged`
#'   (character vector of participant ids), `$summary_df` (a
#'   data.frame with per-participant correlation and a `flag` logical).
#' @seealso [rel_split_half()], [rel_icc()], [run_within()]
#' @export
rel_loo <- function(signal_matrix,
                    flag_threshold_sd = 2.5) {
  validate_signal_matrix(signal_matrix)
  n <- ncol(signal_matrix)
  if (is.null(colnames(signal_matrix))) {
    colnames(signal_matrix) <- sprintf("p%03d", seq_len(n))
  }
  full <- rowMeans(signal_matrix)
  cors <- numeric(n)
  names(cors) <- colnames(signal_matrix)
  for (i in seq_len(n)) {
    held_out <- rowMeans(signal_matrix[, -i, drop = FALSE])
    cors[i]  <- stats::cor(full, held_out)
  }

  mean_r <- mean(cors)
  sd_r   <- stats::sd(cors)
  threshold <- mean_r - flag_threshold_sd * sd_r
  flagged <- names(cors)[cors < threshold]

  summary_df <- data.frame(
    participant_id = names(cors),
    correlation    = unname(cors),
    flag           = cors < threshold,
    stringsAsFactors = FALSE
  )

  new_rcicrely_loo(
    correlations = cors,
    mean_r       = mean_r,
    sd_r         = sd_r,
    threshold    = threshold,
    flagged      = flagged,
    summary_df   = summary_df
  )
}
