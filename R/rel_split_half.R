#' Permuted split-half reliability with Spearman-Brown correction
#'
#' Measures how stable a condition's group-level CI is across random
#' halves of its producer set. For each of `n_permutations`
#' iterations: randomly partition the producers (columns of
#' `signal_matrix`) into two roughly-equal halves, average each half
#' to get two group vectors, and compute the Pearson correlation.
#'
#' For **odd N**, one randomly-chosen producer is dropped per
#' permutation (re-drawn each iteration, not fixed) so both halves
#' have `floor(N/2)` producers.
#'
#' Reports both the mean per-permutation `r` and the Spearman-Brown
#' corrected reliability `r_SB = (2 r) / (1 + r)`. 95% CI is the
#' 2.5 / 97.5 percentile of the permutation distribution.
#'
#' @param signal_matrix Pixels x participants, base-subtracted.
#' @param n_permutations Integer number of random splits. Default
#'   2000  --  cheap and keeps Monte Carlo error on tail probabilities
#'   below 0.01.
#' @param seed Optional integer; if set, results are reproducible.
#'   The caller's global RNG state is restored on exit.
#' @param progress Show a `cli` progress bar.
#' @return An object of class `rcicrely_split_half`. Fields:
#'   `$r_hh`, `$r_sb`, `$ci_95`, `$ci_95_sb`, `$distribution`
#'   (length `n_permutations` vector of per-split `r`),
#'   `$n_participants`, `$n_permutations`.
#' @seealso [rel_loo()], [rel_icc()], [run_within()]
#' @references
#' Brinkman, L., Goffin, S., van de Schoot, R., van Haren, N. E. M.,
#' Dotsch, R., & Aarts, H. (2017). Quantifying the informational
#' value of classification images. *Behavior Research Methods*.
#'
#' Shrout, P. E., & Fleiss, J. L. (1979). Intraclass correlations:
#' uses in assessing rater reliability. *Psychological Bulletin*.
#' @export
#' @examples
#' \dontrun{
#' r <- rel_split_half(signal_matrix, n_permutations = 2000, seed = 1)
#' print(r); plot(r)
#' }
rel_split_half <- function(signal_matrix,
                           n_permutations = 2000L,
                           seed           = NULL,
                           progress       = TRUE) {
  validate_signal_matrix(signal_matrix)
  n_participants <- ncol(signal_matrix)
  half_size <- n_participants %/% 2L
  odd_n <- (n_participants %% 2L) == 1L
  n_permutations <- as.integer(n_permutations)
  if (n_permutations < 10L) {
    cli::cli_warn(
      "{.arg n_permutations} = {n_permutations} is very low; \\
       CI will be unstable."
    )
  }

  pid <- progress_start(n_permutations, "split-half", show = progress)
  on.exit(progress_done(pid), add = TRUE)

  dist <- with_seed(seed, {
    out <- numeric(n_permutations)
    for (i in seq_len(n_permutations)) {
      idx <- sample.int(n_participants)
      if (odd_n) {
        idx <- idx[-1L]  # drop one (rotates each perm)
      }
      h1 <- idx[seq_len(half_size)]
      h2 <- idx[(half_size + 1L):length(idx)]
      v1 <- rowMeans(signal_matrix[, h1, drop = FALSE])
      v2 <- rowMeans(signal_matrix[, h2, drop = FALSE])
      out[i] <- stats::cor(v1, v2)
      progress_tick(pid)
    }
    out
  })

  r_hh <- mean(dist, na.rm = TRUE)
  r_sb <- (2 * r_hh) / (1 + r_hh)
  ci_95     <- stats::quantile(dist, c(0.025, 0.975), na.rm = TRUE,
                               names = FALSE)
  dist_sb   <- (2 * dist) / (1 + dist)
  ci_95_sb  <- stats::quantile(dist_sb, c(0.025, 0.975), na.rm = TRUE,
                               names = FALSE)

  new_rcicrely_split_half(
    r_hh           = r_hh,
    r_sb           = r_sb,
    ci_95          = ci_95,
    ci_95_sb       = ci_95_sb,
    distribution   = dist,
    n_participants = n_participants,
    n_permutations = n_permutations
  )
}
