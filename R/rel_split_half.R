#' Permuted split-half reliability with Spearman-Brown correction
#'
#' @description
#' Measures how stable a condition's group-level CI is across random
#' halves of its producer set. For each of `n_permutations`
#' iterations: randomly partition the producers (columns of
#' `signal_matrix`) into two roughly-equal halves, average each half
#' to get two group vectors, and compute the Pearson correlation.
#'
#' Use this to answer "would I get the same group CI if I'd run a
#' different half of my participants?".
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
#'   2000, cheap and keeps Monte Carlo error on tail probabilities
#'   below 0.01.
#' @param seed Optional integer; if set, results are reproducible.
#'   The caller's global RNG state is restored on exit.
#' @param mask Optional logical vector of length
#'   `nrow(signal_matrix)` restricting computation to a region of
#'   the image (e.g., from [face_mask()] or [load_face_mask()]).
#'   When supplied, the matrix is row-subsetted before any
#'   permutation. See `vignette("tutorial")` §6.6 for the
#'   apply-symmetrically rule.
#' @param progress Show a `cli` progress bar.
#' @section Reading the result:
#' * `$r_hh`, mean per-permutation Pearson r between the two halves.
#' * `$r_sb`, Spearman-Brown projected reliability of the full
#'   sample. This is usually the headline number to report.
#' * `$ci_95`, `$ci_95_sb`, percentile 95% CIs on each.
#' * `$distribution`, the full per-permutation r vector for plotting
#'   or custom CIs.
#' * `$n_participants`, `$n_permutations`, metadata.
#'
#' @section Common mistakes:
#' * Reporting `$r_hh` instead of `$r_sb` when the audience cares
#'   about the full-sample CI's reliability (almost always).
#' * Running with a small `n_permutations` (< 200) to save time, then
#'   trusting the CI bounds, CIs widen with low n_perm.
#'
#' @section Reliability metrics expect raw masks:
#' Operates on the raw mask; results may be distorted if
#' `signal_matrix` was extracted from rendered (scaled) PNGs. See
#' `vignette("tutorial", package = "rcicrely")` chapter 3.
#'
#' @return An object of class `rcicrely_split_half`. Fields described
#'   in **Reading the result** above.
#' @seealso [rel_loo()], [rel_icc()], [run_within()]
#' @references
#' Brinkman, L., Goffin, S., van de Schoot, R., van Haren, N. E. M.,
#' Dotsch, R., & Aarts, H. (2019). Quantifying the informational
#' value of classification images. *Behavior Research Methods*,
#' 51(5), 2059-2073. \doi{10.3758/s13428-019-01232-2}
#'
#' Shrout, P. E., & Fleiss, J. L. (1979). Intraclass correlations:
#' uses in assessing rater reliability. *Psychological Bulletin*,
#' 86(2), 420-428. \doi{10.1037/0033-2909.86.2.420}
#' @export
#' @examples
#' \dontrun{
#' r <- rel_split_half(signal_matrix, n_permutations = 2000, seed = 1)
#' print(r); plot(r)
#' }
rel_split_half <- function(signal_matrix,
                           n_permutations = 2000L,
                           mask           = NULL,
                           seed           = NULL,
                           progress       = TRUE) {
  validate_signal_matrix(signal_matrix)
  signal_matrix <- apply_mask_to_signal(signal_matrix, mask)
  if (looks_scaled(signal_matrix)) warn_looks_scaled("signal_matrix")
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
