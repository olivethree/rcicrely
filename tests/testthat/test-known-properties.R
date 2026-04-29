## Known-property tests for every reliability metric.
##
## For each metric, hand-build inputs whose expected behaviour is
## algebraically determined, then assert the metric matches. These
## are the "if a future refactor breaks the math, this catches it"
## tests. Inputs:
##
##   * `zero` — pure noise, no signal. Reliability and discrim
##     metrics should be near their no-signal floor.
##   * `pure` — signal only, no noise. Reliability should approach 1
##     and dissimilarity between two pure-signal conditions should
##     be large.
##   * `sign_flip` — pure signal multiplied by random ±1 per producer.
##     Group-mean cancels to zero; LOO correlations should be highly
##     variable.
##   * `permuted` — random per-producer permutation of pixels. ICC
##     and split-half should drop to floor.

# Fixtures `raw_pure_signal()`, `raw_zero_signal()`,
# `raw_sign_flipped_signal()` live in `helper-fixtures.R` so other
# test files (e.g. `test-rel_split_half-null.R`) can reuse them.

# ----- ICC -----------------------------------------------------------

test_that("rel_icc returns ~1 on pure-signal data", {
  sig <- raw_pure_signal()
  res <- rel_icc(sig)
  expect_gt(res$icc_3_1, 0.95)
  expect_gt(res$icc_3_k, 0.99)
})

test_that("rel_icc is near zero on zero-signal data", {
  sig <- raw_zero_signal()
  res <- rel_icc(sig)
  # No producer-level structure beyond noise, ICC(3,1) should be small
  # in absolute value (allow tolerance for finite-sample variability).
  expect_lt(abs(res$icc_3_1), 0.10)
})

# ----- Split-half ----------------------------------------------------

test_that("rel_split_half returns ~1 on pure-signal data", {
  sig <- raw_pure_signal()
  res <- suppressWarnings(
    rel_split_half(sig, n_permutations = 200L, seed = 1L,
                   progress = FALSE)
  )
  expect_gt(res$r_hh, 0.95)
  expect_gt(res$r_sb, 0.95)
})

test_that("rel_split_half drops on noisy sign-flipped data", {
  # Sign-flip on top of per-producer noise: in expectation the
  # signal-bearing half cancels across producers, so split-half r
  # collapses to the noise floor. Pure-signal-with-noise stays high.
  set.seed(7L)
  n_pix <- 1024L
  n_p   <- 24L
  signal_vec <- c(rep(1, n_pix / 2L), rep(0, n_pix / 2L))

  noise_pure <- matrix(stats::rnorm(n_pix * n_p, sd = 0.05),
                       n_pix, n_p)
  sig_pure   <- noise_pure +
    outer(signal_vec, stats::runif(n_p, 0.15, 0.25))
  attr(sig_pure, "img_dims") <- c(32L, 32L)
  attr(sig_pure, "source")   <- "raw"

  noise_flip <- matrix(stats::rnorm(n_pix * n_p, sd = 0.05),
                       n_pix, n_p)
  signs      <- rep(c(-1, 1), each = n_p / 2L)  # balanced
  sig_flip   <- noise_flip +
    outer(signal_vec, signs * stats::runif(n_p, 0.15, 0.25))
  attr(sig_flip, "img_dims") <- c(32L, 32L)
  attr(sig_flip, "source")   <- "raw"

  res_pure <- suppressWarnings(
    rel_split_half(sig_pure, n_permutations = 300L, seed = 1L,
                   progress = FALSE)
  )
  res_flip <- suppressWarnings(
    rel_split_half(sig_flip, n_permutations = 300L, seed = 1L,
                   progress = FALSE)
  )
  expect_gt(res_pure$r_hh, 0.7)
  expect_lt(res_flip$r_hh, 0.3)
})

# ----- Pixel t-test --------------------------------------------------

test_that("pixel_t_test returns ~0 t-values on identical conditions", {
  sig <- raw_pure_signal(n_p = 16L)
  half_a <- sig[, 1:8]
  half_b <- sig[, 9:16]
  attr(half_a, "source") <- "raw"
  attr(half_b, "source") <- "raw"
  t_vec <- pixel_t_test(half_a, half_b)
  # Both halves carry the same underlying signal; per-pixel t should
  # be small in expectation.
  expect_lt(stats::median(abs(t_vec)), 1.0)
})

test_that("pixel_t_test recovers a planted difference", {
  sig_a <- raw_pure_signal(n_p = 12L, seed = 10L)
  sig_b <- raw_pure_signal(n_p = 12L, seed = 11L)
  # Inject a strong signal in B's first half of pixels
  n_pix <- nrow(sig_a)
  shift <- c(rep(2, n_pix / 2L), rep(0, n_pix / 2L))
  sig_b <- sig_b + shift
  attr(sig_a, "source") <- "raw"
  attr(sig_b, "source") <- "raw"
  t_vec <- pixel_t_test(sig_a, sig_b)
  # Pixels in B's signal region should have strong negative t (B > A)
  signal_region <- seq_len(n_pix / 2L)
  expect_lt(stats::median(t_vec[signal_region]), -1)
})

# ----- Dissimilarity -------------------------------------------------

test_that("rel_dissimilarity Euclidean is small for identical conditions", {
  sig <- raw_pure_signal(n_p = 16L)
  half_a <- sig[, 1:8]
  half_b <- sig[, 9:16]
  attr(half_a, "source") <- "raw"
  attr(half_b, "source") <- "raw"
  res <- suppressWarnings(
    rel_dissimilarity(half_a, half_b, n_boot = 100L,
                      seed = 1L, progress = FALSE)
  )
  # Same underlying signal -> small Euclidean distance, normalised
  # close to 0.
  expect_lt(res$euclidean_normalised, 0.1)
})

test_that("rel_dissimilarity Euclidean is large for different conditions", {
  sig_a <- raw_pure_signal(seed = 10L)
  sig_b <- raw_pure_signal(seed = 11L)
  # Flip B so its mean signal is the opposite pattern
  sig_b <- -sig_b
  attr(sig_a, "source") <- "raw"
  attr(sig_b, "source") <- "raw"
  res <- suppressWarnings(
    rel_dissimilarity(sig_a, sig_b, n_boot = 100L,
                      seed = 1L, progress = FALSE)
  )
  expect_gt(res$euclidean_normalised, 0.5)
})

# ----- LOO -----------------------------------------------------------

test_that("rel_loo correlations are near 1 on pure-signal data", {
  sig <- raw_pure_signal()
  res <- suppressWarnings(rel_loo(sig))
  # Every leave-one-out CI is highly correlated with the full CI.
  expect_gt(min(res$correlations), 0.99)
})

test_that("rel_loo flags an injected outlier", {
  sig <- raw_pure_signal()
  outlier <- matrix(stats::runif(nrow(sig), -1, 1), nrow(sig), 1L)
  attr(outlier, "source") <- "raw"
  sig_plus <- cbind(sig, outlier)
  colnames(sig_plus) <- c(sprintf("p%02d", seq_len(ncol(sig))),
                          "outlier")
  attr(sig_plus, "source") <- "raw"
  res <- suppressWarnings(rel_loo(sig_plus))
  expect_true("outlier" %in% res$flagged)
})
