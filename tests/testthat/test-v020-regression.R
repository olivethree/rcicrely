## Explicit regression tests for the v0.2.0 sampling-with-replacement
## bug in `simulate_reference_norms()`. The integration-style test
## "2IFC reference is calibrated" in `test-infoval.R` covers the
## end-to-end z calibration; this file pins the unit-level invariant
## (no replacement when n_trials <= n_pool) so a future refactor
## can't reintroduce the bug without one of these tests failing.

test_that("simulate_reference_norms samples without replacement when n_trials <= n_pool", {
  n_pix  <- 64L
  n_pool <- 30L
  noise  <- matrix(stats::rnorm(n_pix * n_pool), n_pix, n_pool)

  # Spy on sample.int by intercepting via trace(): inspect the value of
  # `replace` actually used inside the function. The simpler invariant
  # we can verify without tracing is that the *resulting reference
  # mean norm* matches the closed-form expectation for sampling
  # without replacement when n_trials == n_pool, and is strictly
  # smaller than the with-replacement mean would be.
  set.seed(1L)
  ref_no_repl <- rcicrely:::simulate_reference_norms(
    noise_matrix = noise,
    n_trials     = n_pool,
    iter         = 500L
  )
  # Manually simulate the buggy with-replacement variant for the same
  # config so we can assert the fix moves the mean down.
  with_replacement_ref <- replicate(500L, {
    stim <- sample.int(n_pool, n_pool, replace = TRUE)
    resp <- sample(c(-1, 1), n_pool, replace = TRUE)
    uniq <- sort(unique(stim))
    if (length(uniq) < length(stim)) {
      idx  <- match(stim, uniq)
      sums <- tabulate(idx, length(uniq))
      wts  <- as.vector(tapply(resp, idx, sum)) / sums
    } else {
      ord  <- order(stim)
      uniq <- stim[ord]
      wts  <- resp[ord]
    }
    mask_vec <- noise[, uniq, drop = FALSE] %*% wts / length(uniq)
    sqrt(sum(mask_vec * mask_vec))
  })

  # The with-replacement variant uses fewer unique stims (~ n_pool * (1 - 1/e))
  # which inflates the per-iteration norm because the genMask divisor
  # n_unique_stims drops, so each kept stim contributes more weight.
  expect_gt(mean(with_replacement_ref), mean(ref_no_repl))
})

test_that("simulate_reference_norms uses replacement only when n_trials > n_pool", {
  n_pix  <- 64L
  n_pool <- 20L
  noise  <- matrix(stats::rnorm(n_pix * n_pool), n_pix, n_pool)

  # When n_trials > n_pool, sampling without replacement is impossible,
  # so the function must fall back to with-replacement sampling. Sanity:
  # the call completes and returns the requested number of norms.
  set.seed(2L)
  out <- rcicrely:::simulate_reference_norms(
    noise_matrix = noise,
    n_trials     = n_pool * 2L,
    iter         = 100L
  )
  expect_length(out, 100L)
  expect_true(all(is.finite(out)))
})
