## §1.5: tests for agreement_map_test()

test_that("agreement_map_test returns the expected shape", {
  sig <- raw_pure_signal(n_pix = 256L, n_p = 16L)
  res <- suppressWarnings(
    agreement_map_test(sig, n_permutations = 200L, seed = 1L,
                       progress = FALSE)
  )
  expect_s3_class(res, "rcicrely_agreement_map_test")
  expect_length(res$observed_t, 256L)
  expect_length(res$pmap, 256L)
  expect_length(res$significant_mask, 256L)
  expect_length(res$null_distribution, 200L)
})

test_that("agreement_map_test detects planted signal", {
  sig <- raw_pure_signal(n_pix = 256L, n_p = 20L)
  res <- suppressWarnings(
    agreement_map_test(sig, n_permutations = 500L, seed = 1L,
                       progress = FALSE)
  )
  # Half the pixels carry constant >0 signal across producers, so a
  # large fraction must be significant
  signal_region <- seq_len(128L)
  expect_gt(mean(res$significant_mask[signal_region]), 0.7)
  # Off-signal half should have ~0 significant pixels (FWER controlled)
  off_region <- (129L):256L
  expect_lt(mean(res$significant_mask[off_region]), 0.05)
})

test_that("agreement_map_test under H0 controls FWER at alpha", {
  # Generate many independent zero-signal datasets, ask: does the
  # max-|t| null produce significant pixels at the nominal alpha rate?
  set.seed(99L)
  n_runs <- 30L
  any_sig <- logical(n_runs)
  for (i in seq_len(n_runs)) {
    sig <- raw_zero_signal(n_pix = 64L, n_p = 12L,
                           seed = 100L + i)
    res <- suppressWarnings(
      agreement_map_test(sig, n_permutations = 200L,
                         seed = i, progress = FALSE)
    )
    any_sig[i] <- any(res$significant_mask)
  }
  # Empirical "any significant" rate should be near alpha = 0.05;
  # with 30 runs we have very wide confidence band, just check that
  # it is not catastrophically inflated.
  expect_lt(mean(any_sig), 0.20)
})

test_that("mask restricts the test to a region", {
  sig <- raw_pure_signal(n_pix = 256L, n_p = 16L)
  mask <- c(rep(TRUE, 128L), rep(FALSE, 128L))
  res <- suppressWarnings(
    agreement_map_test(sig, n_permutations = 200L,
                       mask = mask, seed = 1L, progress = FALSE)
  )
  # Outside-mask pixels return NA in observed/pmap and FALSE in mask
  expect_true(all(is.na(res$observed_t[!mask])))
  expect_true(all(is.na(res$pmap[!mask])))
  expect_true(all(!res$significant_mask[!mask]))
})

test_that("agreement_map_test errors on rendered without acknowledge_scaling", {
  sig <- raw_pure_signal()
  attr(sig, "source") <- "rendered"
  expect_error(
    agreement_map_test(sig, n_permutations = 50L, progress = FALSE),
    "rendered"
  )
})
