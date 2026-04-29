## §1.3: empirical-null tests for rel_dissimilarity()
## Covers `null = "permutation"` in both between-subjects
## (condition-label permutation) and paired (sign-flip) modes.

test_that("null = 'none' (default) keeps backward compat: null fields are NA / NULL", {
  s <- make_sig_pair()
  attr(s$a, "source") <- "raw"
  attr(s$b, "source") <- "raw"
  res <- suppressWarnings(
    rel_dissimilarity(s$a, s$b, n_boot = 50L, seed = 1L,
                      progress = FALSE)
  )
  expect_identical(res$null, "none")
  expect_null(res$null_distribution)
  expect_true(is.na(res$d_null_p95))
  expect_true(is.na(res$d_z))
  expect_true(is.na(res$d_ratio))
  expect_false(isTRUE(res$paired))
})

test_that("null = 'permutation' between-subjects builds a chance baseline", {
  s <- make_sig_pair()
  attr(s$a, "source") <- "raw"
  attr(s$b, "source") <- "raw"
  res <- suppressWarnings(
    rel_dissimilarity(s$a, s$b, n_boot = 50L,
                      null = "permutation", n_permutations = 200L,
                      seed = 1L, progress = FALSE)
  )
  expect_identical(res$null, "permutation")
  expect_length(res$null_distribution, 200L)
  expect_true(all(is.finite(res$null_distribution)))
  # make_sig_pair has different signal patterns in A vs B, so the
  # observed distance should sit above the null distribution and
  # produce a meaningfully positive z.
  expect_gt(res$d_z, 1)
  expect_gt(res$d_ratio, 1)
})

test_that("null = 'permutation' paired uses sign-flip", {
  s <- make_sig_pair()
  colnames(s$a) <- sprintf("p%02d", seq_len(ncol(s$a)))
  colnames(s$b) <- colnames(s$a)
  attr(s$a, "source") <- "raw"
  attr(s$b, "source") <- "raw"
  res <- suppressWarnings(
    rel_dissimilarity(s$a, s$b, paired = TRUE, n_boot = 50L,
                      null = "permutation", n_permutations = 200L,
                      seed = 1L, progress = FALSE)
  )
  expect_identical(res$null, "permutation")
  expect_true(isTRUE(res$paired))
  expect_length(res$null_distribution, 200L)
  expect_true(is.finite(res$d_z))
  expect_gt(res$d_z, 0)
})

test_that("null distribution distinguishes signal from no-signal", {
  # Two zero-signal conditions: observed distance should sit close
  # to the null median and z should be small.
  s_a <- raw_zero_signal(n_p = 20L, seed = 11L)
  s_b <- raw_zero_signal(n_p = 20L, seed = 12L)
  res <- suppressWarnings(
    rel_dissimilarity(s_a, s_b, n_boot = 50L,
                      null = "permutation", n_permutations = 200L,
                      seed = 1L, progress = FALSE)
  )
  expect_lt(abs(res$d_z), 2)   # should hover near zero
})
