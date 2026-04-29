## §1.1: empirical-null tests for rel_split_half()
## Covers `null = "permutation"` and `null = "random_responders"`,
## including the noise_matrix gate, the result-shape changes, and the
## standalone `rel_split_half_null()`.

test_that("null = 'none' (default) keeps backward compat: null fields are NA / NULL", {
  sig <- raw_pure_signal(n_p = 12L)
  res <- suppressWarnings(
    rel_split_half(sig, n_permutations = 50L, seed = 1L,
                   progress = FALSE)
  )
  expect_identical(res$null, "none")
  expect_null(res$null_distribution)
  expect_true(is.na(res$r_hh_null_p95))
  expect_true(is.na(res$r_hh_excess))
  expect_true(is.na(res$r_sb_excess))
})

test_that("null = 'permutation' returns a chance-baseline distribution", {
  sig <- raw_pure_signal(n_p = 12L)
  res <- suppressWarnings(
    rel_split_half(sig, n_permutations = 100L, null = "permutation",
                   seed = 1L, progress = FALSE)
  )
  expect_identical(res$null, "permutation")
  expect_length(res$null_distribution, 100L)
  expect_true(all(is.finite(res$null_distribution)))
  # Pure-signal r_hh ~ 1; permutation null r_hh ~ 0; so excess large.
  expect_gt(res$r_hh_excess, 0.5)
  expect_lt(abs(stats::median(res$null_distribution)), 0.2)
})

test_that("null = 'random_responders' aborts without noise_matrix", {
  sig <- raw_pure_signal(n_p = 12L)
  expect_error(
    suppressWarnings(
      rel_split_half(sig, n_permutations = 20L,
                     null = "random_responders",
                     progress = FALSE)
    ),
    "noise_matrix"
  )
})

test_that("null = 'random_responders' returns a finite distribution", {
  sig    <- raw_pure_signal(n_pix = 256L, n_p = 10L)
  noise  <- matrix(stats::rnorm(256L * 30L), 256L, 30L)
  res <- suppressWarnings(
    rel_split_half(sig, n_permutations = 30L,
                   null = "random_responders",
                   noise_matrix = noise,
                   seed = 1L, progress = FALSE)
  )
  expect_identical(res$null, "random_responders")
  expect_length(res$null_distribution, 30L)
  expect_true(all(is.finite(res$null_distribution)))
  expect_true(is.finite(res$r_hh_null_p95))
  expect_true(is.finite(res$r_hh_excess))
})

test_that("rel_split_half_null computes a standalone null", {
  noise <- matrix(stats::rnorm(256L * 30L), 256L, 30L)
  null_d <- suppressWarnings(
    rel_split_half_null(
      n_participants = 12L,
      n_pixels       = 256L,
      null           = "random_responders",
      noise_matrix   = noise,
      n_permutations = 50L,
      seed           = 1L
    )
  )
  expect_length(null_d, 50L)
  expect_true(all(is.finite(null_d)))
})

test_that("rel_split_half_null aborts on random_responders without noise_matrix", {
  expect_error(
    rel_split_half_null(
      n_participants = 12L, n_pixels = 256L,
      null = "random_responders", n_permutations = 10L
    ),
    "noise_matrix"
  )
})

test_that("noise_matrix shape mismatch errors", {
  sig   <- raw_pure_signal(n_pix = 256L, n_p = 10L)
  noise <- matrix(stats::rnorm(64L * 30L), 64L, 30L)
  expect_error(
    suppressWarnings(
      rel_split_half(sig, n_permutations = 10L,
                     null = "random_responders",
                     noise_matrix = noise, progress = FALSE)
    ),
    "Row counts"
  )
})

test_that("run_within forwards null + noise_matrix to rel_split_half", {
  sig   <- raw_pure_signal(n_pix = 256L, n_p = 10L)
  noise <- matrix(stats::rnorm(256L * 30L), 256L, 30L)
  rep <- suppressWarnings(
    run_within(sig, n_permutations = 30L,
               null = "random_responders",
               noise_matrix = noise,
               seed = 1L, progress = FALSE)
  )
  expect_s3_class(rep, "rcicrely_report")
  expect_identical(rep$results$split_half$null, "random_responders")
  expect_length(rep$results$split_half$null_distribution, 30L)
})
