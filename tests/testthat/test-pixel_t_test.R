test_that("pixel_t_test matches per-pixel Welch t-test via stats::t.test", {
  set.seed(5)
  n_pix <- 20L
  n_a <- 8L
  n_b <- 10L
  a <- matrix(rnorm(n_pix * n_a), n_pix, n_a)
  b <- matrix(rnorm(n_pix * n_b, mean = 0.3), n_pix, n_b)

  vec <- suppressWarnings(pixel_t_test(a, b))
  expect_equal(length(vec), n_pix)

  # reference via stats::t.test row by row
  ref <- vapply(
    seq_len(n_pix),
    function(i) stats::t.test(a[i, ], b[i, ], var.equal = FALSE)$statistic,
    numeric(1L)
  )
  expect_equal(vec, unname(ref), tolerance = 1e-10)
})

test_that("pixel_t_test handles zero-variance pixels without NaN spam", {
  # both matrices constant -> zero variance per row
  a <- matrix(1, 10L, 5L)
  b <- matrix(1, 10L, 5L)
  vec <- suppressWarnings(pixel_t_test(a, b))
  expect_equal(vec, rep(0, 10L))
})

test_that("pixel_t_test errors on pixel mismatch", {
  a <- matrix(0, 10L, 5L)
  b <- matrix(0, 20L, 5L)
  expect_error(suppressWarnings(pixel_t_test(a, b)), "differ")
})
