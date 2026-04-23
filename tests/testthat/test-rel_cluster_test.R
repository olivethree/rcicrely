test_that("rel_cluster_test finds signal clusters", {
  pair <- make_sig_pair(n_pix = 256L, n_p = 10L, seed = 1L)
  r <- suppressWarnings(
    rel_cluster_test(
      pair$a, pair$b,
      img_dims = pair$img_dims,
      n_permutations = 100L,
      cluster_threshold = 2.0,
      alpha = 0.05,
      seed = 1L, progress = FALSE
    )
  )
  expect_s3_class(r, "rcicrely_cluster_test")
  expect_true(any(r$clusters$significant))
  expect_true(all(abs(r$clusters$mass) > 0))
  expect_equal(
    sort(unique(r$clusters$direction)),
    c("neg", "pos")
  )
})

test_that("A-vs-A split-halves rarely produce significant clusters", {
  sig <- make_sig(256L, 20L, signal_strength = 0.3, seed = 4L)
  # split sig into two halves of 10 each; under H0 of no difference,
  # significant clusters should be rare
  half1 <- sig[, 1:10]
  half2 <- sig[, 11:20]
  attr(half1, "img_dims") <- c(16L, 16L)
  attr(half2, "img_dims") <- c(16L, 16L)
  r <- suppressWarnings(
    rel_cluster_test(
      half1, half2,
      img_dims = c(16L, 16L),
      n_permutations = 200L, seed = 2L, progress = FALSE
    )
  )
  # loose sanity bound — should be 0 or very small
  expect_lte(sum(r$clusters$significant), 2L)
})

test_that("rel_cluster_test is reproducible under seed", {
  pair <- make_sig_pair(64L, 8L, seed = 11L)
  a <- suppressWarnings(
    rel_cluster_test(pair$a, pair$b, img_dims = pair$img_dims,
                     n_permutations = 50L, seed = 7L, progress = FALSE)
  )
  b <- suppressWarnings(
    rel_cluster_test(pair$a, pair$b, img_dims = pair$img_dims,
                     n_permutations = 50L, seed = 7L, progress = FALSE)
  )
  expect_identical(a$observed_t,        b$observed_t)
  expect_identical(a$null_distribution, b$null_distribution)
})

test_that("img_dims validation catches pixel-count mismatch", {
  pair <- make_sig_pair(64L, 8L)
  expect_error(
    suppressWarnings(
      rel_cluster_test(pair$a, pair$b, img_dims = c(10L, 10L),
                       n_permutations = 10L, progress = FALSE)
    ),
    "inconsistent"
  )
})
