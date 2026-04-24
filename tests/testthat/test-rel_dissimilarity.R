test_that("rel_dissimilarity returns well-formed result", {
  pair <- make_sig_pair(256L, 10L, seed = 1L)
  d <- suppressWarnings(
    rel_dissimilarity(pair$a, pair$b, n_boot = 200L,
                      seed = 1L, progress = FALSE)
  )
  expect_s3_class(d, "rcicrely_dissim")
  expect_true(is.finite(d$correlation))
  expect_true(is.finite(d$euclidean))
  expect_length(d$boot_cor, 200L)
  expect_length(d$boot_dist, 200L)
  expect_true(d$ci_cor[1] < d$correlation + 0.01)
  expect_true(d$ci_cor[2] > d$correlation - 0.01)
})

test_that("A-vs-A more similar than A-vs-B", {
  pair <- make_sig_pair(256L, 10L, seed = 1L)
  d_ab <- suppressWarnings(rel_dissimilarity(pair$a, pair$b,
                                             n_boot = 200L, seed = 1L,
                                             progress = FALSE))
  # A vs A across split halves
  half1 <- pair$a[, 1:5]; attr(half1, "img_dims") <- pair$img_dims
  half2 <- pair$a[, 6:10]; attr(half2, "img_dims") <- pair$img_dims
  d_aa <- suppressWarnings(rel_dissimilarity(half1, half2,
                                             n_boot = 200L, seed = 2L,
                                             progress = FALSE))
  expect_gt(d_aa$correlation, d_ab$correlation)
  expect_lt(d_aa$euclidean,   d_ab$euclidean)
})

test_that("rel_dissimilarity is reproducible under seed", {
  pair <- make_sig_pair(64L, 8L)
  a <- suppressWarnings(rel_dissimilarity(pair$a, pair$b, n_boot = 100L,
                                          seed = 3L, progress = FALSE))
  b <- suppressWarnings(rel_dissimilarity(pair$a, pair$b, n_boot = 100L,
                                          seed = 3L, progress = FALSE))
  expect_identical(a$boot_cor,  b$boot_cor)
  expect_identical(a$boot_dist, b$boot_dist)
})

test_that("euclidean_normalised equals euclidean / sqrt(n_pixels)", {
  pair <- make_sig_pair(256L, 10L, seed = 4L)
  d <- suppressWarnings(
    rel_dissimilarity(pair$a, pair$b, n_boot = 100L, seed = 4L,
                      progress = FALSE)
  )
  expect_equal(d$euclidean_normalised, d$euclidean / sqrt(256))
  expect_equal(d$n_pixels, 256L)
})

test_that("rel_dissimilarity recovers at realistic scale (64x64, N=30)", {
  skip_on_cran()
  pair <- make_realistic_pair(n_side = 64L, n_p = 30L, snr = 0.3, seed = 21L)
  d <- suppressWarnings(
    rel_dissimilarity(pair$a, pair$b, n_boot = 300L,
                      seed = 21L, progress = FALSE)
  )
  expect_s3_class(d, "rcicrely_dissim")
  # Disjoint spatial signals at medium SNR should produce a clearly
  # non-zero Euclidean distance with a CI that excludes 0.
  expect_gt(d$euclidean, 0)
  expect_gt(d$ci_dist[1], 0)
  expect_equal(d$n_pixels, 64L * 64L)
  expect_equal(d$euclidean_normalised, d$euclidean / sqrt(64 * 64))
})

test_that("rel_dissimilarity CI matches quantile() on the boot distribution", {
  # Consistency check: the reported percentile CI must equal
  # quantile(boot_dist, c(tail, 1-tail)) by construction. This
  # validates the CI-extraction step without a second engine.
  skip_on_cran()
  pair <- make_realistic_pair(n_side = 32L, n_p = 30L, snr = 0.3, seed = 31L)
  d <- suppressWarnings(
    rel_dissimilarity(pair$a, pair$b, n_boot = 500L,
                      seed = 31L, progress = FALSE)
  )
  manual_ci <- unname(stats::quantile(d$boot_dist, c(0.025, 0.975)))
  expect_equal(d$ci_dist, manual_ci)
  manual_ci_cor <- unname(stats::quantile(d$boot_cor, c(0.025, 0.975)))
  expect_equal(d$ci_cor, manual_ci_cor)
})

test_that("rel_dissimilarity agrees with an independent bootstrap reference", {
  # External-ish cross-check: implement a second, independent
  # bootstrap loop from scratch with the same seed and verify the
  # bootstrap distribution matches. This validates both the resample
  # scheme and the statistic.
  skip_on_cran()
  pair <- make_realistic_pair(n_side = 32L, n_p = 20L, snr = 0.3, seed = 41L)
  d <- suppressWarnings(
    rel_dissimilarity(pair$a, pair$b, n_boot = 200L,
                      seed = 41L, progress = FALSE)
  )
  ref_dist <- local({
    set.seed(41L)
    n_a <- ncol(pair$a); n_b <- ncol(pair$b)
    out <- numeric(200L)
    for (i in seq_len(200L)) {
      ia <- sample.int(n_a, n_a, replace = TRUE)
      ib <- sample.int(n_b, n_b, replace = TRUE)
      ma <- rowMeans(pair$a[, ia, drop = FALSE])
      mb <- rowMeans(pair$b[, ib, drop = FALSE])
      out[i] <- sqrt(sum((ma - mb)^2))
    }
    out
  })
  expect_equal(d$boot_dist, ref_dist, tolerance = 1e-10)
})
