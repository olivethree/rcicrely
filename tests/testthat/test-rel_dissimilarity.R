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
