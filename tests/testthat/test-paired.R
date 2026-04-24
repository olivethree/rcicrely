make_paired_fixture <- function(n_pix = 256L, n_p = 12L,
                                signal_amp = 0.3, pair_cor = 0.7,
                                seed = 1L) {
  # Same producers in A and B; B = A + producer-specific shift + noise.
  set.seed(seed)
  base <- matrix(stats::rnorm(n_pix * n_p, sd = 0.05), n_pix, n_p)
  # pair_cor controls how much of B is carried from A (within-producer).
  a <- base + matrix(stats::rnorm(n_pix * n_p, sd = 0.05), n_pix, n_p)
  b <- pair_cor * base + matrix(stats::rnorm(n_pix * n_p, sd = 0.05),
                                n_pix, n_p)
  # add a disjoint-spatial signal in A only (so A-B has real signal)
  mask <- c(rep(1, n_pix / 2L), rep(0, n_pix / 2L))
  a <- a + signal_amp * outer(mask, stats::runif(n_p, 0.8, 1.2))

  side <- as.integer(sqrt(n_pix))
  colnames(a) <- sprintf("p%02d", seq_len(n_p))
  colnames(b) <- colnames(a)
  attr(a, "img_dims") <- c(side, side)
  attr(b, "img_dims") <- c(side, side)
  list(a = a, b = b, img_dims = c(side, side))
}

test_that("pixel_t_test(paired = TRUE) matches stats::t.test row-by-row", {
  pair <- make_paired_fixture(n_pix = 64L, n_p = 12L, seed = 1L)
  t_pkg <- pixel_t_test(pair$a, pair$b, paired = TRUE)
  t_ref <- vapply(
    seq_len(nrow(pair$a)),
    function(i) {
      tr <- stats::t.test(pair$a[i, ], pair$b[i, ], paired = TRUE)
      unname(tr$statistic)
    },
    numeric(1L)
  )
  expect_equal(t_pkg, t_ref, tolerance = 1e-10)
})

test_that("pixel_t_test(paired = TRUE) aborts on mismatched column names", {
  pair <- make_paired_fixture(n_pix = 64L, n_p = 12L)
  b_wrong <- pair$b
  colnames(b_wrong) <- rev(colnames(b_wrong))
  expect_error(
    pixel_t_test(pair$a, b_wrong, paired = TRUE),
    regexp = "[Cc]olumn names"
  )
})

test_that("pixel_t_test(paired = TRUE) aborts on mismatched column counts", {
  pair <- make_paired_fixture(n_pix = 64L, n_p = 12L)
  b_short <- pair$b[, 1:8, drop = FALSE]
  expect_error(
    pixel_t_test(pair$a, b_short, paired = TRUE),
    "different"
  )
})

test_that("rel_cluster_test(paired = TRUE) uses sign-flip null", {
  pair <- make_paired_fixture(n_pix = 256L, n_p = 15L, seed = 2L)
  r <- suppressWarnings(
    rel_cluster_test(
      pair$a, pair$b,
      img_dims       = pair$img_dims,
      paired         = TRUE,
      n_permutations = 100L,
      seed           = 2L, progress = FALSE
    )
  )
  expect_s3_class(r, "rcicrely_cluster_test")
  expect_true(r$paired)
  expect_equal(r$method, "threshold")
  # injected signal in A only -> positive cluster expected somewhere
  expect_true(any(r$clusters$significant))
})

test_that("rel_cluster_test(paired = TRUE) + TFCE path runs", {
  pair <- make_paired_fixture(n_pix = 64L, n_p = 12L, seed = 3L)
  r <- suppressWarnings(
    rel_cluster_test(
      pair$a, pair$b,
      img_dims       = pair$img_dims,
      method         = "tfce",
      paired         = TRUE,
      n_permutations = 50L,
      tfce_n_steps   = 30L,
      seed           = 3L, progress = FALSE
    )
  )
  expect_equal(r$method, "tfce")
  expect_true(r$paired)
  expect_equal(length(r$tfce_pmap), 64L)
})

test_that("rel_dissimilarity(paired = TRUE) uses shared resample index", {
  pair <- make_paired_fixture(n_pix = 256L, n_p = 15L, seed = 4L)
  d <- suppressWarnings(
    rel_dissimilarity(pair$a, pair$b,
                      paired = TRUE, n_boot = 200L,
                      seed = 4L, progress = FALSE)
  )
  expect_s3_class(d, "rcicrely_dissim")
  expect_equal(length(d$boot_dist), 200L)

  # Hand-rolled reference: same seed, shared resample index
  ref <- local({
    set.seed(4L)
    n <- ncol(pair$a)
    out <- numeric(200L)
    for (i in seq_len(200L)) {
      idx <- sample.int(n, n, replace = TRUE)
      m_a <- rowMeans(pair$a[, idx, drop = FALSE])
      m_b <- rowMeans(pair$b[, idx, drop = FALSE])
      out[i] <- sqrt(sum((m_a - m_b)^2))
    }
    out
  })
  expect_equal(d$boot_dist, ref, tolerance = 1e-10)
})

test_that("paired dissimilarity CI narrower than unpaired when matched", {
  # When producers are matched and there's correlated producer-level
  # variance, paired bootstrap should produce a tighter CI on
  # Euclidean distance than unpaired.
  pair <- make_paired_fixture(n_pix = 1024L, n_p = 30L,
                              signal_amp = 0.3, pair_cor = 0.9,
                              seed = 5L)
  d_p <- suppressWarnings(
    rel_dissimilarity(pair$a, pair$b, paired = TRUE,
                      n_boot = 400L, seed = 5L, progress = FALSE)
  )
  d_u <- suppressWarnings(
    rel_dissimilarity(pair$a, pair$b, paired = FALSE,
                      n_boot = 400L, seed = 5L, progress = FALSE)
  )
  paired_width   <- diff(d_p$ci_dist)
  unpaired_width <- diff(d_u$ci_dist)
  expect_lt(paired_width, unpaired_width)
})
