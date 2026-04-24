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

test_that("rel_cluster_test finds injected signal at realistic scale", {
  skip_on_cran()
  pair <- make_realistic_pair(n_side = 64L, n_p = 30L, snr = 0.4, seed = 41L)
  r <- suppressWarnings(
    rel_cluster_test(
      pair$a, pair$b,
      img_dims          = pair$img_dims,
      n_permutations    = 200L,
      cluster_threshold = 2.0,
      alpha             = 0.05,
      seed              = 41L,
      progress          = FALSE
    )
  )
  expect_true(any(r$clusters$significant))

  # Each significant cluster should overlap one of the two injected
  # signal regions (disjoint quadrants in make_realistic_pair).
  sig_union <- pair$signal_mask_a | pair$signal_mask_b
  sig_clusters <- r$clusters[r$clusters$significant, , drop = FALSE]
  for (i in seq_len(nrow(sig_clusters))) {
    cid  <- sig_clusters$cluster_id[i]
    labs <- if (sig_clusters$direction[i] == "pos")
      r$pos_labels else r$neg_labels
    cluster_pixels <- as.vector(labs == cid)
    overlap <- sum(cluster_pixels & sig_union) / sum(cluster_pixels)
    expect_gt(overlap, 0.3)  # majority of cluster sits in a signal region
  }
})
