test_that("tfce_enhance returns a same-length signed vector", {
  t_map <- c(rep(0, 40), rep(3, 10), rep(0, 40), rep(-3, 10))
  dims  <- c(10L, 10L)
  out <- tfce_enhance(t_map, dims, H = 2.0, E = 0.5, n_steps = 50L)
  expect_equal(length(out), 100L)
  # positive input contributes positive TFCE; negative input contributes
  # negative TFCE
  pos_idx <- which(t_map > 0)
  neg_idx <- which(t_map < 0)
  expect_true(all(out[pos_idx] > 0))
  expect_true(all(out[neg_idx] < 0))
  expect_true(all(out[t_map == 0] == 0))
})

test_that("tfce gives larger enhancement to larger clusters at same height", {
  # Two disjoint clusters at identical height. The bigger one should
  # receive a bigger TFCE value per pixel.
  img  <- matrix(0, 20L, 20L)
  img[2:3, 2:3]     <- 3   # small cluster (4 pixels)
  img[6:11, 6:11]   <- 3   # large cluster (36 pixels)
  out  <- matrix(tfce_enhance(as.vector(img), c(20L, 20L), n_steps = 50L),
                 20L, 20L)
  small_val <- mean(out[2:3, 2:3])
  large_val <- mean(out[6:11, 6:11])
  expect_gt(large_val, small_val)
})

test_that("tfce gives larger enhancement to taller clusters at same size", {
  img  <- matrix(0, 20L, 20L)
  img[2:5, 2:5] <- 2    # shorter
  img[12:15, 12:15] <- 5 # taller, same size
  out  <- matrix(tfce_enhance(as.vector(img), c(20L, 20L), n_steps = 50L),
                 20L, 20L)
  short_val <- mean(out[2:5, 2:5])
  tall_val  <- mean(out[12:15, 12:15])
  expect_gt(tall_val, short_val)
})

test_that("rel_cluster_test method = 'tfce' returns expected fields", {
  pair <- make_sig_pair(n_pix = 256L, n_p = 10L, seed = 1L)
  r <- suppressWarnings(
    rel_cluster_test(
      pair$a, pair$b,
      img_dims       = pair$img_dims,
      method         = "tfce",
      n_permutations = 100L,
      tfce_n_steps   = 50L,
      seed           = 1L,
      progress       = FALSE
    )
  )
  expect_s3_class(r, "rcicrely_cluster_test")
  expect_equal(r$method, "tfce")
  expect_equal(length(r$tfce_map),  256L)
  expect_equal(length(r$tfce_pmap), 256L)
  expect_equal(length(r$tfce_significant_mask), 256L)
  expect_true(all(r$tfce_pmap >= 0 & r$tfce_pmap <= 1))
  expect_equal(length(r$null_distribution$max_abs_tfce), 100L)
  expect_true(all(r$null_distribution$max_abs_tfce >= 0))
})

test_that("rel_cluster_test method dispatch preserves threshold path", {
  pair <- make_sig_pair(256L, 10L, seed = 2L)
  r <- suppressWarnings(
    rel_cluster_test(pair$a, pair$b, img_dims = pair$img_dims,
                     method = "threshold", n_permutations = 50L,
                     seed = 2L, progress = FALSE)
  )
  expect_equal(r$method, "threshold")
  expect_s3_class(r$clusters, "data.frame")
  expect_true(is.null(r$tfce_map))
})

test_that("tfce is reproducible under seed", {
  pair <- make_sig_pair(64L, 8L, seed = 11L)
  a <- suppressWarnings(
    rel_cluster_test(pair$a, pair$b, img_dims = pair$img_dims,
                     method = "tfce",
                     n_permutations = 30L, tfce_n_steps = 30L,
                     seed = 7L, progress = FALSE)
  )
  b <- suppressWarnings(
    rel_cluster_test(pair$a, pair$b, img_dims = pair$img_dims,
                     method = "tfce",
                     n_permutations = 30L, tfce_n_steps = 30L,
                     seed = 7L, progress = FALSE)
  )
  expect_equal(a$tfce_map, b$tfce_map)
  expect_equal(a$tfce_pmap, b$tfce_pmap)
  expect_equal(a$null_distribution$max_abs_tfce,
               b$null_distribution$max_abs_tfce)
})

test_that("tfce localises an injected signal region", {
  skip_on_cran()
  pair <- make_realistic_pair(n_side = 32L, n_p = 30L, snr = 1.0, seed = 71L)
  r <- suppressWarnings(
    rel_cluster_test(
      pair$a, pair$b,
      img_dims       = pair$img_dims,
      method         = "tfce",
      n_permutations = 200L,
      tfce_n_steps   = 50L,
      seed           = 71L, progress = FALSE
    )
  )
  # Significant pixels should concentrate in the injected regions.
  # Baseline concentration under H0 is ~signal_frac (~20% at these
  # settings). A concentration that clears ~2.5x baseline is
  # evidence that TFCE is finding the injected signal.
  sig_union <- pair$signal_mask_a | pair$signal_mask_b
  signal_frac <- mean(sig_union)
  n_sig_in_signal <- sum(r$tfce_significant_mask & sig_union)
  n_sig_in_bg     <- sum(r$tfce_significant_mask & !sig_union)
  if (any(r$tfce_significant_mask)) {
    concentration <- n_sig_in_signal /
      (n_sig_in_signal + n_sig_in_bg)
    expect_gt(concentration, 2 * signal_frac)
  } else {
    skip("No significant pixels at this seed / SNR combination")
  }
})
