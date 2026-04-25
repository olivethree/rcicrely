## Tests for the v0.2.x `mask = NULL` argument on the rel_*() family.
## Verifies that:
##  - Passing a mask is equivalent to manually subsetting the signal
##    matrix beforehand (drop pattern: rel_split_half, rel_loo, rel_icc,
##    pixel_t_test, rel_dissimilarity).
##  - For rel_cluster_test (zero-out pattern), masked pixels never join
##    a cluster.
##  - mask validation errors are clean (wrong type, wrong length, too few
##    pixels).

suppressMessages(reset_session_warnings())

# Use a 16x16 signal matrix with a known central signal patch
make_test_signal <- function(n_p = 12L, n_side = 16L, snr = 0.3,
                             seed = 1L) {
  set.seed(seed)
  n_pix <- n_side * n_side
  rr <- row(matrix(0, n_side, n_side)) / n_side - 0.5
  cc <- col(matrix(0, n_side, n_side)) / n_side - 0.5
  d  <- sqrt(rr^2 + cc^2)
  m  <- as.vector(pmax(0, 1 - d / 0.25))
  out <- snr * outer(m, runif(n_p, 0.7, 1.3)) +
         matrix(rnorm(n_pix * n_p, sd = 0.05), n_pix, n_p)
  attr(out, "img_dims") <- c(n_side, n_side)
  out
}

test_that("rel_split_half mask = full equals direct subset", {
  sig  <- make_test_signal(seed = 1L)
  mask <- face_mask(c(16L, 16L), region = "full")
  r_mask <- suppressWarnings(rel_split_half(sig, n_permutations = 50L,
                                            mask = mask, seed = 1L,
                                            progress = FALSE))
  r_man  <- suppressWarnings(rel_split_half(sig[mask, , drop = FALSE],
                                            n_permutations = 50L,
                                            seed = 1L, progress = FALSE))
  expect_equal(r_mask$r_hh, r_man$r_hh)
  expect_equal(r_mask$distribution, r_man$distribution)
})

test_that("rel_loo mask works and matches direct subset", {
  sig  <- make_test_signal(seed = 2L)
  mask <- face_mask(c(16L, 16L), region = "full")
  l_mask <- suppressWarnings(rel_loo(sig, mask = mask))
  l_man  <- suppressWarnings(rel_loo(sig[mask, , drop = FALSE]))
  expect_equal(l_mask$correlations, l_man$correlations)
})

test_that("rel_icc mask works and matches direct subset", {
  sig  <- make_test_signal(seed = 3L)
  mask <- face_mask(c(16L, 16L), region = "full")
  i_mask <- suppressWarnings(rel_icc(sig, mask = mask))
  i_man  <- suppressWarnings(rel_icc(sig[mask, , drop = FALSE]))
  expect_equal(i_mask$icc_3_1, i_man$icc_3_1)
  expect_equal(i_mask$icc_3_k, i_man$icc_3_k)
})

test_that("pixel_t_test mask returns length sum(mask)", {
  a <- make_test_signal(seed = 4L)
  b <- make_test_signal(seed = 5L)
  mask <- face_mask(c(16L, 16L), region = "eyes")
  t_mask <- pixel_t_test(a, b, mask = mask)
  expect_length(t_mask, sum(mask))
  # Should equal the manual subset version
  t_man <- pixel_t_test(a[mask, , drop = FALSE],
                        b[mask, , drop = FALSE])
  expect_equal(t_mask, t_man)
})

test_that("rel_dissimilarity mask updates n_pixels and matches direct", {
  a <- make_test_signal(seed = 6L)
  b <- make_test_signal(seed = 7L)
  mask <- face_mask(c(16L, 16L), region = "full")
  d_mask <- suppressWarnings(rel_dissimilarity(
    a, b, n_boot = 100L, mask = mask, seed = 1L, progress = FALSE
  ))
  d_man  <- suppressWarnings(rel_dissimilarity(
    a[mask, , drop = FALSE], b[mask, , drop = FALSE],
    n_boot = 100L, seed = 1L, progress = FALSE
  ))
  expect_equal(d_mask$euclidean, d_man$euclidean)
  expect_equal(d_mask$n_pixels, sum(mask))
  expect_equal(d_mask$euclidean_normalised,
               d_mask$euclidean / sqrt(sum(mask)))
})

test_that("rel_cluster_test mask zeroes outside-mask pixels", {
  a <- make_test_signal(snr = 0.6, seed = 8L)
  b <- make_test_signal(snr = 0.0, seed = 9L)
  mask <- face_mask(c(16L, 16L), region = "full")
  cl <- suppressWarnings(rel_cluster_test(
    a, b, img_dims = c(16L, 16L),
    n_permutations = 100L, cluster_threshold = 1.5,
    mask = mask, seed = 1L, progress = FALSE
  ))
  expect_s3_class(cl, "rcicrely_cluster_test")
  # Pixels outside the mask should have observed_t = 0
  expect_true(all(cl$observed_t[!mask] == 0))
  # No cluster pixels outside the mask
  if (any(cl$pos_labels != 0L)) {
    expect_true(all(cl$pos_labels[!mask] == 0L))
  }
  if (any(cl$neg_labels != 0L)) {
    expect_true(all(cl$neg_labels[!mask] == 0L))
  }
})

test_that("mask argument errors clearly on bad input", {
  sig <- make_test_signal(seed = 10L)
  expect_error(
    suppressWarnings(rel_split_half(sig, mask = c("a", "b"),
                                    n_permutations = 10L,
                                    progress = FALSE)),
    "logical vector"
  )
  expect_error(
    suppressWarnings(rel_split_half(sig, mask = rep(TRUE, 5L),
                                    n_permutations = 10L,
                                    progress = FALSE)),
    "length must match"
  )
  expect_error(
    suppressWarnings(rel_split_half(sig,
                                    mask = c(rep(TRUE, 2L),
                                             rep(FALSE, 254L)),
                                    n_permutations = 10L,
                                    progress = FALSE)),
    "too few pixels"
  )
})

test_that("run_within and run_between thread mask through correctly", {
  sig <- make_test_signal(seed = 11L)
  pair_b <- make_test_signal(seed = 12L)
  mask <- face_mask(c(16L, 16L), region = "full")
  rep_w <- suppressWarnings(run_within(sig, n_permutations = 50L,
                                       mask = mask, seed = 1L,
                                       progress = FALSE))
  expect_s3_class(rep_w$results$split_half, "rcicrely_split_half")
  # Should match the direct rel_split_half call with same mask
  r_direct <- suppressWarnings(rel_split_half(sig, mask = mask,
                                              n_permutations = 50L,
                                              seed = 1L,
                                              progress = FALSE))
  expect_equal(rep_w$results$split_half$r_hh, r_direct$r_hh)
})
