test_that("split_half returns a well-formed result", {
  sig <- make_sig(n_pix = 256L, n_p = 12L, signal_strength = 0.3, seed = 2L)
  r <- suppressWarnings(
    rel_split_half(sig, n_permutations = 200L, seed = 1L, progress = FALSE)
  )
  expect_s3_class(r, "rcicrely_split_half")
  expect_equal(length(r$distribution), 200L)
  expect_true(r$r_hh > 0)            # synthetic signal is aligned
  expect_gte(r$r_sb, r$r_hh)         # Spearman-Brown >= raw when r > 0
  expect_true(r$ci_95[1] <= r$r_hh && r$ci_95[2] >= r$r_hh - 0.01)
})

test_that("split_half is reproducible under seed", {
  sig <- make_sig(256L, 12L, seed = 5L)
  a <- suppressWarnings(
    rel_split_half(sig, n_permutations = 50L, seed = 42L, progress = FALSE)
  )
  b <- suppressWarnings(
    rel_split_half(sig, n_permutations = 50L, seed = 42L, progress = FALSE)
  )
  expect_identical(a$r_hh, b$r_hh)
  expect_identical(a$distribution, b$distribution)
})

test_that("split_half handles odd N via per-perm drop", {
  sig <- make_sig(256L, 11L, seed = 7L)   # odd N
  r <- suppressWarnings(
    rel_split_half(sig, n_permutations = 50L, seed = 1L, progress = FALSE)
  )
  expect_true(is.finite(r$r_hh))
  expect_equal(r$n_participants, 11L)
})

test_that("split_half aborts below minimum participant count", {
  sig <- matrix(rnorm(100L * 3L), 100L, 3L)
  expect_error(
    rel_split_half(sig, n_permutations = 50L, progress = FALSE),
    "at least"
  )
})

test_that("split_half runs at realistic scale (128x128, N=30)", {
  skip_on_cran()
  sig <- make_realistic_sig(n_side = 128L, n_p = 30L, snr = 1.0, seed = 11L)
  r   <- suppressWarnings(
    rel_split_half(sig, n_permutations = 200L, seed = 11L, progress = FALSE)
  )
  expect_s3_class(r, "rcicrely_split_half")
  expect_equal(r$n_participants, 30L)
  expect_equal(length(r$distribution), 200L)
  # Sanity bounds only: precise recovery curves live in the
  # simulation centrepiece (data-raw/simulations/). At SNR = 1.0
  # with 30 producers and a ~13%-area signal region, we expect a
  # clearly non-zero r_hh; loose bound protects against Monte Carlo
  # noise and against signal-coverage fluctuations.
  expect_gt(r$r_hh, 0)
  expect_gte(r$r_sb, r$r_hh)
})

test_that("split_half matches a hand-rolled permutation at realistic scale", {
  skip_on_cran()
  sig <- make_realistic_sig(n_side = 64L, n_p = 30L, snr = 0.3, seed = 12L)

  # rcicrely run
  r <- suppressWarnings(
    rel_split_half(sig, n_permutations = 100L, seed = 7L, progress = FALSE)
  )

  # Hand-rolled reference: same seed semantics (set.seed then loop).
  ref <- local({
    set.seed(7L)
    N  <- ncol(sig)
    h  <- N %/% 2L
    odd <- (N %% 2L) == 1L
    out <- numeric(100L)
    for (i in seq_len(100L)) {
      idx <- sample.int(N)
      if (odd) idx <- idx[-1L]
      h1 <- idx[seq_len(h)]
      h2 <- idx[(h + 1L):length(idx)]
      v1 <- rowMeans(sig[, h1, drop = FALSE])
      v2 <- rowMeans(sig[, h2, drop = FALSE])
      out[i] <- stats::cor(v1, v2)
    }
    out
  })

  expect_equal(r$distribution, ref, tolerance = 1e-10)
  expect_equal(r$r_hh,         mean(ref), tolerance = 1e-10)
})
