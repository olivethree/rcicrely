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
