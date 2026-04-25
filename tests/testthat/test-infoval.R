make_infoval_fixture <- function(n_side = 16L, n_p = 8L, n_pool = 40L,
                                 trials = 100L, snr = 0.3, seed = 1L) {
  n_pix <- n_side * n_side
  set.seed(seed)
  # synthetic noise pool: unit-variance gaussian noise patterns
  noise <- matrix(stats::rnorm(n_pix * n_pool, sd = 1.0),
                  n_pix, n_pool)

  # A circular signal region, anti-aliased
  rr <- row(matrix(0, n_side, n_side)) / n_side - 0.5
  cc <- col(matrix(0, n_side, n_side)) / n_side - 0.5
  d  <- sqrt(rr^2 + cc^2)
  target <- as.vector(pmax(0, 1 - d / 0.2))

  # Producers' signal matrices: simulate responses that partly align
  # with `target`. Each producer does `trials` trials.
  sig_mat <- matrix(NA_real_, n_pix, n_p,
                    dimnames = list(NULL, sprintf("p%02d", seq_len(n_p))))
  tc <- stats::setNames(rep(trials, n_p), colnames(sig_mat))
  for (j in seq_len(n_p)) {
    stim <- sample.int(n_pool, trials, replace = TRUE)
    # preference strength = snr for signal pixels, 0 elsewhere
    proj <- crossprod(noise[, stim, drop = FALSE], target)
    p <- stats::plogis(snr * proj)   # chance response leans toward target
    resp <- ifelse(stats::runif(trials) < p, 1, -1)
    # genMask: mean response by stim, mask build
    uniq <- sort(unique(stim))
    idx  <- match(stim, uniq)
    wts  <- as.vector(tapply(resp, idx, mean))
    sig_mat[, j] <- noise[, uniq, drop = FALSE] %*% wts / length(uniq)
  }

  list(signal_matrix = sig_mat,
       noise_matrix  = noise,
       trial_counts  = tc,
       n_side        = n_side)
}

test_that("infoval returns well-formed per-producer z scores", {
  fx <- make_infoval_fixture(n_side = 16L, n_p = 8L, n_pool = 40L,
                             trials = 100L, snr = 0.4, seed = 1L)
  iv <- suppressWarnings(
    infoval(fx$signal_matrix, fx$noise_matrix, fx$trial_counts,
            iter = 300L, seed = 1L, progress = FALSE)
  )
  expect_s3_class(iv, "rcicrely_infoval")
  expect_equal(length(iv$infoval), 8L)
  expect_equal(length(iv$norms),   8L)
  expect_true(all(iv$norms > 0))
  expect_setequal(names(iv$reference), as.character(unique(fx$trial_counts)))
})

test_that("infoval recovers positive z for signal-aligned producers", {
  # Higher SNR producers should have z > 0 on average; random (SNR=0)
  # producers should have z near 0.
  fx_sig <- make_infoval_fixture(n_side = 16L, n_p = 10L, n_pool = 40L,
                                 trials = 120L, snr = 0.6, seed = 11L)
  iv_sig <- suppressWarnings(
    infoval(fx_sig$signal_matrix, fx_sig$noise_matrix, fx_sig$trial_counts,
            iter = 500L, seed = 11L, progress = FALSE)
  )
  expect_gt(mean(iv_sig$infoval), 0)

  fx_null <- make_infoval_fixture(n_side = 16L, n_p = 10L, n_pool = 40L,
                                  trials = 120L, snr = 0, seed = 22L)
  iv_null <- suppressWarnings(
    infoval(fx_null$signal_matrix, fx_null$noise_matrix,
            fx_null$trial_counts,
            iter = 500L, seed = 22L, progress = FALSE)
  )
  # null should have a roughly zero-centred z distribution
  expect_lt(abs(stats::median(iv_null$infoval)), 2)
})

test_that("infoval keys reference distributions on trial count", {
  # Two groups: one with trials = 80, another with trials = 200;
  # the function should simulate two distinct references.
  fx <- make_infoval_fixture(n_side = 16L, n_p = 6L, n_pool = 40L,
                             trials = 80L, snr = 0.3, seed = 3L)
  fx$trial_counts[1:3] <- 80L
  fx$trial_counts[4:6] <- 200L
  iv <- suppressWarnings(
    infoval(fx$signal_matrix, fx$noise_matrix, fx$trial_counts,
            iter = 200L, seed = 3L, progress = FALSE)
  )
  expect_setequal(names(iv$reference), c("80", "200"))
  # reference norms should differ: more trials -> different (typically
  # smaller) expected norm
  expect_false(identical(iv$reference$`80`, iv$reference$`200`))
})

test_that("infoval mask is applied symmetrically to observed and reference", {
  fx <- make_infoval_fixture(n_side = 16L, n_p = 5L, n_pool = 40L,
                             trials = 100L, snr = 0.3, seed = 4L)
  m  <- face_mask(c(16L, 16L))
  iv_masked <- suppressWarnings(
    infoval(fx$signal_matrix, fx$noise_matrix, fx$trial_counts,
            iter = 200L, mask = m, seed = 4L, progress = FALSE)
  )
  iv_full <- suppressWarnings(
    infoval(fx$signal_matrix, fx$noise_matrix, fx$trial_counts,
            iter = 200L, seed = 4L, progress = FALSE)
  )
  # masked norms must be no larger than full-image norms
  expect_true(all(iv_masked$norms <= iv_full$norms + 1e-10))
})

test_that("infoval cache round-trips on disk", {
  fx <- make_infoval_fixture(n_side = 16L, n_p = 5L, n_pool = 40L,
                             trials = 100L, snr = 0.3, seed = 5L)
  cache <- tempfile(fileext = ".rds")
  on.exit(unlink(cache), add = TRUE)

  iv1 <- suppressWarnings(
    infoval(fx$signal_matrix, fx$noise_matrix, fx$trial_counts,
            iter = 200L, cache_path = cache,
            seed = 5L, progress = FALSE)
  )
  # second call with cache should reuse the reference distributions
  iv2 <- suppressWarnings(
    infoval(fx$signal_matrix, fx$noise_matrix, fx$trial_counts,
            iter = 200L, cache_path = cache,
            seed = 999L, progress = FALSE)
  )
  # norms reproduce (observed data is deterministic)
  expect_equal(iv1$norms, iv2$norms)
  expect_equal(iv1$reference, iv2$reference)
})

test_that("infoval aborts on shape / name mismatch", {
  fx <- make_infoval_fixture(n_side = 16L, n_p = 5L, n_pool = 40L)
  expect_error(
    suppressWarnings(infoval(fx$signal_matrix, fx$noise_matrix[1:50, ],
                             fx$trial_counts, iter = 50L,
                             progress = FALSE)),
    "Row counts"
  )
  tc_bad <- fx$trial_counts
  names(tc_bad) <- NULL
  expect_error(
    suppressWarnings(infoval(fx$signal_matrix, fx$noise_matrix, tc_bad,
                             iter = 50L, progress = FALSE)),
    "named"
  )
})

test_that("2IFC reference is calibrated: random producer z is centred ~ 0", {
  # Regression test for the v0.2.0 bug where the reference distribution
  # sampled stim ids with replacement, inflating the reference Frobenius
  # norm in 2IFC designs (every producer sees every pool item once) and
  # systematically pushing per-producer z below zero. Discovered while
  # replicating Oliveira et al. (2019) Study 1.
  set.seed(123L)
  n_pix  <- 256L
  n_pool <- 60L
  noise  <- matrix(stats::rnorm(n_pix * n_pool), n_pix, n_pool)

  # Many random producers, each one gets a fresh random ±1 over all stims
  n_producers <- 80L
  random_signals <- sapply(seq_len(n_producers), function(j) {
    resp <- sample(c(-1, 1), n_pool, replace = TRUE)
    (noise %*% resp) / n_pool
  })
  colnames(random_signals) <- sprintf("p%02d", seq_len(n_producers))
  tc <- stats::setNames(rep(n_pool, n_producers), colnames(random_signals))

  iv <- suppressWarnings(infoval(
    random_signals, noise, tc,
    iter = 1500L, seed = 1L, progress = FALSE
  ))
  # Under H0 the median should be near 0 (within Monte Carlo).
  # The buggy v0.2.0 implementation produced a median around -1, which
  # this bound catches.
  expect_lt(abs(median(iv$infoval)), 0.5)
})

test_that("face_mask produces a roughly elliptical logical vector", {
  m <- face_mask(c(64L, 64L))
  expect_length(m, 64L * 64L)
  expect_true(is.logical(m))
  # Fraction should be near pi * 0.35 * 0.45 ~ 0.49
  frac <- sum(m) / length(m)
  expect_gt(frac, 0.35)
  expect_lt(frac, 0.55)
})

test_that("face_mask region argument selects sub-regions correctly", {
  full   <- face_mask(c(128L, 128L), region = "full")
  eyes   <- face_mask(c(128L, 128L), region = "eyes")
  mouth  <- face_mask(c(128L, 128L), region = "mouth")
  nose   <- face_mask(c(128L, 128L), region = "nose")
  upper  <- face_mask(c(128L, 128L), region = "upper_face")
  lower  <- face_mask(c(128L, 128L), region = "lower_face")
  # Each sub-region is a strict subset of the full face
  expect_true(all(eyes  <= full))
  expect_true(all(mouth <= full))
  expect_true(all(nose  <= full))
  # Sub-regions are smaller than the full face
  expect_lt(sum(eyes), sum(full))
  expect_lt(sum(mouth), sum(full))
  expect_lt(sum(nose), sum(full))
  # Upper / lower partition (approximately) the full face
  expect_equal(sum(upper) + sum(lower), sum(full),
               tolerance = sum(full) * 0.02)
  # Eyes, mouth, nose are mutually disjoint
  expect_false(any(eyes  & mouth))
  expect_false(any(eyes  & nose))
  expect_false(any(mouth & nose))
})

test_that("face_mask integrates with infoval as a region restrictor", {
  fx <- make_infoval_fixture(n_side = 32L, n_p = 5L, n_pool = 40L,
                             trials = 100L, snr = 0.4, seed = 1L)
  eyes_mask <- face_mask(c(32L, 32L), region = "eyes")
  expect_no_error(
    suppressWarnings(infoval(fx$signal_matrix, fx$noise_matrix,
                             fx$trial_counts,
                             iter = 100L, mask = eyes_mask,
                             seed = 1L, progress = FALSE))
  )
})
