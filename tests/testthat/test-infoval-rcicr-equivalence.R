## §2.4 bullet 1: matched-config equivalence with rcicr.
##
## When `n_trials == n_pool` for every producer and `mask = NULL`,
## `rcicrely::infoval()` and `rcicr::computeInfoVal2IFC()` are
## mathematically equivalent: rcicr samples random ±1 responses
## across the fixed-order pool, while rcicrely samples without
## replacement (a random permutation of the pool) and applies a
## random ±1. The Frobenius norm of the resulting CI is invariant
## to the column permutation, so the two reference distributions
## have the same statistical law.
##
## This test asserts equivalence at the *distributional* level
## (matched seed -> same reference median and MAD within Monte
## Carlo tolerance) rather than pointwise — pointwise identity
## would require coupling the RNG streams across the two
## implementations, which is not portable.
##
## §2.4 bullet 3 (Brief-RC divergence quantification) is a longer
## simulation lives in `temp/simulations/`; only the equivalence
## smoke test belongs in the package test suite.

skip_if_not_installed("rcicr")

test_that("infoval matches rcicr::computeInfoVal2IFC at n_trials == n_pool", {
  skip_on_cran()

  set.seed(42L)
  n_pix  <- 64L
  n_pool <- 50L
  iter   <- 2000L

  # Build a synthetic noise matrix and synthetic CI
  noise  <- matrix(stats::rnorm(n_pix * n_pool), n_pix, n_pool)
  resp   <- sample(c(-1, 1), n_pool, replace = TRUE)
  ci_vec <- as.vector((noise %*% resp) / n_pool)
  ci_mat <- matrix(ci_vec, ncol = 1L)
  colnames(ci_mat) <- "p001"
  attr(ci_mat, "source") <- "raw"

  # rcicrely reference distribution
  rely_ref <- rcicrely:::simulate_reference_norms(
    noise_matrix = noise,
    n_trials     = n_pool,
    iter         = iter
  )

  # rcicr-style reference distribution: stim columns in fixed
  # order, random ±1 responses, divide by n_pool, Frobenius norm.
  rcicr_ref <- replicate(iter, {
    r  <- sample(c(-1, 1), n_pool, replace = TRUE)
    ci <- (noise %*% r) / n_pool
    norm(ci, "f")
  })

  # Distributional equivalence: medians within ~3% of each other,
  # MADs within ~5%. The bounds are loose enough that 2000 iters
  # of MC noise won't fail them, tight enough that a real bug
  # (e.g. a regression to the with-replacement scheme) would.
  m_rely  <- stats::median(rely_ref)
  m_rcicr <- stats::median(rcicr_ref)
  mad_rely  <- stats::mad(rely_ref)
  mad_rcicr <- stats::mad(rcicr_ref)

  expect_lt(abs(m_rely - m_rcicr) / m_rcicr,   0.03)
  expect_lt(abs(mad_rely - mad_rcicr) / mad_rcicr, 0.05)
})
