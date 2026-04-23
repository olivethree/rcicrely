test_that("loo returns a correlation per producer", {
  sig <- make_sig(256L, 10L, seed = 3L)
  lo <- suppressWarnings(rel_loo(sig))
  expect_s3_class(lo, "rcicrely_loo")
  expect_equal(length(lo$correlations), 10L)
  expect_true(all(is.finite(lo$correlations)))
})

test_that("loo flags an injected outlier", {
  set.seed(10)
  sig <- make_sig(256L, 10L, seed = 3L)
  outlier <- matrix(runif(256L, -1, 1), 256L, 1L)
  sig_plus <- cbind(sig, outlier)
  colnames(sig_plus) <- c(sprintf("p%02d", 1:10), "outlier")
  lo <- suppressWarnings(rel_loo(sig_plus, flag_threshold_sd = 2.5))
  expect_true("outlier" %in% lo$flagged)
})

test_that("loo summary_df has the expected columns", {
  sig <- make_sig(256L, 10L)
  lo <- suppressWarnings(rel_loo(sig))
  expect_setequal(
    colnames(lo$summary_df),
    c("participant_id", "correlation", "flag")
  )
  expect_equal(nrow(lo$summary_df), 10L)
})
