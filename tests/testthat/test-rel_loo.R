test_that("loo returns a correlation per producer", {
  sig <- make_sig(256L, 10L, seed = 3L)
  lo <- suppressWarnings(rel_loo(sig))
  expect_s3_class(lo, "rcicrely_loo")
  expect_equal(length(lo$correlations), 10L)
  expect_true(all(is.finite(lo$correlations)))
})

test_that("loo flags an injected outlier under SD rule", {
  set.seed(10)
  sig <- make_sig(256L, 10L, seed = 3L)
  outlier <- matrix(runif(256L, -1, 1), 256L, 1L)
  sig_plus <- cbind(sig, outlier)
  colnames(sig_plus) <- c(sprintf("p%02d", 1:10), "outlier")
  lo <- suppressWarnings(
    rel_loo(sig_plus, flag_threshold = 2.5, flag_method = "sd")
  )
  expect_true("outlier" %in% lo$flagged)
  expect_equal(lo$flag_method, "sd")
})

test_that("loo flags an injected outlier under MAD rule", {
  set.seed(11)
  sig <- make_sig(256L, 10L, seed = 3L)
  outlier <- matrix(runif(256L, -1, 1), 256L, 1L)
  sig_plus <- cbind(sig, outlier)
  colnames(sig_plus) <- c(sprintf("p%02d", 1:10), "outlier")
  lo_mad <- suppressWarnings(
    rel_loo(sig_plus, flag_threshold = 2.5, flag_method = "mad")
  )
  expect_true("outlier" %in% lo_mad$flagged)
  expect_equal(lo_mad$flag_method, "mad")
  expect_true(is.finite(lo_mad$median_r))
  expect_true(is.finite(lo_mad$mad_r))
})

test_that("flag_threshold_sd is accepted as a deprecated alias", {
  sig <- make_sig(256L, 10L, seed = 3L)
  lo_old <- suppressWarnings(rel_loo(sig, flag_threshold_sd = 2.5))
  lo_new <- suppressWarnings(rel_loo(sig, flag_threshold    = 2.5))
  expect_equal(lo_old$threshold, lo_new$threshold)
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

test_that("loo result carries flag_method and flag_threshold metadata", {
  sig <- make_sig(256L, 10L)
  lo <- suppressWarnings(rel_loo(sig, flag_threshold = 3, flag_method = "mad"))
  expect_equal(lo$flag_method, "mad")
  expect_equal(lo$flag_threshold, 3)
})
