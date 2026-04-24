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

test_that("loo summary_df has the expected columns including z_score", {
  sig <- make_sig(256L, 10L)
  lo <- suppressWarnings(rel_loo(sig))
  expect_setequal(
    colnames(lo$summary_df),
    c("participant_id", "correlation", "z_score", "flag")
  )
  expect_equal(nrow(lo$summary_df), 10L)
})

test_that("z_scores are returned and match flag_method", {
  sig <- make_sig(256L, 10L)
  lo_sd  <- suppressWarnings(rel_loo(sig, flag_method = "sd"))
  lo_mad <- suppressWarnings(rel_loo(sig, flag_method = "mad"))
  expect_equal(length(lo_sd$z_scores), 10L)
  expect_equal(names(lo_sd$z_scores), names(lo_sd$correlations))
  # Under SD rule, z_scores have mean ~0 and sd ~1.
  expect_lt(abs(mean(lo_sd$z_scores)), 1e-8)
  expect_lt(abs(stats::sd(lo_sd$z_scores) - 1), 1e-8)
  # Under MAD rule, z_scores have median ~0 and mad ~1.
  expect_lt(abs(stats::median(lo_mad$z_scores)), 1e-8)
  expect_lt(abs(stats::mad(lo_mad$z_scores) - 1), 1e-8)
})

test_that("rel_loo_z accepts either a signal matrix or a result", {
  sig <- make_sig(256L, 10L, seed = 3L)
  df_from_mat <- suppressWarnings(rel_loo_z(sig))
  lo          <- suppressWarnings(rel_loo(sig))
  df_from_obj <- rel_loo_z(lo)
  expect_equal(nrow(df_from_mat), 10L)
  expect_setequal(
    colnames(df_from_mat),
    c("participant_id", "correlation", "z_score", "flag")
  )
  # sorted ascending by z_score
  expect_equal(df_from_mat$z_score, sort(df_from_mat$z_score))
  expect_equal(df_from_obj, lo$summary_df)
})

test_that("rel_loo_z errors on invalid input", {
  expect_error(rel_loo_z("not a matrix"), class = "rlang_error")
  expect_error(rel_loo_z(list(a = 1)),    class = "rlang_error")
})

test_that("loo result carries flag_method and flag_threshold metadata", {
  sig <- make_sig(256L, 10L)
  lo <- suppressWarnings(rel_loo(sig, flag_threshold = 3, flag_method = "mad"))
  expect_equal(lo$flag_method, "mad")
  expect_equal(lo$flag_threshold, 3)
})
