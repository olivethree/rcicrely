test_that("run_within wraps split_half + icc (reliability metrics only)", {
  sig <- make_sig(256L, 10L, seed = 1L)
  rep <- suppressWarnings(
    run_within(sig, n_permutations = 50L, seed = 1L, progress = FALSE)
  )
  expect_s3_class(rep, "rcicrely_report")
  expect_setequal(names(rep$results), c("split_half", "icc"))
  expect_s3_class(rep$results$split_half, "rcicrely_split_half")
  expect_s3_class(rep$results$icc,        "rcicrely_icc")
  expect_equal(rep$method, "within")
})

test_that("run_within no longer includes loo; use rel_loo() separately", {
  sig <- make_sig(256L, 10L, seed = 2L)
  rep <- suppressWarnings(
    run_within(sig, n_permutations = 30L, seed = 2L, progress = FALSE)
  )
  expect_false("loo" %in% names(rep$results))
  # but rel_loo() still works standalone
  loo <- suppressWarnings(rel_loo(sig))
  expect_s3_class(loo, "rcicrely_loo")
})

test_that("report round-trips through saveRDS", {
  sig <- make_sig(256L, 10L)
  rep <- suppressWarnings(
    run_within(sig, n_permutations = 30L, seed = 1L, progress = FALSE)
  )
  tmp <- tempfile(fileext = ".rds")
  saveRDS(rep, tmp)
  back <- readRDS(tmp)
  expect_identical(rep$results$split_half$r_hh,
                   back$results$split_half$r_hh)
  expect_identical(rep$results$icc$icc_3_1,
                   back$results$icc$icc_3_1)
})
