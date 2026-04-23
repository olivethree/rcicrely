test_that("rel_icc returns ICC(3,1) and ICC(3,k) by default", {
  sig <- make_sig(256L, 10L, seed = 2L)
  ic <- suppressWarnings(rel_icc(sig))
  expect_s3_class(ic, "rcicrely_icc")
  expect_true(is.finite(ic$icc_3_1))
  expect_true(is.finite(ic$icc_3_k))
  expect_true(is.na(ic$icc_2_1))   # not requested
  expect_true(is.na(ic$icc_2_k))
  # ICC(3,k) >= ICC(3,1) for positive signal
  expect_gte(ic$icc_3_k, ic$icc_3_1)
  expect_match(ic$model, "ICC\\(3,\\*\\)")
})

test_that("variants argument exposes ICC(2,*)", {
  sig <- make_sig(256L, 10L)
  ic <- suppressWarnings(rel_icc(sig, variants = c("3_1", "3_k",
                                                   "2_1", "2_k")))
  expect_true(all(is.finite(c(ic$icc_3_1, ic$icc_3_k,
                              ic$icc_2_1, ic$icc_2_k))))
})

test_that("rel_icc matches psych::ICC on a small matrix", {
  skip_if_not_installed("psych")
  skip_if_not_installed("lme4")  # psych::ICC needs lme4 for some calls
  set.seed(10)
  m <- matrix(rnorm(50L * 8L), 50L, 8L)
  pkg <- suppressWarnings(rel_icc(m, variants = c("3_1", "3_k",
                                                  "2_1", "2_k")))
  ref <- suppressWarnings(psych::ICC(m)$results)
  get_ref <- function(type) {
    ref$ICC[ref$type == type]
  }
  expect_equal(pkg$icc_3_1, get_ref("ICC3"),  tolerance = 1e-6)
  expect_equal(pkg$icc_3_k, get_ref("ICC3k"), tolerance = 1e-6)
  expect_equal(pkg$icc_2_1, get_ref("ICC2"),  tolerance = 1e-6)
  expect_equal(pkg$icc_2_k, get_ref("ICC2k"), tolerance = 1e-6)
})
