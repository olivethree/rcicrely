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

test_that("rel_icc recovers at realistic scale (64x64, N=30)", {
  skip_on_cran()
  sig <- make_realistic_sig(n_side = 64L, n_p = 30L, snr = 0.3, seed = 51L)
  ic  <- suppressWarnings(rel_icc(sig))
  expect_s3_class(ic, "rcicrely_icc")
  expect_gt(ic$icc_3_1, 0)
  expect_gte(ic$icc_3_k, ic$icc_3_1)
})

test_that("rel_icc matches psych::ICC at realistic scale (2048 pixels)", {
  skip_on_cran()
  skip_if_not_installed("psych")
  skip_if_not_installed("lme4")
  # 2048 pixels is a compromise: enough to stress numerical stability,
  # small enough that psych::ICC() does not blow memory. Inject a
  # small per-target signal so both estimators sit clearly above the
  # no-signal edge: psych::ICC uses lme4-backed REML and clamps
  # near-zero variance components at 0, while rcicrely's mean-squares
  # formula (Shrout & Fleiss 1979) can return small negative values
  # at pure noise. Both implementations are correct in their domain;
  # they only diverge when ms_rows < ms_error (no producer signal).
  set.seed(52L)
  signal <- stats::rnorm(2048L)
  m <- matrix(stats::rnorm(2048L * 15L, sd = 0.5),
              2048L, 15L) + signal
  pkg <- suppressWarnings(rel_icc(m, variants = c("3_1", "3_k",
                                                  "2_1", "2_k")))
  ref <- suppressWarnings(psych::ICC(m)$results)
  get_ref <- function(type) ref$ICC[ref$type == type]
  # 1e-4 tolerance accommodates the small floating-point drift between
  # two different ICC implementations on 30,720 cells; the small
  # 50x8 test above keeps the bit-tight 1e-6 cross-check.
  expect_equal(pkg$icc_3_1, get_ref("ICC3"),  tolerance = 1e-4)
  expect_equal(pkg$icc_3_k, get_ref("ICC3k"), tolerance = 1e-4)
  expect_equal(pkg$icc_2_1, get_ref("ICC2"),  tolerance = 1e-4)
  expect_equal(pkg$icc_2_k, get_ref("ICC2k"), tolerance = 1e-4)
})

test_that("rel_icc emits large-image warning above 50k pixels", {
  skip_on_cran()
  reset_session_warnings()
  # Use raw-mask-scale noise (SD ~ 0.03) and 30 raters so the
  # warnings we do NOT want firing (small-N, looks_scaled) are
  # silent; only the ICC large-image advisory should fire.
  m <- matrix(stats::rnorm(51200L * 30L, sd = 0.03), 51200L, 30L)

  # Warning fires when ICC(3,k) is requested.
  reset_session_warnings()
  expect_warning(
    rel_icc(m, variants = c("3_1", "3_k")),
    regexp = "ICC\\(3,k\\) tends toward 1"
  )

  # Silenced via option.
  reset_session_warnings()
  op <- options(rcicrely.silence_icc_warning = TRUE)
  on.exit(options(op), add = TRUE)
  expect_no_warning(
    rel_icc(m, variants = c("3_1", "3_k")),
    message = "ICC\\(3,k\\) tends toward 1"
  )
})
