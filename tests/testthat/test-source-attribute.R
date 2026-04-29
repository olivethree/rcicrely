## Tests for §1.2 source-attribute enforcement.
##
## The shared `assert_raw_signal()` helper is wired into every
## variance-based metric. A signal matrix tagged
## `attr(., "source") = "rendered"` must abort by default and pass
## under `acknowledge_scaling = TRUE`. A matrix tagged `"raw"` (or
## untagged but with raw-looking values) must pass silently.

raw_sig <- function(...) {
  s <- make_sig_pair(...)
  attr(s$a, "source") <- "raw"
  attr(s$b, "source") <- "raw"
  s
}

rendered_sig <- function(...) {
  s <- make_sig_pair(...)
  attr(s$a, "source") <- "rendered"
  attr(s$b, "source") <- "rendered"
  s
}

test_that("constructors set source attribute", {
  s_raw <- raw_sig()
  expect_identical(attr(s_raw$a, "source"), "raw")
  expect_identical(attr(s_raw$b, "source"), "raw")

  s_rend <- rendered_sig()
  expect_identical(attr(s_rend$a, "source"), "rendered")
})

test_that("apply_mask_to_signal preserves source attribute", {
  s <- raw_sig()
  mask <- rep(c(TRUE, FALSE), length.out = nrow(s$a))
  out <- apply_mask_to_signal(s$a, mask)
  expect_identical(attr(out, "source"), "raw")
})

test_that("pixel_t_test errors on rendered without acknowledge_scaling", {
  s <- rendered_sig()
  expect_error(
    pixel_t_test(s$a, s$b),
    "rendered"
  )
  out <- suppressWarnings(
    pixel_t_test(s$a, s$b, acknowledge_scaling = TRUE)
  )
  expect_true(is.numeric(out))
  expect_length(out, nrow(s$a))
})

test_that("rel_icc errors on rendered without acknowledge_scaling", {
  s <- rendered_sig()
  expect_error(
    rel_icc(s$a),
    "rendered"
  )
  res <- suppressWarnings(
    rel_icc(s$a, acknowledge_scaling = TRUE)
  )
  expect_s3_class(res, "rcicrely_icc")
})

test_that("rel_dissimilarity errors on rendered without acknowledge_scaling", {
  s <- rendered_sig()
  expect_error(
    rel_dissimilarity(s$a, s$b, n_boot = 10L, progress = FALSE),
    "rendered"
  )
  res <- suppressWarnings(
    rel_dissimilarity(s$a, s$b, n_boot = 10L, progress = FALSE,
                      acknowledge_scaling = TRUE)
  )
  expect_s3_class(res, "rcicrely_dissim")
})

test_that("rel_cluster_test inherits enforcement via pixel_t_test", {
  s <- rendered_sig()
  expect_error(
    rel_cluster_test(s$a, s$b, n_permutations = 10L,
                     progress = FALSE),
    "rendered"
  )
})

test_that("rel_cluster_test acknowledge_scaling cascades", {
  s <- rendered_sig()
  res <- suppressWarnings(
    rel_cluster_test(s$a, s$b, n_permutations = 10L,
                     progress = FALSE,
                     acknowledge_scaling = TRUE)
  )
  expect_s3_class(res, "rcicrely_cluster_test")
})

test_that("raw-tagged matrices pass silently", {
  s <- raw_sig()
  expect_silent(
    suppressMessages(suppressWarnings(
      pixel_t_test(s$a, s$b)
    ))
  )
})

test_that("untagged + heuristic-clean falls through without error", {
  s <- make_sig_pair()
  # No source attribute set; values are raw-like (overall SD < 0.15).
  expect_no_error(
    suppressMessages(suppressWarnings(
      pixel_t_test(s$a, s$b)
    ))
  )
})

test_that("orchestrators forward acknowledge_scaling", {
  s <- rendered_sig()
  expect_error(
    run_between(s$a, s$b,
                n_permutations = 10L, n_boot = 10L,
                progress = FALSE),
    "rendered"
  )
  res <- suppressWarnings(
    run_between(s$a, s$b,
                n_permutations = 10L, n_boot = 10L,
                progress = FALSE, acknowledge_scaling = TRUE)
  )
  expect_s3_class(res, "rcicrely_report")

  expect_error(
    run_within(s$a, n_permutations = 10L, progress = FALSE),
    "rendered"
  )
})
