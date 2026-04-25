## Tests for the new plot helpers added in v0.2.x:
##   - plot_agreement_map()
##   - compare_dissimilarity()
##
## These don't assert on appearance — only that the functions run
## without errors and return the expected structure.

test_that("plot_agreement_map runs and returns t_map of right length", {
  set.seed(1)
  n_side <- 32L
  n_pix  <- n_side * n_side
  rr <- row(matrix(0, n_side, n_side)) / n_side - 0.5
  cc <- col(matrix(0, n_side, n_side)) / n_side - 0.5
  d  <- sqrt(rr^2 + cc^2)
  mask_vec <- as.vector(pmax(0, 1 - d / 0.2))
  signal <- 0.5 * outer(mask_vec, runif(20L, 0.7, 1.3)) +
            matrix(rnorm(n_pix * 20L), n_pix, 20L)
  attr(signal, "img_dims") <- c(n_side, n_side)

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  out <- suppressWarnings(plot_agreement_map(signal))
  expect_type(out, "list")
  expect_named(out, c("t_map", "n", "img_dims", "mask"))
  expect_length(out$t_map, n_pix)
  expect_equal(out$n, 20L)
  # Signal pixels (centre) should have higher t than periphery
  centre_idx     <- which(mask_vec > 0.5)
  periphery_idx  <- which(mask_vec == 0)
  expect_gt(mean(out$t_map[centre_idx]),
            mean(out$t_map[periphery_idx]))
})

test_that("plot_agreement_map respects mask and threshold args", {
  set.seed(2)
  signal <- matrix(rnorm(64L * 10L), 64L, 10L)
  attr(signal, "img_dims") <- c(8L, 8L)

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  m <- rep(c(TRUE, FALSE), each = 32L)
  expect_no_error(
    plot_agreement_map(signal, mask = m, threshold = 0.5)
  )
})

test_that("plot_agreement_map errors on bad input", {
  expect_error(plot_agreement_map("not a matrix"),
               "numeric matrix")
  signal <- matrix(rnorm(50L * 5L), 50L, 5L)
  expect_error(plot_agreement_map(signal),
               "Cannot infer")  # not a perfect square, no img_dims
})

test_that("compare_dissimilarity returns a tidy table and runs", {
  pair <- make_sig_pair(n_pix = 256L, n_p = 8L, seed = 1L)
  d1 <- suppressWarnings(rel_dissimilarity(pair$a, pair$b,
    n_boot = 100L, seed = 1L, progress = FALSE))
  d2 <- suppressWarnings(rel_dissimilarity(pair$a, pair$a,
    n_boot = 100L, seed = 2L, progress = FALSE))

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  out <- compare_dissimilarity("A vs B" = d1, "A vs A" = d2)
  expect_s3_class(out, "data.frame")
  expect_named(out,
               c("label", "observed", "ci_low", "ci_high"))
  expect_equal(nrow(out), 2L)
  # A-vs-A should have observed distance lower than A-vs-B
  expect_lt(out$observed[out$label == "A vs A"],
            out$observed[out$label == "A vs B"])
})

test_that("compare_dissimilarity errors on missing names or wrong class", {
  pair <- make_sig_pair(n_pix = 64L, n_p = 5L, seed = 1L)
  d <- suppressWarnings(rel_dissimilarity(pair$a, pair$b,
    n_boot = 50L, seed = 1L, progress = FALSE))
  expect_error(compare_dissimilarity(d), "named")
  expect_error(compare_dissimilarity("A" = list(foo = 1)),
               "rcicrely_dissim")
  expect_error(compare_dissimilarity(), "at least one")
})
