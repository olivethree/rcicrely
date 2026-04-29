## §1.6: smoke + composition tests for plot_ci_overlay()

test_that("plot_ci_overlay returns a 3D raster of the right size", {
  sig <- raw_pure_signal(n_pix = 256L, n_p = 12L)
  base <- matrix(stats::runif(256L), 16L, 16L)
  grDevices::pdf(file = tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  out <- plot_ci_overlay(sig, base, img_dims = c(16L, 16L))
  expect_equal(dim(out), c(16L, 16L, 3L))
  expect_true(all(out >= 0 & out <= 1))
})

test_that("plot_ci_overlay accepts a numeric vector signal", {
  sig_vec <- stats::rnorm(256L)
  base    <- matrix(stats::runif(256L), 16L, 16L)
  grDevices::pdf(file = tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  out <- plot_ci_overlay(sig_vec, base, img_dims = c(16L, 16L))
  expect_equal(dim(out), c(16L, 16L, 3L))
})

test_that("plot_ci_overlay errors on dim mismatch", {
  sig  <- raw_pure_signal(n_pix = 256L, n_p = 12L)
  base <- matrix(stats::runif(64L), 8L, 8L)
  grDevices::pdf(file = tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_error(
    plot_ci_overlay(sig, base, img_dims = c(8L, 8L)),
    "Pixel count"
  )
})

test_that("plot_ci_overlay accepts an agreement_map_test result for contours", {
  sig <- raw_pure_signal(n_pix = 256L, n_p = 16L)
  base <- matrix(stats::runif(256L), 16L, 16L)
  test <- suppressWarnings(
    agreement_map_test(sig, n_permutations = 100L,
                       seed = 1L, progress = FALSE)
  )
  grDevices::pdf(file = tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(
    plot_ci_overlay(sig, base, img_dims = c(16L, 16L), test = test)
  )
})

test_that("plot_ci_overlay rejects a wrong-class test argument", {
  sig <- raw_pure_signal(n_pix = 256L, n_p = 12L)
  base <- matrix(stats::runif(256L), 16L, 16L)
  grDevices::pdf(file = tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_error(
    plot_ci_overlay(sig, base, img_dims = c(16L, 16L),
                    test = list(a = 1)),
    "agreement_map_test"
  )
})

test_that("threshold zeroes out small-magnitude pixels", {
  sig_vec <- c(rep(0.01, 128L), rep(0.5, 128L))
  base    <- matrix(0.5, 16L, 16L)
  grDevices::pdf(file = tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  out_no_thr <- plot_ci_overlay(sig_vec, base, img_dims = c(16L, 16L))
  out_thr    <- plot_ci_overlay(sig_vec, base, img_dims = c(16L, 16L),
                                threshold = 0.1)
  # Pixels that fall under threshold should equal the base
  # (alpha = 0 -> base shows through unchanged).
  base_pixels <- which(abs(sig_vec) < 0.1)
  expect_true(all(
    abs(out_thr[, , 1][base_pixels] - 0.5) < 1e-12
  ))
})
