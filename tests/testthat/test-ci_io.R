test_that("read_cis reads a directory of PNGs into pixels x participants", {
  skip_if_not_installed("png")
  tmp <- withr::local_tempdir()
  side <- 16L
  set.seed(1)
  for (i in seq_len(4L)) {
    img <- matrix(runif(side * side), side, side)
    write_tiny_png(file.path(tmp, sprintf("p%02d.png", i)), img)
  }
  cis <- read_cis(tmp)
  expect_true(is.matrix(cis))
  expect_equal(nrow(cis), side * side)
  expect_equal(ncol(cis), 4L)
  expect_setequal(colnames(cis), sprintf("p%02d", 1:4))
  expect_equal(attr(cis, "img_dims"), c(side, side))
})

test_that("read_cis errors cleanly when directory is missing or empty", {
  expect_error(read_cis("/this/does/not/exist"), "Directory not found")
  tmp <- withr::local_tempdir()
  expect_error(read_cis(tmp), "No image files")
})

test_that("read_cis errors when image dimensions differ across files", {
  skip_if_not_installed("png")
  tmp <- withr::local_tempdir()
  write_tiny_png(file.path(tmp, "a.png"), matrix(runif(64L), 8L, 8L))
  write_tiny_png(file.path(tmp, "b.png"), matrix(runif(100L), 10L, 10L))
  expect_error(read_cis(tmp), "Image dimensions differ")
})

test_that("load_signal_matrix composes read_cis + extract_signal", {
  skip_if_not_installed("png")
  tmp <- withr::local_tempdir()
  side <- 8L
  base <- matrix(0.5, side, side)
  write_tiny_png(file.path(tmp, "base.png"), base)  # base in its own dir OK
  base_path <- file.path(tmp, "base.png")
  ci_dir <- file.path(tmp, "cis")
  dir.create(ci_dir)
  for (i in 1:3) {
    noise <- matrix(runif(side * side, 0, 0.1), side, side)
    write_tiny_png(file.path(ci_dir, sprintf("p%d.png", i)), base + noise)
  }
  sig <- load_signal_matrix(ci_dir, base_image_path = base_path)
  expect_equal(dim(sig), c(side * side, 3L))
  # every pixel should be approximately noise (0 to ~0.1), base subtracted
  expect_true(all(sig >= -0.01 & sig <= 0.12))
})
