test_that("extract_signal subtracts the base image column-wise", {
  skip_if_not_installed("png")
  side <- 8L
  base <- matrix(0.4, side, side)
  base_path <- tempfile(fileext = ".png")
  png::writePNG(base, base_path)

  cis <- cbind(
    as.vector(base + 0.1),
    as.vector(base + 0.2),
    as.vector(base)
  )
  sig <- suppressWarnings(extract_signal(cis, base_image_path = base_path))
  expect_equal(sig[, 1], rep(0.1, side * side), tolerance = 1e-6)
  expect_equal(sig[, 2], rep(0.2, side * side), tolerance = 1e-6)
  expect_equal(sig[, 3], rep(0.0, side * side), tolerance = 1e-6)
})

test_that("extract_signal errors on pixel-count mismatch", {
  skip_if_not_installed("png")
  base_path <- tempfile(fileext = ".png")
  png::writePNG(matrix(0.5, 4L, 4L), base_path)
  cis <- matrix(0.5, nrow = 100L, ncol = 2L)   # wrong pixel count
  expect_error(
    extract_signal(cis, base_image_path = base_path),
    "does not match"
  )
})

test_that("extract_signal emits the Mode-1 scaling warning", {
  skip_if_not_installed("png")
  side <- 8L
  base_path <- tempfile(fileext = ".png")
  png::writePNG(matrix(0.4, side, side), base_path)
  cis <- matrix(0.5, nrow = side * side, ncol = 2L)
  rcicrely:::reset_session_warnings()
  expect_warning(extract_signal(cis, base_image_path = base_path),
                 "rendered")
})

test_that("acknowledge_scaling = TRUE silences the warning", {
  skip_if_not_installed("png")
  side <- 8L
  base_path <- tempfile(fileext = ".png")
  png::writePNG(matrix(0.4, side, side), base_path)
  cis <- matrix(0.5, nrow = side * side, ncol = 2L)
  rcicrely:::reset_session_warnings()
  expect_no_warning(extract_signal(cis, base_image_path = base_path,
                                   acknowledge_scaling = TRUE))
})
