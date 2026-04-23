test_that("read_noise_matrix roundtrips an .Rdata containing `stimuli`", {
  tmp <- tempfile(fileext = ".Rdata")
  stimuli <- matrix(rnorm(200L), nrow = 50L, ncol = 4L)
  save(stimuli, file = tmp)
  nm <- read_noise_matrix(tmp)
  expect_true(is.matrix(nm))
  expect_equal(dim(nm), c(50L, 4L))
  expect_equal(nm, stimuli, tolerance = 1e-10)
})

test_that("read_noise_matrix accepts whitespace-delimited text files", {
  tmp <- tempfile(fileext = ".txt")
  m <- matrix(seq_len(12L), nrow = 4L, ncol = 3L)
  utils::write.table(m, tmp, row.names = FALSE, col.names = FALSE)
  nm <- read_noise_matrix(tmp)
  expect_equal(dim(nm), c(4L, 3L))
  expect_equal(nm, m, tolerance = 1e-10)
})

test_that("read_noise_matrix errors when the named object is missing", {
  tmp <- tempfile(fileext = ".Rdata")
  other_name <- matrix(1:10, 5L, 2L)
  save(other_name, file = tmp)
  expect_error(read_noise_matrix(tmp), "not found")
})
