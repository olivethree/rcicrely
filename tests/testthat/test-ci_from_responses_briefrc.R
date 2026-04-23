test_that("Brief-RC mask equals Schmitz genMask by hand", {
  skip_if_not_installed("png")
  set.seed(1)
  n_pix <- 16L * 16L
  n_pool <- 100L
  noise_matrix <- matrix(rnorm(n_pix * n_pool), n_pix, n_pool)

  base_path <- tempfile(fileext = ".png")
  png::writePNG(matrix(0.5, 16L, 16L), base_path)

  responses <- data.frame(
    participant_id = rep(c("p1", "p2"), each = 30L),
    trial          = rep(1:30, 2L),
    stimulus       = sample.int(n_pool, 60L, replace = TRUE),
    response       = sample(c(-1L, 1L), 60L, replace = TRUE)
  )

  res <- ci_from_responses_briefrc(
    responses       = responses,
    noise_matrix    = noise_matrix,
    base_image_path = base_path
  )
  expect_equal(dim(res$signal_matrix), c(n_pix, 2L))

  # hand computation for p1
  sub <- responses[responses$participant_id == "p1", ]
  X <- data.table::data.table(response = sub$response, stim = sub$stimulus)
  X <- X[, list(response = mean(response)), by = "stim"]
  hand_mask <- as.numeric(noise_matrix[, X$stim] %*% X$response) / nrow(X)
  expect_equal(res$signal_matrix[, "p1"], hand_mask, tolerance = 1e-10)
})

test_that("duplicate-stimulus responses average (genMask rule)", {
  skip_if_not_installed("png")
  n_pix <- 64L
  noise_matrix <- matrix(c(rep(1, n_pix), rep(2, n_pix), rep(3, n_pix)),
                         nrow = n_pix, ncol = 3L)
  base_path <- tempfile(fileext = ".png")
  png::writePNG(matrix(0.5, 8L, 8L), base_path)

  # p1 chose stim=1 twice with opposite responses (should cancel)
  # and stim=2 once with response=1 (unique).
  # Expected: mean(response) by stim -> stim=1 weight 0, stim=2 weight 1
  # mask = noise[,1]*0 + noise[,2]*1, divided by 2 unique stims
  responses <- data.frame(
    participant_id = "p1",
    trial          = 1:3,
    stimulus       = c(1L, 1L, 2L),
    response       = c(1L, -1L, 1L)
  )
  res <- ci_from_responses_briefrc(
    responses       = responses,
    noise_matrix    = noise_matrix,
    base_image_path = base_path
  )
  expected <- (noise_matrix[, 1] * 0 + noise_matrix[, 2] * 1) / 2
  expect_equal(as.numeric(res$signal_matrix[, "p1"]), expected,
               tolerance = 1e-10)
})

test_that("Brief-RC rejects bad response coding", {
  skip_if_not_installed("png")
  base_path <- tempfile(fileext = ".png")
  png::writePNG(matrix(0.5, 4L, 4L), base_path)
  noise_matrix <- matrix(0, nrow = 16L, ncol = 5L)
  responses <- data.frame(
    participant_id = "p1", trial = 1:3,
    stimulus = c(1L, 2L, 3L),
    response = c(0L, 1L, -1L)   # 0 is invalid
  )
  expect_error(
    ci_from_responses_briefrc(responses, noise_matrix = noise_matrix,
                              base_image_path = base_path),
    "must contain only values"
  )
})

test_that("briefrc20 method aborts with a clear message", {
  skip_if_not_installed("png")
  base_path <- tempfile(fileext = ".png")
  png::writePNG(matrix(0.5, 4L, 4L), base_path)
  responses <- data.frame(
    participant_id = "p1", trial = 1L, stimulus = 1L, response = 1L
  )
  expect_error(
    ci_from_responses_briefrc(
      responses, noise_matrix = matrix(0, 16L, 2L),
      base_image_path = base_path, method = "briefrc20"
    ),
    "not supported"
  )
})

test_that("default scaling = 'none' returns no rendered_ci field", {
  skip_if_not_installed("png")
  set.seed(2)
  n_pix <- 8L * 8L
  noise_matrix <- matrix(rnorm(n_pix * 20L), n_pix, 20L)
  base_path <- tempfile(fileext = ".png")
  png::writePNG(matrix(0.5, 8L, 8L), base_path)
  responses <- data.frame(
    participant_id = rep("p1", 10L),
    trial          = 1:10,
    stimulus       = sample.int(20L, 10L),
    response       = sample(c(-1L, 1L), 10L, replace = TRUE)
  )
  res <- ci_from_responses_briefrc(
    responses       = responses,
    noise_matrix    = noise_matrix,
    base_image_path = base_path
  )
  expect_null(res$rendered_ci)
  expect_equal(res$scaling, "none")
})

test_that("scaling = 'matched' adds rendered_ci, signal_matrix unchanged", {
  skip_if_not_installed("png")
  set.seed(3)
  n_pix <- 8L * 8L
  noise_matrix <- matrix(rnorm(n_pix * 20L), n_pix, 20L)
  # base must have within-image variance for "matched" to do anything
  base_img <- matrix(seq(0, 1, length.out = n_pix), 8L, 8L)
  base_path <- tempfile(fileext = ".png")
  png::writePNG(base_img, base_path)
  responses <- data.frame(
    participant_id = rep(c("p1", "p2"), each = 8L),
    trial          = c(1:8, 1:8),
    stimulus       = sample.int(20L, 16L, replace = TRUE),
    response       = sample(c(-1L, 1L), 16L, replace = TRUE)
  )
  res_none <- ci_from_responses_briefrc(
    responses, noise_matrix = noise_matrix,
    base_image_path = base_path, scaling = "none"
  )
  res_match <- ci_from_responses_briefrc(
    responses, noise_matrix = noise_matrix,
    base_image_path = base_path, scaling = "matched"
  )
  expect_equal(res_match$signal_matrix, res_none$signal_matrix)
  expect_false(is.null(res_match$rendered_ci))
  expect_equal(dim(res_match$rendered_ci), dim(res_match$signal_matrix))
  # matched stretches mask range to base range, then adds base.
  # Per-column rendered range should be at least the base range
  # (since rendered = base + matched(mask), and matched(mask) has
  # the same range as base).
  base_rng <- diff(range(base_img))
  rendered_rngs <- vapply(seq_len(ncol(res_match$rendered_ci)),
                          function(j) diff(range(res_match$rendered_ci[, j])),
                          numeric(1L))
  expect_true(all(rendered_rngs >= base_rng - 1e-3))
})

test_that("scaling = 'constant' applies the requested multiplier", {
  skip_if_not_installed("png")
  set.seed(4)
  n_pix <- 8L * 8L
  noise_matrix <- matrix(rnorm(n_pix * 20L), n_pix, 20L)
  # Use a base value that round-trips PNG losslessly: png stores
  # 8-bit integers so 0.5 becomes 127/255 != 0.5.
  base_img <- matrix(127 / 255, 8L, 8L)
  base_path <- tempfile(fileext = ".png")
  png::writePNG(base_img, base_path)
  responses <- data.frame(
    participant_id = rep("p1", 10L),
    trial          = 1:10,
    stimulus       = sample.int(20L, 10L),
    response       = sample(c(-1L, 1L), 10L, replace = TRUE)
  )
  res <- ci_from_responses_briefrc(
    responses, noise_matrix = noise_matrix,
    base_image_path  = base_path,
    scaling          = "constant",
    scaling_constant = 3
  )
  # rendered = base_read_back + 3 * mask
  base_read <- as.vector(png::readPNG(base_path))
  expect_equal(res$rendered_ci[, 1L],
               base_read + 3 * res$signal_matrix[, 1L],
               tolerance = 1e-6)
})

test_that("scaling = 'constant' without scaling_constant aborts", {
  skip_if_not_installed("png")
  base_path <- tempfile(fileext = ".png")
  png::writePNG(matrix(0.5, 4L, 4L), base_path)
  responses <- data.frame(
    participant_id = "p1", trial = 1L, stimulus = 1L, response = 1L
  )
  expect_error(
    ci_from_responses_briefrc(
      responses, noise_matrix = matrix(0, 16L, 2L),
      base_image_path = base_path,
      scaling = "constant"
    ),
    "scaling_constant"
  )
})
