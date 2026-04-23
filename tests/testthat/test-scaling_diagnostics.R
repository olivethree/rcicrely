test_that("looks_scaled returns FALSE on raw-mask-shaped matrices", {
  set.seed(1)
  raw <- matrix(rnorm(100L * 5L, sd = 0.02), 100L, 5L)
  expect_false(rcicrely:::looks_scaled(raw))
})

test_that("looks_scaled returns TRUE on rendered-CI-shaped matrices", {
  set.seed(2)
  # rendered = base + scaling(noise); typical range ~1
  rendered <- matrix(0.5 + rnorm(100L * 5L, sd = 0.3), 100L, 5L)
  expect_true(rcicrely:::looks_scaled(rendered))
})

test_that("looks_scaled is robust to NAs and empty matrices", {
  expect_false(rcicrely:::looks_scaled(matrix(numeric(0L), 0L, 0L)))
  expect_false(rcicrely:::looks_scaled(NULL))
  expect_false(rcicrely:::looks_scaled(matrix(NA_real_, 4L, 2L)))
})

test_that("warn_looks_scaled flips session-state and respects the option", {
  # direct: call the helper, check state
  rcicrely:::reset_session_warnings()
  expect_warning(rcicrely:::warn_looks_scaled("x"), "looks like")
  state <- rcicrely:::.rcicrely_session_state
  expect_true(state$looks_scaled_warning_emitted)

  # second call: silent (once per session)
  expect_no_warning(rcicrely:::warn_looks_scaled("x"))

  # global option silences even after reset
  rcicrely:::reset_session_warnings()
  withr::with_options(
    list(rcicrely.silence_scaling_warning = TRUE),
    expect_no_warning(rcicrely:::warn_looks_scaled("x"))
  )
})
