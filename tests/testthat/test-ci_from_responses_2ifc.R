test_that("2IFC wrapper smokes against the original rcicr if available", {
  skip_if_not_installed("rcicr")
  skip_if_not_installed("jpeg")
  skip_if_not_installed("foreach")
  skip_if_not_installed("tibble")
  skip_if_not_installed("dplyr")
  # generate a tiny stimulus set on the fly. img_size has to be even.
  tmp <- withr::local_tempdir()
  base_path <- file.path(tmp, "base.jpg")
  jpeg::writeJPEG(matrix(0.5, 32L, 32L), base_path, quality = 1)

  withr::with_dir(tmp, {
    set.seed(1)
    suppressMessages(suppressWarnings(
      rcicr::generateStimuli2IFC(
        base_face_files     = list(base = base_path),
        n_trials            = 40L,
        img_size            = 32L,
        stimulus_path       = file.path(tmp, "stim"),
        label               = "test",
        ncores              = 1L,
        return_as_dataframe = TRUE,
        save_as_png         = FALSE,
        seed                = 1L
      )
    ))
  })

  rdata_path <- list.files(
    file.path(tmp, "stim"),
    pattern = "\\.rdata$", ignore.case = TRUE, full.names = TRUE
  )[1]
  skip_if(is.na(rdata_path) || !file.exists(rdata_path),
          "rcicr did not produce rdata (environment-specific issue)")

  responses <- data.frame(
    participant_id = rep(c("p1", "p2", "p3", "p4"), each = 20L),
    stimulus       = rep(1:40, 2L),
    response       = sample(c(-1L, 1L), 80L, replace = TRUE)
  )

  res <- ci_from_responses_2ifc(
    responses       = responses,
    rdata_path      = rdata_path,
    baseimage       = "base",
    save_as_png     = FALSE,
    targetpath      = file.path(tmp, "cis")
  )
  expect_equal(dim(res$signal_matrix), c(32L * 32L, 4L))
  expect_setequal(res$participants, c("p1", "p2", "p3", "p4"))
})

test_that("2IFC wrapper rejects bad response coding", {
  skip_if_not_installed("rcicr")
  # Just check the pre-rcicr validation path; no need for real rdata.
  tmp <- tempfile(fileext = ".Rdata")
  base_faces <- list(base = matrix(0.5, 8L, 8L))
  save(base_faces, file = tmp)
  responses <- data.frame(
    participant_id = "p1",
    stimulus       = 1L,
    response       = 5L    # not in {-1, 1}
  )
  expect_error(
    ci_from_responses_2ifc(responses = responses, rdata_path = tmp,
                           baseimage = "base"),
    "must contain only values"
  )
})
