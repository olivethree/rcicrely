## These tests just call every S3 method to make sure they don't error.
## We don't assert on plot appearance.

test_that("print/summary/plot methods exist and do not error", {
  sig <- make_sig(64L, 10L)
  r <- suppressWarnings(
    rel_split_half(sig, n_permutations = 50L, progress = FALSE)
  )
  expect_invisible(print(r))
  expect_invisible(summary(r))
  expect_silent({ grDevices::pdf(NULL); plot(r); grDevices::dev.off() })

  lo <- suppressWarnings(rel_loo(sig))
  expect_invisible(print(lo))
  expect_invisible(summary(lo))
  expect_silent({ grDevices::pdf(NULL); plot(lo); grDevices::dev.off() })

  ic <- suppressWarnings(rel_icc(sig))
  expect_invisible(print(ic))
  expect_invisible(summary(ic))
  expect_silent({ grDevices::pdf(NULL); plot(ic); grDevices::dev.off() })

  pair <- make_sig_pair(64L, 8L)
  ct <- suppressWarnings(
    rel_cluster_test(pair$a, pair$b, img_dims = pair$img_dims,
                     n_permutations = 30L, progress = FALSE)
  )
  expect_invisible(print(ct))
  expect_invisible(summary(ct))
  expect_silent({ grDevices::pdf(NULL); plot(ct); grDevices::dev.off() })

  d <- suppressWarnings(
    rel_dissimilarity(pair$a, pair$b, n_boot = 50L, progress = FALSE)
  )
  expect_invisible(print(d))
  expect_invisible(summary(d))
  expect_silent({ grDevices::pdf(NULL); plot(d); grDevices::dev.off() })

  rep <- suppressWarnings(
    run_within(sig, n_permutations = 30L, progress = FALSE)
  )
  expect_invisible(print(rep))
  expect_invisible(summary(rep))
  expect_silent({ grDevices::pdf(NULL); plot(rep); grDevices::dev.off() })
})
