test_that("run_between wraps cluster_test + dissimilarity in a report", {
  pair <- make_sig_pair(256L, 10L, seed = 1L)
  rep <- suppressWarnings(
    run_between(pair$a, pair$b, img_dims = pair$img_dims,
                n_permutations = 50L, n_boot = 100L,
                seed = 1L, progress = FALSE)
  )
  expect_s3_class(rep, "rcicrely_report")
  expect_setequal(names(rep$results),
                  c("cluster_test", "dissimilarity"))
  expect_s3_class(rep$results$cluster_test,  "rcicrely_cluster_test")
  expect_s3_class(rep$results$dissimilarity, "rcicrely_dissim")
  expect_equal(rep$method, "between")
  expect_equal(rep$img_dims, pair$img_dims)
})
