## Tests target the non-exported label_components_4conn and
## find_clusters helpers. Use `:::` to reach them.

test_that("4-connectivity labels a simple connected blob as one component", {
  m <- matrix(FALSE, 5L, 5L)
  m[2:4, 2:4] <- TRUE
  labels <- rcicrely:::label_components_4conn(m)
  expect_equal(length(unique(as.vector(labels[labels != 0L]))), 1L)
  expect_equal(sum(labels != 0L), 9L)
})

test_that("diagonal neighbours are NOT connected under 4-conn", {
  m <- matrix(FALSE, 4L, 4L)
  m[1L, 1L] <- TRUE
  m[2L, 2L] <- TRUE   # diagonal
  m[3L, 3L] <- TRUE
  labels <- rcicrely:::label_components_4conn(m)
  expect_equal(length(unique(as.vector(labels[labels != 0L]))), 3L)
})

test_that("find_clusters returns matching pos and neg masses", {
  # construct a small t-map with one obvious pos and one obvious neg cluster
  tmat <- matrix(0, 10L, 10L)
  tmat[2:4, 2:4] <-  3.0    # pos cluster, mass = 27
  tmat[7:9, 7:9] <- -3.5    # neg cluster, mass = -31.5
  out <- rcicrely:::find_clusters(as.vector(tmat), c(10L, 10L), threshold = 2)
  expect_equal(length(out$pos_masses), 1L)
  expect_equal(length(out$neg_masses), 1L)
  expect_equal(out$pos_masses, 27)
  expect_equal(out$neg_masses, -31.5)
})

test_that("find_clusters handles an all-below-threshold map", {
  tmat <- matrix(0.5, 8L, 8L)
  out <- rcicrely:::find_clusters(as.vector(tmat), c(8L, 8L), threshold = 2)
  expect_equal(length(out$pos_masses), 0L)
  expect_equal(length(out$neg_masses), 0L)
})
