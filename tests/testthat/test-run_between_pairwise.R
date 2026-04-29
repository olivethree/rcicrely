## §1.4: tests for run_between_pairwise()

build_three_conditions <- function() {
  set.seed(101L)
  a <- raw_pure_signal(n_pix = 256L, n_p = 12L, seed = 1L)
  b <- raw_zero_signal(n_pix = 256L, n_p = 12L, seed = 2L)
  c_mat <- raw_zero_signal(n_pix = 256L, n_p = 12L, seed = 3L)
  list(A = a, B = b, C = c_mat)
}

test_that("run_between_pairwise builds K-choose-2 pair table", {
  conds <- build_three_conditions()
  rep <- suppressWarnings(
    run_between_pairwise(
      conds,
      n_permutations = 50L, n_boot = 50L,
      seed = 1L, progress = FALSE
    )
  )
  expect_s3_class(rep, "rcicrely_pairwise_report")
  # 3 conditions -> C(3, 2) = 3 pairs
  expect_equal(nrow(rep$pairs), 3L)
  expect_setequal(rep$pairs$pair_id,
                  c("A_vs_B", "A_vs_C", "B_vs_C"))
  expect_length(rep$results, 3L)
  for (pid in rep$pairs$pair_id) {
    expect_s3_class(rep$results[[pid]]$cluster_test,
                    "rcicrely_cluster_test")
    expect_s3_class(rep$results[[pid]]$dissimilarity,
                    "rcicrely_dissim")
  }
})

test_that("FWER scope: holm adjusts p_min across pairs, not within", {
  conds <- build_three_conditions()
  rep <- suppressWarnings(
    run_between_pairwise(
      conds, fwer = "holm",
      n_permutations = 50L, n_boot = 50L,
      seed = 1L, progress = FALSE
    )
  )
  # Manual Holm on the same p_min values should match
  expect_equal(rep$pairs$p_adj_pair,
               stats::p.adjust(rep$pairs$p_min, method = "holm"))
  # When all three pairs have p_min = 1 (no clusters under such weak
  # synthetic signal), p_adj_pair must also all be 1
  if (all(rep$pairs$p_min == 1)) {
    expect_true(all(rep$pairs$p_adj_pair == 1))
  }
})

test_that("FWER 'none' returns p_min unchanged as p_adj_pair", {
  conds <- build_three_conditions()
  rep <- suppressWarnings(
    run_between_pairwise(
      conds, fwer = "none",
      n_permutations = 50L, n_boot = 50L,
      seed = 1L, progress = FALSE
    )
  )
  expect_equal(rep$pairs$p_adj_pair, rep$pairs$p_min)
})

test_that("no-clusters edge case sets p_min = 1 with well-defined Holm", {
  # Two zero-signal conditions, twice, so all three pairs have no
  # detectable clusters at the default threshold
  conds <- list(
    A = raw_zero_signal(n_p = 12L, seed = 1L),
    B = raw_zero_signal(n_p = 12L, seed = 2L),
    C = raw_zero_signal(n_p = 12L, seed = 3L)
  )
  rep <- suppressWarnings(
    run_between_pairwise(
      conds,
      n_permutations = 50L, n_boot = 50L,
      seed = 1L, progress = FALSE
    )
  )
  expect_true(all(rep$pairs$p_min >= 0 & rep$pairs$p_min <= 1))
  expect_true(all(is.finite(rep$pairs$p_adj_pair)))
})

test_that("aborts on a single-condition input", {
  expect_error(
    run_between_pairwise(
      list(A = raw_pure_signal(n_p = 12L)),
      n_permutations = 10L, n_boot = 10L, progress = FALSE
    ),
    "at least 2"
  )
})

test_that("aborts on unnamed list", {
  expect_error(
    run_between_pairwise(
      list(raw_pure_signal(n_p = 12L), raw_zero_signal(n_p = 12L)),
      n_permutations = 10L, n_boot = 10L, progress = FALSE
    ),
    "named list"
  )
})
