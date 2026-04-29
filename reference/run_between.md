# Run every between-condition discriminability metric

Convenience orchestrator that runs
[`rel_cluster_test()`](https://olivethree.github.io/rcicrely/reference/rel_cluster_test.md)
and
[`rel_dissimilarity()`](https://olivethree.github.io/rcicrely/reference/rel_dissimilarity.md)
on two condition signal matrices and wraps both results in an
`rcicrely_report` for joint printing / plotting.

Use this when you want both the spatial-pattern test and the overall
magnitude test in one call.

## Usage

``` r
run_between(
  signal_matrix_a,
  signal_matrix_b,
  img_dims = NULL,
  n_permutations = 2000L,
  n_boot = 2000L,
  cluster_threshold = 2,
  alpha = 0.05,
  ci_level = 0.95,
  mask = NULL,
  seed = NULL,
  progress = TRUE,
  acknowledge_scaling = FALSE
)
```

## Arguments

- signal_matrix_a, signal_matrix_b:

  Pixels x participants, base-subtracted. Row counts must match.

- img_dims:

  Integer `c(nrow, ncol)`. If `NULL`, inferred from the `img_dims`
  attribute on `signal_matrix_a`, or from `sqrt(n_pixels)` if that's a
  whole number.

- n_permutations:

  Passed to
  [`rel_cluster_test()`](https://olivethree.github.io/rcicrely/reference/rel_cluster_test.md).
  Default 2000.

- n_boot:

  Passed to
  [`rel_dissimilarity()`](https://olivethree.github.io/rcicrely/reference/rel_dissimilarity.md).
  Default 2000.

- cluster_threshold:

  Passed to
  [`rel_cluster_test()`](https://olivethree.github.io/rcicrely/reference/rel_cluster_test.md).
  Default 2.0.

- alpha:

  Passed to
  [`rel_cluster_test()`](https://olivethree.github.io/rcicrely/reference/rel_cluster_test.md).
  Default 0.05.

- ci_level:

  Passed to
  [`rel_dissimilarity()`](https://olivethree.github.io/rcicrely/reference/rel_dissimilarity.md).
  Default 0.95.

- mask:

  Optional logical vector of length `nrow(signal_matrix_a)`. Threaded
  through to both
  [`rel_cluster_test()`](https://olivethree.github.io/rcicrely/reference/rel_cluster_test.md)
  (zero-out pattern, preserving 2D structure) and
  [`rel_dissimilarity()`](https://olivethree.github.io/rcicrely/reference/rel_dissimilarity.md)
  (drop pattern, both matrices subsetted identically). See
  [`vignette("tutorial")`](https://olivethree.github.io/rcicrely/articles/tutorial.md)
  §6.6.

- seed:

  Optional integer.

- progress:

  Show `cli` progress bars.

- acknowledge_scaling:

  Logical. Forwarded to both
  [`rel_cluster_test()`](https://olivethree.github.io/rcicrely/reference/rel_cluster_test.md)
  and
  [`rel_dissimilarity()`](https://olivethree.github.io/rcicrely/reference/rel_dissimilarity.md);
  when `FALSE` (default), known-rendered inputs error.

## Value

Object of class `rcicrely_report` with `$results` = named list of two
`rcicrely_*` objects (`cluster_test`, `dissimilarity`) and
`$method = "between"`.

## Reading the result

`$results$cluster_test` and `$results$dissimilarity`, one result object
each, fields as in the standalone functions. `$method = "between"`.

## Reliability metrics expect raw masks

Both downstream metrics are scale-sensitive: the cluster test uses
variance-based Welch t, and Euclidean distance in
[`rel_dissimilarity()`](https://olivethree.github.io/rcicrely/reference/rel_dissimilarity.md)
is sensitive to any scaling. See
[`vignette("tutorial", package = "rcicrely")`](https://olivethree.github.io/rcicrely/articles/tutorial.md)
chapter 3.

## See also

[`rel_cluster_test()`](https://olivethree.github.io/rcicrely/reference/rel_cluster_test.md),
[`rel_dissimilarity()`](https://olivethree.github.io/rcicrely/reference/rel_dissimilarity.md),
[`run_within()`](https://olivethree.github.io/rcicrely/reference/run_within.md)
