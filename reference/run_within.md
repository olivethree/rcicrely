# Run every within-condition reliability metric

Convenience orchestrator that runs
[`rel_split_half()`](https://olivethree.github.io/rcicrely/reference/rel_split_half.md),
[`rel_loo()`](https://olivethree.github.io/rcicrely/reference/rel_loo.md),
and
[`rel_icc()`](https://olivethree.github.io/rcicrely/reference/rel_icc.md)
on a single condition's signal matrix and wraps the three results in an
`rcicrely_report` for joint printing / plotting.

Use this when you want the full within-condition reliability report in
one call. Pass each metric individually if you need to tune arguments
per metric.

## Usage

``` r
run_within(
  signal_matrix,
  n_permutations = 2000L,
  flag_threshold = 2.5,
  flag_method = c("sd", "mad"),
  flag_threshold_sd = NULL,
  icc_variants = c("3_1", "3_k"),
  seed = NULL,
  progress = TRUE
)
```

## Arguments

- signal_matrix:

  Pixels x participants, base-subtracted.

- n_permutations:

  Passed to
  [`rel_split_half()`](https://olivethree.github.io/rcicrely/reference/rel_split_half.md).
  Default 2000.

- flag_threshold:

  Passed to
  [`rel_loo()`](https://olivethree.github.io/rcicrely/reference/rel_loo.md).
  Default 2.5.

- flag_method:

  Passed to
  [`rel_loo()`](https://olivethree.github.io/rcicrely/reference/rel_loo.md).
  Default `"sd"`.

- flag_threshold_sd:

  Deprecated alias for `flag_threshold`.

- icc_variants:

  Passed to
  [`rel_icc()`](https://olivethree.github.io/rcicrely/reference/rel_icc.md).

- seed:

  Optional integer; used for the split-half permutations.

- progress:

  Show `cli` progress bars.

## Value

Object of class `rcicrely_report` with `$results` = named list of three
`rcicrely_*` objects (`split_half`, `loo`, `icc`) and
`$method = "within"`.

## Reading the result

`$results$split_half`, `$results$loo`, `$results$icc`, one result object
each, with the same fields as the standalone functions.
`$method = "within"`.

## Reliability metrics expect raw masks

All three downstream metrics expect the raw mask. See
[`vignette("tutorial", package = "rcicrely")`](https://olivethree.github.io/rcicrely/articles/tutorial.md)
chapter 3.

## See also

[`rel_split_half()`](https://olivethree.github.io/rcicrely/reference/rel_split_half.md),
[`rel_loo()`](https://olivethree.github.io/rcicrely/reference/rel_loo.md),
[`rel_icc()`](https://olivethree.github.io/rcicrely/reference/rel_icc.md),
[`run_between()`](https://olivethree.github.io/rcicrely/reference/run_between.md)
