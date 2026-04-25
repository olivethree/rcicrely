# Run every within-condition reliability metric

Convenience orchestrator that runs the reliability metrics proper
(split-half with Spearman-Brown correction and ICC) on a single
condition's signal matrix and wraps the two results in an
`rcicrely_report` for joint printing and plotting.

Use this when you want the full within-condition reliability report in
one call. Pass each metric individually if you need to tune arguments
per metric.

## Usage

``` r
run_within(
  signal_matrix,
  n_permutations = 2000L,
  icc_variants = c("3_1", "3_k"),
  mask = NULL,
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

- icc_variants:

  Passed to
  [`rel_icc()`](https://olivethree.github.io/rcicrely/reference/rel_icc.md).

- mask:

  Optional logical vector of length `nrow(signal_matrix)` restricting
  all metrics to a region (e.g., from
  [`face_mask()`](https://olivethree.github.io/rcicrely/reference/face_mask.md)).
  Threaded through to both
  [`rel_split_half()`](https://olivethree.github.io/rcicrely/reference/rel_split_half.md)
  and
  [`rel_icc()`](https://olivethree.github.io/rcicrely/reference/rel_icc.md).

- seed:

  Optional integer; used for the split-half permutations.

- progress:

  Show `cli` progress bars.

## Value

Object of class `rcicrely_report` with `$results` = named list of two
`rcicrely_*` objects (`split_half`, `icc`) and `$method = "within"`.

## What is included (and what is not)

`run_within()` returns the two metrics that quantify the reliability of
the group-level classification image proper: **split-half** (a
permutation-based estimate of group-CI stability with Spearman-Brown
projection to the full sample) and **ICC(3,\*)** (the psychometric
variance decomposition). These are non-redundant — split-half is
resolution-stable and assumption-light; ICC decomposes variance and
reports both single-producer and average-measures reliability.

Leave-one-out influence screening lives in
[`rel_loo()`](https://olivethree.github.io/rcicrely/reference/rel_loo.md)
and is **not** bundled here. Its output is an influence diagnostic, not
a reliability statistic, and mixing it into the reliability report
invites mis-reading `r_loo` values (which are near 1 by construction) as
reliability.

## Reading the result

`$results$split_half`, `$results$icc`, one result object each, with the
same fields as the standalone functions. `$method = "within"`.

## Reliability metrics expect raw masks

Both downstream metrics expect the raw mask. See
[`vignette("tutorial", package = "rcicrely")`](https://olivethree.github.io/rcicrely/articles/tutorial.md)
chapter 3.

## See also

[`rel_split_half()`](https://olivethree.github.io/rcicrely/reference/rel_split_half.md),
[`rel_icc()`](https://olivethree.github.io/rcicrely/reference/rel_icc.md),
[`rel_loo()`](https://olivethree.github.io/rcicrely/reference/rel_loo.md)
for the influence diagnostic,
[`run_between()`](https://olivethree.github.io/rcicrely/reference/run_between.md).
