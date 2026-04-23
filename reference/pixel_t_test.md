# Vectorised pixel-wise Welch t-test

Computes Welch's t (unequal variances, Student's is **not** appropriate
– conditions can differ in both N and variance) per pixel between two
condition signal matrices. Fully vectorised: no per-row
[`apply()`](https://rdrr.io/r/base/apply.html). Returns a
length-`n_pixels` numeric vector.

Use this when you want the raw t-map (e.g. for custom plotting or
threshold sensitivity analysis) without the cluster-level inference
machinery of
[`rel_cluster_test()`](https://olivethree.github.io/rcicrely/reference/rel_cluster_test.md).

## Usage

``` r
pixel_t_test(signal_matrix_a, signal_matrix_b)
```

## Arguments

- signal_matrix_a, signal_matrix_b:

  Pixels x participants, base-subtracted. Row counts must match.

## Value

Numeric vector of length `nrow(signal_matrix_a)`.

## Reading the result

Numeric vector, one t-value per pixel. Pixels with zero variance in both
conditions get `0` rather than `NaN` so cluster utilities don't have to
special-case them.

## Reliability metrics expect raw masks

Welch t is variance-based and sensitive to scaling. See
[`vignette("tutorial", package = "rcicrely")`](https://olivethree.github.io/rcicrely/articles/tutorial.md)
chapter 3.

## See also

[`rel_cluster_test()`](https://olivethree.github.io/rcicrely/reference/rel_cluster_test.md)
