# Vectorised pixel-wise t-test (independent or paired)

At every pixel, tests whether condition A's mean signal differs from
condition B's. Two modes:

- `paired = FALSE` (default): independent-samples Welch t per pixel.
  Correct when producers in A and B are different people
  (between-subjects design).

- `paired = TRUE`: paired t per pixel on the per-producer difference
  `A - B`. Correct when the same producers contributed to both
  conditions (within-subjects design). Paired t is strictly more
  powerful than Welch t whenever producer-level variance is present.

Returns a numeric vector of t-values, length `n_pixels`. Used as an
intermediate by
[`rel_cluster_test()`](https://olivethree.github.io/rcicrely/reference/rel_cluster_test.md);
not intended as a standalone inferential test.

## Usage

``` r
pixel_t_test(signal_matrix_a, signal_matrix_b, paired = FALSE, mask = NULL)
```

## Arguments

- signal_matrix_a, signal_matrix_b:

  Pixels x participants, base-subtracted. Row counts must match. When
  `paired = TRUE` the column counts must also match, and column names
  must correspond to the same producer across matrices.

- paired:

  Logical. `FALSE` (default) uses independent Welch t; `TRUE` uses
  paired t.

- mask:

  Optional logical vector of length `nrow(signal_matrix_a)` restricting
  the per-pixel t computation to a region. **Both** matrices are
  subsetted with the same mask before computing t — see
  [`vignette("tutorial")`](https://olivethree.github.io/rcicrely/articles/tutorial.md)
  §6.6 for the apply-symmetrically rule. The returned vector is then of
  length `sum(mask)`, not `nrow(signal_matrix_a)`.

## Value

Numeric vector of length `nrow(signal_matrix_a)` (or `sum(mask)` if
`mask` is supplied).

## Reading the result

Numeric vector, one t-value per pixel. Pixels with zero variance get `0`
rather than `NaN` so cluster utilities don't have to special-case them.

## Reliability metrics expect raw masks

Welch t and paired t are variance-based and sensitive to scaling. See
[`vignette("tutorial", package = "rcicrely")`](https://olivethree.github.io/rcicrely/articles/tutorial.md)
chapter 3.

## See also

[`rel_cluster_test()`](https://olivethree.github.io/rcicrely/reference/rel_cluster_test.md)
