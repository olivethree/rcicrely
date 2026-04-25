# Compare Euclidean dissimilarity intervals across multiple contrasts on a single plot

Side-by-side comparison of
[`rel_dissimilarity()`](https://olivethree.github.io/rcicrely/reference/rel_dissimilarity.md)
results: one point per contrast, with its 95% percentile bootstrap CI as
a horizontal range. Useful for paper figures showing whether two
contrasts (e.g., Trust vs Friendly, Dominant vs Competent) have
overlapping or non-overlapping CIs without forcing the reader to read
four numbers from a table.

## Usage

``` r
compare_dissimilarity(
  ...,
  metric = c("euclidean", "euclidean_normalised"),
  main = "Between-condition Euclidean distance"
)
```

## Arguments

- ...:

  Named `rcicrely_dissim` objects (use `name = obj` syntax). Names
  become the y-axis labels.

- metric:

  `"euclidean"` (default) or `"euclidean_normalised"`.

- main:

  Plot title.

## Value

Invisibly a data frame of (label, observed, ci_low, ci_high) for further
use.

## See also

[`rel_dissimilarity()`](https://olivethree.github.io/rcicrely/reference/rel_dissimilarity.md)

## Examples

``` r
if (FALSE) { # \dontrun{
d_AB <- rel_dissimilarity(sig_A, sig_B, n_boot = 500L, seed = 1L)
d_AC <- rel_dissimilarity(sig_A, sig_C, n_boot = 500L, seed = 1L)
compare_dissimilarity("A vs B" = d_AB, "A vs C" = d_AC)
} # }
```
