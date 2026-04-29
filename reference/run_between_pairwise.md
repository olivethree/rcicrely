# Run all pairwise between-condition comparisons across K conditions

Generalises
[`run_between()`](https://olivethree.github.io/rcicrely/reference/run_between.md)
from a 2-condition comparison to K conditions: runs
[`rel_cluster_test()`](https://olivethree.github.io/rcicrely/reference/rel_cluster_test.md)
and
[`rel_dissimilarity()`](https://olivethree.github.io/rcicrely/reference/rel_dissimilarity.md)
on every K-choose-2 pair and applies a family-wise error correction
across pairs.

## Usage

``` r
run_between_pairwise(
  signal_matrices,
  fwer = c("holm", "bonferroni", "none"),
  img_dims = NULL,
  paired = FALSE,
  method = c("threshold", "tfce"),
  n_permutations = 2000L,
  n_boot = 2000L,
  cluster_threshold = 2,
  alpha = 0.05,
  ci_level = 0.95,
  dissim_null = c("none", "permutation"),
  mask = NULL,
  seed = NULL,
  progress = TRUE,
  acknowledge_scaling = FALSE
)
```

## Arguments

- signal_matrices:

  Named list of pixels x participants signal matrices, one per
  condition. Names become condition labels in the output. Must have at
  least 2 elements.

- fwer:

  One of `"holm"` (default), `"bonferroni"`, or `"none"`. The
  across-pairs family-wise correction applied to the per-pair minimum
  cluster *p*-value.

- img_dims:

  Integer `c(nrow, ncol)`. If `NULL`, inferred from the `img_dims`
  attribute on the first matrix.

- paired:

  Logical. When `TRUE`, all pairs use the paired variant. All matrices
  must have identical column names; rows correspond to the same producer
  across matrices.

- method:

  Cluster-test method. Forwarded to
  [`rel_cluster_test()`](https://olivethree.github.io/rcicrely/reference/rel_cluster_test.md).
  Default `"threshold"`.

- n_permutations, n_boot, cluster_threshold, alpha, ci_level:

  Forwarded to the per-pair
  [`rel_cluster_test()`](https://olivethree.github.io/rcicrely/reference/rel_cluster_test.md)
  and
  [`rel_dissimilarity()`](https://olivethree.github.io/rcicrely/reference/rel_dissimilarity.md)
  calls.

- dissim_null:

  Forwarded to
  [`rel_dissimilarity()`](https://olivethree.github.io/rcicrely/reference/rel_dissimilarity.md)
  as `null`. Default `"none"` (skip the per-pair dissimilarity null to
  keep wall time bounded).

- mask:

  Optional logical vector forwarded to both per-pair calls.

- seed:

  Optional integer.

- progress:

  Show `cli` progress bars per per-pair call.

- acknowledge_scaling:

  Logical. Forwarded.

## Value

Object of class `rcicrely_pairwise_report` with:

- `$pairs`, data.frame with one row per K-choose-2 pair: columns
  `pair_id`, `cond_a`, `cond_b`, `n_clusters`, `p_min` (minimum
  within-pair cluster *p*-value), `p_adj_pair` (FWER-adjusted across
  pairs under `fwer`), `significant` (logical: `p_adj_pair < alpha`),
  `euclidean`, `euclidean_normalised`.

- `$results`, named list of per-pair `rel_cluster_test` and
  `rel_dissimilarity` results (keyed by `pair_id`).

- `$conditions`, the input names.

- `$fwer`, the correction method used.

- `$alpha`, the across-pairs alpha.

## FWER scope

Cluster-level *p*-values within each pair are already max-statistic
FWER-controlled by
[`rel_cluster_test()`](https://olivethree.github.io/rcicrely/reference/rel_cluster_test.md).
The Holm/Bonferroni adjustment in `run_between_pairwise()` controls
family-wise error across the K-choose-2 *pair comparisons* (a second
layer above the cluster test's internal control), not over individual
pixels or clusters.

For each pair, the statistic carried into the across-pairs adjustment is
the **minimum cluster-level *p*-value within that pair** (i.e., the
most-significant cluster's *p*). Within-pair cluster *p*-values are not
re-adjusted: they remain the max-statistic FWER-controlled values from
[`rel_cluster_test()`](https://olivethree.github.io/rcicrely/reference/rel_cluster_test.md).
A pair with no clusters at all contributes `p_min = 1.0` so the Holm
ordering is well-defined.

## See also

[`run_between()`](https://olivethree.github.io/rcicrely/reference/run_between.md),
[`rel_cluster_test()`](https://olivethree.github.io/rcicrely/reference/rel_cluster_test.md),
[`rel_dissimilarity()`](https://olivethree.github.io/rcicrely/reference/rel_dissimilarity.md)
