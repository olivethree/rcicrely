# Cluster-based permutation test with max-statistic FWER control

Pixel-level discriminability test between two condition signal matrices.
Use this to identify **where** in the image two conditions' CIs differ,
with family-wise error control across pixels. Pair with
[`rel_dissimilarity()`](https://olivethree.github.io/rcicrely/reference/rel_dissimilarity.md)
for an overall magnitude. Implements the Maris & Oostenveld (2007)
approach with Nichols & Holmes (2002) max-statistic family-wise error
control:

1.  Compute observed pixel-wise Welch t.

2.  Threshold at `|t| > cluster_threshold` (separately for positive and
    negative tails); label connected components with **4-connectivity**.

3.  For each observed cluster, cluster mass = sum of t-values in the
    cluster (not pixel count).

4.  Build the null distribution via `n_permutations` **stratified**
    label permutations (each preserves `(N_A, N_B)`). For each, find
    clusters and record max positive mass and max negative mass.

5.  Observed cluster's p-value = fraction of null max-masses (matching
    sign) that equal or exceed the observed mass in absolute value.

Stratification (sec.7.5): unstratified permutation lets the two group
sizes drift, which changes Welch degrees of freedom under the null and
biases the max-mass distribution. We preserve `(N_A, N_B)` exactly.

## Usage

``` r
rel_cluster_test(
  signal_matrix_a,
  signal_matrix_b,
  img_dims = NULL,
  n_permutations = 2000L,
  cluster_threshold = 2,
  alpha = 0.05,
  seed = NULL,
  progress = TRUE
)
```

## Arguments

- signal_matrix_a, signal_matrix_b:

  Pixels x participants, base-subtracted.

- img_dims:

  Integer `c(nrow, ncol)`. Can be inferred if the signal matrices carry
  an `img_dims` attribute.

- n_permutations:

  Number of label permutations. 2000 default.

- cluster_threshold:

  Absolute t-value threshold for forming clusters. Default 2.0 (~
  two-tailed p \< 0.05 at large df).

- alpha:

  Significance level for flagging clusters. Default 0.05.

- seed:

  Optional integer; RNG state restored on exit.

- progress:

  Show a `cli` progress bar.

## Value

Object of class `rcicrely_cluster_test`. Fields described in **Reading
the result** above.

## Reading the result

- `$observed_t` – per-pixel Welch t vector for the observed condition
  labelling.

- `$clusters` – data.frame with `cluster_id`, `direction`
  (`"pos"`/`"neg"`), `mass`, `size`, `p_value`, `significant`. One row
  per supra-threshold cluster, sorted by direction then mass.

- `$null_distribution` – list with `$pos` and `$neg` vectors of
  per-permutation max masses, useful for plotting the null.

- `$pos_labels`, `$neg_labels` – integer matrices the same shape as the
  image, where each non-zero value labels a cluster. Drives the contour
  overlay in the result's
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method.

- `$cluster_threshold`, `$alpha`, `$n_permutations`,
  `$n_participants_a`, `$n_participants_b` – metadata.

## Common mistakes

- Reading `cluster_threshold` (default 2.0) as a p-value cutoff. It is
  the t-value above which pixels join a cluster; significance is decided
  afterwards by comparing cluster mass to the null.

- Trusting cluster significance with `n_permutations < 1000` – tail
  probabilities are noisy.

- Permuting **pixels** instead of condition labels (the function does
  the right thing internally; this is just a warning if you reimplement
  it yourself).

## Reliability metrics expect raw masks

Welch t and cluster mass are variance-based and sensitive to any
scaling. If `signal_matrix_a` / `_b` came from rendered PNGs, results
are distorted. See
[`vignette("tutorial", package = "rcicrely")`](https://olivethree.github.io/rcicrely/articles/tutorial.md)
chapter 3.

## References

Maris, E., & Oostenveld, R. (2007). Nonparametric statistical testing of
EEG- and MEG-data. *Journal of Neuroscience Methods*.

Nichols, T. E., & Holmes, A. P. (2002). Nonparametric permutation tests
for functional neuroimaging: a primer with examples. *Human Brain
Mapping*.

## See also

[`rel_dissimilarity()`](https://olivethree.github.io/rcicrely/reference/rel_dissimilarity.md),
[`run_between()`](https://olivethree.github.io/rcicrely/reference/run_between.md)
