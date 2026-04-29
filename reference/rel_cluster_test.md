# Cluster-based permutation test with family-wise error control

Pixel-level discriminability test between two condition signal matrices.
Returns either (a) the classical threshold-based cluster test (Maris &
Oostenveld, 2007) with cluster mass as the statistic, or (b)
threshold-free cluster enhancement (TFCE; Smith & Nichols, 2009), which
integrates across thresholds and removes the arbitrary cluster-forming
cutoff. Both paths use the same stratified permutation scheme and the
maximum-statistic approach to family- wise error control (Nichols &
Holmes, 2002).

Pair with
[`rel_dissimilarity()`](https://olivethree.github.io/rcicrely/reference/rel_dissimilarity.md)
when you also want an overall magnitude summary.

## Usage

``` r
rel_cluster_test(
  signal_matrix_a,
  signal_matrix_b,
  img_dims = NULL,
  method = c("threshold", "tfce"),
  paired = FALSE,
  n_permutations = 2000L,
  cluster_threshold = 2,
  tfce_H = 2,
  tfce_E = 0.5,
  tfce_n_steps = 100L,
  alpha = 0.05,
  mask = NULL,
  seed = NULL,
  progress = TRUE,
  acknowledge_scaling = FALSE
)
```

## Arguments

- signal_matrix_a, signal_matrix_b:

  Pixels x participants, base-subtracted.

- img_dims:

  Integer `c(nrow, ncol)`. Can be inferred if the signal matrices carry
  an `img_dims` attribute.

- method:

  One of `"threshold"` (default) or `"tfce"`. The threshold-based method
  requires `cluster_threshold`; TFCE does not.

- paired:

  Logical. `FALSE` (default) for a between-subjects design: Welch t per
  pixel and stratified label permutation for the null. `TRUE` for a
  within-subjects design: paired t per pixel on `A - B` and random
  sign-flip of matched pairs for the null. When `paired = TRUE`, the two
  matrices must have the same number of columns and, if named, identical
  column names.

- n_permutations:

  Number of stratified label permutations (between-subjects) or
  sign-flip permutations (paired). Default 2000.

- cluster_threshold:

  Absolute t-value threshold for forming clusters under
  `method = "threshold"`. Default 2.0. Ignored under `method = "tfce"`.

- tfce_H, tfce_E:

  TFCE height and extent exponents. Defaults match Smith & Nichols
  (2009): `H = 2.0`, `E = 0.5`. Ignored under `method = "threshold"`.

- tfce_n_steps:

  Number of thresholds in the TFCE integration grid. Default 100. Finer
  grids cost more per permutation. Ignored under `method = "threshold"`.

- alpha:

  Significance level. Default 0.05.

- mask:

  Optional logical vector of length `n_pixels` restricting cluster /
  TFCE inference to a region (e.g., from
  [`face_mask()`](https://olivethree.github.io/rcicrely/reference/face_mask.md)
  or
  [`load_face_mask()`](https://olivethree.github.io/rcicrely/reference/load_face_mask.md)).
  Implementation uses the **zero-out pattern**, not row-subsetting:
  pixels where `!mask` have their per-pixel t set to 0 (in both observed
  and permutation t-maps) before cluster identification. This way the 2D
  image structure is preserved for 4-connectivity logic, but masked-out
  pixels can never join a cluster (`|t| = 0 < cluster_threshold`) and
  contribute 0 to TFCE. Apply the same mask in any companion analyses
  for consistency.

- seed:

  Optional integer; RNG state restored on exit.

- progress:

  Show a `cli` progress bar.

- acknowledge_scaling:

  Logical. When `FALSE` (default), the internal
  [`pixel_t_test()`](https://olivethree.github.io/rcicrely/reference/pixel_t_test.md)
  errors on a known-rendered matrix. Set `TRUE` to override; cascades
  down to that call.

## Value

Object of class `rcicrely_cluster_test`. Fields described in **Reading
the result**.

## Details

Common scaffold (both methods):

1.  Compute observed pixel-wise Welch t (see
    [`pixel_t_test()`](https://olivethree.github.io/rcicrely/reference/pixel_t_test.md)).

2.  Build the null via `n_permutations` **stratified** label
    permutations, each preserving `(N_A, N_B)` exactly.

3.  Max-statistic over the null controls FWER in the strong sense.

Method-specific step:

- `method = "threshold"` (default): threshold at
  `|t| > cluster_threshold`, label connected components with
  4-connectivity, cluster mass = sum of t-values within the cluster. Per
  observed cluster, p-value = fraction of null max-masses (matching
  sign) that equal or exceed the observed mass in absolute value.

- `method = "tfce"`: enhance the observed t-map into a TFCE map
  (integral over thresholds of `size^E * h^H dh`, pos and neg tails
  enhanced separately and recombined with sign preserved; H = 2.0, E =
  0.5 default). Compute the same enhancement for each permuted t-map and
  record `max(|TFCE|)`. Per-pixel p-value = fraction of null max-TFCE
  values that equal or exceed the observed `|TFCE|` at that pixel. No
  free threshold parameter.

## Reading the result (threshold method)

- `$observed_t`, per-pixel Welch t vector.

- `$clusters`, data.frame with `cluster_id`, `direction`
  (`"pos"`/`"neg"`), `mass`, `size`, `p_value`, `significant`.

- `$null_distribution$pos`, `$neg`, per-permutation max masses.

- `$pos_labels`, `$neg_labels`, integer matrices for plotting.

- `$cluster_threshold`, `$alpha`, `$n_permutations`,
  `$n_participants_a`, `$n_participants_b`, `$method = "threshold"`.

## Reading the result (TFCE method)

- `$observed_t`, per-pixel Welch t vector (before enhancement).

- `$tfce_map`, per-pixel signed TFCE values.

- `$tfce_pmap`, per-pixel p-values against the max-TFCE null.

- `$tfce_significant_mask`, logical vector flagging pixels with
  `p < alpha`.

- `$null_distribution$max_abs_tfce`, per-permutation max\|TFCE\|.

- `$tfce_H`, `$tfce_E`, `$tfce_n_steps`, `$alpha`, `$n_permutations`,
  `$n_participants_a`, `$n_participants_b`, `$method = "tfce"`.

## Common mistakes

- Reading `cluster_threshold` as a p-value cutoff. It is the t-value
  above which pixels join a cluster; significance is decided afterwards
  by comparing cluster mass to the null.

- Trusting cluster significance with `n_permutations < 1000`.

- Permuting **pixels** instead of condition labels (the function does
  the right thing internally; this is a warning if you reimplement it
  yourself).

## Reliability metrics expect raw masks

Welch t and cluster mass / TFCE are variance-based and sensitive to any
scaling. Inputs with `attr(., "source") == "rendered"` (set
automatically by Mode 1 readers like
[`extract_signal()`](https://olivethree.github.io/rcicrely/reference/extract_signal.md))
error unless `acknowledge_scaling = TRUE` cascades through the internal
[`pixel_t_test()`](https://olivethree.github.io/rcicrely/reference/pixel_t_test.md)
call. See
[`vignette("tutorial", package = "rcicrely")`](https://olivethree.github.io/rcicrely/articles/tutorial.md)
chapter 3.

## References

Maris, E., & Oostenveld, R. (2007). Nonparametric statistical testing of
EEG- and MEG-data. *Journal of Neuroscience Methods*, 164(1), 177-190.
[doi:10.1016/j.jneumeth.2007.03.024](https://doi.org/10.1016/j.jneumeth.2007.03.024)

Nichols, T. E., & Holmes, A. P. (2002). Nonparametric permutation tests
for functional neuroimaging: a primer with examples. *Human Brain
Mapping*, 15(1), 1-25.
[doi:10.1002/hbm.1058](https://doi.org/10.1002/hbm.1058)

Smith, S. M., & Nichols, T. E. (2009). Threshold-free cluster
enhancement: addressing problems of smoothing, threshold dependence and
localisation in cluster inference. *NeuroImage*, 44(1), 83-98.
[doi:10.1016/j.neuroimage.2008.03.061](https://doi.org/10.1016/j.neuroimage.2008.03.061)

## See also

[`rel_dissimilarity()`](https://olivethree.github.io/rcicrely/reference/rel_dissimilarity.md),
[`run_between()`](https://olivethree.github.io/rcicrely/reference/run_between.md)
