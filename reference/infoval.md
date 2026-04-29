# Per-producer informational value with trial-count-matched reference

Computes a z-scored informational value (infoVal) for each producer's
classification image, using a reference distribution matched to that
producer's trial count. Handles both 2IFC and Brief-RC paradigms with a
single function: the difference is entirely in what the user passes as
`noise_matrix`, not in how the statistic is computed.

## Usage

``` r
infoval(
  signal_matrix,
  noise_matrix,
  trial_counts,
  iter = 10000L,
  mask = NULL,
  cache_path = NULL,
  seed = NULL,
  progress = TRUE
)
```

## Arguments

- signal_matrix:

  Pixels x participants numeric matrix of raw masks (as returned by
  [`ci_from_responses_2ifc()`](https://olivethree.github.io/rcicrely/reference/ci_from_responses_2ifc.md)
  or
  [`ci_from_responses_briefrc()`](https://olivethree.github.io/rcicrely/reference/ci_from_responses_briefrc.md)).

- noise_matrix:

  Pixels x pool-size numeric matrix of noise patterns (for 2IFC: the
  `stimuli` object returned by
  `rcicr::generateStimuli2IFC(..., return_as_dataframe = TRUE)`; for
  Brief-RC: the same pool). Row count must match `signal_matrix`.

- trial_counts:

  Named integer vector of trial counts per producer. Names must match
  `colnames(signal_matrix)`.

- iter:

  Reference-distribution Monte Carlo size. Default 10000 (matches rcicr
  convention).

- mask:

  Optional logical vector of length `nrow(signal_matrix)`. When
  supplied, both observed and reference norms are computed on the masked
  pixel subset.

- cache_path:

  Optional path to an `.rds` file. When set and the file exists, the
  reference distributions are loaded from it (keyed on trial count);
  when set but the file does not exist, the computed distributions are
  written there after running. This accelerates repeated runs on the
  same `(noise_matrix, iter, mask, seed)` configuration. The cache
  stores only reference norms, never observed data.

- seed:

  Optional integer; RNG state restored on exit.

- progress:

  Show a `cli` progress bar.

## Value

Object of class `rcicrely_infoval` with fields:

- `$infoval`, per-producer z-scores (named numeric vector).

- `$norms`, per-producer observed Frobenius norms.

- `$reference`, named list of reference norm vectors, keyed by trial
  count.

- `$ref_median`, `$ref_mad`, named numeric vectors (centre / spread) per
  trial count.

- `$trial_counts`, `$mask`, `$iter`, `$n_pool`, `$seed`, metadata.

## Details

The observed statistic is the Frobenius norm of producer j's
classification image mask, optionally restricted to a logical `mask`
(e.g. an oval face region via
[`face_mask()`](https://olivethree.github.io/rcicrely/reference/face_mask.md)):

    norm_j = sqrt(sum( signal_matrix[mask, j]^2 ))

The reference distribution for producer j is built by simulating `iter`
random masks at the same trial count `trial_counts[j]`, mirroring the
`genMask()` construction exactly (Schmitz, Rougier, & Yzerbyt, 2024):

1.  Sample `trial_counts[j]` stimulus ids uniformly from `1:n_pool`
    (without replacement when `trial_counts[j] <= n_pool`, the typical
    case for both 2IFC and standard Brief-RC; with replacement only when
    the producer recorded more contributions than the pool size). See
    the v0.2.1 calibration note below.

2.  Sample `trial_counts[j]` responses uniformly from `{-1, +1}` (always
    with replacement).

3.  Collapse `response` by `stim` via
    [`mean()`](https://rdrr.io/r/base/mean.html) (this is the
    duplicate-stim rule used by Brief-RC).

4.  Build the mask:
    `noise_matrix[, unique_stims] %*% mean_response / n_unique_stims`.

5.  Apply `mask` (if supplied) and compute the Frobenius norm.

The per-producer z-score is `(norm_j - median(ref)) / mad(ref)`, where
the reference distribution is the one matched to producer j's trial
count. Producers sharing a trial count share a reference; this makes the
simulation efficient when several producers have the same number of
trials.

Choice of statistic: the **Frobenius norm** is the correct form per
Schmitz, Rougier, & Yzerbyt (2020, comment) and the 2020 erratum. The
original
[`rcicr::computeInfoVal2IFC()`](https://rdrr.io/pkg/rcicr/man/computeInfoVal2IFC.html)
already uses it (`norm(matrix(target_ci[["ci"]]), "f")`); this function
matches that convention for both 2IFC and Brief-RC data.

Trial-count matching: the original
[`rcicr::generateReferenceDistribution2IFC()`](https://rdrr.io/pkg/rcicr/man/generateReferenceDistribution2IFC.html)
builds the reference using the full pool size
(`n_trials = ncol(noise_matrix)`). For 2IFC this is appropriate because
every producer typically responds to every pool item, so per-producer
trial count equals the pool size. For Brief-RC the producer's *cognitive
exposure* per trial is richer than 2IFC (12 noisy faces vs 2), but only
one chosen `(stim, +/-1)` is recorded per trial, so the number of mask
contributions is `n_trials` and is typically smaller than `n_pool`. A
pool-size reference is therefore a mismatch for Brief-RC and biases
infoVal downward. `infoval()` closes this gap by keying the reference on
the actual per-producer recorded trial count.

What infoVal measures (and what it does not). The Frobenius norm is a
**magnitude** statistic, not a target-alignment one: it answers "is the
mask larger than chance?" but not "is it pointing at the right
pattern?". Two consequences worth flagging when interpreting
per-producer z-scores:

- **Cross-paradigm comparisons need care.** Brief-RC and 2IFC are placed
  on the same z-scale by `infoval()` (per-producer trial- count-matched
  reference, identical mask construction), but the *cognitive* processes
  generating those masks differ. A producer who benefits from Brief-RC's
  richer per-trial context might produce a more accurately localized but
  not necessarily larger mask, and the magnitude metric will not reward
  that. Treat absolute Brief-RC vs 2IFC z-score differences with
  caution.

- **Stability and discriminability are different questions.**
  [`rel_split_half()`](https://olivethree.github.io/rcicrely/reference/rel_split_half.md)
  /
  [`rel_icc()`](https://olivethree.github.io/rcicrely/reference/rel_icc.md)
  ask whether the signal is stable;
  [`rel_cluster_test()`](https://olivethree.github.io/rcicrely/reference/rel_cluster_test.md)
  /
  [`rel_dissimilarity()`](https://olivethree.github.io/rcicrely/reference/rel_dissimilarity.md)
  ask whether conditions are separable. These are the right complements
  to a magnitude-based infoVal.

## References

Brinkman, L., Goffin, S., van de Schoot, R., van Haren, N. E. M.,
Dotsch, R., & Aarts, H. (2019). Quantifying the informational value of
classification images. *Behavior Research Methods*, 51(5), 2059-2073.
[doi:10.3758/s13428-019-01232-2](https://doi.org/10.3758/s13428-019-01232-2)

Schmitz, M., Rougier, M., & Yzerbyt, V. (2020). Comment on "Quantifying
the informational value of classification images": A miscomputation of
the infoVal metric. *Behavior Research Methods*, 52(3), 1383-1386.
[doi:10.3758/s13428-019-01295-1](https://doi.org/10.3758/s13428-019-01295-1)

Schmitz, M., Rougier, M., Yzerbyt, V., Brinkman, L., & Dotsch, R.
(2020). Erratum to: Comment on "Quantifying the informational value of
classification images": Miscomputation of infoVal metric was a minor
issue and is now corrected. *Behavior Research Methods*, 52(4),
1800-1801.
[doi:10.3758/s13428-020-01367-7](https://doi.org/10.3758/s13428-020-01367-7)

Schmitz, M., Rougier, M., & Yzerbyt, V. (2024). Introducing the brief
reverse correlation: an improved tool to assess visual representations.
*European Journal of Social Psychology*.
[doi:10.1002/ejsp.3100](https://doi.org/10.1002/ejsp.3100)

## See also

[`face_mask()`](https://olivethree.github.io/rcicrely/reference/face_mask.md),
[`ci_from_responses_2ifc()`](https://olivethree.github.io/rcicrely/reference/ci_from_responses_2ifc.md),
[`ci_from_responses_briefrc()`](https://olivethree.github.io/rcicrely/reference/ci_from_responses_briefrc.md).
