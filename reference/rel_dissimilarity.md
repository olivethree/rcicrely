# Between-condition dissimilarity with bootstrap confidence intervals

Quantifies overall dissimilarity between two conditions' group-level
classification images. The primary statistic is **Euclidean distance**
between the two group-mean CIs, reported both raw and normalised by
`sqrt(n_pixels)` for cross-resolution comparability. Percentile
bootstrap 95% confidence intervals are computed by resampling
participants with replacement within each condition.

Pair with
[`rel_cluster_test()`](https://olivethree.github.io/rcicrely/reference/rel_cluster_test.md)
when you also want to know **where** the two conditions differ.

## Usage

``` r
rel_dissimilarity(
  signal_matrix_a,
  signal_matrix_b,
  paired = FALSE,
  n_boot = 2000L,
  ci_level = 0.95,
  mask = NULL,
  seed = NULL,
  progress = TRUE
)
```

## Arguments

- signal_matrix_a, signal_matrix_b:

  Pixels x participants, base-subtracted. Row counts must match.

- paired:

  Logical. `FALSE` (default) for a between-subjects design: participants
  are resampled within A and B independently. `TRUE` for a
  within-subjects design: A and B share a single resample index per
  replicate so the paired covariance structure is preserved. Requires
  `ncol(a) == ncol(b)` and, if named, identical column names.

- n_boot:

  Bootstrap replicates. Default 2000.

- ci_level:

  Confidence level. Default 0.95.

- mask:

  Optional logical vector of length `nrow(signal_matrix_a)` restricting
  the Euclidean / correlation computation to a region. Both matrices are
  subsetted with the same mask. The reported `n_pixels` and
  `euclidean_normalised` (= `euclidean / sqrt(n_pixels)`) reflect the
  masked count.

- seed:

  Optional integer; RNG state restored on exit.

- progress:

  Show a `cli` progress bar.

## Value

Object of class `rcicrely_dissim`. Fields described in **Reading the
result** above.

## Details

For the observed statistics and each bootstrap replicate `i`:

    mean_a      = rowMeans(signal_matrix_a)
    mean_b      = rowMeans(signal_matrix_b)
    observed_dist            = sqrt(sum((mean_a - mean_b)^2))
    observed_dist_normalised = observed_dist / sqrt(n_pixels)

Percentile CI via base R
[`quantile()`](https://rdrr.io/r/stats/quantile.html); no `boot`
dependency. BCa intervals are on the roadmap for a future release.

## Why Euclidean and not Pearson correlation as the primary

An earlier release emphasised Pearson correlation between the two
group-mean CIs as a co-primary metric. The current release retains those
fields **for backwards compatibility only**; they are scheduled for
removal in v0.3 and **should not be reported as primary** in new work.
The reason is a positive baseline: two base-subtracted CIs share
systematic image-domain spatial structure (face shape, oval signal
support, low-frequency Gaussian-noise smoothness) that pushes their
correlation above zero even when the underlying mental representations
are unrelated. Absolute correlation values therefore do not cleanly mean
"these conditions are similar"; they conflate real similarity with
shared image-domain structure.

Euclidean distance does not share this failure mode. It is an honest
magnitude summary of how far the two group-mean CIs sit from each other
in pixel-signal units. For cross-study comparisons, use
`$euclidean_normalised` (divided by `sqrt(n_pixels)` so the metric is
resolution-invariant).

## Reading the result

- `$euclidean`, observed Euclidean distance between group means on the
  full sample (primary statistic).

- `$euclidean_normalised`, `$euclidean / sqrt(n_pixels)`. Use this for
  cross-resolution comparisons.

- `$boot_dist`, `$ci_dist`, `$boot_se_dist`, bootstrap distribution,
  percentile CI, and SE of the Euclidean distance.

- `$correlation`, `$boot_cor`, `$ci_cor`, `$boot_se_cor`: Pearson
  correlation of the group means and its bootstrap summary.
  **Deprecated.** Retained for backwards compatibility with v0.1.x code
  and scheduled for removal in v0.3. Treat these as within-study
  descriptive statistics only; see **Why Euclidean and not Pearson
  correlation**.

- `$n_boot`, `$ci_level`, metadata.

## Common mistakes

- Reporting `$correlation` near 1 as "the conditions look the same". Two
  CIs with identical spatial pattern but different magnitudes produce a
  correlation near 1; magnitude differences only show up in
  `$euclidean`.

- Comparing `$euclidean` across studies with different image
  resolutions. Use `$euclidean_normalised` for that.

- Reading `$correlation` as a standalone "similarity score". The
  baseline is above zero for independent CIs (see above).

## Reliability metrics expect raw masks

Euclidean distance is sensitive to **any** scaling; Pearson r survives a
single uniform scaling but breaks under per-CI `"matched"`-style
scaling. If `signal_matrix_a` / `_b` came from
[`read_cis()`](https://olivethree.github.io/rcicrely/reference/read_cis.md)
/
[`extract_signal()`](https://olivethree.github.io/rcicrely/reference/extract_signal.md)
on rendered PNGs, treat results as approximate. See
[`vignette("tutorial", package = "rcicrely")`](https://olivethree.github.io/rcicrely/articles/tutorial.md)
chapter 3.

## References

Efron, B., & Tibshirani, R. J. (1994). *An introduction to the
bootstrap*. Chapman & Hall / CRC.

## See also

[`rel_cluster_test()`](https://olivethree.github.io/rcicrely/reference/rel_cluster_test.md),
[`run_between()`](https://olivethree.github.io/rcicrely/reference/run_between.md)
