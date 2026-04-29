# Leave-one-out influence screening

Screens for influential producers by asking: "if this producer is
removed from the sample, how much does the group classification image
change?" Producers whose removal moves the group pattern
disproportionately are flagged for inspection. This is an **influence /
outlier screening** tool, not a reliability statistic — see below.

## Usage

``` r
rel_loo(
  signal_matrix,
  flag_threshold = 2.5,
  flag_method = c("mad", "sd"),
  flag_threshold_sd = NULL,
  mask = NULL
)
```

## Arguments

- signal_matrix:

  Pixels x participants, base-subtracted.

- flag_threshold:

  Numeric multiplier on `sd` (or `mad`) below the centre, defining the
  outlier cutoff. Default 2.5. Was `flag_threshold_sd` before v0.1.1;
  the old name still works as an alias.

- flag_method:

  One of `"mad"` (default since v0.3) or `"sd"`. The MAD/median rule is
  robust to the very influential producers it is meant to detect (one
  outlier inflates `sd_r` and pulls `mean_r`, masking itself); SD/mean
  is retained for backwards compatibility and emits a one-time
  per-session deprecation message. SD is scheduled for removal in v0.4.
  See **Details** for the mathematical justification.

- flag_threshold_sd:

  Deprecated alias for `flag_threshold`. Kept for backwards
  compatibility with v0.1.0.

- mask:

  Optional logical vector of length `nrow(signal_matrix)` restricting
  computation to a region (e.g., from
  [`face_mask()`](https://olivethree.github.io/rcicrely/reference/face_mask.md)
  or
  [`load_face_mask()`](https://olivethree.github.io/rcicrely/reference/load_face_mask.md)).

## Value

Object of class `rcicrely_loo`. Fields described in **Reading the
result** above.

## Details

For each producer `i`, compute the Pearson correlation between the
full-sample group CI and the group CI recomputed without that producer:

    full        <- rowMeans(signal_matrix)
    r_loo[i]    <- cor(full, rowMeans(signal_matrix[, -i]))

Because the full-sample mean and the leave-one-out mean share
`(N - 1) / N` of their data, `r_loo` values are near 1 by construction
even on noisy data, typically `[0.95, 0.999]` at `N = 30`. **The
absolute level of `r_loo` is not informative.** What is informative is
the *relative ordering*: producers whose `r_loo` sits clearly below the
pack are candidates for inspection. For this reason the function also
returns a **z-scored** version of `r_loo` in `$z_scores`, using the same
centre / spread estimators as the flagging rule. `$z_scores` is the
recommended quantity to plot or report.

Two flagging rules:

- `"sd"` (default): flag producers with
  `r_loo < mean(r) - flag_threshold * sd(r)`, equivalently
  `z_scores < -flag_threshold`. Standard convention; sensitive to the
  very outliers it is trying to detect, so the flagging threshold can be
  pulled into the outliers' region.

- `"mad"`: flag producers with
  `r_loo < median(r) - flag_threshold * mad(r)`, equivalently
  `z_scores < -flag_threshold` with a robust centre/spread. Robust to
  the few atypical producers that typically show up in RC datasets and
  to skewed correlation distributions.

The default `flag_threshold = 2.5` is calibrated so that a 30-producer
dataset flags roughly 0.3 producers by chance under `"sd"`, rather than
the ~1.5 a 2-SD rule would produce. Under `"mad"` it is roughly
comparable thanks to MAD's 1.4826 consistency factor.

## What this function is, and is not

`rel_loo()` is an **influence-screening diagnostic**. It answers "which
producers disproportionately shape the group CI?", not "how reliable is
the group CI?". For reliability, use
[`rel_split_half()`](https://olivethree.github.io/rcicrely/reference/rel_split_half.md)
or
[`rel_icc()`](https://olivethree.github.io/rcicrely/reference/rel_icc.md).
A flag does not mean the producer is "bad"; it means the producer's
individual CI sits far enough from the group pattern that the data
deserve a second look (response coding, fatigue, task misunderstanding,
or a genuinely atypical mental representation).

## Reading the result

- `$z_scores`, named numeric vector, per-producer standardised
  influence. This is the recommended quantity to plot, report, or
  threshold. Values near 0 = typical producer; values below
  `-flag_threshold` = flagged.

- `$correlations`, named numeric vector, raw per-producer `r_loo`
  values. Included for transparency; see the note above about why the
  raw level is not informative.

- `$mean_r`, `$sd_r`, `$median_r`, `$mad_r`, centre / spread of the
  correlation distribution under each rule.

- `$threshold`, the raw cutoff value on `r_loo` under the chosen
  `flag_method`.

- `$flagged`, character vector of producer ids below threshold.

- `$summary_df`, one row per producer with `correlation`, `z_score`, and
  `flag`, sorted by `z_score`.

- `$flag_method`, `$flag_threshold`, what was used.

## Common mistakes

- Reading `r_loo` as a reliability. An `r_loo` of .98 does not mean the
  CI is 98% reliable; it means a single producer's removal changed the
  group mean by 2%, which is the expected scale at `N = 30`. Report
  reliability via
  [`rel_split_half()`](https://olivethree.github.io/rcicrely/reference/rel_split_half.md)
  or
  [`rel_icc()`](https://olivethree.github.io/rcicrely/reference/rel_icc.md).

- Treating `$flagged` as "drop these producers". Investigate first
  (response coding, fatigue, atypical strategy). Cross-check with the
  `rcicrdiagnostics` companion package to rule out response-coding
  errors.

- Lowering `flag_threshold` below 2 to flag more producers. That trades
  real signal for noise; use `flag_method = "mad"` instead if the SD
  rule is dominated by the outliers it would otherwise catch.

## Reliability metrics expect raw masks

Operates on the raw mask; results may be distorted if `signal_matrix`
was extracted from rendered (scaled) PNGs. See
[`vignette("tutorial", package = "rcicrely")`](https://olivethree.github.io/rcicrely/articles/tutorial.md)
chapter 3.

## See also

[`rel_loo_z()`](https://olivethree.github.io/rcicrely/reference/rel_loo_z.md)
for a tidy z-score accessor;
[`rel_split_half()`](https://olivethree.github.io/rcicrely/reference/rel_split_half.md),
[`rel_icc()`](https://olivethree.github.io/rcicrely/reference/rel_icc.md)
for reliability metrics proper;
[`run_within()`](https://olivethree.github.io/rcicrely/reference/run_within.md).
