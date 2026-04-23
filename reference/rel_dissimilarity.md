# Representational dissimilarity with bootstrap CIs

Compute two scalar dissimilarity metrics between two conditions'
group-level signals – Pearson correlation and Euclidean distance – with
percentile bootstrap confidence intervals.

Use this when you want a single number (or two) summarising "how
different are these two conditions' CIs?", with uncertainty. Pair with
[`rel_cluster_test()`](https://olivethree.github.io/rcicrely/reference/rel_cluster_test.md)
when you also want to know **where** they differ.

Bootstrap: within each condition, resample participants with
replacement; recompute the group mean and the two metrics. `n_boot`
replicates give percentile CIs via
[`quantile()`](https://rdrr.io/r/stats/quantile.html). No `boot`
dependency. BCa is future work.

## Usage

``` r
rel_dissimilarity(
  signal_matrix_a,
  signal_matrix_b,
  n_boot = 2000L,
  ci_level = 0.95,
  seed = NULL,
  progress = TRUE
)
```

## Arguments

- signal_matrix_a, signal_matrix_b:

  Pixels x participants, base-subtracted. Row counts must match.

- n_boot:

  Bootstrap replicates. Default 2000.

- ci_level:

  Confidence level. Default 0.95.

- seed:

  Optional integer; RNG state restored on exit.

- progress:

  Show a `cli` progress bar.

## Value

Object of class `rcicrely_dissim`. Fields described in **Reading the
result** above.

## Reading the result

- `$correlation`, `$euclidean` – observed values on the full sample.
  Pearson r is scale-tolerant for a single uniform scaling; Euclidean
  distance is sensitive to any scaling.

- `$boot_cor`, `$boot_dist` – full bootstrap distributions.

- `$ci_cor`, `$ci_dist` – percentile CIs at `ci_level`.

- `$boot_se_cor`, `$boot_se_dist` – bootstrap SEs.

## Common mistakes

- Reading `$correlation` near 1 as "the conditions look the same". Two
  CIs with identical spatial pattern but different magnitudes produce r
  near 1; the magnitude difference shows up in `$euclidean`.

- Treating `$euclidean` as comparable across pixel resolutions – it
  scales with sqrt(n_pixels). Use it as a within-study metric.

## Reliability metrics expect raw masks

Pearson r survives a single uniform linear scaling but breaks under
per-CI `"matched"`-style scaling; Euclidean distance is sensitive to
**any** scaling. If `signal_matrix_a` / `_b` came from
[`read_cis()`](https://olivethree.github.io/rcicrely/reference/read_cis.md)
/
[`extract_signal()`](https://olivethree.github.io/rcicrely/reference/extract_signal.md)
on rendered PNGs, treat results as approximate. See
[`vignette("tutorial", package = "rcicrely")`](https://olivethree.github.io/rcicrely/articles/tutorial.md)
chapter 3.

## See also

[`rel_cluster_test()`](https://olivethree.github.io/rcicrely/reference/rel_cluster_test.md),
[`run_between()`](https://olivethree.github.io/rcicrely/reference/run_between.md)
