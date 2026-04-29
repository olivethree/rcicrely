# Per-pixel inferential agreement map with FWER control

Within a single condition, tests at each pixel whether the
producer-level signal differs from zero (one-sample *t*). The
permutation null is built by random sign-flip per producer (exact under
the assumption that, under H0, the producer's signal contribution is
symmetric around zero). Family-wise error is controlled across pixels by
the **maximum \|t\|** statistic.

Use this when you want a per-pixel inferential overlay on a descriptive
agreement-map plot — typically paired with
[`plot_ci_overlay()`](https://olivethree.github.io/rcicrely/reference/plot_ci_overlay.md)
(§1.6) so the significance contours are rendered on top of the observed
group CI.

## Usage

``` r
agreement_map_test(
  signal_matrix,
  n_permutations = 5000L,
  alpha = 0.05,
  mask = NULL,
  seed = NULL,
  progress = TRUE,
  acknowledge_scaling = FALSE
)
```

## Arguments

- signal_matrix:

  Pixels x participants, base-subtracted.

- n_permutations:

  Integer. Number of sign-flip iterations. Default 5000.

- alpha:

  Numeric. Significance level. Default 0.05.

- mask:

  Optional logical vector of length `nrow(signal_matrix)`; the test is
  computed on the masked pixel subset. Pixels outside the mask are
  returned as `NA_real_` per-pixel and `FALSE` in the significant mask.

- seed:

  Optional integer.

- progress:

  Show a `cli` progress bar.

- acknowledge_scaling:

  Logical. Forwarded to `assert_raw_signal()`.

## Value

Object of class `rcicrely_agreement_map_test` with:

- `$observed_t` — per-pixel one-sample *t*.

- `$pmap` — per-pixel *p*-value under the max-\|t\| null.

- `$significant_mask` — logical, `pmap < alpha`.

- `$null_distribution` — numeric vector of `max_abs_t` per permutation.

- `$alpha`, `$n_permutations`, `$n_participants`, `$mask`.

## See also

[`plot_ci_overlay()`](https://olivethree.github.io/rcicrely/reference/plot_ci_overlay.md),
[`rel_cluster_test()`](https://olivethree.github.io/rcicrely/reference/rel_cluster_test.md)
