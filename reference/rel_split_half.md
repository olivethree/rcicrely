# Permuted split-half reliability with Spearman-Brown correction

Measures how stable a condition's group-level CI is across random halves
of its producer set. For each of `n_permutations` iterations: randomly
partition the producers (columns of `signal_matrix`) into two
roughly-equal halves, average each half to get two group vectors, and
compute the Pearson correlation.

Use this to answer "would I get the same group CI if I'd run a different
half of my participants?".

For **odd N**, one randomly-chosen producer is dropped per permutation
(re-drawn each iteration, not fixed) so both halves have `floor(N/2)`
producers.

Reports both the mean per-permutation `r` and the Spearman-Brown
corrected reliability `r_SB = (2 r) / (1 + r)`. 95% CI is the 2.5 / 97.5
percentile of the permutation distribution.

## Usage

``` r
rel_split_half(
  signal_matrix,
  n_permutations = 2000L,
  mask = NULL,
  seed = NULL,
  progress = TRUE
)
```

## Arguments

- signal_matrix:

  Pixels x participants, base-subtracted.

- n_permutations:

  Integer number of random splits. Default 2000, cheap and keeps Monte
  Carlo error on tail probabilities below 0.01.

- mask:

  Optional logical vector of length `nrow(signal_matrix)` restricting
  computation to a region of the image (e.g., from
  [`face_mask()`](https://olivethree.github.io/rcicrely/reference/face_mask.md)
  or
  [`load_face_mask()`](https://olivethree.github.io/rcicrely/reference/load_face_mask.md)).
  When supplied, the matrix is row-subsetted before any permutation. See
  [`vignette("tutorial")`](https://olivethree.github.io/rcicrely/articles/tutorial.md)
  §6.6 for the apply-symmetrically rule.

- seed:

  Optional integer; if set, results are reproducible. The caller's
  global RNG state is restored on exit.

- progress:

  Show a `cli` progress bar.

## Value

An object of class `rcicrely_split_half`. Fields described in **Reading
the result** above.

## Reading the result

- `$r_hh`, mean per-permutation Pearson r between the two halves.

- `$r_sb`, Spearman-Brown projected reliability of the full sample. This
  is usually the headline number to report.

- `$ci_95`, `$ci_95_sb`, percentile 95% CIs on each.

- `$distribution`, the full per-permutation r vector for plotting or
  custom CIs.

- `$n_participants`, `$n_permutations`, metadata.

## Common mistakes

- Reporting `$r_hh` instead of `$r_sb` when the audience cares about the
  full-sample CI's reliability (almost always).

- Running with a small `n_permutations` (\< 200) to save time, then
  trusting the CI bounds, CIs widen with low n_perm.

## Reliability metrics expect raw masks

Operates on the raw mask; results may be distorted if `signal_matrix`
was extracted from rendered (scaled) PNGs. See
[`vignette("tutorial", package = "rcicrely")`](https://olivethree.github.io/rcicrely/articles/tutorial.md)
chapter 3.

## References

Brinkman, L., Goffin, S., van de Schoot, R., van Haren, N. E. M.,
Dotsch, R., & Aarts, H. (2019). Quantifying the informational value of
classification images. *Behavior Research Methods*, 51(5), 2059-2073.
[doi:10.3758/s13428-019-01232-2](https://doi.org/10.3758/s13428-019-01232-2)

Shrout, P. E., & Fleiss, J. L. (1979). Intraclass correlations: uses in
assessing rater reliability. *Psychological Bulletin*, 86(2), 420-428.
[doi:10.1037/0033-2909.86.2.420](https://doi.org/10.1037/0033-2909.86.2.420)

## See also

[`rel_loo()`](https://olivethree.github.io/rcicrely/reference/rel_loo.md),
[`rel_icc()`](https://olivethree.github.io/rcicrely/reference/rel_icc.md),
[`run_within()`](https://olivethree.github.io/rcicrely/reference/run_within.md)

## Examples

``` r
if (FALSE) { # \dontrun{
r <- rel_split_half(signal_matrix, n_permutations = 2000, seed = 1)
print(r); plot(r)
} # }
```
