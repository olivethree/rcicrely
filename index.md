# rcicrely

> **Reliability assessment for reverse-correlation classification
> images.** Within-condition consistency + between-condition
> discriminability, without a second-phase trait-rating study.

`rcicrely` assesses the reliability of classification images (CIs)
produced by reverse-correlation experiments. It answers the questions
`infoVal` doesn’t: is the signal in your CI *stable*, and is it
*distinguishable* from a comparison condition?

- **Within-condition**: permuted split-half with Spearman–Brown
  correction, leave-one-out sensitivity, ICC(3,1) and ICC(3,k) via
  direct mean-squares.
- **Between-condition**: vectorised Welch pixel *t*, cluster-based
  permutation with max-statistic FWER control, bootstrap
  representational dissimilarity.

Supports both standard **2IFC** (via the canonical
[`rcicr`](https://github.com/rdotsch/rcicr)) and **Brief-RC 12**
(Schmitz, Rougier & Yzerbyt, 2024; implemented natively).

Companion to
[`rcicrdiagnostics`](https://github.com/olivethree/rcicrdiagnostics):
diagnostics catches silent data-processing errors *before* CI
computation; `rcicrely` asks whether the CIs *themselves* are reliable
and discriminable *after* they have been cleanly computed.

## Installation

`rcicrely` is GitHub-only. `rcicr` (required only for the 2IFC path) is
also GitHub-only:

``` r
# install.packages("remotes")
remotes::install_github("rdotsch/rcicr")      # 2IFC path only
remotes::install_github("olivethree/rcicrely")
```

If you only use Brief-RC, you can skip `rcicr`; the Brief-RC CI
implementation is fully native.

## Quick start

``` r
library(rcicrely)

# Option A: from a directory of pre-computed CI images
signal_A <- load_signal_matrix(
  dir             = "data/cis_condition_A/",
  base_image_path = "data/base.jpg"
)
signal_B <- load_signal_matrix("data/cis_condition_B/", "data/base.jpg")

# Option B: from raw trial-level responses (2IFC)
res_A <- ci_from_responses_2ifc(
  responses  = my_responses_A,
  rdata_path = "data/rcicr_stimuli.Rdata",
  baseimage  = "base"
)
signal_A <- res_A$signal_matrix

# Reliability reports
within  <- run_within(signal_A, n_permutations = 2000, seed = 1L)
between <- run_between(signal_A, signal_B,
                       n_permutations = 2000, n_boot = 2000, seed = 1L)

print(within);  plot(within)
print(between); plot(between)
```

See the [user guide
vignette](https://olivethree.github.io/rcicrely/articles/tutorial.html)
for a function-by-function walkthrough covering both paradigms,
cluster-map interpretation, sample-size warnings, and the ICC-variant
decision.

## Important: raw vs. rendered CIs

Reliability metrics in this package operate on the **raw mask** (the
participant’s noise contribution before any display transformation). The
two ways to get there:

- **Mode 2 (raw responses)** -
  [`ci_from_responses_2ifc()`](https://olivethree.github.io/rcicrely/reference/ci_from_responses_2ifc.md)
  and
  [`ci_from_responses_briefrc()`](https://olivethree.github.io/rcicrely/reference/ci_from_responses_briefrc.md)
  return the raw mask directly. Safe.
- **Mode 1 (PNGs on disk)** -
  [`read_cis()`](https://olivethree.github.io/rcicrely/reference/read_cis.md)
  /
  [`extract_signal()`](https://olivethree.github.io/rcicrely/reference/extract_signal.md)
  /
  [`load_signal_matrix()`](https://olivethree.github.io/rcicrely/reference/load_signal_matrix.md)
  load the *rendered* CI (`base + scaling(mask)`) because that’s what
  was saved to disk; subtracting the base then yields `scaling(mask)`,
  **not** the raw mask.

Variance-based reliability metrics
([`rel_icc()`](https://olivethree.github.io/rcicrely/reference/rel_icc.md),
Euclidean half of
[`rel_dissimilarity()`](https://olivethree.github.io/rcicrely/reference/rel_dissimilarity.md),
[`pixel_t_test()`](https://olivethree.github.io/rcicrely/reference/pixel_t_test.md),
[`rel_cluster_test()`](https://olivethree.github.io/rcicrely/reference/rel_cluster_test.md))
are sensitive to the scaling transform. The canonical 2IFC `infoVal`
path
([`rcicr::computeInfoVal2IFC()`](https://rdrr.io/pkg/rcicr/man/computeInfoVal2IFC.html))
extracts the raw `$ci` from the rcicr CI-list internally and is
unaffected; hand-rolled `infoVal` implementations (Brief-RC, custom
code) need the raw mask and should not be fed PNG-derived data. The
package emits a one-time per-session warning when you use Mode 1;
silence with `acknowledge_scaling = TRUE` or
`options(rcicrely.silence_scaling_warning = TRUE)` once you’ve read
[chapter 3 of the user
guide](https://olivethree.github.io/rcicrely/articles/tutorial.html).

## What’s different from the trait-rating ICCs in the RC literature

Most ICCs published in reverse-correlation papers are **trait-rating**
reliability (phase-2 naive raters scoring CIs on dimensions like
“trustworthy” or “competent”). `rcicrely`’s ICC is a structurally
different object: it operates on the pixel-level signal produced by the
original producers, so no phase-2 rating study is involved. This
sidesteps the two-phase design Cone, Brown-Iannuzzi, Lei, & Dotsch
(2021) showed inflates Type I error.

The package reports **ICC(3,1)** and **ICC(3,k)** by positive argument,
not convention: pixels are a fixed `img_size × img_size` grid (two-way
mixed, pixels fixed, participants random), not a random sample from a
pixel population. ICC(2,\*) is available via the `variants` argument if
your reviewer explicitly asks.

## Citation

    Oliveira, M. (2026). rcicrely: Reliability assessment for
    reverse-correlation classification images. R package v0.1.0.
    https://github.com/olivethree/rcicrely

## License

MIT © Manuel Oliveira.

## Acknowledgements

This package owes its test harness and code-review discipline to work
done in collaboration with Anthropic’s Claude.
