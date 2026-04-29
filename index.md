# rcicrely

## Why this package?

You ran a reverse correlation study, each participant saw pairs of noisy
faces and picked the one that looked more *trustworthy*. After averaging
across participants, you have a plausible looking group-level
classification image (CI). Before reporting it, two questions deserve an
honest answer though:

1.  **Did the participants actually agree?** If you split your sample in
    half, do the two halves produce similar CIs, or did you just average
    a pile of noise?
2.  **Is the “trustworthy” CI really different from your comparison
    condition** (say, “untrustworthy”, or “healthy”), or could the
    apparent difference be mere chance?

`rcicrely` answers both by working directly on the pixel-level signal
produced by your participants. No second-phase study where independent
raters subjectively rate CIs on rating scales, and no inheriting the
Type I error inflation of that two-phase design (Cone, Brown-Iannuzzi,
Lei, & Dotsch, 2021).

It works whether you ran a standard **2IFC** task (the classic
reverse-correlation paradigm) or a **Brief-RC 12** task (Schmitz,
Rougier & Yzerbyt, 2024 — 12 noisy faces per trial instead of 2).

## Installation

``` r
# install.packages("remotes")
remotes::install_github("olivethree/rcicrely")
```

If you ran a 2IFC study, you’ll also need `rcicr` (used to compute the
individual CIs):

``` r
remotes::install_github("rdotsch/rcicr")
```

If you only run Brief-RC, you can skip `rcicr` — the Brief-RC code is
fully native to `rcicrely`.

## User’s guide

The **[full user
guide](https://olivethree.github.io/rcicrely/articles/tutorial.html)**
walks through every exported function with worked examples, including
real-data examples from a published 2IFC reverse- correlation study
([Oliveira et al., 2019](https://doi.org/10.1002/ejsp.2569)).

## Quick start

Suppose you ran a study with two conditions, *trustworthy* and
*untrustworthy*. Here’s the whole workflow.

### Starting from raw trial-level data

``` r
library(rcicrely)

# 1. Load your trial-level experiment data (one row per trial; needs
#    columns for participant id, stimulus number, and +/- 1 response)
trust_responses   <- read.csv("data/trust_condition_response_data.csv")    # experiment data from the trustworthy condition
untrust_responses <- read.csv("data/untrust_condition_response_data.csv")  # experiment data from the untrustworthy condition

# 2. Build per-participant CIs from the trial data
trustworthy <- ci_from_responses_2ifc(
  responses  = trust_responses,
  rdata_path = "data/rcicr_stimuli.Rdata",   # from rcicr stimulus gen
  baseimage  = "base"
)
untrustworthy <- ci_from_responses_2ifc(
  responses  = untrust_responses,
  rdata_path = "data/rcicr_stimuli.Rdata",
  baseimage  = "base"
)

# 3. Did participants in each condition agree with each other?
print(run_within(trustworthy$signal_matrix,   seed = 1))
print(run_within(untrustworthy$signal_matrix, seed = 1))

# 4. Are the two conditions actually distinguishable?
result <- run_between(
  trustworthy$signal_matrix,
  untrustworthy$signal_matrix,
  seed = 1
)
print(result)
plot(result)   # cluster map of pixels where conditions differ
```

### Starting from CI images already saved as PNG/JPEG

If you’ve already generated your CIs and have them as image files on
disk:

``` r
library(rcicrely)

trustworthy   <- load_signal_matrix("data/cis_trustworthy/",   "data/base.jpg")
untrustworthy <- load_signal_matrix("data/cis_untrustworthy/", "data/base.jpg")

print(run_within(trustworthy,   seed = 1))
print(run_within(untrustworthy, seed = 1))
print(run_between(trustworthy, untrustworthy, seed = 1))
```

> **Heads-up.** Results from PNG/JPEG inputs may differ slightly from
> raw-response inputs because PNGs encode the *display-scaled* CI, not
> the raw signal. The package warns you once per session when this
> matters; see [the
> tutorial](https://olivethree.github.io/rcicrely/articles/tutorial.html)
> for the full story. When possible always use CIdata instead of
> deriving it from existing images.

For a function-by-function walkthrough, including how to interpret
cluster maps, sample-size warnings, when to use ICC(3,1) vs ICC(3,k),
and Brief-RC end-to-end; see the **[full user
guide](https://olivethree.github.io/rcicrely/articles/tutorial.html)**.

## Behind the scenes

If you want the technical details:

- **Within-condition reliability**: permuted split-half with
  Spearman–Brown correction, leave-one-out influence diagnostic,
  ICC(3,1) and ICC(3,k) computed via direct mean-squares.
- **Between-condition discriminability**: vectorised Welch pixel *t*,
  cluster-based permutation with max-statistic FWER control (or
  threshold-free cluster enhancement), bootstrap representational
  dissimilarity (Euclidean distance with percentile CIs).
- **Per-producer informational value
  ([`infoval()`](https://olivethree.github.io/rcicrely/reference/infoval.md))**:
  Frobenius-norm z-score with a reference distribution matched to each
  producer’s trial count — for both 2IFC and Brief-RC.

Companion to
[`rcicrdiagnostics`](https://github.com/olivethree/rcicrdiagnostics):
that one catches silent data-processing errors *before* CI computation;
`rcicrely` asks whether the CIs *themselves* are reliable and
discriminable *after* they’ve been cleanly computed.

## Citation

If `rcicrely` helps your research, please cite it:

    Oliveira, M. (2026). rcicrely: Toolkit for reliability analysis of
    classification images from reverse correlation studies in social
    psychology. R package v0.2.2. Zenodo.
    https://doi.org/10.5281/zenodo.19772888

Or run `citation("rcicrely")` in R for a BibTeX entry.

## License

MIT © Manuel Oliveira.

## Credits

Manuel Oliveira — <https://www.manueloliveira.nl>. Development was
assisted by Claude (Anthropic); the author is responsible for all
content and decisions.
