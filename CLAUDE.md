# rcicrely

An R package for reliability assessment of classification images (CIs)
produced by reverse correlation experiments. Covers within-condition
consistency (permuted split-half with Spearman-Brown, leave-one-out,
ICC) and between-condition discriminability (cluster-based permutation
testing, representational dissimilarity). Works with both standard 2IFC
(via canonical `rcicr`, Dotsch v1.0.1) and Brief-RC (Schmitz, Rougier, &
Yzerbyt, 2024 — owned natively by this package).

Companion package to `rcicrdiagnostics`. Diagnostics catches silent
data-processing errors *before* CI computation; rcicrely asks whether
the CIs *themselves* are reliable and discriminable *after* they have
been cleanly computed. Independent: neither imports the other.

This file is the design reference for rcicrely. It is complete enough
that a competent R developer — or an AI assistant — can rebuild the
package from scratch using nothing but this document. Every non-obvious
implementation choice, every rcicr integration gotcha that cost real
time during the companion rcicrdiagnostics build, and every statistical
commitment is captured here.

CLAUDE_rcicrely.md is gitignored and developer-facing only.

## 1. Purpose and audience

### Purpose

`infoVal` (Brinkman et al., 2019) tells you whether a CI contains more
signal than chance. It does not tell you whether that signal is stable,
replicable, or distinguishable from another condition’s CI. This package
fills that gap with a reliability assessment toolkit that operates at
the pixel level, directly on the CI noise (signal) components, without
requiring the problematic two-phase rating design (Cone et al., 2020).
The checks are:

- **Within-condition** (how consistent is a condition’s CIs across
  participants?): permuted split-half r with Spearman-Brown correction,
  leave-one-out sensitivity, ICC(3,1) and ICC(3,k) via direct
  mean-squares (two-way mixed: pixels fixed, participants random; see
  §7.3 for why not ICC(2,\*)).
- **Between-condition** (is condition A’s CI distinguishable from
  condition B’s?): Welch pixel-wise t-test, cluster-based permutation
  with 4-connectivity and max-statistic FWER control, bootstrap
  correlation and Euclidean distance on full-image signals.

### Audience

RC researchers who have already computed individual CIs (either via
`rcicr` for 2IFC, or via rcicrely’s own Brief-RC path) and need to
evaluate CI quality before publication. Users know R at an intermediate
level and understand what a CI is; they are not expected to have a
psychometrics or permutation-inference background, nor to know
`data.table`, S3 class internals, or advanced package machinery.

## 2. Stack

- R \>= 4.1.

- **Imports** (CRAN, installed automatically): `data.table`, `cli`. Keep
  this list minimal — only declare what the shipped R code actually
  calls. Bootstrap CIs in
  [`rel_dissimilarity()`](https://olivethree.github.io/rcicrely/reference/rel_dissimilarity.md)
  use base R [`quantile()`](https://rdrr.io/r/stats/quantile.html)
  (percentile method), so `boot` is **not** an Import. Do **not**
  pre-declare `png`, `jpeg`, `pbapply`, or `psych` — those are either
  Suggests or never needed; unused Imports trigger
  `R CMD check --as-cran` NOTEs.

- **Suggests**: `rcicr`, `foreach`, `tibble`, `dplyr`, `testthat`,
  `knitr`, `rmarkdown`, `withr`. `rcicr` is only needed for the 2IFC
  path (users who only work on Brief-RC don’t install it).
  `foreach`/`tibble`/`dplyr` are attached at runtime by the `rcicr`
  wrappers because `rcicr` uses `%dopar%`, `tribble`, `%>%` during
  evaluation — see §8.

- **Remotes**: `rdotsch/rcicr` (default branch). **Do not** write
  `rdotsch/rcicr@development` — that branch does not exist on
  `rdotsch/rcicr` and breaks CI with `pkgdepends resolution error`.
  Older READMEs of rcicr may claim otherwise; trust `git ls-remote`.
  Canonical install instructions are Dotsch’s own:

  ``` r
  install.packages("rcicr")                   # stable from CRAN
  devtools::install_github("rdotsch/rcicr")   # dev from GitHub
  ```

  As of 2026 canonical rcicr is v1.0.1 (2023-01-13) by Ron Dotsch,
  GPL-2, with exactly 17 exported functions — **none** of them
  Brief-RC-specific.

- **License**: MIT (same as rcicrdiagnostics). Set up via
  `usethis::use_mit_license("Manuel Oliveira")`, which produces:

  - `LICENSE` — two-line stub: `YEAR: 2026`,
    `COPYRIGHT HOLDER: Manuel Oliveira`. Ships with the package.
  - `LICENSE.md` — full MIT legal text. `usethis` adds `^LICENSE\.md$`
    to `.Rbuildignore` automatically, so the full text is GitHub-visible
    but doesn’t duplicate in the build.
  - DESCRIPTION declares `License: MIT + file LICENSE`. The
    `+ file LICENSE` clause attaches the copyright-holder stub.
  - **Do not** `.Rbuildignore` the `LICENSE` stub — MIT requires it to
    ship.

- **Distribution**: GitHub-only (`olivethree/rcicrely`). No CRAN
  submission planned. Users install via
  `remotes::install_github("olivethree/rcicrely")`.

- No compiled code (no Rcpp/C++). Pure R.

### Provenance warning — canonical rcicr vs forks

During development of `rcicrdiagnostics`, a local package also named
`rcicr v0.1.0` (authored by Manuel Oliveira) was inadvertently installed
and treated as canonical for a while. That fork had Brief-RC additions —
`batchGenerateCI_brief`, `computeInfoVal_brief`,
`generateStimuli_brief`, `generateReferenceDistribution_brief`,
`generate_briefrc_data`, `prep_briefrc_data`, `process_briefrc_cis`,
`simulate_briefrc_responses`, plus `applyMask`, `applyScaling`,
`combine`, `img_pos`, `saveToImage`, `generate_and_assess_cis`,
`generate_and_assess_pixel_importance`, and several `_data` objects.
**None of these exist in canonical rcicr v1.0.1.** The fork was a local
derivative only — collaborators and CI runners see Dotsch’s v1.0.1 and
any code that depends on `_brief` functions fails with
`could not find function "batchGenerateCI_brief"`. Lesson: never assume
any `_brief`/`briefrc` function is canonical rcicr. Verify via
`ls("package:rcicr")` or by inspecting `temp/rcicr-main/NAMESPACE`.

## 3. Package structure

    rcicrely/
    ├── CLAUDE_rcicrely.md                     # this file (gitignored)
    ├── DESCRIPTION
    ├── NAMESPACE                              # roxygen-generated
    ├── LICENSE                                # MIT stub (YEAR + COPYRIGHT HOLDER)
    ├── LICENSE.md                             # full MIT legal text (in .Rbuildignore)
    ├── .Rbuildignore
    ├── .gitignore                             # ignores CLAUDE_rcicrely.md, data-raw/generated/, docs/, .claude/
    ├── README.md
    ├── _pkgdown.yml
    │
    ├── .github/
    │   └── workflows/
    │       ├── R-CMD-check.yaml               # standard r-lib/actions matrix
    │       └── pkgdown.yaml                   # deploys to gh-pages
    │
    ├── R/
    │   ├── rcicrely-package.R                 # _PACKAGE sentinel
    │   ├── utils.R                            # validators, ensure_attached helper
    │   ├── globals.R                          # utils::globalVariables + .datatable.aware <- TRUE
    │   ├── rcicrely_result.R                    # S3 class for metric results
    │   ├── report.R                           # print/summary/plot methods for rcicrely_report
    │   │
    │   ├── ci_io.R                            # read_cis(): read CI PNG/JPEG directories into pixel matrix
    │   ├── noise_extraction.R                 # extract_signal(): subtract base image
    │   ├── ci_from_responses_2ifc.R           # thin wrapper over rcicr::generateCI2IFC()
    │   ├── ci_from_responses_briefrc.R        # native Brief-RC CI (Schmitz et al., 2024)
    │   │
    │   ├── rel_split_half.R                   # permuted split-half + Spearman-Brown
    │   ├── rel_loo.R                          # leave-one-out sensitivity
    │   ├── rel_icc.R                          # ICC(3,1) and ICC(3,k) from mean squares
    │   │
    │   ├── pixel_t_test.R                     # vectorised Welch t per pixel
    │   ├── cluster_utils.R                    # 2D connected components (4-connectivity, BFS)
    │   ├── rel_cluster_test.R                 # cluster-based permutation (max-stat FWER)
    │   ├── rel_dissimilarity.R                # bootstrap correlation + Euclidean distance
    │   │
    │   ├── run_within.R                       # orchestrator: all within-condition metrics
    │   ├── run_between.R                      # orchestrator: all between-condition metrics
    │   ├── plot_methods.R                     # t-maps, cluster maps, difference maps, diagnostics
    │   └── internal_helpers.R                 # small utilities (image dims, progress)
    │
    ├── data-raw/                              # .Rbuildignore-d; data-raw/generated/ also gitignored
    │   ├── README.md                          # how to run the generators
    │   ├── 01_generate_stimuli.R              # rcicr + native Brief-RC stimulus generation
    │   ├── 02_generate_bogus_2ifc.R           # 30-participant simulated 2IFC with known signal patterns
    │   ├── 02_generate_bogus_briefrc12.R      # 30-participant simulated Brief-RC 12
    │   ├── 03_generate_example_cis.R          # runs CI computation, saves small PNGs to inst/extdata/ci_examples/
    │   ├── stimuli/                           # base face goes here (user-supplied)
    │   └── generated/                         # script outputs (gitignored entirely)
    │
    ├── inst/
    │   └── extdata/
    │       └── ci_examples/                   # small (64x64 or 128x128) illustrative CIs only
    │                                          # Do NOT ship 512x512 full-size CIs (too big)
    │                                          # Do NOT ship the 16384 x N noise matrix (too big)
    │
    ├── tests/
    │   ├── testthat.R
    │   └── testthat/
    │       ├── test-ci_io.R
    │       ├── test-noise_extraction.R
    │       ├── test-ci_briefrc.R              # native Brief-RC CI computation
    │       ├── test-rel_split_half.R
    │       ├── test-rel_loo.R
    │       ├── test-rel_icc.R
    │       ├── test-cluster_utils.R           # connected component labeling edge cases
    │       ├── test-rel_cluster_test.R
    │       ├── test-rel_dissimilarity.R
    │       ├── test-run_within.R
    │       ├── test-run_between.R
    │       └── test-2ifc-integration.R        # gated on rcicr availability
    │
    ├── vignettes/
    │   └── tutorial.Rmd                       # one unified article — reliability workflow
    │
    └── man/                                   # roxygen-generated

### What’s *not* in this layout (deliberately)

- **No three separate vignettes.** An earlier draft planned
  `within-condition-reliability.Rmd`, `between-condition-inference.Rmd`,
  `briefrc-workflow.Rmd`. A single `tutorial.Rmd` with chapter-sized
  sections is easier to build on CI (one vignette, one build), easier to
  maintain (no cross-vignette duplication), and maps cleanly to a future
  Quarto ebook. Split later if the tutorial grows past ~1500 lines.
- **No large data in `inst/extdata/`.** A 512×512 CI for 30 participants
  is ~60 MB; a 16384 × 1000 noise matrix is ~130 MB. Shipping those
  inside the installed package is unreasonable. Keep full-size generated
  outputs in `data-raw/generated/` (gitignored).
  `inst/extdata/ci_examples/` ships only small illustrative CIs (64×64
  or 128×128), enough for examples and tests.
- **No tidyverse.** `data.table` for tabular manipulation; `cli` for
  messages. `dplyr`/`tibble` appear only as soft deps attached at
  runtime because `rcicr` needs them (§8).

## 4. Two input modes

The package accepts data at two levels of pre-processing, so users can
skip the CI-computation step if they already have CIs.

### 4.1. Mode 1: pre-computed CI images

User has already generated individual CIs as PNG/JPEG files (one per
participant per condition), stored in directories. Simplest path.

**Required inputs**:

- Directory path(s) to individual CI images.
- Path to the base image used at stimulus-generation time (for signal
  extraction).

**Flow**:
[`read_cis()`](https://olivethree.github.io/rcicrely/reference/read_cis.md)
reads all images from a directory via `png` / `jpeg` (loaded only when
actually called), converts to grayscale using ITU-R BT.709 luminance
weights (`0.2126 R + 0.7152 G + 0.0722 B`), returns a pixels ×
participants numeric matrix.
[`extract_signal()`](https://olivethree.github.io/rcicrely/reference/extract_signal.md)
subtracts the base image vector from each column to isolate the
noise/signal component.

**Mode 1 caveat — PNG pixels are rendered, not raw.** Since v0.1.1,
[`read_cis()`](https://olivethree.github.io/rcicrely/reference/read_cis.md)
/
[`extract_signal()`](https://olivethree.github.io/rcicrely/reference/extract_signal.md)
/
[`load_signal_matrix()`](https://olivethree.github.io/rcicrely/reference/load_signal_matrix.md)
emit a once-per-session `cli` warning explaining that PNGs encode
`base + scaling(mask)`, so the extracted signal is `scaling(mask)`, not
the raw mask. Variance-based reliability metrics (`rel_icc`, Euclidean
half of `rel_dissimilarity`, `pixel_t_test`, `rel_cluster_test`) are
sensitive; only Pearson-based metrics survive a single uniform linear
scaling unmodified, and even those break under per-CI “matched” scaling.
**The canonical 2IFC `infoVal` path is NOT affected**:
[`rcicr::computeInfoVal2IFC()`](https://rdrr.io/pkg/rcicr/man/computeInfoVal2IFC.html)
does `norm(matrix(target_ci[["ci"]]), "f")` internally, always reading
the raw `$ci` from the rcicr CI-list regardless of scaling settings. The
infoVal pitfall only applies to hand-rolled implementations (Brief-RC
has no canonical infoVal; custom 2IFC scripts that re-implement the
norm). Silence with the `acknowledge_scaling = TRUE` arg or
`options(rcicrely.silence_scaling_warning = TRUE)`. Implementation:
session-private env in `R/utils.R` (`warn_mode1_scaling()`,
`reset_session_warnings()`); heuristic backstop in every `rel_*()` /
`run_*()` via `looks_scaled()` in `R/scaling_diagnostics.R`.

A convenience wrapper `load_signal_matrix(dir, base_image_path)` exports
the common path in one call
([`read_cis()`](https://olivethree.github.io/rcicrely/reference/read_cis.md) +
[`extract_signal()`](https://olivethree.github.io/rcicrely/reference/extract_signal.md)
composed);
[`read_cis()`](https://olivethree.github.io/rcicrely/reference/read_cis.md)
and
[`extract_signal()`](https://olivethree.github.io/rcicrely/reference/extract_signal.md)
remain exported for power users who want to intervene between the two
steps (e.g. to mask pixels, crop, or swap the base).

The **signal matrix** (pixels × participants, base-subtracted) is the
canonical internal representation. Every downstream function takes this
matrix.

### 4.2. Mode 2: from raw response data

User has trial-level response data and wants the package to compute
individual CIs internally. Dispatches by `method`:

- **`method = "2ifc"`**: delegate to
  [`rcicr::generateCI2IFC()`](https://rdrr.io/pkg/rcicr/man/generateCI2IFC.html)
  (or
  [`rcicr::batchGenerateCI2IFC()`](https://rdrr.io/pkg/rcicr/man/batchGenerateCI2IFC.html)
  for per-participant) from canonical rcicr v1.0.1. Required inputs:
  response data frame with columns `participant_id`, `stimulus`,
  `response ∈ {-1, +1}`; the `.RData` file from
  [`rcicr::generateStimuli2IFC()`](https://rdrr.io/pkg/rcicr/man/generateStimuli2IFC.html);
  base image label. See §8 for rcicr integration gotchas.

- **`method = "briefrc"`**: compute natively inside rcicrely. Canonical
  rcicr does not ship Brief-RC machinery (§5). Required inputs: response
  data frame with columns `participant_id`, `trial`, `stimulus`,
  `response ∈ {-1, +1}`; the noise matrix; base image path.

  **Noise-matrix sourcing.** Crucially,
  [`rcicr::generateStimuli2IFC()`](https://rdrr.io/pkg/rcicr/man/generateStimuli2IFC.html)
  does **not** save the noise matrix to the `.Rdata` it writes. The
  rdata contains `stimuli_params` (per-trial 4092-wide parameter rows)
  and `p` (the sinusoid / gabor basis), from which individual noise
  patterns are reconstructable via
  [`rcicr::generateNoiseImage()`](https://rdrr.io/pkg/rcicr/man/generateNoiseImage.html).
  The pixels × n_trials matrix proper is only available as the **return
  value** of `generateStimuli2IFC(return_as_dataframe = TRUE)`. The
  rcicrely helper
  [`read_noise_matrix()`](https://olivethree.github.io/rcicrely/reference/read_noise_matrix.md)
  accepts any of: an `.rds` file of a previously-captured return value,
  an rcicr `.Rdata` (it reconstructs via `generateNoiseImage` internally
  — needs rcicr installed), a pre-saved `stimuli` object inside an
  `.Rdata`, or a whitespace-delimited text file. See §8.11.

### 4.3. Unified interface

All analysis functions accept either:

- a pre-built signal matrix (pixels × participants) — already
  base-subtracted;
- a directory path + base image path (Mode 1);
- raw response data + pipeline-specific inputs (Mode 2).

Signal matrix is the canonical internal representation — all roads lead
there.

## 5. Brief-RC: rcicrely owns this implementation

Canonical `rcicr` v1.0.1 has **no Brief-RC-specific functions**. Schmitz
et al. (2024) describe their Brief-RC method by **adapting** rcicr’s
`generateStimuli2IFC()` and `genCI()` — their adapted scripts are on
their OSF repository and in a small companion package (`schmitz`). Those
adaptations were never folded into canonical rcicr.

Consequence: rcicrely **implements Brief-RC natively**, citing Schmitz
et al. (2024) and Brinkman et al. (2019) in each function’s roxygen
`@references`. This is not an accident of architecture — it’s the right
call:

- Canonical rcicr is stable. Relying on non-canonical forks is fragile
  (branches rename, functions get retracted).
- The Brief-RC math is straightforward — we do not need a black box.
  Implementing directly is auditable.
- Under CRAN rules, a Suggests dep must be CRAN-installable (with narrow
  `Remotes:` exceptions). `rcicr::*_brief` cannot be assumed to exist on
  a CRAN reviewer’s machine.

### 5.1. Schmitz et al. (2024) Brief-RC canonical format

This release targets **Brief-RC 12** (12 faces per trial). Other
face-counts described by Schmitz et al. are not supported in v0.1 and
are listed as future work.

Per the paper (§3.1.2, p. 6) and Schmitz’s own `genMask()`
implementation:

- Each Brief-RC 12 trial presents **12 noisy faces: 6 oriented
  (`base + noise_i`) + 6 inverted (`base − noise_i`)**, drawn from a
  shared pool of noise patterns.

- Each noise pattern is one of a pair; each pair has a unique id in
  `1..n_pool`.

- Participant picks **one** face per trial.

- Record format — **one row per trial**:

  | Column           | Value                                             |
  |------------------|---------------------------------------------------|
  | `participant_id` | producer id                                       |
  | `trial`          | 1..n_trials (per participant)                     |
  | `stimulus`       | pool id of the chosen noise pattern (1..n_pool)   |
  | `response`       | `+1` if oriented version chosen, `−1` if inverted |
  | `rt` (optional)  | response time in ms                               |

Unchosen faces are not recorded. **Do not** use the “expanded” format
(one row per shown face with `response = 1` for chosen / `−1/11` for
unchosen) — that’s an artefact of a misread of the paper and not what
Schmitz’s adapted `genCI` expects.

### 5.2. Native Brief-RC CI formula (rcicrely’s implementation)

**rcicrely implements Schmitz’s `genMask()` exactly, including the
`mean(response) by stim` collapse.** This is the shipped reference
implementation:

``` r
# rcicrely's Brief-RC mask (identical to Schmitz's genMask in RC.R):
X    <- data.table::data.table(response = response, stim = stim)
X    <- X[, .(response = mean(response)), stim]
mask <- (noiseMatrix[, X$stim] %*% X$response) / length(X$response)
ci   <- base + scaling(mask)
```

Where `noiseMatrix` has one column per pool id, each column a sinusoidal
(or equivalent) noise pattern of the same dimensions as the base image.

**Duplicate-stimulus handling.** If a participant chooses the same pool
id on more than one trial, `mean(response) by stim` averages the
responses for that id, and the denominator is `length(X$response)` —
i.e. the number of **unique** pool ids chosen by that participant, not
the raw trial count. This matches Schmitz’s own implementation. In
practice duplicates are rare (pool sizes are usually much larger than
`n_trials`), but the tie-breaking rule is part of the spec and unit
tests must cover it.

**Naive intuition (not the shipped formula).** For exposition only, when
there are no duplicate stimulus choices the mask reduces to:

    mask = (1 / n_trials) × Σ_t sign_t × noise_{stim_t}
         = noise_matrix[, chosen_pool_ids] %*% response_signs / n_trials

This is only equivalent to `genMask()` in the no-duplicate case; the two
formulas differ in both numerator and denominator when a participant
chooses the same pool id twice with different responses. **Do not ship
this version.**

Critical implementation details:

- `noise_matrix` can be obtained from an rcicr stimulus generation run
  via `rcicr::generateStimuli2IFC(..., return_as_dataframe = TRUE)` (a
  v1.0.1 canonical parameter). That **returns** a pixels × n_trials data
  frame where each column is one trial’s noise pattern. The return value
  must be captured — rcicr does **not** save it to the `.Rdata` (§8.11).
  The canonical pattern in the rcicrely test harness is: capture the
  return value, save to `noise_matrix.rds`, and pass that to downstream
  scripts. If only the `.Rdata` is available, rcicrely’s
  [`read_noise_matrix()`](https://olivethree.github.io/rcicrely/reference/read_noise_matrix.md)
  reconstructs the noise matrix by calling
  `rcicr::generateNoiseImage(stimuli_params[i, ], p)` for each row —
  correct but slow (~4s per 120 trials at 64², scales linearly).
- Response coding must be exactly `±1`. Fraction of ±1 responses
  represents the signed-selection; no other weighting is applied.
- Scaling of `mask` before combining with the base image follows
  Schmitz’s options: `"matched"` (stretches mask range to base range —
  suboptimal but common), `"constant"` (multiply by a numeric), or
  `"none"` (default in rcicrely). **The shipped `signal_matrix` is
  always the raw mask, regardless of `scaling`.** Since v0.1.1,
  [`ci_from_responses_briefrc()`](https://olivethree.github.io/rcicrely/reference/ci_from_responses_briefrc.md)
  exposes a `scaling` argument and produces a parallel `$rendered_ci`
  field (`base + scaling(mask)`) when `scaling != "none"`. That field is
  for visualisation / saving to disk only — never feed it to `rel_*()`
  or to any external `infoVal` computation. The 2IFC wrapper has a
  symmetric `keep_rendered = TRUE` argument that extracts rcicr’s
  `$combined`.

### 5.3. Brief-RC infoVal: deferred, not faked

InfoVal for Brief-RC needs a **reference distribution matched to each
participant’s trial count** (not the rdata pool size). That is
non-trivial and is out of scope for this package’s initial release.
`rel_*` reliability metrics do **not** depend on infoVal — they operate
directly on the signal matrix. If infoVal is needed, either:

- compute it externally using rcicrdiagnostics’ 2IFC path (works if the
  user is willing to treat their Brief-RC as a pool of signed choices —
  with the trial-count-scale caveat acknowledged); or
- implement it natively here as a future enhancement with a custom
  permutation-based reference matched to per-participant trial count.

rcicrdiagnostics currently returns `"skip"` for `method = "briefrc"` in
its `compute_infoval_summary()`, `check_response_inversion()`, and
`cross_validate_rt_infoval()`. When rcicrely eventually implements
Brief-RC infoVal, rcicrdiagnostics can delegate to it.

## 6. API surface and conventions

### 6.1. Naming

- snake_case, lowercase everywhere.
- Exported analysis functions are **prefixed with `rel_`**:
  [`rel_split_half()`](https://olivethree.github.io/rcicrely/reference/rel_split_half.md),
  [`rel_loo()`](https://olivethree.github.io/rcicrely/reference/rel_loo.md),
  [`rel_icc()`](https://olivethree.github.io/rcicrely/reference/rel_icc.md),
  [`rel_cluster_test()`](https://olivethree.github.io/rcicrely/reference/rel_cluster_test.md),
  [`rel_dissimilarity()`](https://olivethree.github.io/rcicrely/reference/rel_dissimilarity.md).
- Orchestrators:
  [`run_within()`](https://olivethree.github.io/rcicrely/reference/run_within.md),
  [`run_between()`](https://olivethree.github.io/rcicrely/reference/run_between.md).
- I/O and helpers:
  [`read_cis()`](https://olivethree.github.io/rcicrely/reference/read_cis.md),
  [`read_noise_matrix()`](https://olivethree.github.io/rcicrely/reference/read_noise_matrix.md),
  [`extract_signal()`](https://olivethree.github.io/rcicrely/reference/extract_signal.md).
- Plot functions: `plot_tmap()`, `plot_clusters()`, `plot_difference()`,
  `plot_split_half_dist()` — though most are reached through S3
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods on
  result objects.
- Internal helpers: not exported, not prefixed, roxygen
  `@keywords internal` + `@noRd`.

### 6.2. Signature pattern (within-condition)

``` r
rel_split_half(
  signal_matrix,          # pixels × participants, base-subtracted
  n_permutations = 2000,
  seed = NULL,
  progress = TRUE
)

rel_loo(
  signal_matrix,
  flag_threshold_sd = 2.5
)

rel_icc(
  signal_matrix,
  variants = c("3_1", "3_k")   # "2_1", "2_k" also permitted
)
```

`flag_threshold_sd` defaults to 2.5 rather than 2 because a 2-SD rule
flags roughly one participant per dataset by chance at N≈30, which
trains users to ignore the signal. A robust `flag_method = "mad"`
alternative is on the future-work list.

### 6.3. Signature pattern (between-condition)

``` r
rel_cluster_test(
  signal_matrix_a,         # pixels × participants for condition A
  signal_matrix_b,         # pixels × participants for condition B
  img_dims,                # integer c(nrow, ncol)
  n_permutations = 2000,
  cluster_threshold = 2.0,
  alpha = 0.05,
  seed = NULL,
  progress = TRUE
)

rel_dissimilarity(
  signal_matrix_a,
  signal_matrix_b,
  n_boot = 2000,
  ci_level = 0.95,
  seed = NULL,
  progress = TRUE
)
```

**Permutation/bootstrap default is 2000 across the package.** 1000 is on
the low side for stable 95% percentile CIs; 2000 is cheap and gives
Monte Carlo error below 0.01 on tail probabilities. Pick one number and
use it everywhere so users don’t have to remember which metric defaults
to what.

### 6.4. Return type: `rcicrely_*` S3 classes

Every result class gets three methods: `print`, `summary`, `plot`.
Classes:

- `rcicrely_split_half` — `$r_hh` (mean split-half r), `$r_sb`
  (Spearman-Brown corrected), `$ci_95`, `$ci_95_sb`, `$distribution`
  (permutation vector), `$n_participants`, `$n_permutations`.
- `rcicrely_loo` — `$correlations` (named vector), `$mean_r`, `$sd_r`,
  `$threshold`, `$flagged` (character vector), `$summary_df`.
- `rcicrely_icc` — `$icc_3_1`, `$icc_3_k`, `$ms_rows`, `$ms_cols`,
  `$ms_error`, `$n_raters`, `$n_targets`, `$model` (character:
  model-specification description shown by the print method). If the
  user requests ICC(2,\*) via `variants`, `$icc_2_1` and/or `$icc_2_k`
  are added alongside. See §7.3 for the model-selection argument.
- `rcicrely_cluster_test` — `$observed_t` (pixel-wise t vector),
  `$clusters` (data frame with
  `cluster_id, direction, mass, p_value, significant`),
  `$null_distribution` (list with `pos` and `neg` max masses),
  `$img_dims`, `$pos_labels`, `$neg_labels` (for plotting).
- `rcicrely_dissim` — `$correlation`, `$euclidean`, `$boot_cor`,
  `$boot_dist`, `$ci_cor`, `$ci_dist`, `$boot_se_cor`, `$boot_se_dist`.

All orchestrator outputs wrapped in `rcicrely_report` with `$results`
(named list of `rcicrely_*` objects) and `$method`.

### 6.5. Error handling

- [`cli::cli_abort()`](https://cli.r-lib.org/reference/cli_abort.html)
  only for programming errors (wrong arg types, missing files).
- Validate that between-condition signal matrices have the same number
  of rows (pixels). Abort if not.
- Validate that at least 4 participants exist per condition (minimum for
  meaningful split-half and ICC). Abort below that.
- **Warn** (do not abort) if fewer than 30 participants — metrics will
  be noisy but still computable. Cone et al. (2020) recommend N ≥ 60 for
  stable reliability assessment; at N \< 30 the reliability estimates
  themselves become unreliable. The warning message should reference
  Cone et al. (2020) so users know where the threshold came from.
- All functions with `seed` call
  [`set.seed()`](https://rdrr.io/r/base/Random.html) internally and
  restore the RNG state on exit via
  [`withr::with_seed()`](https://withr.r-lib.org/reference/with_seed.html)
  (or manual save and `on.exit`). Do not leak a seed into the caller’s
  global RNG.

## 7. Reliability metric reference

### 7.1. Permuted split-half with Spearman-Brown

`rel_split_half(signal_matrix, n_permutations = 2000, seed = NULL)`.

1.  For each permutation: randomly partition the `N` participants into
    two roughly equal halves. For **odd N**, one participant is randomly
    dropped **per permutation** (not the same participant every time —
    re-drawn each iteration) so both halves contain `floor(N/2)`
    participants. Document this in the help page. Average each half’s
    columns to get two group-level signal vectors. Compute Pearson
    correlation.
2.  Mean `r_hh` = average of per-permutation `r`.
3.  Spearman-Brown correction: `r_SB = (2 × r_hh) / (1 + r_hh)`.
4.  95% CI: 2.5th and 97.5th percentile of the permutation distribution
    (percentile method). Separately for `r` and `r_SB`.

Note: permutation is over **participants**, not trials. Participant
halves are averaged before correlation — the averaging matters
mathematically and is the standard approach in the psychometrics
literature.

### 7.2. Leave-one-out sensitivity

`rel_loo(signal_matrix, flag_threshold_sd = 2.5)`.

For each participant `i`: correlate the full-sample group CI
(`rowMeans(signal_matrix)`) with the group CI computed without
participant `i` (`rowMeans(signal_matrix[, -i])`). Returns per-
participant correlations. Participants whose correlation is more than
`flag_threshold_sd` standard deviations below the mean are flagged as
influential outliers.

### 7.3. ICC(3,1) and ICC(3,k) from mean squares

`rel_icc(signal_matrix, variants = c("3_1", "3_k"))`.

**How this ICC differs from most ICCs in the RC literature.** Most ICCs
published in reverse-correlation papers are **trait-rating reliability**
statistics: phase-2 naive raters score group-level or individual CIs on
trait dimensions (trustworthy, competent, dominant) and the ICC reports
inter-rater agreement on those Likert scores. That design is exactly the
two-phase pipeline Cone et al. (2020) criticised — the reported
reliability inherits any biases introduced by the rating phase (rater
pool, scale anchors, contextual framing).

rcicrely’s ICC is a **different statistical object**. It operates
directly on the pixel-level signal produced by the original producers,
so it sidesteps phase-2 rating entirely. Treat this as a feature, not a
deviation: the whole package is about quantifying reliability without
relying on downstream raters.

\*\*Model choice — ICC(3,\*) by positive argument, not convention.\*\*
Map the signal matrix: pixels = targets, participants = raters. Then:

- \*\*ICC(1,\*) **(one-way random, no crossed design):** rejected.\*\*
  Our design is crossed — every participant’s signal is evaluated at
  every pixel. Using ICC(1,\*) ignores the crossed variance structure.
- \*\*ICC(2,\*) **(two-way random; pixels and raters both random):**
  rejected.\*\* Pixels are a fixed `img_size × img_size` grid, not a
  random sample from a population of pixels. Users are not generalising
  to “other pixels they might have measured”; they are using every pixel
  in the image. The random-targets assumption is a category error, even
  when ICC(2,*) and ICC(3,*) happen to give similar numbers (see note
  below).
- **ICC(3,1)**, **ICC(3,k)** (two-way mixed; pixels fixed, participants
  random): **correct specification.** Participants *are* a random sample
  from a plausible-producer population. Pixels *are* fixed by the
  experimental image.

**Why ship both ICC(3,1) and ICC(3,k).** They answer different
questions:

- **ICC(3,1)** — reliability of a **single participant’s** CI signal as
  an estimate of the group-level pattern. Answers *how informative is
  one individual CI?*
- **ICC(3,k)** — reliability of the **group-mean CI** across `k`
  participants. Answers *how stable is the group CI this experiment
  produced?* This is usually the headline number in a paper.

**Numerical proximity caveat.** At `n_pixels` ≈ 262,144 and
`n_participants` ≈ 30, ICC(2,*) and ICC(3,*) typically differ by \< 0.01
because pixel and rater main-effect variance components are tiny
relative to the interaction term. The model choice still matters for
interpretation: ICC(3,*) reports a defensible population-generalisation
statement; ICC(2,*) reports a nonsensical one that happens to give a
similar number.

**Compatibility option.** `variants` permits `"2_1"` and `"2_k"` for
users whose reviewers explicitly ask for the ICC(2,\*) form. Default is
`c("3_1", "3_k")` and the print method states the shipped model.

**Implementation — direct mean squares, not
[`psych::ICC()`](https://rdrr.io/pkg/psych/man/ICC.html).**
Non-negotiable: [`psych::ICC()`](https://rdrr.io/pkg/psych/man/ICC.html)
on a 262,144 × 30 matrix allocates an intermediate that blows memory.

``` r
# signal_matrix is pixels (targets) × participants (raters)
n <- nrow(signal_matrix)   # n targets (pixels) — fixed
k <- ncol(signal_matrix)   # k raters (participants) — random

row_means  <- rowMeans(signal_matrix)
col_means  <- colMeans(signal_matrix)
grand_mean <- mean(signal_matrix)

ss_rows  <- k * sum((row_means - grand_mean)^2)
ss_cols  <- n * sum((col_means - grand_mean)^2)
ss_total <- sum((signal_matrix - grand_mean)^2)
ss_error <- ss_total - ss_rows - ss_cols

ms_rows  <- ss_rows  / (n - 1)
ms_cols  <- ss_cols  / (k - 1)
ms_error <- ss_error / ((n - 1) * (k - 1))

# ICC(3,1), ICC(3,k) — two-way mixed, pixels fixed, raters random
icc_3_1  <- (ms_rows - ms_error) /
            (ms_rows + (k - 1) * ms_error)

icc_3_k  <- (ms_rows - ms_error) / ms_rows

# ICC(2,1), ICC(2,k) — exposed only via variants = c("2_1", "2_k")
icc_2_1  <- (ms_rows - ms_error) /
            (ms_rows + (k - 1) * ms_error +
             (k / n) * (ms_cols - ms_error))

icc_2_k  <- (ms_rows - ms_error) /
            (ms_rows + (ms_cols - ms_error) / n)
```

**Return shape**: `$icc_3_1`, `$icc_3_k`, `$ms_rows`, `$ms_cols`,
`$ms_error`, `$n_raters`, `$n_targets`, `$model` (character scalar:
e.g. `"ICC(3,*) / two-way mixed; pixels fixed, participants random"`).
If `variants` includes `"2_1"` or `"2_k"`, add `$icc_2_1` / `$icc_2_k`.

**References.** Shrout & Fleiss (1979) for the original six-ICC
taxonomy; McGraw & Wong (1996) for the modern re-derivation and the
random/fixed-effects clarification. Cite both in `@references` of
[`rel_icc()`](https://olivethree.github.io/rcicrely/reference/rel_icc.md).
Test against [`psych::ICC()`](https://rdrr.io/pkg/psych/man/ICC.html) on
a small-matrix case (≤ 50 × 10) during development — numbers must match
to 6 decimals.

### 7.4. Pixel-wise Welch t-test

`pixel_t_test(signal_matrix_a, signal_matrix_b)`.

Vectorised Welch’s t (unequal variances) per pixel. **Not** Student’s t
— conditions may have different N and different variances. Returns a
vector of t-values, length = n_pixels.

### 7.5. Cluster-based permutation test

`rel_cluster_test(signal_matrix_a, signal_matrix_b, img_dims, ...)`.

Implementation:

1.  Compute observed pixel-wise Welch t (§7.4).
2.  Build the observed cluster map: threshold `|t| > cluster_threshold`
    (default 2.0), separately for positive and negative clusters. The
    2.0 default approximates p \< 0.05 for large df, but **for small
    samples** (N_A + N_B \< 20) it under-rejects; users should inspect
    the observed-t histogram before trusting the result, and consider
    raising the threshold or increasing N.
3.  Find connected components with **4-connectivity** (not
    8-connectivity — the more conservative choice; document this).
4.  For each cluster, **cluster mass = sum of t-values within the
    cluster** (not count of pixels). This captures both extent and
    magnitude.
5.  Null distribution: for each of `n_permutations` permutations,
    randomly shuffle condition labels across **participants** (not
    pixels — never pixels; preserving spatial correlation in each
    participant’s CI is the point of the test). Permutation is
    **stratified**: every permutation preserves `(N_A, N_B)` exactly.
    Unstratified permutation would let the two groups’ sizes drift,
    changing the Welch degrees of freedom under the null and biasing the
    max-mass distribution. Recompute pixel-wise Welch t on shuffled
    labels, find clusters, record the **maximum positive and maximum
    negative cluster mass**.
6.  For each observed cluster: p-value is the fraction of permutation
    null-distribution masses (matching direction) that exceed its mass.
    Maximum-statistic approach gives FWER control.

Reference: Maris & Oostenveld (2007) for the cluster-based permutation
approach; Nichols & Holmes (2002) for max-stat FWER.

Implementation notes:

- Connected component labeling: BFS over a thresholded binary matrix.
  Pure R, no compiled code. Pre-allocate a label matrix and a queue as
  integer vectors with head/tail pointers. Process positive and negative
  masks separately.
- **Memory management**: never copy the full signal matrices inside the
  permutation loop — only permute integer column indices and recompute t
  on the shuffled view.

### 7.6. Representational dissimilarity with bootstrap CIs

`rel_dissimilarity(signal_matrix_a, signal_matrix_b, n_boot = 2000)`.

Compute two scalar dissimilarity metrics between the two conditions’
group-level signals:

- Pearson correlation of group CIs.
- Euclidean distance of group CIs.

Bootstrap both: resample participants with replacement within each
condition, recompute, get a bootstrap distribution. 95% CI via the
**percentile method** — `quantile(boot_stat, c(0.025, 0.975))`, no
external dependency. BCa is on the future-work list; it would require
either the `boot` package or a hand-rolled jackknife-plus-z0
computation, neither of which is justified for the first release.

## 8. rcicr integration — gotcha section

These are the issues that cost real time during the rcicrdiagnostics
build. They apply equally to rcicrely’s 2IFC path.

### 8.1. No `@development` branch on `rdotsch/rcicr`

Older rcicr READMEs referred to a `development` branch. It does not
exist on `rdotsch/rcicr`. Verify via `git ls-remote`. If DESCRIPTION has
`Remotes: rdotsch/rcicr@development`, CI fails with
`pkgdepends resolution error: Can't find reference @development`. Use
`Remotes: rdotsch/rcicr` (no ref) so the default branch resolves.

### 8.2. `base_face_files` must be a `list`, not a named vector

`rcicr::generateStimuli2IFC(base_face_files = c(base = path))` errors
with `Please provide base face file name as named list`. Use
`list(base = path)`.

### 8.3. `%dopar%` is not re-exported by rcicr

[`rcicr::generateStimuli2IFC()`](https://rdrr.io/pkg/rcicr/man/generateStimuli2IFC.html)
uses `foreach::foreach(...) %dopar%`. `%dopar%` is an infix that R’s
lookup resolves via the search path, not via `foreach::`. rcicr’s
NAMESPACE does not re-export it. Attach `foreach` before calling rcicr
from a package or data-raw script:

``` r
# Inside package code (not data-raw): use attachNamespace via ensure_attached helper
ensure_attached(c("foreach", "tibble", "dplyr"))
```

Error symptom: `could not find function "%dopar%"`.

### 8.4. `ncores = 1L` as the cross-platform safe default

Default
`rcicr::generateStimuli2IFC(ncores = parallel::detectCores() - 1)`
spawns 7+ PSOCK workers. On macOS this races with socket teardown badly
enough that workers die with
`Error in unserialize(node$con): error reading from connection`, and the
function aborts before `parallel::stopCluster(cl)` and the post-loop
[`save()`](https://rdrr.io/r/base/save.html), leaving PNGs written but
no `.RData`.

**Fix**: rcicrely’s 2IFC wrapper defaults `ncores = 1L` for
cross-platform safety. The macOS PSOCK race is the specific motivator;
Linux workers are generally fine at higher counts but we pick the safe
default so users don’t hit the macOS failure mode unexpectedly and so
stimulus-generation runs are reproducible across machines. The arg is
**exposed** — Linux users who want parallelism can pass e.g.
`ncores = 4L` explicitly and accept the macOS caveat.

**Scope**: this race applies only to
[`rcicr::generateStimuli2IFC()`](https://rdrr.io/pkg/rcicr/man/generateStimuli2IFC.html).
The **CI-generation** path (`batchGenerateCI2IFC()`) does not spawn a
cluster unless the `participants` argument is supplied (which rcicrely’s
wrapper does not do — we pass the data unpartitioned and rely on rcicr’s
internal per-unit loop). So no `ncores` plumbing is needed on the
CI-generation side.

### 8.5. `computeInfoVal2IFC()` expects the whole CI list, not `$ci`

[`rcicr::batchGenerateCI2IFC()`](https://rdrr.io/pkg/rcicr/man/batchGenerateCI2IFC.html)
returns a named list where each element is a list with `$ci`, `$scaled`,
`$base`, `$combined`. Pass `cis[[nm]]` (the whole thing) to
[`rcicr::computeInfoVal2IFC()`](https://rdrr.io/pkg/rcicr/man/computeInfoVal2IFC.html),
not `cis[[nm]]$ci`. The function internally does `target_ci[["ci"]]`;
passing `$ci` directly errors with `subscript out of bounds`.

### 8.6. `computeInfoVal2IFC()` mutates the rdata file — rcicrely copies defensively

First call generates a 10 000-iteration reference distribution and
**writes it back into the supplied `.RData` file** (adds a
`reference_norms` object). By rcicr’s design, but surprising.

**rcicrely’s 2IFC wrapper copies the user-supplied `.RData` to a temp
path (`tempfile(fileext = ".RData")`) before calling
[`rcicr::computeInfoVal2IFC()`](https://rdrr.io/pkg/rcicr/man/computeInfoVal2IFC.html)**,
so the user’s file is never mutated. The `reference_norms` addition
lives on the temp copy; if the user wants to persist it, they save the
temp path output to a location of their choosing. Document this in the
wrapper’s help page: “the original rdata file is not modified; rcicr’s
mutation-in-place is intercepted by operating on a temp copy.”

### 8.7. rcicr’s soft deps must be attached, not just loaded

`computeInfoVal2IFC()` uses
[`tibble::tribble`](https://tibble.tidyverse.org/reference/tribble.html),
[`dplyr::filter`](https://dplyr.tidyverse.org/reference/filter.html),
`%>%` at evaluation time. Without these attached you hit
`could not find function "tribble"` or
`match.arg(method): 'arg' must be NULL or a character vector` (the
second is dplyr’s `filter.default` failing because
[`dplyr::filter`](https://dplyr.tidyverse.org/reference/filter.html)
isn’t dispatched). Package code must use
[`attachNamespace()`](https://rdrr.io/r/base/ns-load.html) via an
`ensure_attached()` helper —
[`library()`](https://rdrr.io/r/base/library.html) in package code is
flagged by R CMD check.

### 8.8. `.Rdata` extension case matters on case-sensitive filesystems

[`rcicr::generateStimuli2IFC()`](https://rdrr.io/pkg/rcicr/man/generateStimuli2IFC.html)
writes its timestamped file with lowercase `.Rdata` (one cap, `Rdata`).
macOS `find -name "*.RData"` and R’s `list.files(pattern = "\\.RData$")`
are case-sensitive by default. Always use `ignore.case = TRUE` when
searching.

### 8.9. Base image: square, grayscale, `img_size × img_size`

[`rcicr::generateStimuli2IFC()`](https://rdrr.io/pkg/rcicr/man/generateStimuli2IFC.html)
does not auto-resize or convert. Passing an RGB 400×500 JPEG produces
`Error in { : task 1 failed - "non-conformable arrays"`. Add a pure-R
`normalise_image()` helper in any data-raw script that calls rcicr: read
via
[`jpeg::readJPEG`](https://rdrr.io/pkg/jpeg/man/readJPEG.html)/[`png::readPNG`](https://rdrr.io/pkg/png/man/readPNG.html),
collapse RGB to luminance (`0.2126 R + 0.7152 G + 0.0722 B` or the
similar BT.601 weights), center-crop to square, nearest-neighbour resize
to `img_size`, write out as grayscale JPEG to a temp file, pass that
temp file to rcicr. Leave the user’s original untouched.

### 8.10. For rcicrely: do **not** call `rcicr::*_brief` functions

They do not exist in canonical rcicr v1.0.1. Any code that tries to will
fail with `could not find function "batchGenerateCI_brief"` on a clean
`remotes::install_github("rdotsch/rcicr")`. Own Brief-RC natively (§5).

### 8.11. rcicr’s `.Rdata` does **not** contain the noise matrix

A trap that cost real time. An earlier draft of §4.2 / §5.2 implied that
the `.Rdata` from
[`rcicr::generateStimuli2IFC()`](https://rdrr.io/pkg/rcicr/man/generateStimuli2IFC.html)
stored a pixels × n_trials `stimuli` object. **It does not.** Inspecting
the rdata via [`ls()`](https://rdrr.io/r/base/ls.html) shows:

    base_face_files, base_faces, generator_version, img_size, label,
    n_trials, noise_type, p, seed, stimuli_params, stimulus_path, trial,
    use_same_parameters

No `stimuli` object. `stimuli_params[[baseimage]]` is an n_trials × 4092
matrix of noise parameters; `p` is the sinusoid / gabor basis. The
actual noise matrix has to come from one of three places:

1.  **Captured return value** of
    `generateStimuli2IFC(return_as_dataframe = TRUE)` (what the rcicrely
    test harness does — save it to `noise_matrix.rds` on disk).
2.  **Reconstructed on demand** via
    `rcicr::generateNoiseImage(stimuli_params[i, ], p)` per row.
    [`rcicrely::read_noise_matrix()`](https://olivethree.github.io/rcicrely/reference/read_noise_matrix.md)
    does this automatically when given an rcicr rdata — but requires
    rcicr to be installed, and is slow (linear in n_trials × n_pixels).
3.  **User-saved** `stimuli` object (rare, but supported if present in
    the rdata).

Consequence for test-harness scripts: always capture
`generateStimuli2IFC()`’s return value; don’t rely on the rdata to
provide the noise matrix unless you’re prepared to reconstruct. And
always cache the captured return value to an `.rds` once, so downstream
scripts can load it in ms rather than seconds.

Secondary trap: **stale `.Rdata` accumulates** if you run the stimulus
script more than once, because rcicr timestamps its filename.
`list.files(pattern = "\\.rdata$")` returns all of them, and
`load(length_2_vector, envir = env)` errors with
`invalid 'description' argument`. The fix is either to wipe previous-run
outputs on entry, or to always take the most recent by mtime
(`which.max(file.info(paths)$mtime)`).

## 9. data.table inside a package

### 9.1. `.datatable.aware <- TRUE` is required

A package that uses `[.data.table` but only `data.table::`-prefixes (no
`import(data.table)` in NAMESPACE) is not recognised as
“data.table-aware” by `cedta`. `[.data.table` falls back to
`[.data.frame`, and symbols like `.N`, `.SD` resolve as regular
variables — producing `object '.N' not found` at runtime.

Fix: create `R/globals.R` containing both:

``` r
utils::globalVariables(c(".N", ".SD", "response"))
.datatable.aware <- TRUE
```

**Column names used inside data.table aggregations** need to be in the
[`globalVariables()`](https://rdrr.io/r/utils/globalVariables.html) list
too — R CMD check parses
`X[, list(response = mean(response)), by = "stim"]` and can’t tell that
`response` is a column being aggregated, not an undeclared global. Every
such name goes in. rcicrely currently declares `.N, .SD, response`; if a
new aggregation introduces a new column name, extend the vector.

**Don’t put `utils` in Imports.**
[`utils::globalVariables()`](https://rdrr.io/r/utils/globalVariables.html)
works at load-time without declaring `utils` in Imports — and R CMD
check will emit a NOTE (“Namespace in Imports field not imported from”)
if you do declare it, because there’s no runtime call to a `utils::`
function. The shipped rcicrely DESCRIPTION lists only `cli`,
`data.table`, `grDevices`, `graphics`, `stats`, `tools` under Imports.

**Don’t declare `"n"` (or other common short names) globally.** Any
aggregation that would produce a column called `n` should be renamed at
the site (`n_obs`, `n_trials`, `n_pixels`) — declaring `"n"` in
`globalVariables` silences R CMD check but also silences real typos.
`.N` and `.SD` are specific `data.table` sentinels; there’s no collision
risk.

`globalVariables` silences R CMD check NOTEs about unbound globals
(static analysis); `.datatable.aware <- TRUE` fixes the runtime
dispatch. Both are needed.

### 9.2. `fread` returns integer for whole-number columns

[`data.table::fread()`](https://rdrr.io/pkg/data.table/man/fread.html)
reads a column of `{-1, 1}` as integer, not numeric. Strict type
comparisons like `identical(uniq, c(-1, 1))` silently return FALSE
because `c(-1L, 1L)` is integer and `c(-1, 1)` is double. Coerce to
numeric before identity checks:

``` r
finite <- as.numeric(vals[is.finite(vals)])
uniq <- sort(unique(finite))
```

Include regression tests with `sample(c(-1L, 1L), ...)` alongside
`sample(c(-1, 1), ...)` to catch this.

### 9.3. Use `data.table`, not dplyr

Project convention.
[`data.table::fread`](https://rdrr.io/pkg/data.table/man/fread.html) for
CSVs, `dt[, list(...), by = ...]` for grouping, `.SD`/`.SDcols` for
dynamic columns. dplyr is only pulled in as a runtime dep of rcicr
(§8.7), never by rcicrely’s own code.

## 10. R CMD check hygiene

Target: **0 errors, 0 warnings, 0 notes** with `--as-cran`.

- **MIT license**: ship both `LICENSE` (stub) and `LICENSE.md` (full
  text). `License: MIT + file LICENSE` in DESCRIPTION. Add
  `^LICENSE\.md$` to `.Rbuildignore` (usethis does this automatically);
  do **not** `.Rbuildignore` the `LICENSE` stub.
- **Unused Imports**: only declare what the shipped code calls. Dropping
  any Import that appears unused to R CMD check clears a NOTE. In
  particular, **`boot` is not an Import** —
  [`rel_dissimilarity()`](https://olivethree.github.io/rcicrely/reference/rel_dissimilarity.md)
  uses base R [`quantile()`](https://rdrr.io/r/stats/quantile.html) for
  percentile CIs; there is no `boot::` call in the shipped code.
- **`.N`/`.SD` unbound globals**: declare in `R/globals.R` (§9.1).
- **Cross-references in roxygen**: on a fresh scaffold, roxygen
  `@seealso` / `[fn()]` links can point at functions that haven’t yet
  been `@export`ed into NAMESPACE, and the first `devtools::document()`
  warns about them. Two fixes: ensure every target of a `[fn()]` link
  carries its `@export` tag on the first pass, or run
  `devtools::document()` twice so the second pass sees the
  fully-populated NAMESPACE. The first fix is the real one; the second
  is a reliable workaround.
- **Vignette build during check**: requires pandoc. On macOS point the
  env var at Quarto’s bundled pandoc:
  `RSTUDIO_PANDOC=/Applications/quarto/bin/tools R -e 'devtools::check()'`.
- **No [`library()`](https://rdrr.io/r/base/library.html) in package
  code**. Use `requireNamespace` + `attachNamespace` via
  `ensure_attached()` for the rcicr soft-dep attach.
  [`library()`](https://rdrr.io/r/base/library.html) inside a package
  triggers `Dependence on R services not declared`.
- **Future file timestamps NOTE**: a benign NOTE when the check runner
  can’t reach worldclockapi. Ignore; it vanishes on CI runners with
  normal network access.
- **ASCII-only R source**. `R CMD check --as-cran` raises a WARNING for
  any non-ASCII character in `R/*.R` or NAMESPACE (comments included —
  the exception in the WARNING text is theoretical, in practice it
  fires). During the rcicrely build, em-dashes (—), the section sign
  (§), times (×), approximately (≈), unicode minus (−), and right arrow
  (→) all showed up from paste-ins of this CLAUDE doc. Remedy: run
  [`tools::showNonASCIIfile()`](https://rdrr.io/r/tools/showNonASCII.html)
  on each `R/*.R` to locate, then replace with ASCII substitutes. Under
  locale-constrained environments
  [`gsub()`](https://rdrr.io/r/base/grep.html) on multibyte literals can
  itself fail (`'pattern' is invalid`); wrap the cleanup script in
  `LC_ALL=en_US.UTF-8`.
- **Roxygen `[x, y]` parsed as a link target**. Writing
  `"... values in [0, 1]"` inside a `#'` block yields a warning about an
  unresolved link topic named `"0, 1"`. Write “values in the 0-1 range”
  (or escape the brackets) instead.
- **pkgdown picks up `CLAUDE.md`** as a top-level HTML page unless it’s
  in `.Rbuildignore` / ignored by pkgdown. `CLAUDE.md` is
  developer-facing and gitignored; `.Rbuildignore` it explicitly
  (`^CLAUDE\.md$`) so the pkgdown site doesn’t publish it.

## 11. Statistical commitments

Each of these is a concrete choice. If the paper deviates from one,
update both the code and this section; do not drift silently.

- **Signal extraction**: `CI_signal = CI_pixels - base_image_pixels`.
  All reliability computations operate on the signal, never on the raw
  CI pixels (which include the shared base image and inflate
  correlations). **Caveat (raw vs rendered, v0.1.1):** when `CI_pixels`
  comes from a rendered PNG (Mode 1), it is `base + scaling(mask)`, so
  `CI_signal` is `scaling(mask)`, not the raw mask. Variance-based
  reliability metrics (ICC, Euclidean dissimilarity, pixel t, cluster
  mass) are sensitive; only Pearson-based metrics survive a single
  uniform linear scaling, and even those break under per-CI “matched”
  scaling. The canonical 2IFC
  [`rcicr::computeInfoVal2IFC()`](https://rdrr.io/pkg/rcicr/man/computeInfoVal2IFC.html)
  is unaffected — its source does
  `norm(matrix(target_ci[["ci"]]), "f")`, always reading the raw `$ci`
  from the rcicr CI-list. Only hand-rolled `infoVal` implementations
  (Brief-RC; custom 2IFC code that re-implements the norm) need the raw
  mask explicitly. Mode 2 (`ci_from_responses_*`) returns the raw mask
  directly. The package emits a one-time per-session warning at the
  Mode-1 input boundary and a heuristic backstop at every `rel_*` entry.
- **Brief-RC mask**: Schmitz’s `genMask()` exactly — `mean(response)`
  per unique chosen `stim`, denominator = number of unique chosen pool
  ids (not raw trial count). See §5.2 for the duplicate-stimulus
  tie-breaking rule. Schmitz et al. (2024) §3.1.2.
- **Split-half reliability**: participants are permuted; halves are
  averaged before correlation; for odd N one participant is randomly
  dropped per permutation (re-drawn each iteration). Spearman-Brown
  correction `r_SB = (2 r_hh) / (1 + r_hh)`. Brinkman et al. (2017) for
  RC-specific practice; Shrout & Fleiss (1979) for psychometric
  foundation.
- **ICC**: ICC(3,1) and ICC(3,k) via mean-squares — two-way mixed model,
  pixels fixed, participants random. Not ICC(2,\*) (would require pixels
  to be a random sample from a pixel population, which they are not) and
  not [`psych::ICC()`](https://rdrr.io/pkg/psych/man/ICC.html) (memory).
  See §7.3 for the full model-selection argument. Shrout & Fleiss
  (1979); McGraw & Wong (1996).
- **Pixel-wise t**: Welch’s (unequal variances), not Student’s.
- **Cluster test permutation**: condition labels, across participants,
  **stratified** — every permutation preserves `(N_A, N_B)`. Never pixel
  labels.
- **Cluster connectivity**: 4-connectivity (more conservative than
  8-connectivity). Document this choice in the help pages of
  [`rel_cluster_test()`](https://olivethree.github.io/rcicrely/reference/rel_cluster_test.md)
  and
  [`run_between()`](https://olivethree.github.io/rcicrely/reference/run_between.md)
  (the only two exposures).
- **Cluster mass**: sum of t-values within the cluster, not pixel count.
- **FWER control**: max-statistic across permutations. Maris &
  Oostenveld (2007); Nichols & Holmes (2002).
- **Bootstrap CIs**: percentile method (2.5th / 97.5th) via base R
  [`quantile()`](https://rdrr.io/r/stats/quantile.html). No `boot`
  dependency. BCa is future work.

When in doubt, test the implementation on synthetic data with known
ground truth before trusting the output on real data.

## 12. Bogus data generation

Scripts in `data-raw/` produce test datasets that exercise every
reliability metric across pass/warn/fail scenarios. All outputs go to
`data-raw/generated/` (gitignored).

### 12.1. `01_generate_stimuli.R`

- Normalises the user’s base face (pure R, §8.9).
- Calls
  [`rcicr::generateStimuli2IFC()`](https://rdrr.io/pkg/rcicr/man/generateStimuli2IFC.html)
  with `ncores = 1L` (§8.4), `return_as_dataframe = TRUE`, and the
  user’s normalised image. Saves the `.RData` under a canonical name.
- Uses the same noise matrix for both 2IFC bogus data and Brief-RC 12
  bogus data (they draw from a shared pool).

### 12.2. `02_generate_bogus_2ifc.R`

**Design** (canonical spec — keep one version, don’t drift):

- 30 producers, each contributing sessions under **both** target-
  pattern conditions (A vs B, different latent signals — e.g. A
  emphasises eyes, B emphasises mouth). Total sessions = 30 × 2 = 60.
- Producers are split into three **trial-count levels** (300 / 500 /
  1000 trials per session), 10 producers per level. Each of those 10
  producers runs at that trial count in both A and B.
- Within each (trial-count × condition) cell of 10 producers, the
  **archetype mix** is 7 good / 2 weak / 1 random.
- Archetypes:
  - **“good signal”**: moderate bias toward the condition’s target
    latent pattern, SNR ≈ 0.3.
  - **“weak signal”**: lower SNR (≈ 0.1); borderline metrics expected.
  - **“random”**: coin-flip responses; noise only.

Expected within-condition results (at trial count ≈ 500, per condition):
split-half r ≈ 0.4–0.6; LOO without dramatic outliers at default 2.5-SD
threshold; ICC(3,1) ≈ 0.1–0.3; ICC(3,k) well above 0.5 (averaging boosts
reliability).

Expected between-condition results:
[`rel_cluster_test()`](https://olivethree.github.io/rcicrely/reference/rel_cluster_test.md)
finds significant clusters in the face regions that differ between the
two target patterns.
[`rel_dissimilarity()`](https://olivethree.github.io/rcicrely/reference/rel_dissimilarity.md)
returns a moderate Euclidean distance and a moderate positive
correlation.

### 12.3. `02_generate_bogus_briefrc12.R`

Same producer layout, but with Brief-RC 12 trial structure. Each trial
samples 12 distinct pool ids, randomly assigns 6 as oriented and 6 as
inverted, participant picks one according to archetype. Output is the
per-trial format (§5.1).

### 12.4. `03_generate_example_cis.R`

Loads the bogus response data, computes individual CIs (via rcicr for
2IFC, natively for Brief-RC), writes small (64×64 or 128×128)
illustrative PNGs to `inst/extdata/ci_examples/` for use in vignette
examples and tests. **Do not** write 512×512 CIs here — too big to ship.

## 13. Vignette approach

### 13.1. Rmd, not Qmd

`knitr` is already available to anyone building an R package. Quarto
would require the Quarto CLI on every CI runner — extra dependency.
Convert to `.qmd` later if an ebook version is wanted; Rmd→Qmd is
near-mechanical (same chunk syntax, same YAML).

### 13.2. Structure

One comprehensive `tutorial.Rmd` with chapter-sized sections:

1.  What the package is and is not.
2.  Installation and requirements.
3.  Data modes (pre-computed CIs vs raw responses).
4.  Within-condition reliability (split-half, LOO, ICC).
5.  Between-condition inference (cluster permutation, dissimilarity).
6.  Brief-RC workflow end-to-end.
7.  Interpreting results and common pitfalls.
8.  References.

### 13.3. `eval = FALSE` for rcicr-dependent chunks, cached outputs read back in

Chunks that call
[`rcicr::generateStimuli2IFC`](https://rdrr.io/pkg/rcicr/man/generateStimuli2IFC.html)
or
[`rcicr::batchGenerateCI2IFC`](https://rdrr.io/pkg/rcicr/man/batchGenerateCI2IFC.html):
set `eval = FALSE`. Requires rcicr + real `.RData` + minutes of compute,
which doesn’t work reliably on CI. Show the pattern without executing.

**But the vignette still shows results.** Small CIs, summary tables and
example plots are pre-computed by `data-raw/03_generate_example_cis.R`
and cached in `inst/extdata/ci_examples/`. The vignette reads from that
cache in separate chunks (with `eval = TRUE`) and renders the outputs,
so readers see a complete end-to-end example even though the expensive
compute chunk didn’t run at knit time. Pair every `eval = FALSE` rcicr
chunk with a following `eval = TRUE` chunk that loads and displays the
cached result.

### 13.4. `VignetteBuilder: knitr` in DESCRIPTION

Required. Without it the vignette is not installed by
`R CMD INSTALL --build-vignettes`.

## 14. GitHub deployment

### 14.1. Source-of-truth repo

`github.com/olivethree/rcicrely`. Default branch `main`. Users install
via:

``` r
remotes::install_github("olivethree/rcicrely")
# or:
pak::pak("olivethree/rcicrely")
```

### 14.2. CI workflows

Two workflows in `.github/workflows/`:

1.  **`R-CMD-check.yaml`** — standard r-lib/actions matrix (macOS,
    Windows, Ubuntu; R release). Use
    `usethis::use_github_action_check_standard()` to generate.
2.  **`pkgdown.yaml`** — builds the site and deploys to `gh-pages` via
    `JamesIves/github-pages-deploy-action`. Use
    `usethis::use_pkgdown_github_pages()` to set up.

### 14.3. Repo settings that matter

Two settings the workflow files themselves can’t set:

- **Settings → Actions → General → Workflow permissions** = “Read and
  write permissions”. Otherwise the pkgdown job 403s on push.
- **Settings → Pages → Source** = “Deploy from a branch”, branch
  `gh-pages`, folder `/ (root)`. `usethis::use_pkgdown_github_pages()`
  sets this via API if your PAT has `repo` scope.

Both can be set either via the web UI or scripted with `gh`:

``` bash
# Workflow permissions
gh api repos/olivethree/rcicrely \
  -X PATCH \
  -F 'default_workflow_permissions=write' \
  -F 'can_approve_pull_request_reviews=true'

# Pages source (once the gh-pages branch exists)
gh api repos/olivethree/rcicrely/pages \
  -X POST \
  -F 'source[branch]=gh-pages' \
  -F 'source[path]=/'
```

Use whichever fits the setup workflow; the `gh` commands are
reproducible, which matters if you end up rebuilding the repo.

First pkgdown build takes 7–10 minutes (cold cache, compiles rcicr from
source). Subsequent ~2 minutes.

### 14.4. PAT setup

`usethis::use_pkgdown_github_pages()`, `devtools::install()`,
`remotes::install_github()` all hit the GitHub API and need a PAT. Store
it in the gitcreds keychain, **not** as `GITHUB_PAT` / `GITHUB_TOKEN` /
`GH_TOKEN` env vars (those override the keychain and are the usual
source of 401s):

``` r
usethis::create_github_token()
gitcreds::gitcreds_set()
gh::gh_whoami()                # verify
```

GitHub Apps are not the right answer — the R toolchain is built around
classic PATs.

## 15. Build order for a from-scratch rebuild

If starting from empty directory + this file, build in this order so
each step has what it needs:

1.  **Metadata**: `DESCRIPTION`, `NAMESPACE` (empty initially),
    `.Rbuildignore` (minimum viable set: `^rcicrely\.Rproj$`,
    `^\.Rproj\.user$`, `^LICENSE\.md$`, `^temp$`, `^data-raw$`,
    `^CLAUDE\.md$`, `^CLAUDE_rcicrely\.md$`, `^\.claude$`,
    `^_pkgdown\.yml$`, `^docs$`, `^pkgdown$`, `^\.github$`,
    `^cran-comments\.md$` — the `_pkgdown`/`docs`/`.github`/`CLAUDE.md`
    entries matter because pkgdown auto-discovers top-level markdown and
    R CMD check treats the generated pkgdown site as non-standard files
    otherwise), `.gitignore` (with `CLAUDE.md`, `.claude/`,
    `data-raw/generated/`, `docs/`, `inst/doc/`, `.Rproj.user`, `temp/`,
    `.vscode/`, `.Rhistory`, `.RData`, `.Ruserdata`, `.DS_Store`),
    `LICENSE` (MIT two-line stub), `LICENSE.md` (full MIT text).
    Fastest: scaffold then run
    `usethis::use_mit_license("Manuel Oliveira")`.
2.  **Core infrastructure**: `R/globals.R` (`globalVariables` +
    `.datatable.aware <- TRUE` — §9.1), `R/utils.R` (validators,
    `ensure_attached` helper — §8), `R/rcicrely-package.R` (`"_PACKAGE"`
    sentinel), `R/rcicrely_result.R` + `R/report.R` (S3 classes for
    metric results).
3.  **I/O**: `R/ci_io.R` (PNG/JPEG directory reader),
    `R/noise_extraction.R` (base-image subtraction). Write tests against
    synthetic PNGs.
4.  **CI-from-responses**:
    - `R/ci_from_responses_2ifc.R` — thin wrapper around
      [`rcicr::generateCI2IFC()`](https://rdrr.io/pkg/rcicr/man/generateCI2IFC.html)
      / `batchGenerateCI2IFC()`. All §8 gotchas apply.
    - `R/ci_from_responses_briefrc.R` — native implementation per §5. No
      rcicr `_brief` calls. Test with synthetic noise matrix and
      known-signal responses.
5.  **Within-condition metrics**: `R/rel_split_half.R`, `R/rel_loo.R`,
    `R/rel_icc.R`. Write small-N cross-checks against
    [`psych::ICC()`](https://rdrr.io/pkg/psych/man/ICC.html) on tiny
    matrices to verify the mean-squares formula.
6.  **Between-condition primitives**: `R/pixel_t_test.R` (vectorised
    Welch’s t), `R/cluster_utils.R` (BFS 4-connectivity). Test the BFS
    on adversarial small cases (rings, corners, single pixel).
7.  **Between-condition metrics**: `R/rel_cluster_test.R`,
    `R/rel_dissimilarity.R`.
8.  **Orchestrators**: `R/run_within.R`, `R/run_between.R`.
9.  **S3 methods**: `R/plot_methods.R` + print/summary methods for each
    `rcicrely_*` class.
10. **Data-raw scripts**: §12. `01_generate_stimuli.R` first (slow,
    needs rcicr with `ncores = 1L`); then the two `02_*` scripts; then
    `03_generate_example_cis.R` which populates
    `inst/extdata/ci_examples/`.
11. **Tests**: one `test-*.R` per R file in `tests/testthat/`. Gate
    rcicr-dependent tests with `skip_if_not_installed("rcicr")`.
12. **Vignette**: `vignettes/tutorial.Rmd` — Rmd, one unified article,
    rcicr-dependent chunks `eval = FALSE`. Add `VignetteBuilder: knitr`
    to DESCRIPTION.
13. **README.md**: tagline, install command, minimal quick-start, link
    to pkgdown site, citation block, license, author, Claude-assistance
    acknowledgement if using.
14. **`_pkgdown.yml`**: via `usethis::use_pkgdown()`; accept default
    layout.
15. **CI workflows**: `usethis::use_github_action_check_standard()` and
    `usethis::use_pkgdown_github_pages()`.
16. Push; toggle repo settings (§14.3); watch first pkgdown build
    complete.

Target after step 15: `devtools::check(args = "--as-cran")` reports **0
/ 0 / 0** (or 0/0/1 with the benign “unable to verify current time” NOTE
on offline runs).

## 16. Lessons learned from the companion package

Distilled from the rcicrdiagnostics build — every one cost real time:

- **PNG-derived signals are scaled, not raw.**
  [`read_cis()`](https://olivethree.github.io/rcicrely/reference/read_cis.md)
  reads what was *rendered* to disk: `base + scaling(mask)`.
  [`extract_signal()`](https://olivethree.github.io/rcicrely/reference/extract_signal.md)
  therefore yields `scaling(mask)`, not the raw mask. Variance-based
  reliability metrics (ICC, Euclidean dissimilarity, t-tests, cluster
  mass) are wrong on this. Pearson-based metrics survive a single
  uniform scaling but break under per-CI “matched”.
  [`rcicr::computeInfoVal2IFC()`](https://rdrr.io/pkg/rcicr/man/computeInfoVal2IFC.html)
  is NOT affected — verified against rcicr v1.0.1 source: it does
  `norm(matrix(target_ci[["ci"]]), "f")`, always pulling the raw `$ci`
  from the CI-list element. The infoVal pitfall only applies to
  hand-rolled implementations (Brief-RC, custom code). Mode 2
  (`ci_from_responses_*`) is the safe path for everything. Lesson
  realised mid-v0.1.1; remedy in that release is a once-per-session
  warning at the input boundary plus a `looks_scaled()` heuristic
  backstop wired into every `rel_*` / `run_*` entry. Documentation
  chapter 3 of the user-guide vignette is the single source of truth
  users are pointed at.

- **Don’t invent data formats.** An earlier rcicrdiagnostics attempt
  used an “expanded” Brief-RC format (12 rows per trial, weight `1` for
  chosen / `-1/11` for unchosen). That’s not Schmitz. Users got confused
  by the `-0.0909…` values. Read the paper. Match the paper. The
  per-trial signed-selection format is what Schmitz et al.

  2024. §3.1.2 describes and what Schmitz’s own `genMask()` expects.

- **Don’t guess at statistics.** rcicrdiagnostics skips Brief-RC infoVal
  rather than ship a subtly wrong reference distribution. rcicrely
  implements statistics it can audit — anything beyond that gets
  implemented only once it can be validated against synthetic ground
  truth.

- **Don’t hallucinate authors in citations.** I once extended an author
  list for Brinkman et al. (2019) using names I “remembered”. Three of
  six were wrong. Use citations exactly as supplied by the author, or
  verify against DOI.

- **Don’t trust external-package READMEs for branch names.** rcicr’s
  `@development` branch doesn’t exist; an older README said it did.
  `git ls-remote` is authoritative.

- **Rename aggressively pre-1.0.** `rcicrdiagnostics` had
  `run_all_checks()` until it became obvious that only some checks run
  by default (the rest are gated on optional inputs). Renamed to
  `run_diagnostics()`. For rcicrely, consider upfront whether any name
  oversells its behaviour — fix before users depend on it.

- **Rcicr v1.0.1 is the authoritative rcicr.** Canonical. 17 functions.
  No Brief-RC. Anything else is someone’s fork. Verify identity of the
  installed rcicr with `packageVersion("rcicr")`,
  `packageDescription("rcicr")$Author`, `length(ls("package:rcicr"))`.

- **PSOCK clusters on macOS are racy.** `ncores = 1L` is the safe
  default.

- **Case matters on filesystems you thought were case-insensitive.**
  macOS `find -name "*.RData"` is case-sensitive. `ignore.case = TRUE`
  always.

- **Vignette-eval chunks load the *installed* package.** After any
  function rename, reinstall before re-rendering, or
  [`rmarkdown::render()`](https://pkgs.rstudio.com/rmarkdown/reference/render.html)
  executes stale code.

- **Temp folders pollute R CMD check.** If you drop a third-party source
  tree in `temp/` for reference, `.Rbuildignore` it (`^temp$`) and
  `.gitignore` it.

### 16.1. Additional lessons from the rcicrely build itself

These are the findings that showed up during the first pass of the
rcicrely implementation — distinct from the rcicrdiagnostics lessons
above. Every one cost real time.

- **rcicr’s `.Rdata` does not save the noise matrix.** Worth saying
  twice because the earlier draft of this doc implied otherwise. The
  rdata contains `stimuli_params` + `p` + `base_faces` + scalars — no
  `stimuli` object. The pixels × n_trials matrix is either the captured
  return value of `generateStimuli2IFC(return_as_dataframe = TRUE)`, or
  reconstructed via `rcicr::generateNoiseImage(stimuli_params[i, ], p)`
  per row. See §8.11 for the full story.

- **rcicr v1.0.1’s `batchGenerateCI2IFC` signature is
  `data, by, stimuli, responses, baseimage, rdata, ...`** where `by`,
  `stimuli`, `responses` are **column names** (strings) in `data`, not
  values. Easy trap: passing the columns themselves as numeric vectors.

- **The `$ci` element from `batchGenerateCI2IFC`’s return list is
  already the base-less noise component.** Don’t subtract the base again
  — it’s stacked straight into the signal matrix. `$combined` is
  base+scaled(ci); `$scaled` is the scaled ci only.

- **Don’t invent a `base_image_path` argument for the 2IFC wrapper.**
  The base is read from `base_faces[[baseimage]]` inside the rdata.
  Users pass a label (e.g. `baseimage = "base"`), not a file path. An
  early test Rmd that passed `base_image_path = base_jpg` cascaded three
  confusing errors before I spotted it. The Brief-RC wrapper *does* take
  a `base_image_path` because its noise-matrix source (the rdata or an
  rds) doesn’t necessarily carry a base — the base is used only to
  validate dimensions.

- **cli pluralization with multi-bullet abort needs an inline count.**
  This errors at runtime with “Cannot pluralize without a quantity”:

  ``` r
  cli::cli_abort(c(
    "Missing column{?s} in {.arg responses}:",   # {?s} here ...
    "*" = "{.val {missing_cols}}"                # ... can't see this
  ))
  ```

  Fix: put a count-providing variable on the same line as `{?s}`:

  ``` r
  n_missing <- length(missing_cols)
  cli::cli_abort(c(
    "Missing {n_missing} column{?s} in {.arg responses}:",
    "*" = "{.val {missing_cols}}"
  ))
  ```

  `cli`’s inflection only consults the line the `{?s}` lives on.

- **Stale output accumulation is a real trap.** If a data-raw script
  runs more than once and its output filenames are timestamped,
  [`list.files()`](https://rdrr.io/r/base/list.files.html) returns
  multiple paths. [`load()`](https://rdrr.io/r/base/load.html) /
  [`gzfile()`](https://rdrr.io/r/base/connections.html) error with
  `invalid 'description' argument` when handed a length-2 vector. Every
  script that writes timestamped files should either wipe previous-run
  outputs on entry (`file.remove(old)`) or always take the newest by
  `file.info(paths)$mtime`.

- **roxygen parses `[x, y]` as a link.** Writing “values in `[0, 1]`” in
  a `#'` comment emits a warning during `devtools::document()` about an
  unresolved link topic named `"0, 1"`. Rephrase to avoid square
  brackets, or escape them.

- **pkgdown auto-discovers top-level markdown.** Without an explicit
  ignore, `CLAUDE.md` is rendered as a page on the published site.
  `.Rbuildignore` it (`^CLAUDE\.md$`) to keep developer notes
  developer-only.

- **The `.datatable.aware <- TRUE` and `.Random.seed` names trip
  linters** that enforce snake_case. Both are mandated by their
  respective runtimes; ignore the lint.

- **Test-fixtures file needs a plain name.** `testthat` auto-sources any
  `tests/testthat/helper-*.R` at the start of the test run, so shared
  fixtures live there (for rcicrely: `helper-fixtures.R` with
  `make_sig()`, `make_sig_pair()`). Don’t name it `test-helpers.R` —
  that gets picked up as a test file and tries to run the helpers as
  assertions.

- **`ncores = 1L` only matters for stimulus generation.** The macOS
  PSOCK race is in
  [`rcicr::generateStimuli2IFC()`](https://rdrr.io/pkg/rcicr/man/generateStimuli2IFC.html).
  The CI-generation path (`batchGenerateCI2IFC`) does not spawn a
  cluster under the rcicrely call pattern, so no plumbing is needed
  there.

- **[`rmarkdown::render()`](https://pkgs.rstudio.com/rmarkdown/reference/render.html)
  uses the installed package, not the working directory.** Renaming a
  function in `R/*.R` and re-rendering a vignette without reinstalling
  produces stale output.
  `devtools::install(".", quick = TRUE, upgrade = "never")` before each
  render pass. On CI, R CMD check installs into a temp library
  automatically so this problem is local-only.

- **`vignettes/tutorial.html` and `vignettes/tutorial_files/` are build
  artefacts.** Either delete them before commit, or add them to
  `.gitignore`. The Rmd is the source of truth.

## 17. Relationship to rcicrdiagnostics

Independent packages, no mutual imports. Shared design:

- Same two-pipeline model (2IFC + Brief-RC).
- Same `data.table` + `cli` infrastructure.
- Same S3-result discipline: every exported function returns an object
  with a consistent shape, dedicated print/summary/plot methods.
- **Shared bogus data convention**: both packages’ data generators use
  the same column names (`participant_id`, `trial`, `stimulus`,
  `response`, `rt`), same archetype layout (20 signal / 6 biased / 2
  constant / 2 inverted for rcicrdiagnostics’ pattern; rcicrely’s is
  slightly different because it needs signal/noise differentiation
  across conditions). Same
  [`set.seed()`](https://rdrr.io/r/base/Random.html) values where
  reproducibility across packages matters for cross-referenced
  vignettes.

**Intended workflow**:

1.  User runs `rcicrdiagnostics::run_diagnostics()` on raw response data
    to catch coding errors, response bias, alignment bugs. Fix anything
    that fails.
2.  User generates CIs (via
    [`rcicr::generateCI2IFC()`](https://rdrr.io/pkg/rcicr/man/generateCI2IFC.html)
    for 2IFC or rcicrely’s native function for Brief-RC).
3.  User runs
    [`rcicrely::run_within()`](https://olivethree.github.io/rcicrely/reference/run_within.md)
    on each condition to check within-condition reliability.
4.  User runs
    [`rcicrely::run_between()`](https://olivethree.github.io/rcicrely/reference/run_between.md)
    on pairs of conditions to test discriminability.

Eventually, when rcicrely has a clean Brief-RC infoVal implementation,
rcicrdiagnostics’ `compute_infoval_summary(method = "briefrc")` can
delegate to it, lifting the current `"skip"` status.

## 18. Current status

**v0.1.1 in development.** Adds the raw-vs-rendered safety net, expanded
user guide, Brief-RC scaling exposure, and
[`rel_loo()`](https://olivethree.github.io/rcicrely/reference/rel_loo.md)
MAD outlier rule. State at time of writing:

- 20 R source files under `R/` (was 19) — adds `R/scaling_diagnostics.R`
  (`looks_scaled()` heuristic).
- **14 exported functions**: `read_cis`, `extract_signal`,
  `load_signal_matrix`, `read_noise_matrix`, `ci_from_responses_2ifc`,
  `ci_from_responses_briefrc`, `rel_split_half`, `rel_loo`, `rel_icc`,
  `pixel_t_test`, `rel_cluster_test`, `rel_dissimilarity`, `run_within`,
  `run_between`.
- **18 S3 methods** (print / summary / plot × 6 result classes —
  `rcicrely_split_half`, `rcicrely_loo`, `rcicrely_icc`,
  `rcicrely_cluster_test`, `rcicrely_dissim`, `rcicrely_report`).
- **Tests**: 15 `test-*.R` files under `tests/testthat/`. **116 passing,
  1 skipped** (psych cross-check, gated on
  `skip_if_not_installed("psych")`).
- **`R CMD check --as-cran`: 0 errors, 0 warnings, 0 notes** (or 0/0/1
  with the benign “unable to verify current time” note on offline runs).
- **Vignette**: `vignettes/tutorial.Rmd` renders to a ~220 KB HTML. 8
  chapters: what the package is / isn’t, install, data modes,
  within-condition, between-condition, Brief-RC end-to-end,
  interpretation, references.
- **pkgdown**: `_pkgdown.yml` organises reference pages by function
  area.
  [`pkgdown::build_site()`](https://pkgdown.r-lib.org/reference/build_site.html)
  produces a site with 14 reference pages + tutorial article. Not yet
  deployed — push and toggle repo settings (§14.3) to activate.
- **CI workflows**: `.github/workflows/R-CMD-check.yaml` (r-lib standard
  matrix, macOS / Windows / Ubuntu × release + devel + oldrel-1) and
  `.github/workflows/pkgdown.yaml` (build + deploy).
- **`temp/` test harness**: four bogus-data generators under
  `temp/bogus/` (stimuli, 2IFC, Brief-RC 12, Brief-RC 20) and five Rmd
  notebooks under `temp/tests/` — gitignored and `.Rbuildignore`-d. The
  companion to `tests/testthat/`: the testthat suite is the CI gate, the
  temp harness is the human sanity check on real data. See
  `temp/README.md` for the run order.

### Deferred to later versions

- **Brief-RC 20** — `ci_from_responses_briefrc(method = "briefrc20")`
  aborts with a clear message in v0.1. A Brief-RC 20 bogus-data fixture
  is generated under `temp/generated/responses/` so the shape is known
  when v0.2 adds support.
- **Brief-RC infoVal** — still deferred (§5.3). Needs a reference
  distribution matched to each participant’s trial count.
- **BCa bootstrap CIs** in
  [`rel_dissimilarity()`](https://olivethree.github.io/rcicrely/reference/rel_dissimilarity.md)
  — percentile only in v0.1. Base R
  [`quantile()`](https://rdrr.io/r/stats/quantile.html); no `boot` dep.

### Housekeeping

Update this section after every milestone. When this document falls out
of sync with the code, update it. It is the entry point for rebuilding
and the first place someone (human or AI) will look for “why did we do
it this way?” answers.
