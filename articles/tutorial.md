# rcicrely user guide

This guide answers two questions side by side: which function to call
for a given question, and why its output is trustworthy. Each function
gets a plain-language description, the signature, a worked example, how
to read the output, and the footguns to avoid.

## 1. What `rcicrely` is - and isn’t

In reverse correlation (RC), participants’ trial-level choices project
onto a pool of noise patterns to reconstruct the mental representation
they used. The output is a per-producer **classification image** (CI)
plus a group-level CI obtained by averaging.

`infoVal` (Brinkman et al., 2019) tells you whether a CI contains more
signal than chance. It does **not** tell you whether that signal is
*stable across producers*, *replicable on a different subset of
producers*, or *distinguishable from another condition’s CI*.

`rcicrely` fills that gap with a pixel-level reliability toolkit:

- **Within-condition** (how consistent is a condition’s CI across
  producers?): permuted split-half with Spearman-Brown, leave-one-out
  sensitivity, ICC(3,1) and ICC(3,k).
- **Between-condition** (is condition A’s CI distinguishable from
  condition B’s?): vectorised Welch t per pixel, cluster-based
  permutation with max-statistic FWER, bootstrap representational
  dissimilarity.

None of this requires a phase-2 trait-rating study. The metrics read
directly from the pixel-level signal produced by the original producers,
sidestepping the two-phase design Cone, Brown-Iannuzzi, Lei, and Dotsch
(2021) showed inflates Type I error.

Since v0.2 `rcicrely` also ships [`infoval()`](#section-infoval),
covering both 2IFC and Brief-RC with a single function.
[`infoval()`](https://olivethree.github.io/rcicrely/reference/infoval.md)
tests presence of signal; `rel_*()` tests reliability and
discriminability of signal; they are complementary.

## 2. Installation

``` r
# install.packages("remotes")
remotes::install_github("rdotsch/rcicr")      # 2IFC path only
remotes::install_github("olivethree/rcicrely")
```

If you only use Brief-RC you can skip `rcicr` - the Brief-RC CI
implementation is fully native (Schmitz, Rougier & Yzerbyt, 2024).

``` r
library(rcicrely)
```

## 3. Mental model: the signal matrix

Every reliability metric in this package takes a single object: a
**signal matrix**, pixels by participants. Get this right and the rest
of the package is a thin shell over it.

> **A note on the word *mask*.** The term has two distinct meanings in
> RC work and in this package:
>
> 1.  **Noise mask** (often shortened to *mask* or *raw mask*): the
>     participant’s *visual noise contribution* (the noise pattern that,
>     when superimposed on the base face, produces their classification
>     image). This is the core statistical object. Every `rel_*()`
>     function operates on a matrix of these.
> 2.  **Region mask** (or *pixel mask*): a logical vector marking which
>     pixels to include in an analysis (oval face region, eyes region,
>     etc.). This is what
>     [`face_mask()`](https://olivethree.github.io/rcicrely/reference/face_mask.md)
>     returns and what
>     [`infoval()`](https://olivethree.github.io/rcicrely/reference/infoval.md)’s
>     `mask` argument consumes.
>
> Throughout this guide, “mask” without qualification means the noise
> mask. Pixel selection is always called “region mask” or referenced
> through
> [`face_mask()`](https://olivethree.github.io/rcicrely/reference/face_mask.md).

### 3.1 Two paths to the signal matrix

There are two ways to obtain the signal matrix, and the choice matters
because they produce mathematically different objects.

    Mode 1: PNGs on disk           Mode 2: raw responses
              |                              |
       read_cis()                    ci_from_responses_2ifc() / _briefrc()
              |                              |
       extract_signal()                      |
              |                              |
       scaling(mask)                       mask                    <-- raw
              |                              |
              +---------> rel_*() <----------+

Mode 2 returns the **raw mask** (the participant’s noise contribution
before any display transformation). Mode 1 cannot, because PNGs on disk
have always been rendered with a display scaling step baked in;
subtracting the base recovers `scaling(mask)`, not `mask`. §3.2 explains
why this distinction is statistically load-bearing.

Prefer Mode 2 when raw responses are available. If you must use Mode 1,
the package emits a one-time warning per session. To silence:

``` r
# acknowledge per-call:
signal <- load_signal_matrix("data/cis/", "data/base.jpg",
                             acknowledge_scaling = TRUE)

# or silence for the whole session:
options(rcicrely.silence_scaling_warning = TRUE)
```

If you generated the PNGs yourself with
`rcicr::generateCI2IFC(scaling = "none")`, the rendered output is
effectively raw and the warning is a false positive.

### 3.2 Raw mask vs. rendered CI

A CI exists in two forms.

A **raw mask** is the participant’s noise contribution before any
display transformation. It is the object every reliability metric
operates on. For Brief-RC it is
`(noise_pool[, chosen_ids] %*% signed_responses) / n_unique` (Schmitz’s
`genMask()`); for 2IFC it is rcicr’s `$ci` element. Values are typically
in a narrow band around zero (roughly +/- 0.05).

The **rendered CI** is what gets saved to PNG. It is
`base + scaling(mask)`, where `scaling()` may stretch the mask to the
base image’s range (“matched”), apply a constant gain, or do nothing.
Scaling is a *display* transformation, not a statistical one.
Subtracting the base from a rendered CI does **not** recover the raw
mask:

    PNG_pixels - base = scaling(mask)   (NOT  mask)

The metrics introduced in §5 and §6 react differently to scaling:

|                                                                                            Metric | Behavior under scaling                                   |
|--------------------------------------------------------------------------------------------------:|:---------------------------------------------------------|
|                                                                       `rel_split_half`, `rel_loo` | Survives uniform scaling. Breaks under per-CI “matched”. |
|                                                              `rel_dissimilarity` (Pearson r half) | Same.                                                    |
|                                                              `rel_dissimilarity` (Euclidean half) | **Distorted** by any scaling.                            |
|                                                                                         `rel_icc` | **Distorted** by any scaling.                            |
|                                                                `pixel_t_test`, `rel_cluster_test` | **Distorted** by per-CI scaling.                         |
| [`rcicr::computeInfoVal2IFC()`](https://rdrr.io/pkg/rcicr/man/computeInfoVal2IFC.html) (standard) | Uses raw `$ci` internally, unaffected.                   |
|                                                          Hand-rolled `infoVal` (Brief-RC, custom) | **Distorted** by any scaling. Use the raw mask.          |

### 3.3 A concrete demonstration

Here is what scaling does to the same underlying mask.

``` r
set.seed(1)
n_pix <- 32L * 32L
n_p   <- 20L
mask  <- matrix(rnorm(n_pix * n_p, sd = 0.02), n_pix, n_p)

# Render: stretch each column to [-0.5, 0.5] (rcicr's "matched" /
# auto-scale option; the rcicr default for visualisation)
rendered <- apply(mask, 2L, function(m) {
  rng <- diff(range(m)); if (rng == 0) m else m * (1 / rng)
})

c(
  raw_max_range      = round(diff(range(mask)),     3L),
  rendered_max_range = round(diff(range(rendered)), 3L)
)
#>      raw_max_range rendered_max_range 
#>              0.162              1.159

# Pearson r between same-column raw and rendered: ~ 1
# (a perfect linear stretch preserves Pearson)
mean(vapply(seq_len(n_p),
            function(j) cor(mask[, j], rendered[, j]),
            numeric(1L)))
#> [1] 1

# But the variance ratios that drive ICC are very different
var(rowMeans(mask))      # variance of the group mean over pixels
#> [1] 1.923774e-05
var(rowMeans(rendered))  # same metric, on the rendered version
#> [1] 0.001134698
```

Pearson r per column is exactly 1 (linear scaling is correlation-
preserving), but the row-mean variance the ICC formula depends on is
~25x larger on the rendered version. ICC(3,\*) computed on rendered data
is therefore not the same ICC computed on the raw mask. The same logic
applies to Euclidean distance and to `infoVal`.

### 3.4 Base image preparation

Both 2IFC and Brief-RC pipelines depend on a base face image. Both the
visible quality of generated CIs and the reliability statistics that
follow are sensitive to how this image is prepared. Requirements:

- **Square** (e.g. 256x256, 512x512).
- **Grayscale** (single channel). RGB images error out of rcicr with
  non-conformable-array messages.
- **Pixel range `[0, 1]`** (the convention
  [`png::readPNG`](https://rdrr.io/pkg/png/man/readPNG.html) and
  [`jpeg::readJPEG`](https://rdrr.io/pkg/jpeg/man/readJPEG.html)
  produce).
- **Centred face** with eye / nose / mouth roughly at the geometry the
  Schmitz oval mask assumes (eyes near the upper third, mouth near the
  lower third). The built-in `face_mask(region = ...)` geometries assume
  this layout; custom or off-centre faces need their own mask geometry.

The conventional base face is a *morph* of several individual faces,
which removes idiosyncratic features and centres the geometry.
[webmorphR](https://github.com/debruine/webmorphR) (DeBruine, 2022) is
the current best-in-class tool for producing one. Typical pipeline:

``` r
# install.packages("webmorphR")
library(webmorphR)

stim <- read_stim("path/to/raw_face_images/") |>
  auto_delin() |>                       # auto landmark delineation
  align(procrustes = TRUE) |>           # Procrustes alignment
  crop(width = 0.85, height = 0.85) |>  # tight crop around the face
  to_size(c(256, 256)) |>               # resize to rcicr-friendly size
  greyscale() |>                        # single channel, in [0, 1]
  avg()                                 # morph into one average face

write_stim(stim, dir = "stimuli/", names = "base", format = "png")
```

`write_stim()` produces `stimuli/base.png`, ready for
`rcicr::generateStimuli2IFC(base_face_files = list(base = "stimuli/base.png"))`.
For non-square crops, oval-region operations, and feature-aligned masks
(which then plug into [`load_face_mask()`](#section-load-face-mask)),
see `?webmorphR::mask_oval` and `?webmorphR::crop_tem`. Cite DeBruine
(2022).

If you do not have a stack to morph (e.g. you are working from a single
base photograph), prepare it by hand: any image editor will do (GIMP /
Krita / Photoshop / PowerPoint export). Convert to grayscale, crop to a
square that fits the face within the central ~70%, resize to the target
dimension, save as PNG.

The function-reference sections (§4–§6) work on small synthetic data so
the vignette knits in seconds. For an end-to-end demonstration on a
published RC dataset, see §7.

## 4. Function reference: I/O

The six functions in this section get your data into the package’s
standard signal-matrix shape (`pixels x participants`). Two of them
([`ci_from_responses_2ifc()`](https://olivethree.github.io/rcicrely/reference/ci_from_responses_2ifc.md)
and
[`ci_from_responses_briefrc()`](https://olivethree.github.io/rcicrely/reference/ci_from_responses_briefrc.md))
take **trial-level response data**; before diving into them the shape
and coding conventions for that data deserve a quick callout.

### Input response-data shape

One row per trial. The format can be a `data.frame`, `data.table`,
`tibble`, or anything that behaves like a data frame; the source (CSV,
RData, Parquet, programmatic construction) does not matter, only that
the columns below are present.

| Column           | Purpose                                                                                   |
|------------------|-------------------------------------------------------------------------------------------|
| `participant_id` | Identifier per producer. Character or integer.                                            |
| `stimulus`       | Stimulus / pool id. Numbering must match what was used at stimulus generation.            |
| `response`       | `+1` or `-1` (see below).                                                                 |
| `rt` (optional)  | Response time in milliseconds. Not consumed by `rcicrely`; useful for `rcicrdiagnostics`. |

**Response coding.**

- **2IFC**: `response = +1` if the participant chose the *oriented*
  variant on that trial, `-1` if the *inverted* variant. Each trial has
  its own unique stimulus pair, so `stimulus` indexes the *pair* and
  ranges `1:n_trials`.
- **Brief-RC 12** (Schmitz, Rougier & Yzerbyt, 2024): on each trial the
  participant sees 12 noisy faces (6 oriented + 6 inverted from 6
  distinct underlying noise patterns) and picks one. Recorded as one row
  per trial, with `stimulus` = pool id of the chosen pattern and
  `response = +1` if the oriented variant was chosen, `-1` if the
  inverted. The same `stimulus` id can repeat across trials.

> **Watch out for the `{0, 1}` miscoding.** This is the single most
> common silent failure in RC pipelines: experiment software sometimes
> records “left” / “right” or “first” / “second” as `0` / `1`, and the
> analyst forgets to recode. Pass `{0, 1}` responses to
> `ci_from_responses_*()` and the resulting CI is approximately blank
> (the math leans heavily on the sign). Both functions error out at the
> input boundary if responses are not in `{-1, +1}`. If you see that
> error, the fix is one line:
> `responses$response <- 2 * responses$response - 1`.

> **Brief-RC: do not record the 11 unselected faces.** Each Brief-RC
> trial shows 12 faces but the participant chooses one. Record only the
> chosen face (one row per trial, as above). Do *not* record the 11
> unselected faces as additional rows with `response = 0`, and do *not*
> use the “expanded” 12-rows-per-trial format that weights chosen as
> `+1` and unchosen as `-1/11`. Neither convention matches the Schmitz
> `genMask()` formula
> (`mask = noise[, chosen] %*% chosen_responses / n_unique_chosen`) that
> [`ci_from_responses_briefrc()`](https://olivethree.github.io/rcicrely/reference/ci_from_responses_briefrc.md)
> implements. Unselected faces’ noise patterns simply do not enter the
> multiplication; their absence *is* the zero.

### 4.1 `load_signal_matrix()`

**What it does.** Reads every PNG/JPEG in a directory, subtracts the
base image, and returns a pixels-by-participants signal matrix in one
call. The end product is `scaling(mask)` (Mode 1; see §3.2).

**When to use it.** You have one CI image per producer on disk.

``` r
# Tiny synthetic dataset on disk to demonstrate against:
demo_dir  <- tempfile("ci_demo_"); dir.create(demo_dir)
base_img  <- matrix(0.5, 32L, 32L)
base_path <- file.path(demo_dir, "_base.png")
png::writePNG(base_img, base_path)
for (i in 1:6) {
  img <- pmin(pmax(base_img + matrix(rnorm(32 * 32, sd = 0.05), 32, 32),
                   0), 1)
  png::writePNG(img, file.path(demo_dir, sprintf("p%02d.png", i)))
}

signal <- load_signal_matrix(
  dir                 = demo_dir,
  base_image_path     = base_path,
  pattern             = "^p[0-9]+\\.png$",  # exclude _base.png
  acknowledge_scaling = TRUE
)
dim(signal); attr(signal, "img_dims")
#> [1] 1024    6
#> [1] 32 32
```

**Reading the result.** `nrow = n_pixels`, `ncol = n_files`; the
`img_dims` attribute carries `c(nrow, ncol)` of the original image so
plot helpers know how to reshape.

**Common mistakes.** Pointing at the wrong base (a generated CI PNG of
the right size will pass the dimension check). Mixing CIs of different
sizes (aborts loudly). Treating this as the raw mask (it isn’t; see
§3.2).

**Granular API.** If you need to intervene between reading and
subtracting (custom masking, swap the base), the two underlying
functions are exported: `read_cis(dir)` returns the raw pixel matrix,
then `extract_signal(cis, base_image_path)` does the subtraction.
[`load_signal_matrix()`](https://olivethree.github.io/rcicrely/reference/load_signal_matrix.md)
is exactly the composition of those two.

### 4.2 `read_noise_matrix()`

**What it does.** Loads the pixels by trials noise pool used at stimulus
generation. Autodetects three formats: whitespace-delimited text
(`.txt`, the format Schmitz et al. distribute on OSF), saved R `.rds`,
or an rcicr `.Rdata` (reconstructs from `stimuli_params` + `p`).

**When to use it.** You are computing Brief-RC CIs and want to load the
noise pool once, then reuse it across multiple
[`ci_from_responses_briefrc()`](#section-cifr-briefrc) calls.

``` r
# Plain text (the Schmitz et al. 2024 OSF format):
noise <- read_noise_matrix("data/noise_matrix.txt")

# Cached .rds saved from rcicr's return-as-dataframe output (fast):
noise <- read_noise_matrix("data/noise_matrix.rds")

# rcicr .Rdata (slow; reconstructs each trial via
# rcicr::generateNoiseImage):
noise <- read_noise_matrix("data/rcicr_stimuli.Rdata")
```

**Common mistakes.** Assuming the rcicr `.Rdata` contains a `stimuli`
object (it does not; the pixels-by-trials matrix is the *return value*
of `generateStimuli2IFC(return_as_dataframe = TRUE)`, never saved to the
file). Cache to `.rds` once and load from there. Multiple stale `.rdata`
files in one directory (rcicr timestamps filenames; pick the most recent
by mtime if you have several).

For a field guide to which objects live in the rcicr `.Rdata` and what
each one means, see
[`?read_noise_matrix`](https://olivethree.github.io/rcicrely/reference/read_noise_matrix.md).

### 4.3 `ci_from_responses_2ifc()`

**What it does.** Wraps
[`rcicr::batchGenerateCI2IFC()`](https://rdrr.io/pkg/rcicr/man/batchGenerateCI2IFC.html)
with all known gotchas pre-handled, returns a uniform result list.

**When to use it.** You have 2IFC trial-level responses and want
per-producer CIs in one call. Requires `rcicr` installed.

``` r
res <- ci_from_responses_2ifc(
  responses       = my_responses,
  rdata_path      = "data/rcicr_stimuli.Rdata",
  baseimage       = "base",
  participant_col = "participant_id",
  keep_rendered   = FALSE   # set TRUE if you want $rendered_ci for plotting
)
signal <- res$signal_matrix
```

**Reading the result.** `$signal_matrix` is the **raw mask** (rcicr’s
`$ci` per producer) - feed this to `rel_*`. `$participants` and
`$img_dims` are convenience metadata. `$rcicr_result` is the raw rcicr
return value if you want anything else. `$rendered_ci` (present only if
`keep_rendered = TRUE`) is `base + scaling(mask)` per producer -
**visualization only**, do not feed to `rel_*`.

**Common mistakes.** Inventing a `base_image_path` argument - the base
is read from the rdata via `baseimage` *label*, not a path. Passing
`$rendered_ci` to `rel_*`. Trusting `scaling = "autoscale"` to be
“raw” - the rendered output is scaled even when the signal_matrix isn’t.

### 4.4 `ci_from_responses_briefrc()`

**What it does.** Native Brief-RC 12 implementation of Schmitz’s
`genMask()` formula. No `rcicr::*_brief` calls (those don’t exist in the
upstream rcicr v1.0.1).

**When to use it.** You have Brief-RC 12 trial-level responses and the
noise pool used at stimulus generation.

``` r
res <- ci_from_responses_briefrc(
  responses       = my_responses,    # one row per trial
  rdata_path      = "data/rcicr_stimuli.Rdata",
  base_image_path = "data/base.jpg",
  method          = "briefrc12",
  scaling         = "none"           # raw mask only - no $rendered_ci
)
signal <- res$signal_matrix
```

To also get the rendered CI for plotting:

``` r
res <- ci_from_responses_briefrc(
  responses       = my_responses,
  noise_matrix    = my_noise,
  base_image_path = "data/base.jpg",
  scaling         = "matched"        # rcicr auto-scale (used for
                                     # visualisation by Schmitz et al.
                                     # 2024 in Experiment 2)
)
# res$signal_matrix is the raw mask, feed to rel_*
# res$rendered_ci is base + matched(mask), for PNG / plotting only
```

**Reading the result.** `$signal_matrix` is **always** the raw mask,
regardless of `scaling`. `$rendered_ci` exists only when
`scaling != "none"`. `$participants`, `$img_dims`, `$scaling` are
metadata.

**Common mistakes.** Passing the “expanded” 12-rows-per-trial format -
Brief-RC 12 is one row per trial (the chosen face), not 12. Asking for
`method = "briefrc20"` (deferred). Using `$rendered_ci` for downstream
stats.

## 5. Function reference: within-condition reliability

For the rest of this guide we’ll work with a small synthetic dataset
where two conditions have different latent signals.

``` r
set.seed(1)
make_condition <- function(region_top, n_pix = 32L * 32L,
                           n_p = 30L, snr = 0.3) {
  noise <- matrix(rnorm(n_pix * n_p, sd = 0.05), n_pix, n_p)
  mask <- if (region_top) {
    c(rep(1, n_pix / 2L), rep(0, n_pix / 2L))
  } else {
    c(rep(0, n_pix / 2L), rep(1, n_pix / 2L))
  }
  out <- noise + snr * outer(mask, runif(n_p, 0.7, 1.3))
  attr(out, "img_dims") <- c(32L, 32L)
  out
}
sig_A <- make_condition(region_top = TRUE,  n_p = 30L)
sig_B <- make_condition(region_top = FALSE, n_p = 30L)
```

### 5.1 `rel_split_half()`

**What it does.** Estimates how stable a condition’s group CI is across
random halves of producers. Reports the per-permutation correlation
between halves and the Spearman-Brown projection to the full sample.

**When to use it.** You want to answer “would I get the same group CI if
I’d run a different half of my participants?”. `r_SB` is usually the
headline number.

``` r
r_sh <- rel_split_half(sig_A, n_permutations = 500L, seed = 1L,
                       progress = FALSE)
print(r_sh)
#> <rcicrely split-half reliability>
#>   N producers:          30
#>   n_permutations:       500
#>   mean per-split r:     0.993  [0.992, 0.993]
#>   Spearman-Brown r_SB:  0.996  [0.996, 0.997]
```

``` r
plot(r_sh)
```

![Permutation distribution of split-half
r.](tutorial_files/figure-html/plot-split-half-1.png)

Permutation distribution of split-half r.

**Reading the result.** `$r_hh` = mean per-split r. `$r_sb` =
Spearman-Brown projected reliability of the full sample (report this
unless asked otherwise). `$ci_95`, `$ci_95_sb` are percentile CIs.
`$distribution` is the full per-permutation r vector.

**Common mistakes.** Reporting `$r_hh` instead of `$r_sb` when the
audience cares about the full-sample CI. Running with
`n_permutations < 200` and trusting the CIs.

### 5.2 `rel_loo()` - influence screening (not a reliability metric)

**What it does.** For each producer, correlates the full-sample group CI
with the group CI computed without that producer, and returns a
**z-scored** version of those correlations. Producers whose z-score sits
clearly below zero have a disproportionate influence on the group
pattern and are flagged for inspection.

**When to use it.** To screen for producers whose individual CI is far
from the group pattern - often a sign of response miscoding, fatigue, or
a genuinely atypical mental representation. This is a diagnostic; it is
*not* a reliability statistic (see the note at the end of this section).

``` r
r_loo <- rel_loo(sig_A, flag_threshold = 2.5, flag_method = "sd")
#> Warning: `rel_loo(flag_method = "sd")` is deprecated since v0.3 and will be removed in
#> v0.4.
#> ℹ The MAD/median rule (the new default) is robust to the influential producers
#>   LOO is designed to flag.
#> • Drop the `flag_method` argument or pass `flag_method = "mad"` explicitly.
#> ℹ Silence: `options(rcicrely.silence_loo_deprecation = TRUE)`.
print(r_loo)
#> <rcicrely leave-one-out influence screening>
#>   N producers:        30
#>   flag rule:          sd (threshold = 2.50)
#> 
#>   z-scored influence (most influential first):
#>     p002                  z =  -2.19  (r_loo = 0.9999)
#>     p029                  z =  -1.65  (r_loo = 0.9999)
#>     p003                  z =  -1.58  (r_loo = 0.9999)
#>     p001                  z =  -1.43  (r_loo = 0.9999)
#>     p026                  z =  -0.88  (r_loo = 0.9999)
#>     ... (25 more)
#> 
#>   flagged producers:  none
#> 
#>   Note: r_loo values are near 1 by construction because the
#>   full-sample and leave-one-out means share (N-1)/N of their
#>   data. Use z_scores (relative ordering), not r_loo levels,
#>   to interpret influence. This is a diagnostic, not a
#>   reliability statistic.
```

The default is mean +/- SD; use `flag_method = "mad"` for a median +/-
MAD rule that doesn’t get pulled around by the very outliers it tries to
detect:

``` r
r_loo_mad <- rel_loo(sig_A, flag_threshold = 2.5, flag_method = "mad")
c(sd_threshold = r_loo$threshold, mad_threshold = r_loo_mad$threshold)
#>  sd_threshold mad_threshold 
#>     0.9999303     0.9999308
```

``` r
plot(r_loo)
```

![z-scored LOO influence per
producer.](tutorial_files/figure-html/plot-loo-1.png)

z-scored LOO influence per producer.

For a tidy ordered table of producers by influence, use
[`rel_loo_z()`](https://olivethree.github.io/rcicrely/reference/rel_loo_z.md):

``` r
head(rel_loo_z(r_loo))
#   participant_id  correlation  z_score  flag
# 1     outlier01        0.9714    -3.12  TRUE
# 2           p07        0.9942    -0.84  FALSE
# ...
```

**Reading the result.** `$z_scores` is the informative quantity; plot or
report it directly. `$correlations` is retained for transparency, but
**do not read it as reliability**: `r_loo` is almost always in
`[0.95, 0.999]` at N = 30 even on very noisy CIs because the full-sample
mean and the leave-one-out mean share `(N-1)/N` of their data. A naive
reader can mistake “LOO r = 0.98” for “the CIs are 98% reliable”; that
is wrong. For reliability use
[`rel_split_half()`](https://olivethree.github.io/rcicrely/reference/rel_split_half.md)
or
[`rel_icc()`](https://olivethree.github.io/rcicrely/reference/rel_icc.md)
(§5.1, §5.3). `$flagged` lists ids below threshold; `$summary_df` is a
tidy table sorted by `z_score`.

**Common mistakes.** Treating `$flagged` as “drop these producers”;
investigate first. A flagged producer may have a genuinely atypical
mental representation rather than bad data. Cross-check with
`rcicrdiagnostics` to rule out response-coding errors.

### 5.3 `rel_icc()`

**What it does.** Two-way mixed ICC via direct mean squares. Default
returns ICC(3,1) and ICC(3,k). Pixels are fixed (the image grid),
participants are random (drawn from a producer population).

**When to use it.** Single-scalar reliability summary. ICC(3,k) is “how
stable is the group-mean CI?”, ICC(3,1) is “how informative is one
producer’s CI?”.

``` r
r_icc <- rel_icc(sig_A, variants = c("3_1", "3_k"))
print(r_icc)
#> <rcicrely ICC>
#>   model:        ICC(3,*) / two-way mixed; pixels fixed, participants random
#>   N targets:    1024 pixels
#>   N raters:     30 participants
#> 
#>   Primary: ICC(3,1) (single-producer reliability)
#>   ICC(3,1):      0.8788
#> 
#>   Secondary: ICC(3,k) (group-mean reliability)
#>   ICC(3,k):      0.9954
```

To also report ICC(2,\*):

``` r
r_icc_full <- rel_icc(sig_A, variants = c("3_1", "3_k", "2_1", "2_k"))
print(r_icc_full)
#> <rcicrely ICC>
#>   model:        ICC(3,*) / two-way mixed; pixels fixed, participants random
#>   N targets:    1024 pixels
#>   N raters:     30 participants
#> 
#>   Primary: ICC(3,1) (single-producer reliability)
#>   ICC(3,1):      0.8788
#> 
#>   Secondary: ICC(3,k) (group-mean reliability)
#>   ICC(3,k):      0.9954
#> 
#>   Two-way-random variants (for reviewer comparability only):
#>   ICC(2,1):      0.8565
#>   ICC(2,k):      0.9944
```

**Reading the result.** `$icc_3_1`, `$icc_3_k`, `$icc_2_1`, `$icc_2_k`
(the latter two NA unless requested). `$ms_rows`, `$ms_cols`,
`$ms_error` are the underlying ANOVA mean squares for transparency.
`$model` is a description string the print method shows.

**Common mistakes.** Reporting ICC(2,*) without need. ICC(2,*) treats
pixels as a random sample from a pixel population; in RC studies pixels
are a fixed image grid, so ICC(3,*) is the correctly-specified model.
The numbers are usually close at high pixel counts, but ICC(3,*) is the
defensible choice. Use `variants = c("3_1", "3_k", "2_1", "2_k")` when
both are needed.

ICC is **strongly** scale-sensitive. Use the raw mask only.

### 5.4 `run_within()`

**What it does.** Calls all three within-condition metrics and bundles
the result for joint printing/plotting.

``` r
within_A <- run_within(
  sig_A,
  n_permutations = 500L,
  seed           = 1L,
  progress       = FALSE
)
print(within_A)
#> <rcicrely report: within>
#> ---
#> * $results$split_half
#> <rcicrely split-half reliability>
#>   N producers:          30
#>   n_permutations:       500
#>   mean per-split r:     0.993  [0.992, 0.993]
#>   Spearman-Brown r_SB:  0.996  [0.996, 0.997]
#> ---
#> * $results$icc
#> <rcicrely ICC>
#>   model:        ICC(3,*) / two-way mixed; pixels fixed, participants random
#>   N targets:    1024 pixels
#>   N raters:     30 participants
#> 
#>   Primary: ICC(3,1) (single-producer reliability)
#>   ICC(3,1):      0.8788
#> 
#>   Secondary: ICC(3,k) (group-mean reliability)
#>   ICC(3,k):      0.9954
```

`$results$split_half` and `$results$icc` are the two reliability result
objects. `plot(within_A)` lays them out side by side. Leave-one-out
influence screening is run separately with
[`rel_loo()`](https://olivethree.github.io/rcicrely/reference/rel_loo.md)
(it is an influence diagnostic, not a reliability statistic; see §5.2).

## 6. Function reference: between-condition inference

The two examples in §5 (`sig_A`, `sig_B`) use a half-image signal, which
produces *very* strong contrasts (useful for demonstrating that
within-condition reliability metrics work, but the between-condition
cluster figures end up looking like two solid blobs rather than
realistic localized clusters). For §6.1–§6.2 we therefore use a small
additional fixture with **localized circular signal**: two small
“patches” at different positions on a 32 x 32 grid, mimicking how real
RC signal tends to cluster around face features rather than fill the
whole image.

``` r
set.seed(7)
make_localized <- function(centre, n_p = 30L, n_side = 32L,
                           snr = 0.20, radius = 0.18, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  rr <- row(matrix(0, n_side, n_side)) / n_side - centre[1]
  cc <- col(matrix(0, n_side, n_side)) / n_side - centre[2]
  d  <- sqrt(rr^2 + cc^2)
  m  <- as.vector(pmax(0, 1 - d / radius))
  out <- snr * outer(m, runif(n_p, 0.7, 1.3)) +
         matrix(rnorm(n_side^2 * n_p, sd = 0.05), n_side^2, n_p)
  attr(out, "img_dims") <- c(n_side, n_side)
  out
}
# Condition A: signal in upper-left; Condition B: signal in lower-right
sig_A_loc <- make_localized(c(0.30, 0.30), n_p = 30L, seed = 11L)
sig_B_loc <- make_localized(c(0.70, 0.70), n_p = 30L, seed = 12L)
```

### 6.1 `pixel_t_test()`: building block, not standalone test

**What it does.** Computes a Welch’s t-statistic per pixel between two
condition signal matrices. Welch (unequal variances) is the correct
choice because conditions can differ in both N and per-pixel variance.
With `paired = TRUE`, computes a paired t on the per-producer difference
instead.

**Why this function exists.** 65,000 simultaneous t-tests are
statistically meaningless on their own; this is a *building block*, not
an inferential test. Three legitimate uses:

1.  **Internal use by
    [`rel_cluster_test()`](https://olivethree.github.io/rcicrely/reference/rel_cluster_test.md)**:
    every permutation of the cluster test recomputes a t-map;
    [`pixel_t_test()`](https://olivethree.github.io/rcicrely/reference/pixel_t_test.md)
    exposes that computation.
2.  **Diagnostic visualisation**: inspecting the raw t-map *before*
    running the cluster test helps you pick a sensible cluster- forming
    threshold and check whether any spatially-coherent structure is even
    there to test.
3.  **Custom inference pipelines**: implementing TFCE by hand, defining
    non-standard cluster geometries, or applying a different
    multiple-comparison correction (FDR, max-T, etc.).

**Do not** threshold the returned vector at `|t| > 1.96` and report
“significant pixels”. At 65,536 pixels you would get ~3,300 false
positives by chance under H0 (any α), and spatial clustering in real CI
data makes that worse. Use
[`rel_cluster_test()`](https://olivethree.github.io/rcicrely/reference/rel_cluster_test.md)
for inference.

``` r
t_vec <- pixel_t_test(sig_A_loc, sig_B_loc)
range(t_vec)
#> [1] -13.16874  12.78800
```

A quick visualization makes the output structure obvious. Positive (red)
regions indicate where condition A’s signal is stronger than B’s;
negative (blue) regions where B’s is stronger than A’s; near- white
pixels carry no consistent difference. **This is descriptive only**
(there is no inference yet):

``` r
local({
  op <- par(mar = c(1, 1, 3, 6) + 0.1)
  on.exit(par(op))
  n_side <- 32L
  tmap   <- matrix(t_vec, n_side, n_side)
  rng    <- max(abs(t_vec))
  graphics::image(
    seq_len(n_side), seq_len(n_side),
    t(tmap[n_side:1L, ]),
    col  = grDevices::hcl.colors(256L, "RdBu", rev = TRUE),
    zlim = c(-rng, rng),
    main = "Raw pixel-wise t-map (descriptive only)",
    axes = FALSE, xlab = "", ylab = "",
    asp  = 1, useRaster = TRUE
  )
  graphics::box(col = "grey80", lwd = 0.5)
})
```

![Raw pixel-wise Welch t-map between two conditions with localized
signal patches. Positive (red) = condition A larger; negative (blue) =
condition B larger. Note the two distinct hot regions corresponding to
the injected signal locations; the rest of the image is dominated by
noise. This map is descriptive only (significance comes from
rel_cluster_test()).](tutorial_files/figure-html/pixel-t-viz-1.png)

Raw pixel-wise Welch t-map between two conditions with localized signal
patches. Positive (red) = condition A larger; negative (blue) =
condition B larger. Note the two distinct hot regions corresponding to
the injected signal locations; the rest of the image is dominated by
noise. This map is descriptive only (significance comes from
rel_cluster_test()).

The hot blue/red regions are exactly where the two conditions differ;
everything else is noise. To turn this descriptive map into a statement
about *which* differences are reliable, we move on to the cluster test.

### 6.2 `rel_cluster_test()`: spatially-localized inference with FWER control

**What it does.** Identifies *spatially contiguous regions* of the image
where two conditions differ, with family-wise error control across the
whole image (so you do not over-claim significance from the kind of
mass-testing problem §6.1 warned about).

**The vocabulary.** Three terms are unavoidable; here is what each one
means in plain language.

- **Cluster-forming threshold** (`cluster_threshold`, default `2.0`). A
  t-value cutoff used at the *first* step: any pixel with `|t| > 2.0` is
  “in”, everything else is “out”. This is *not* a significance
  threshold; it is the cutoff for deciding which pixels are even
  candidates for being part of a cluster. Significance is decided
  afterwards. Lower thresholds catch more diffuse signal; higher
  thresholds isolate sharper peaks.

- **4-connectivity**. Once we have a binary “in / out” image, we group
  adjacent “in” pixels into clusters. *4-connectivity* treats two pixels
  as adjacent only if they share an edge (up / down / left / right),
  i.e. the four orthogonal neighbours. *8-connectivity*, the
  alternative, also counts diagonals. 4-connectivity is the more
  conservative choice: a one-pixel-wide diagonal stripe is split into
  isolated pixels and forms no cluster, whereas 8-connectivity would
  chain them. Suits face-CI signal, which tends to be contiguous in
  orthogonal directions.

      4-connectivity:    8-connectivity:
           . X .              X X X
           X o X              X o X
           . X .              X X X
       (o is the centre, X are neighbours, . are non-neighbours)

- **Cluster mass** (`mass`). For each candidate cluster, mass is the
  *sum of t-values within it*. Not pixel count. This means a small
  cluster with very strong t-values can have the same mass as a large
  cluster with weak t-values; the metric mixes *extent* (how big) and
  *magnitude* (how strong). Cluster mass was introduced by Maris &
  Oostenveld (2007) and is the dominant choice in EEG / MEG / fMRI
  cluster-permutation work.

- **Stratified label permutation**. To build the null distribution, we
  randomly relabel which producers belong to A vs. B many times,
  recompute the t-map, find clusters, and record the *maximum* cluster
  mass per permutation. “Stratified” means we preserve the original
  group sizes `(N_A, N_B)` exactly each time, so the Welch denominator
  is not perturbed under H0. With `paired = TRUE`, the permutation is
  sign-flip on matched pairs instead.

- **Max-statistic FWER**. Family-wise error rate = the probability of
  any false positive across the *entire image*. By comparing observed
  cluster mass to the *maximum* mass under each permutation, we directly
  bound that probability at `alpha`. No per-pixel correction is needed;
  the cluster-level p-value is already FWE-corrected.

**When to use it.** “Where in the image do these two conditions differ,
and how confident am I that the differences are not noise?”

``` r
clust <- rel_cluster_test(
  sig_A_loc, sig_B_loc,
  img_dims          = c(32L, 32L),
  n_permutations    = 500L,
  cluster_threshold = 2.0,
  alpha             = 0.05,
  seed              = 1L,
  progress          = FALSE
)
print(clust)
#> <rcicrely cluster-based permutation test>
#>   N_A = 30, N_B = 30
#>   image dims:        32 x 32
#>   cluster threshold: |t| > 2.00
#>   n_permutations:    500
#>   alpha:             0.05
#>   clusters found:    46 (2 significant)
#> 
#>  cluster_id direction    mass size p_value significant
#>          11       neg -466.77   79   0.000        TRUE
#>          14       neg   -4.95    2   0.452       FALSE
#>          20       neg   -3.55    1   0.828       FALSE
#>           1       neg   -3.11    1   0.930       FALSE
#>          10       neg   -2.95    1   0.966       FALSE
#>           4       neg   -2.87    1   0.984       FALSE
#>           8       neg   -2.82    1   0.986       FALSE
#>          17       neg   -2.59    1   1.000       FALSE
#>          15       neg   -2.59    1   1.000       FALSE
#>          13       neg   -2.51    1   1.000       FALSE
#>          16       neg   -2.42    1   1.000       FALSE
#>           7       neg   -2.42    1   1.000       FALSE
#>           5       neg   -2.31    1   1.000       FALSE
#>          23       neg   -2.28    1   1.000       FALSE
#>           2       neg   -2.26    1   1.000       FALSE
#>           6       neg   -2.21    1   1.000       FALSE
#>          24       neg   -2.20    1   1.000       FALSE
#>          22       neg   -2.17    1   1.000       FALSE
#>          21       neg   -2.13    1   1.000       FALSE
#>          18       neg   -2.12    1   1.000       FALSE
#>          19       neg   -2.08    1   1.000       FALSE
#>           9       neg   -2.07    1   1.000       FALSE
#>          12       neg   -2.05    1   1.000       FALSE
#>           3       neg   -2.03    1   1.000       FALSE
#>           2       pos  481.53   76   0.000        TRUE
#>           3       pos    4.78    2   0.466       FALSE
#>          11       pos    4.67    2   0.534       FALSE
#>          21       pos    3.03    1   0.954       FALSE
#>          15       pos    2.94    1   0.970       FALSE
#>          19       pos    2.65    1   0.998       FALSE
#>          14       pos    2.64    1   0.998       FALSE
#>           4       pos    2.54    1   0.998       FALSE
#>           1       pos    2.51    1   0.998       FALSE
#>          22       pos    2.48    1   1.000       FALSE
#>          10       pos    2.44    1   1.000       FALSE
#>           7       pos    2.38    1   1.000       FALSE
#>           6       pos    2.37    1   1.000       FALSE
#>           5       pos    2.31    1   1.000       FALSE
#>          17       pos    2.27    1   1.000       FALSE
#>           8       pos    2.15    1   1.000       FALSE
#>          20       pos    2.15    1   1.000       FALSE
#>          18       pos    2.08    1   1.000       FALSE
#>          16       pos    2.06    1   1.000       FALSE
#>           9       pos    2.02    1   1.000       FALSE
#>          13       pos    2.00    1   1.000       FALSE
#>          12       pos    2.00    1   1.000       FALSE
```

``` r
plot(clust)
```

![Pixel-wise Welch t-map (red = condition A larger, blue = condition B
larger) with significant cluster contours overlaid in black. The two
clusters correspond to the two injected signal patches, and only those
clusters are FWE-corrected significant at alpha = 0.05; the surrounding
noise does not produce false-positive clusters thanks to the
max-statistic
correction.](tutorial_files/figure-html/plot-cluster-1.png)

Pixel-wise Welch t-map (red = condition A larger, blue = condition B
larger) with significant cluster contours overlaid in black. The two
clusters correspond to the two injected signal patches, and only those
clusters are FWE-corrected significant at alpha = 0.05; the surrounding
noise does not produce false-positive clusters thanks to the
max-statistic correction.

**Reading the result.** `$observed_t` is the raw t-map. `$clusters` is a
data.frame with one row per candidate cluster, columns: `cluster_id`,
`direction` (`"pos"`/`"neg"`), `mass` (sum of t-values), `size` (pixel
count), `p_value` (mass percentile against the null), `significant`
(logical). `$pos_labels` / `$neg_labels` are integer matrices the same
shape as the image, where each non-zero value tags a cluster (useful for
custom plotting overlays). `$null_distribution` exposes the
per-permutation max masses.

**TFCE alternative.** For threshold-free cluster enhancement, which
integrates across many cluster-forming thresholds and removes the
arbitrary cutoff choice, set `method = "tfce"`. See §6.7 below for the
agreement-map relative.

``` r
clust_tfce <- rel_cluster_test(
  sig_A_loc, sig_B_loc,
  img_dims        = c(32L, 32L),
  method          = "tfce",
  n_permutations  = 500L,
  tfce_n_steps    = 50L,
  alpha           = 0.05,
  seed            = 1L,
  progress        = FALSE
)
plot(clust_tfce)
```

![TFCE-enhanced map (no cluster-forming threshold to choose).
Significant pixels are outlined in black. Compared with the
threshold-based version above, TFCE generally has more power for
diffuse, low-amplitude signal because it integrates evidence across many
thresholds rather than committing to
one.](tutorial_files/figure-html/rel-cluster-tfce-1.png)

TFCE-enhanced map (no cluster-forming threshold to choose). Significant
pixels are outlined in black. Compared with the threshold-based version
above, TFCE generally has more power for diffuse, low-amplitude signal
because it integrates evidence across many thresholds rather than
committing to one.

**Common mistakes.**

- Reading `cluster_threshold = 2.0` as “p \< 0.05”. It is the t-cutoff
  for forming candidate clusters, not the cluster significance
  threshold.
- Trusting cluster p-values with `n_permutations < 1000` (the function
  warns below 100). Tail probabilities below 0.01 need more permutations
  than that.
- Using
  [`pixel_t_test()`](https://olivethree.github.io/rcicrely/reference/pixel_t_test.md)
  standalone for inference (§6.1).

### 6.3 `rel_dissimilarity()`

**What it does.** Two scalar dissimilarity metrics between the two
conditions’ group CIs - Pearson correlation and Euclidean distance -
each with percentile bootstrap CIs.

**When to use it.** “How different are these two conditions overall?”
with uncertainty. Pair with
[`rel_cluster_test()`](https://olivethree.github.io/rcicrely/reference/rel_cluster_test.md)
for spatial information.

``` r
diss <- rel_dissimilarity(
  sig_A, sig_B,
  n_boot   = 500L,
  ci_level = 0.95,
  seed     = 1L,
  progress = FALSE
)
print(diss)
#> <rcicrely representational dissimilarity>
#>   n_boot:               500
#>   CI level:             95%
#>   n_pixels:             1024
#> 
#>   Primary: Euclidean distance between group-mean CIs
#>     Euclidean             = 9.668   [9.293, 10.089]   SE = 0.207
#>     Euclidean / sqrt(n)   = 0.3021   (resolution-normalised)
#> 
#>   [Deprecated - will be removed in v0.3]
#>     Pearson r             = -0.996   [-0.994, -0.991]   SE = 0.001
#>     Note: correlation between two base-subtracted CIs has a
#>     positive baseline from shared image-domain structure and is
#>     not a clean similarity score. Prefer Euclidean distance.
```

``` r
plot(diss)
```

![Bootstrap distributions for Pearson r and Euclidean
distance.](tutorial_files/figure-html/plot-dissim-1.png)

Bootstrap distributions for Pearson r and Euclidean distance.

**Reading the result.** `$correlation`, `$euclidean` are the observed
values; `$ci_cor`, `$ci_dist` the percentile CIs; `$boot_cor`,
`$boot_dist` the full bootstrap distributions.

**Common mistakes.** Reading `$correlation` near 1 as “the conditions
look the same” - identical spatial patterns with different magnitudes
give r ~ 1; the magnitude difference shows up in `$euclidean`.
`$euclidean` scales with sqrt(n_pixels) so it is a within-study metric,
not comparable across resolutions.

**Comparing multiple contrasts side-by-side.** When you have several
[`rel_dissimilarity()`](https://olivethree.github.io/rcicrely/reference/rel_dissimilarity.md)
results and want one figure that shows whether their CIs overlap, use
[`compare_dissimilarity()`](https://olivethree.github.io/rcicrely/reference/compare_dissimilarity.md).
It draws each contrast as a horizontal point with its 95% CI and the
underlying bootstrap density, so the reader can read overlap (or
non-overlap) at a glance. Useful for paper figures.

``` r
d_AB <- rel_dissimilarity(sm_A, sm_B, n_boot = 500L, seed = 1L)
d_CD <- rel_dissimilarity(sm_C, sm_D, n_boot = 500L, seed = 2L)
compare_dissimilarity(
  "A vs B" = d_AB,
  "C vs D" = d_CD
)
```

### 6.4 `run_between()`

**What it does.** Calls both between-condition metrics and bundles the
result.

``` r
between_AB <- run_between(
  sig_A, sig_B,
  img_dims        = c(32L, 32L),
  n_permutations  = 500L,
  n_boot          = 500L,
  seed            = 1L,
  progress        = FALSE
)
print(between_AB)
#> <rcicrely report: between>
#> ---
#> * $results$cluster_test
#> <rcicrely cluster-based permutation test>
#>   N_A = 30, N_B = 30
#>   image dims:        32 x 32
#>   cluster threshold: |t| > 2.00
#>   n_permutations:    500
#>   alpha:             0.05
#>   clusters found:    2 (2 significant)
#> 
#>  cluster_id direction  mass size p_value significant
#>           1       neg -9886  512       0        TRUE
#>           1       pos  9807  512       0        TRUE
#> ---
#> * $results$dissimilarity
#> <rcicrely representational dissimilarity>
#>   n_boot:               500
#>   CI level:             95%
#>   n_pixels:             1024
#> 
#>   Primary: Euclidean distance between group-mean CIs
#>     Euclidean             = 9.668   [9.293, 10.089]   SE = 0.207
#>     Euclidean / sqrt(n)   = 0.3021   (resolution-normalised)
#> 
#>   [Deprecated - will be removed in v0.3]
#>     Pearson r             = -0.996   [-0.994, -0.991]   SE = 0.001
#>     Note: correlation between two base-subtracted CIs has a
#>     positive baseline from shared image-domain structure and is
#>     not a clean similarity score. Prefer Euclidean distance.
```

### 6.5 `infoval()` - producer-level signal vs. chance

**What it does.** For each producer, compares the Frobenius norm of
their raw mask to a reference distribution simulated from random
sign-weighted aggregations of the same noise pool, at that producer’s
trial count. Returns a per-producer z-score against the reference’s
median / MAD, plus the full reference distribution and metadata. A
single function handles both 2IFC and Brief-RC; the difference lives
entirely in what you pass as `noise_matrix`.

**When to use it.** As a complement to the `rel_*()` reliability checks.
[`infoval()`](https://olivethree.github.io/rcicrely/reference/infoval.md)
answers *does each producer’s mask carry more signal than chance?*;
`rel_*()` answers *how reliable / how distinguishable is the group-level
pattern?*. They report different things and a study should usually carry
both.

``` r
# Brief-RC route: infoval() takes the raw signal_matrix from
# ci_from_responses_briefrc() and the noise_matrix used to generate
# the trials.
res <- ci_from_responses_briefrc(responses, rdata_path = "rcic.Rdata",
                                 base_image_path = "base.jpg")
trial_counts <- with(responses,
                     setNames(table(participant_id),
                              names(table(participant_id))))
iv <- infoval(res$signal_matrix, noise_matrix,
              trial_counts = trial_counts,
              iter = 10000L, cache_path = "infoval_cache.rds",
              seed = 1L)
print(iv)
```

**Reading the result.** `$infoval` is a per-producer z-score vector;
`z > 1.96` is roughly one-tailed *p* \< .025 against the chance
reference. `$norms` is the raw observed Frobenius norm per producer.
`$reference` holds the simulated reference norm vectors keyed by trial
count. `$ref_median` and `$ref_mad` are the centre / spread used for z
computation.

#### Troubleshooting low or negative `infoval` z-scores

If you compute
[`infoval()`](https://olivethree.github.io/rcicrely/reference/infoval.md)
and find that *most* or *all* per-producer z-scores sit well below 1.96
(sometimes negative) even though spot checks suggest producers are doing
the task seriously, that is a common pattern, not a sign the data is
bad. Five reasons in order of how often they apply:

1.  **Frobenius norm is a global energy statistic.** It sums squared
    pixel deviations across the *entire* image. Real internal
    representations are usually spatially sparse (eyes, mouth, jaw,
    maybe 10-30% of pixels), so the 70-90% of “background” pixels
    contribute noise of similar magnitude to the chance reference and
    dilute the signal-bearing region. A producer with strong,
    visually-obvious signal in the eyes can have a Frobenius norm only
    marginally above the random reference.

2.  **The reference is strict because it lives in the same subspace.**
    Both the observed mask and the reference are projections onto the
    same low-dimensional sinusoidal noise basis. The reference
    distribution has plenty of overall energy by construction, so the
    only way to clear z = 1.96 is to align signs with a *specific
    subset* of patterns more than chance.

3.  **Per-trial signal is small.** Each 2IFC choice contributes a tiny
    signal increment relative to the per-trial noise amplitude. With 300
    trials the SNR gain is √300 ≈ 17×, but if per-trial signal is ~5% of
    per-trial noise, post-aggregation effective SNR is barely visible to
    a global energy measure.

4.  **Without a face mask, infoVal counts background.**
    [`face_mask()`](https://olivethree.github.io/rcicrely/reference/face_mask.md)
    ships an oval mask that approximates the Schmitz 2024 face region.
    Applying it (`infoval(..., mask = face_mask(c(256, 256)))`)
    concentrates the norm on signal-bearing pixels and typically lifts
    z-scores noticeably.

5.  **Group-level CIs have much higher z than individual CIs.**
    Averaging 20 producers’ masks reduces noise by √20 ≈ 4.5×, so the
    group-mean CI’s effective trial count is `300 × 20 = 6000` for a
    20-producer condition. The infoVal of the group-mean CI (with
    `trial_counts = setNames(300L * n_producers, "group")`) is usually
    5-10× the per-producer median; this is a structural √N consequence
    of averaging, not a defect of the per-producer metric. Brinkman et
    al. (2019) themselves only ever computed infoVal on *individual* CIs
    (mean per-producer infoVal in their empirical 2IFC data was 3.9 lab
    / 2.9 online, with 68% / 54% of producers individually exceeding
    1.96). For group-level reporting they recommended inspecting *the
    distribution of per-producer infoVals contributing to the group CI*
    rather than computing one infoVal on the group-mean CI. Either
    choice is defensible; we show the group-mean CI z below as the
    simplest single number, but the distributional view is what the 2019
    paper endorses.

##### Diagnostic recipe

``` r
sm     <- res$signal_matrix
tc     <- setNames(rep(300L, ncol(sm)), colnames(sm))

# 1. Compare the observed and reference norm distributions
iv  <- infoval(sm, noise_matrix, tc, iter = 1000L,
               cache_path = "iv.rds", seed = 1L)
ref <- iv$reference[[as.character(tc[1])]]
cat(sprintf("observed median = %.4f, reference median = %.4f, %% above = %+.1f%%\n",
            median(iv$norms), median(ref),
            100 * (median(iv$norms) - median(ref)) / median(ref)))

# 2. Apply the face mask
fm  <- face_mask(c(256L, 256L))
iv_masked <- infoval(sm, noise_matrix, tc, mask = fm,
                     iter = 1000L, cache_path = "iv_masked.rds",
                     seed = 1L)
median(iv_masked$infoval)        # typically higher than iv

# 3. Group-mean CI (the headline number to report)
group   <- matrix(rowMeans(sm), ncol = 1, dimnames = list(NULL, "group"))
tc_grp  <- setNames(300L * ncol(sm), "group")
iv_grp  <- infoval(group, noise_matrix, tc_grp, iter = 1000L,
                   seed = 1L)
iv_grp$infoval                   # usually large (e.g. > 5)

# 4. Sanity check with a simulated random responder
set.seed(0)
random_resp <- sample(c(-1, 1), 300L, replace = TRUE)
random_mask <- (noise_matrix %*% random_resp) / 300
iv_rand <- infoval(matrix(random_mask, ncol = 1,
                          dimnames = list(NULL, "rand")),
                   noise_matrix,
                   setNames(300L, "rand"),
                   iter = 1000L, seed = 1L)
iv_rand$infoval                  # should be ~ 0 within MAD noise
```

**Negative z values** mean the observed mask has *less* total Frobenius
energy than the chance reference. This is not “negative signal”; it is
consistent with a producer whose responses align with a structured
low-rank template, since aligning consistently with one template can
reduce the effective dimensionality of the mask and shrink its overall
norm. A negative z is therefore informative: a clearly negative z (say,
\< −2) on a producer who otherwise shows good behaviour usually points
to one of:

- The producer responded coherently but to a *different* feature than
  the trait you intended (idiosyncratic representation).
- The reference distribution is mis-calibrated for your design (see the
  v0.2.1 calibration note in the
  [`?infoval`](https://olivethree.github.io/rcicrely/reference/infoval.md)
  help page; early v0.2.0 versions sampled stim ids with replacement and
  inflated reference norms, biasing per-producer z below zero on
  cooperative producers).
- Per-producer trial count is too low to beat chance for any meaningful
  statistic, even with real signal; fall back to the group-mean infoVal
  as the publishable number.

**Reporting recommendation.** Report **(i)** the group-mean infoVal z
and **(ii)** the proportion of producers whose *individual* z clears
1.96 (with the face mask). Do not report only the median per-producer z
and treat values \< 1.96 as evidence of “no signal”; that conflates the
per-producer noisiness of the metric with the group-level question of
whether the CI carries information.

> **Calibration note.** Since v0.2.1, `simulate_reference_norms()`
> samples stim ids *without* replacement when `n_trials <= n_pool`.
> v0.2.0 over-inflated the reference norm and biased per-producer
> z-scores systematically below zero. If you computed infoVal under
> v0.2.0 and got negative z-scores on cooperative producers, re-run
> after upgrading.

### 6.6 Region-by-region analyses with `face_mask()`

[`face_mask()`](https://olivethree.github.io/rcicrely/reference/face_mask.md)
accepts a `region` argument that returns a logical vector for one of:
`"full"` (default), `"eyes"`, `"nose"`, `"mouth"`, `"upper_face"`,
`"lower_face"`. The companion
[`load_face_mask()`](#section-load-face-mask) reads an externally-
generated mask (e.g., from `webmorphR::mask_oval()` or a hand-painted
PNG) into the same logical-vector format. Pass either to
[`infoval()`](https://olivethree.github.io/rcicrely/reference/infoval.md)
directly via the `mask` argument, or apply it via row-subsetting on a
signal matrix to compute any `rel_*()` metric on a single anatomical
region:

> **Apply the mask symmetrically.** This is the single most important
> rule when working with masks in any reliability or discriminability
> analysis: **the same mask must enter both arms of any comparison**.
> For [`infoval()`](#section-infoval), pass the mask via the `mask`
> argument so observed and reference Frobenius norms are restricted to
> the same pixels. For `rel_*()` reliability metrics, subset the signal
> matrix (`signal_matrix[mask, , drop = FALSE]`) and compute *all*
> metrics on the same subsetted matrix. For
> [`rel_cluster_test()`](#section-cluster) /
> [`rel_dissimilarity()`](#section-dissim), subset *both* condition
> matrices identically before passing them in. Mixing a masked observed
> value with an unmasked reference (or two conditions masked
> differently) yields a number with no defensible interpretation.

``` r
img <- c(256L, 256L)
m_full  <- face_mask(img, region = "full")
m_eyes  <- face_mask(img, region = "eyes")
m_mouth <- face_mask(img, region = "mouth")

# Region-restricted ICC: pass the masked rows of the signal matrix
rel_icc(signal_matrix[m_eyes,  , drop = FALSE])
rel_icc(signal_matrix[m_mouth, , drop = FALSE])

# Region-restricted infoVal: pass the mask directly
infoval(signal_matrix, noise_matrix, trial_counts,
        iter = 1000L, mask = m_eyes,  seed = 1L)
infoval(signal_matrix, noise_matrix, trial_counts,
        iter = 1000L, mask = m_mouth, seed = 1L)
```

Region geometries are heuristic approximations matched to a typical
centred face on a square base image. For non-default base images, tune
the full-face oval via `centre`, `half_width`, `half_height`; sub-region
positions scale relative to that ellipse.

**Region-restricted between-condition inference** uses the same mask but
applied symmetrically to both arms. Because
[`rel_cluster_test()`](https://olivethree.github.io/rcicrely/reference/rel_cluster_test.md)
needs every pixel in the image grid to label connected components, the
way to restrict it is to *zero-fill* outside the region in both signal
matrices: clusters then can only form within the kept region.

``` r
zero_fill <- function(sm, keep) {
  out <- sm
  out[!keep, ] <- 0
  out
}

m_eyes <- face_mask(img_dims, region = "eyes")
clust_eyes <- rel_cluster_test(
  zero_fill(sm_A, m_eyes),
  zero_fill(sm_B, m_eyes),
  img_dims          = img_dims,
  cluster_threshold = 2.0,
  n_permutations    = 1000L,
  seed              = 1L
)
```

§7 walks through both region-restricted ICC and region-restricted
cluster contrasts on a real dataset with figures.

**Critical-review note.** Region masks are heuristic. They do not
account for individual differences in face geometry (the base face’s
eyes, nose, mouth occupy different pixels on different photographs), and
they assume the typical centred-face layout shipped with most rcicr
stimulus pools. For studies using non-centred or
differently-proportioned base faces, define a custom mask via
`ellipse_mask()` (internal) or by hand-curating a logical pixel vector.
The geometries shipped here are calibrated to KDEF- style 256x256 male
faces and should reproduce reasonably for most published RC datasets.

#### Loading an external mask: `load_face_mask()`

If your mask comes from outside R (a hand-painted PNG, a
`webmorphR::mask_oval()` output, or any other binary image),
[`load_face_mask()`](https://olivethree.github.io/rcicrely/reference/load_face_mask.md)
converts it to the logical-vector format
[`face_mask()`](https://olivethree.github.io/rcicrely/reference/face_mask.md)
returns:

``` r
fm <- load_face_mask("masks/oval_256.png",
                     threshold     = 0.5,
                     expected_dims = c(256L, 256L))
mean(fm)               # fraction of pixels included
infoval(signal_matrix, noise_matrix, trial_counts,
        mask = fm, iter = 1000L, seed = 1L)
```

White / light pixels (luminance above `threshold`, default 0.5) become
`TRUE`; dark pixels become `FALSE`. RGB images are converted to
luminance via ITU-R BT.709 weights. Pass `expected_dims` to abort on
mask-vs-signal dimension mismatches at the input boundary; if the saved
image uses the opposite convention (black = include),
`!load_face_mask(...)` inverts it.

### 6.7 Visualising producer agreement: `plot_agreement_map()`

A complement to the group-mean CI image: a heatmap of where producers in
a single condition *agree* on the direction of signal. Each pixel’s
colour is proportional to a one-sample t-statistic across producers
(mean / standard error), so positive saturation means agreement on
positive signal and negative saturation means agreement on negative
signal.

``` r
# Continuous t-map (default)
plot_agreement_map(signal_matrix)

# Region-restricted, mouth only
plot_agreement_map(signal_matrix,
                   mask = face_mask(c(256L, 256L), region = "mouth"))

# Thresholded view (only pixels with |t| > 2.0 are coloured)
plot_agreement_map(signal_matrix, threshold = 2.0)
```

Use this for **descriptive** figures of where producers agree within a
condition; pair with
[`rel_cluster_test()`](https://olivethree.github.io/rcicrely/reference/rel_cluster_test.md)
(which is the **inferential** counterpart for between-condition
pixel-level contrasts).

## 7. Demo: applying `rcicrely` to Oliveira et al. (2019), Study 1

This section walks the package end-to-end on a published 2IFC dataset,
in the order a researcher would actually run the analyses. Every output
shown is real (chunks are `eval = FALSE` because the data is not
bundled).

> Oliveira, M., Garcia-Marques, T., Dotsch, R., & Garcia-Marques, L.
> (2019). Dominance and competence face to face: Dissociations obtained
> with a reverse correlation approach. *European Journal of Social
> Psychology*. <https://doi.org/10.1002/ejsp.2569>. Open data:
> <https://doi.org/10.17605/osf.io/hr5pd>.

In Study 1, 200 participants completed a 2IFC reverse-correlation task
with **300 trials each** on a 256x256 male base face, across 10 trait
conditions in a between-subjects design (20 producers per trait):
Dominant, Submissive, Trust, Untrust, Friendly, Unfriendly, Intelligent,
Unintelligent, Competent, Incompetent.

### 7.1 Setup: per-trait signal matrices

Load the response data and the noise pool. The pattern below works for
any 2IFC study with `subject`, `trial`, `stimulus`, `response` columns
(the column names vary between datasets; rename to fit).

``` r
library(rcicrely)
library(dplyr)

# Trial-level data: one row per (subject, trial) pair.
responses <- read.csv2("study1data.csv")
glimpse(responses)
#> $ subject  <int> ...           # producer id
#> $ trait    <chr> "TRUST" ...    # condition label
#> $ trial    <int> 1, 2, 3, ...
#> $ stimulus <int> ...            # noise pattern shown on this trial
#> $ response <int> +1 / -1        # +1 = oriented chosen, -1 = inverted

# Pixels x pool-size matrix of noise patterns used at stimulus generation.
noise_matrix <- read_noise_matrix("rcic_stimuli.Rdata", baseimage = "male")
dim(noise_matrix)
#> 65536 x 300
```

Build the signal matrix for **one trait** (Trust). A producer’s raw mask
is the sign-weighted average of the noise patterns they chose:

> mask = (noise\[, chosen_stimuli\] %\*% chosen_responses) / n_trials

``` r
# Filter to Trust trials, ordered by producer + trial number.
trust_trials <- responses %>%
  filter(trait == "TRUST") %>%
  arrange(subject, trial)

producers <- unique(trust_trials$subject)             # 20 ids
n_pixels  <- nrow(noise_matrix)
sm_trust  <- matrix(NA_real_, nrow = n_pixels, ncol = length(producers))
colnames(sm_trust) <- as.character(producers)

# One column per producer: their raw mask.
for (i in seq_along(producers)) {
  this_producer <- trust_trials %>% filter(subject == producers[i])
  sm_trust[, i] <- (noise_matrix[, this_producer$stimulus] %*%
                      this_producer$response) / nrow(this_producer)
}
attr(sm_trust, "img_dims") <- c(256L, 256L)
attr(sm_trust, "source")   <- "raw"
```

Repeat the same block (lines 4 onwards) for `"DOMINANT"`, `"COMPETENT"`,
and `"FRIENDLY"` to get `sm_dominant`, `sm_competent`, `sm_friendly`. To
do all 10 traits in one pass, wrap the loop in a function:

``` r
build_sm <- function(trait_label) {
  sub <- responses %>%
    filter(trait == trait_label) %>%
    arrange(subject, trial)
  ids <- unique(sub$subject)
  m   <- matrix(NA_real_, nrow = n_pixels, ncol = length(ids))
  colnames(m) <- as.character(ids)
  for (i in seq_along(ids)) {
    pr      <- sub %>% filter(subject == ids[i])
    m[, i]  <- (noise_matrix[, pr$stimulus] %*% pr$response) / nrow(pr)
  }
  attr(m, "img_dims") <- c(256L, 256L)
  attr(m, "source")   <- "raw"
  m
}

signal_by_trait <- lapply(unique(responses$trait), build_sm)
names(signal_by_trait) <- unique(responses$trait)

sm_trust     <- signal_by_trait[["TRUST"]]
sm_dominant  <- signal_by_trait[["DOMINANT"]]
sm_competent <- signal_by_trait[["COMPETENT"]]
sm_friendly  <- signal_by_trait[["FRIENDLY"]]
```

### 7.2 Within-condition reliability

Did the 20 producers in each condition agree?

``` r
rel_split_half(sm_trust,    n_permutations = 500L, seed = 1L)
#> Trust:    r_hh = 0.405, r_sb = 0.577, 95% CI on r_sb = [0.539, 0.616]
rel_split_half(sm_dominant, n_permutations = 500L, seed = 1L)
#> Dominant: r_hh = 0.262, r_sb = 0.415, 95% CI on r_sb = [0.380, 0.450]

rel_icc(sm_trust)
#> ICC(3,1) = 0.1008,  ICC(3,k) = 0.6917
rel_icc(sm_dominant)
#> ICC(3,1) = 0.0489,  ICC(3,k) = 0.5068
```

Trust replicates with clearly higher group-level reliability than
Dominant on both metrics. ICC(3,1) (the resolution-comparable
single-producer reliability) and Spearman-Brown-projected `r_SB` agree
on the ordering. ICC(3,k) is high in absolute terms because of the
resolution asymptote at large pixel counts; ICC(3,1) is the cross-study
comparable number.

LOO inspection on Trust:

``` r
loo <- rel_loo(sm_trust, flag_method = "mad")
head(rel_loo_z(loo), 3)
#>   participant_id  correlation  z_score  flag
#> 1          18043       0.9905    -1.68  FALSE
#> 2          18051       0.9911    -0.99  FALSE
#> 3          18131       0.9913    -0.79  FALSE
```

Median raw `r_loo` is .992 (MAD = .001), near 1 by construction. Zero
producers cleared the MAD flag threshold, but the z-scored ordering
still ranks producers by influence on the group CI: the top of the table
is the inspection priority list.

### 7.3 Per-producer signal vs. chance

[`infoval()`](https://olivethree.github.io/rcicrely/reference/infoval.md)
answers “does each producer’s mask carry more than chance?”. With 300
trials per producer the trial-count-matched reference is simulated at
300:

``` r
# Trial count per producer, computed from the data (named numeric).
trial_counts <- responses %>%
  filter(trait == "TRUST") %>%
  count(subject) %>%
  pull(n, name = subject)

infoval(sm_trust, noise_matrix, trial_counts,
        iter = 1000L, seed = 1L)
#> per-producer: median z = 0.35, range [-1.69, 2.30], n above 1.96 = 2/20

# Group-mean CI (effective n_trials = 300 x 20 = 6000):
group_ci <- matrix(rowMeans(sm_trust), ncol = 1)
colnames(group_ci) <- "TRUST_group"
group_n  <- c(TRUST_group = 6000L)

infoval(group_ci, noise_matrix, group_n, iter = 1000L, seed = 1L)
#> group-mean CI z = 11.39

# Face-mask lift on per-producer z:
infoval(sm_trust, noise_matrix, trial_counts,
        mask = face_mask(c(256L, 256L)),
        iter = 1000L, seed = 1L)
#> masked median = 0.50  (vs unmasked median = 0.35; lift = +0.15)
```

This is exactly the pattern §6.5 describes: per-producer median z is
0.35 with only 2/20 producers individually above 1.96, *yet* the
group-mean z is 11.39, unmistakable signal. The published Trust/Friendly
correlation of .69 (Table 2 of the paper) and the trait dissociations
visible in Figure 1 are consistent with the group-level z, not with
reading the per-producer median as “weak signal”. Recommended report:
“the group-mean classification image showed strong informational value
(z = 11.39, n_eff = 6000); 2 of 20 producers individually exceeded z =
1.96.”

### 7.4 Where in the image do conditions disagree?

The Dominant vs Competent contrast yields spatially-localized,
FWE-controlled clusters:

``` r
rel_cluster_test(sm_dominant, sm_competent,
                 img_dims          = c(256L, 256L),
                 cluster_threshold = 2.0,
                 n_permutations    = 1000L,
                 seed              = 1L)
```

![Cluster test, Dominant vs Competent group CIs (Oliveira et al., 2019,
Study 1). Pixel-wise Welch t-map rendered as a translucent diverging
overlay on the male base face (positive/red = Dominant larger,
negative/blue = Competent larger). Significant clusters at FWE-corrected
alpha = .05 are outlined in
black.](figures/oliveira_2019/cluster_dominant_vs_competent.png)

Cluster test, Dominant vs Competent group CIs (Oliveira et al., 2019,
Study 1). Pixel-wise Welch t-map rendered as a translucent diverging
overlay on the male base face (positive/red = Dominant larger,
negative/blue = Competent larger). Significant clusters at FWE-corrected
alpha = .05 are outlined in black.

### 7.5 How far apart are conditions overall?

[`rel_dissimilarity()`](https://olivethree.github.io/rcicrely/reference/rel_dissimilarity.md)
complements the spatial picture with a single-number magnitude.
Comparing a “similar” trait pair against a “dissimilar” one:

``` r
rel_dissimilarity(sm_trust,    sm_friendly,  n_boot = 500L, seed = 1L)
#> Euclidean = 16.81, 95% CI [19.26, 24.86], normalised = 0.0657
rel_dissimilarity(sm_dominant, sm_competent, n_boot = 500L, seed = 2L)
#> Euclidean = 21.08, 95% CI [21.53, 30.75], normalised = 0.0823
```

Dominant-vs-Competent yields a substantially larger Euclidean distance
than Trust-vs-Friendly, the same direction as the published group-CI
Pearson correlations (.69 for Trust × Friendly, −.22 for Dominant ×
Competent in Table 2 of the paper). Pearson is the within-shape
question; Euclidean is the how-far-apart question.

### 7.6 Where on the face do producers agree, by region?

Region-restricted ICC(3,1) on the same Trust signal matrix shows where
producer agreement concentrates:

``` r
# rel_icc(sm_trust[face_mask(c(256L,256L), region = R), , drop = FALSE])
#> region        n_pixels   ICC(3,1)   ICC(3,k)
#> full           32,192     0.109      0.710
#> eyes            1,022     0.069      0.595
#> nose              770     0.028      0.363
#> mouth             966     0.092      0.670
#> upper_face     16,096     0.087      0.655
#> lower_face     16,096     0.130      0.749
```

Agreement is strongest in the lower-face region (mouth + jaw), weakest
in the nose. The pattern aligns with the visible features of the
published Trust group CI (Figure 1 of the paper). Region- restricted
infoVal works as a sensitivity check on the same data:

``` r
# infoval(sm_trust, noise_matrix, tc, mask = face_mask(c(256L,256L), region = R), ...)
#> region   median z   range          n above 1.96
#> full     +0.51      [-1.81, +3.13]      4/20
#> eyes     +0.55      [-0.98, +1.56]      0/20
#> mouth    +0.52      [-0.76, +4.67]      3/20
```

Eye and mouth regions both replicate the full-face median; the eye
region’s per-producer range is narrower, suggesting producers agree on
direction in the eye region but the energy is modest.

Region-restricted **between-condition** clusters (Trustworthy vs
Dominant), zero-filling outside each region so clusters can only form
within it:

``` r
img_dims <- c(256L, 256L)
m_eyes   <- face_mask(img_dims, region = "eyes")
m_mouth  <- face_mask(img_dims, region = "mouth")

# rel_cluster_test() needs the full image grid to label connected
# components, so we restrict by "blanking out" pixels outside the
# region: set them to zero. Clusters can then only form inside the
# kept region.
zero_fill <- function(sm, keep_pixels) {
  out <- sm
  out[!keep_pixels, ] <- 0
  out
}

# Eyes region.
sm_trust_eyes    <- zero_fill(sm_trust,    m_eyes)
sm_dominant_eyes <- zero_fill(sm_dominant, m_eyes)
rel_cluster_test(sm_trust_eyes, sm_dominant_eyes,
                 img_dims          = img_dims,
                 cluster_threshold = 2.0,
                 n_permutations    = 1000L,
                 seed              = 1L)

# Mouth region.
sm_trust_mouth    <- zero_fill(sm_trust,    m_mouth)
sm_dominant_mouth <- zero_fill(sm_dominant, m_mouth)
rel_cluster_test(sm_trust_mouth, sm_dominant_mouth,
                 img_dims          = img_dims,
                 cluster_threshold = 2.0,
                 n_permutations    = 1000L,
                 seed              = 2L)
```

![Region-restricted cluster test, Trustworthy vs Dominant group CIs
(Oliveira et al., 2019, Study 1). Pixel-wise Welch t-map rendered as a
translucent diverging overlay on the male base face; positive/red =
Trustworthy larger, negative/blue = Dominant larger. Eyes region (left
panel) and mouth region (right panel). Significant clusters at
FWE-corrected alpha = .05 are outlined in
black.](figures/oliveira_2019/cluster_trust_vs_dominant_eyes_mouth.png)

Region-restricted cluster test, Trustworthy vs Dominant group CIs
(Oliveira et al., 2019, Study 1). Pixel-wise Welch t-map rendered as a
translucent diverging overlay on the male base face; positive/red =
Trustworthy larger, negative/blue = Dominant larger. Eyes region (left
panel) and mouth region (right panel). Significant clusters at
FWE-corrected alpha = .05 are outlined in black.

## 8. Brief-RC end-to-end

Schmitz, Rougier & Yzerbyt (2024) introduced Brief-RC 12: each trial
shows 12 noisy faces (6 oriented, 6 inverted, drawn from a shared noise
pool); participant picks one. The recorded format is one row per trial,
columns `participant_id`, `trial`, `stimulus` (chosen pool id),
`response` (`+1` if oriented, `-1` if inverted). `rcicrely`’s
implementation is native (no rcicr `_brief` functions, which don’t exist
in the upstream rcicr v1.0.1).

The Schmitz materials distribute the noise pool as a plain-text
`noise_matrix.txt` (pixels × pool size, whitespace-delimited), so that’s
the typical entry point if you replicate from their OSF repo. If you
generate the pool yourself with `rcicr`, the captured return value can
be cached as `.rds` for fast reloads (see step 1 below).

``` r
# 1a. Replicating Schmitz et al.: read the published noise_matrix.txt
#     directly. read_noise_matrix() autodetects the format.
noise <- read_noise_matrix("data/noise_matrix.txt")

# 1b. Generating your own pool with rcicr (one-off, slow):
rcicr::generateStimuli2IFC(
  base_face_files     = list(base = "data/base.jpg"),
  n_trials            = 600L,        # becomes the Brief-RC pool size
  img_size            = 256L,
  stimulus_path       = "data/stim",
  label               = "mystudy",
  ncores              = 1L,          # macOS-safe default
  return_as_dataframe = TRUE,        # capture the noise matrix
  seed                = 1L
) |> saveRDS("data/noise_matrix.rds")  # cache for fast reloads

# 2. Collect Brief-RC responses in the per-trial format.

# 3. Compute individual CIs.
res <- ci_from_responses_briefrc(
  responses       = my_brief_rc_responses,
  noise_matrix    = noise,           # or readRDS("data/noise_matrix.rds")
  base_image_path = "data/base.jpg",
  method          = "briefrc12",
  scaling         = "matched"        # also rendered for plotting
)

# 4. Reliability assessment.
within  <- run_within(res$signal_matrix, seed = 1L)
between <- run_between(res_A$signal_matrix, res_B$signal_matrix,
                       seed = 1L)

# 5. Save rendered CIs to PNG (visualisation only; do not feed to rel_*).
if (requireNamespace("png", quietly = TRUE)) {
  for (j in seq_along(res$participants)) {
    img <- matrix(res$rendered_ci[, j], res$img_dims[1], res$img_dims[2])
    img <- pmin(pmax(img, 0), 1)
    png::writePNG(img, sprintf("plots/p%02d.png", j))
  }
}
```

**Brief-RC `infoVal`** ships in v0.2 as
[`infoval()`](https://olivethree.github.io/rcicrely/reference/infoval.md)
(see §6.5). The reference distribution is matched to each producer’s
recorded trial count (number of `(stim, +/-1)` contributions entering
the mask), not to the pool size. For 2IFC the two are equal so the
original rcicr default is fine; for Brief-RC each producer typically
records fewer choices than there are pool items, so a pool-size
reference biases infoVal downward.
[`infoval()`](https://olivethree.github.io/rcicrely/reference/infoval.md)
closes that calibration gap. Pass the **raw mask**
(`res$signal_matrix`), not `res$rendered_ci` or a PNG.

Two notes on interpretation. First, “subset” here means the number of
**recorded** choices that enter the mask (one `(stim, +/-1)` per trial),
not the producer’s perceptual exposure: Brief-RC trials are
*cognitively* richer than 2IFC trials (12 noisy faces of context vs. 2),
so producers do not see less of the pool. The trial-count match is the
right calibration for infoVal’s magnitude statistic, not a claim about
cognitive exposure. Second, infoVal measures the **magnitude** of the
mask (Frobenius norm) against a chance reference, so it does not
distinguish a more accurate mask from a larger one. Treat absolute
Brief-RC vs 2IFC z-score differences with caution;
[`rel_split_half()`](https://olivethree.github.io/rcicrely/reference/rel_split_half.md)/[`rel_icc()`](https://olivethree.github.io/rcicrely/reference/rel_icc.md)
(signal stability) and
[`rel_cluster_test()`](https://olivethree.github.io/rcicrely/reference/rel_cluster_test.md)/[`rel_dissimilarity()`](https://olivethree.github.io/rcicrely/reference/rel_dissimilarity.md)
(between-condition separability) are the right complements when pattern
accuracy matters more than magnitude.

## 9. Interpreting results and common pitfalls

**Sample size.** `rcicrely` aborts below 4 producers per condition and
warns below 30; a sample of N \>= 60 producers is recommended for stable
reliability assessment. At N = 10 the metrics are computable but their
CIs are wide.

**Split-half.** Think of `r_hh` ~ 0.3 as weak, 0.5 as moderate, 0.7+ as
strong producer-level consistency. `r_SB` is what to report when asked
about the full-sample CI.

**LOO flagging.** A flag is a starting point for investigation, not a
verdict. Use `flag_method = "mad"` if the SD rule is dominated by the
very outliers it tries to find.

**ICC(3,*) vs ICC(2,*).** Numerically close at high pixel counts;
ICC(3,\*) is the correctly-specified model (pixels fixed, participants
random). Use `variants = c("3_1", "3_k", "2_1", "2_k")` when both are
needed.

**Cluster threshold != p-value.** `cluster_threshold = 2.0` is the
t-cutoff for forming clusters. Significance is decided by mass relative
to the null. The default is roughly two-tailed p \< 0.05 at large df,
but at small N it under-rejects - inspect the t histogram.

For raw-vs-rendered considerations and the `infoVal` scaling caveat, see
§3.2.

## 10. References

- Brinkman, L., Todorov, A., & Dotsch, R. (2017). Visualising mental
  representations: A primer on noise-based reverse correlation in social
  psychology. *European Review of Social Psychology*, 28(1), 333-361.
  <https://doi.org/10.1080/10463283.2017.1381469>
- Brinkman, L., Goffin, S., van de Schoot, R., van Haren, N. E. M.,
  Dotsch, R., & Aarts, H. (2019). Quantifying the informational value of
  classification images. *Behavior Research Methods*, 51(5), 2059-2073.
  <https://doi.org/10.3758/s13428-019-01232-2>
- Cone, J., Brown-Iannuzzi, J. L., Lei, R., & Dotsch, R. (2021). Type I
  Error Is Inflated in the Two-Phase Reverse Correlation Procedure.
  *Social Psychological and Personality Science*, 12(5), 760-768.
  <https://doi.org/10.1177/1948550620938616>
- Dotsch, R. (2023). *rcicr*: Reverse-Correlation Image-Classification
  Toolbox. R package v1.0.1. <https://github.com/rdotsch/rcicr>
- Maris, E., & Oostenveld, R. (2007). Nonparametric statistical testing
  of EEG- and MEG-data. *Journal of Neuroscience Methods*, 164(1),
  177-190. <https://doi.org/10.1016/j.jneumeth.2007.03.024>
- McGraw, K. O., & Wong, S. P. (1996). Forming inferences about some
  intraclass correlation coefficients. *Psychological Methods*, 1(1),
  30-46. <https://doi.org/10.1037/1082-989X.1.1.30>
- Nichols, T. E., & Holmes, A. P. (2002). Nonparametric permutation
  tests for functional neuroimaging: A primer with examples. *Human
  Brain Mapping*, 15(1), 1-25. <https://doi.org/10.1002/hbm.1058>
- Oliveira, M., Garcia-Marques, T., Dotsch, R., & Garcia-Marques, L.
  (2019). Dominance and competence face to face: Dissociations obtained
  with a reverse correlation approach. *European Journal of Social
  Psychology*. <https://doi.org/10.1002/ejsp.2569>. Open data and
  analysis scripts: <https://doi.org/10.17605/osf.io/hr5pd>.
- Schmitz, M., Rougier, M., & Yzerbyt, V. (2024). Introducing the brief
  reverse correlation: An improved tool to assess visual
  representations. *European Journal of Social Psychology*.
  <https://doi.org/10.1002/ejsp.3100>
- Shrout, P. E., & Fleiss, J. L. (1979). Intraclass correlations: Uses
  in assessing rater reliability. *Psychological Bulletin*, 86(2),
  420-428. <https://doi.org/10.1037/0033-2909.86.2.420>
