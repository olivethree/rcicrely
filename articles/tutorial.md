# rcicrely user guide

This guide is written for two readers in parallel: the researcher who
wants to know “which function do I call for the question I have?”, and
the methodologically-curious user who also wants to know “and why should
I trust what it returns?”. Each function gets a short plain-language
description, the actual signature, a worked example, how to read the
output, and the specific footguns to avoid.

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

Critically, none of this requires a phase-2 trait-rating study. The
metrics read directly from the pixel-level signal produced by the
original producers, sidestepping the two-phase design Cone et al. (2021)
flagged.

### `infoVal` and `rcicrely` together

`rcicrely` does not compute `infoVal` itself (Brief-RC `infoVal` is
deferred to a later release; for 2IFC, use
[`rcicr::computeInfoVal2IFC()`](https://rdrr.io/pkg/rcicr/man/computeInfoVal2IFC.html)
directly). The two are complementary: `infoVal` tests presence of
signal; `rcicrely` tests reliability and discriminability of signal.

If you use
[`rcicr::computeInfoVal2IFC()`](https://rdrr.io/pkg/rcicr/man/computeInfoVal2IFC.html)
on the output of
[`rcicr::batchGenerateCI2IFC()`](https://rdrr.io/pkg/rcicr/man/batchGenerateCI2IFC.html),
you do not need to worry about scaling: the function pulls the raw `$ci`
element out of the CI-list internally
(`norm(matrix(target_ci[["ci"]]), "f")` in its source). The scaling
settings affect the rendered CI but not the infoVal computation.

The scaling caveat applies when you compute `infoVal` outside that path:

- **Brief-RC** has no canonical `infoVal` function. Anyone computing it
  for Brief-RC is rolling their own, feed the raw mask
  (`res$signal_matrix` from
  [`ci_from_responses_briefrc()`](https://olivethree.github.io/rcicrely/reference/ci_from_responses_briefrc.md)),
  not `res$rendered_ci` or anything PNG-derived.
- **PNG-derived data**: you cannot pass
  [`read_cis()`](https://olivethree.github.io/rcicrely/reference/read_cis.md)
  output to `computeInfoVal2IFC()` directly (it expects an rcicr
  CI-list, not a bare matrix). If you reconstruct one or hand-roll the
  norm calculation, use the raw mask, not the PNG pixel values minus
  base, which are `scaling(mask)`.

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

### 3.1 Raw mask vs. rendered CI

A CI exists in two forms.

The **raw mask** is the participant’s noise contribution before any
display transformation. It is what the math is about. For Brief-RC, it
is `(noise_pool[, chosen_ids] %*% signed_responses) / n_unique`
(Schmitz’s `genMask()`). For 2IFC, it is rcicr’s `$ci` element. Values
are typically in a narrow band around zero (~ +/- 0.05).

The **rendered CI** is what gets saved to PNG. It is
`base + scaling(mask)`, where `scaling()` may stretch the mask to the
base image’s range (“matched”), apply a constant gain, or do nothing.
Scaling is a *display* transformation, not a statistical one.

Subtracting the base from a rendered CI does **not** recover the raw
mask:

    PNG_pixels - base = scaling(mask)   (NOT  mask)

|                                                                                             Metric | Behavior under scaling                                   |
|---------------------------------------------------------------------------------------------------:|:---------------------------------------------------------|
|                                                                        `rel_split_half`, `rel_loo` | Survives uniform scaling. Breaks under per-CI “matched”. |
|                                                               `rel_dissimilarity` (Pearson r half) | Same.                                                    |
|                                                               `rel_dissimilarity` (Euclidean half) | **Distorted** by any scaling.                            |
|                                                                                          `rel_icc` | **Distorted** by any scaling.                            |
|                                                                 `pixel_t_test`, `rel_cluster_test` | **Distorted** by per-CI scaling.                         |
| [`rcicr::computeInfoVal2IFC()`](https://rdrr.io/pkg/rcicr/man/computeInfoVal2IFC.html) (canonical) | Uses raw `$ci` internally, unaffected.                   |
|                                                           Hand-rolled `infoVal` (Brief-RC, custom) | **Distorted** by any scaling. Use the raw mask.          |

### 3.2 Two paths to the signal matrix

    Mode 1: PNGs on disk           Mode 2: raw responses
              |                              |
       read_cis()                    ci_from_responses_2ifc() / _briefrc()
              |                              |
       extract_signal()                      |
              |                              |
       scaling(mask)                       mask                    <-- raw
              |                              |
              +---------> rel_*() <----------+

**Mode 2 returns the raw mask**; Mode 1 cannot, because PNGs on disk
have always been rendered. Prefer Mode 2 when raw responses are
available. If you must use Mode 1, the package emits a one-time warning
per session reminding you of the caveat. To silence:

``` r
# acknowledge per-call:
signal <- load_signal_matrix("data/cis/", "data/base.jpg",
                             acknowledge_scaling = TRUE)

# or silence for the whole session:
options(rcicrely.silence_scaling_warning = TRUE)
```

If you generated the PNGs yourself with
`rcicr::generateCI2IFC(scaling = "none")`, the rendered output is
effectively raw, and the warning is genuinely a false positive.

### 3.3 A concrete demonstration

Here is what scaling does to the same underlying mask.

``` r
set.seed(1)
n_pix <- 32L * 32L
n_p   <- 20L
mask  <- matrix(rnorm(n_pix * n_p, sd = 0.02), n_pix, n_p)

# Render: stretch each column to [-0.5, 0.5] (Schmitz "matched")
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

## 4. Function reference: I/O

### 4.1 `read_cis()`

**What it does.** Reads every PNG/JPEG in a directory, converts each to
grayscale, returns a pixels by participants numeric matrix.

**When to use it.** You have one CI image per producer already saved on
disk and you want to load them into R as a numeric matrix. Pair with
[`extract_signal()`](#section-extract-signal) to get a base-subtracted
signal matrix.

``` r
# Make a tiny synthetic dataset on disk to demonstrate against.
demo_dir <- tempfile("ci_demo_"); dir.create(demo_dir)
base_img <- matrix(0.5, 32L, 32L)
png_path <- function(i) file.path(demo_dir, sprintf("p%02d.png", i))
if (requireNamespace("png", quietly = TRUE)) {
  for (i in 1:6) {
    img <- pmin(pmax(base_img + matrix(rnorm(32 * 32, sd = 0.05), 32, 32),
                     0), 1)
    png::writePNG(img, png_path(i))
  }
}
```

``` r
cis <- read_cis(demo_dir, acknowledge_scaling = TRUE)
dim(cis)
#> [1] 1024    6
attr(cis, "img_dims")
#> [1] 32 32
```

**Reading the result.** `nrow = n_pixels`, `ncol = n_files`. The
`img_dims` attribute carries `c(nrow, ncol)` of the original image so
plotting helpers know how to reshape.

**Common mistakes.** Forgetting to subtract the base before passing to a
`rel_*` function (use [`extract_signal()`](#section-extract-signal) or
the convenience wrapper
[`load_signal_matrix()`](#section-load-signal-matrix)). Mixing CIs of
different sizes - aborts loudly. Trusting reliability numbers without
reading the raw vs. rendered note (chapter 3).

### 4.2 `extract_signal()`

**What it does.** Subtracts the base image from each column of a CI
matrix, so reliability metrics see only the participant-specific
contribution.

**When to use it.** You have the output of [`read_cis()`](#read_cis) and
you want a base-subtracted signal matrix.

``` r
base_path <- file.path(demo_dir, "_base.png")
png::writePNG(base_img, base_path)
signal <- extract_signal(cis, base_image_path = base_path,
                         acknowledge_scaling = TRUE)
dim(signal); attr(signal, "img_dims")
#> [1] 1024    6
#> [1] 32 32
```

**Reading the result.** Same shape as `cis`, with the `img_dims`
attribute preserved. Column names propagate.

**Common mistakes.** Pointing at the wrong base (a generated CI PNG of
the right size will pass the dimension check). Treating the result as
the raw mask; it is `scaling(mask)` if PNGs were rendered with any
scaling. See chapter 3.

### 4.3 `load_signal_matrix()`

**What it does.** Convenience wrapper that calls
[`read_cis()`](#read_cis) and
[`extract_signal()`](#section-extract-signal) in one go.

**When to use it.** Most of the time, when you don’t need to do anything
between reading and subtracting (e.g. masking pixels, swapping the
base).

``` r
signal <- load_signal_matrix(
  dir                 = demo_dir,
  base_image_path     = base_path,
  pattern             = "^p[0-9]+\\.png$",  # exclude _base.png
  acknowledge_scaling = TRUE
)
dim(signal)
#> [1] 1024    6
```

### 4.4 `read_noise_matrix()`

**What it does.** Loads the pixels by trials noise pool used at stimulus
generation. Three formats: rcicr `.Rdata` (reconstructs from
`stimuli_params` + `p`), saved `.rds`, or whitespace-delimited text.

**When to use it.** You are computing Brief-RC CIs and want to load the
noise pool once, then reuse it across multiple
[`ci_from_responses_briefrc()`](#section-cifr-briefrc) calls.

``` r
# from the rcicr rdata (slow - reconstructs each trial via
# rcicr::generateNoiseImage):
noise <- read_noise_matrix("data/rcicr_stimuli.Rdata")

# from a cached rds you saved from the rcicr return value (fast):
noise <- read_noise_matrix("data/noise_matrix.rds")
```

**Common mistakes.** Assuming the rcicr `.Rdata` contains a `stimuli`
object (it does not - the pixels-by-trials matrix is the *return value*
of `generateStimuli2IFC(return_as_dataframe = TRUE)`, never saved to the
file). Cache to `.rds` once and load from there. Multiple stale `.rdata`
files in one directory: rcicr timestamps filenames; pick the most recent
by mtime if you have several.

### 4.5 `ci_from_responses_2ifc()`

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

### 4.6 `ci_from_responses_briefrc()`

**What it does.** Native Brief-RC 12 implementation of Schmitz’s
`genMask()` formula. No `rcicr::*_brief` calls (those don’t exist in
canonical rcicr v1.0.1).

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
  scaling         = "matched"        # Schmitz default for visualisation
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

### 5.2 `rel_loo()`

**What it does.** For each producer, correlates the full-sample group CI
with the group CI computed without that producer. Producers with
unusually low LOO correlation are flagged.

**When to use it.** Spotting producers whose CI sits far from the group
pattern.

``` r
r_loo <- rel_loo(sig_A, flag_threshold = 2.5, flag_method = "sd")
print(r_loo)
#> <rcicrely leave-one-out sensitivity>
#>   N producers:        30
#>   mean LOO r:         1.000
#>   SD:                 0.000
#>   median LOO r:       1.000
#>   MAD:                0.000
#>   flag rule:          r < mean - 2.50 * SD   (=> r < 1.000)
#>   flagged producers:  none
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

![Per-producer LOO
correlation.](tutorial_files/figure-html/plot-loo-1.png)

Per-producer LOO correlation.

**Reading the result.** `$correlations` is the per-producer vector.
`$flagged` lists ids below threshold. `$mean_r` / `$sd_r` / `$median_r`
/ `$mad_r` are the centre/spread under both rules so you can compare.
`$flag_method`, `$flag_threshold` record what was used.

**Common mistakes.** Treating `$flagged` as “drop these producers” -
investigate first (Cone et al. 2021 caveat: a flagged producer may have
a genuinely atypical mental representation). Cross-check with
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
#>   ICC(3,1):      0.8788
#>   ICC(3,k):      0.9954
```

To also report ICC(2,\*) for a reviewer who expects it:

``` r
r_icc_full <- rel_icc(sig_A, variants = c("3_1", "3_k", "2_1", "2_k"))
print(r_icc_full)
#> <rcicrely ICC>
#>   model:        ICC(3,*) / two-way mixed; pixels fixed, participants random
#>   N targets:    1024 pixels
#>   N raters:     30 participants
#> 
#>   ICC(3,1):      0.8788
#>   ICC(3,k):      0.9954
#>   ICC(2,1):      0.8565
#>   ICC(2,k):      0.9944
```

**Reading the result.** `$icc_3_1`, `$icc_3_k`, `$icc_2_1`, `$icc_2_k`
(the latter two NA unless requested). `$ms_rows`, `$ms_cols`,
`$ms_error` are the underlying ANOVA mean squares for transparency.
`$model` is a description string the print method shows.

**Common mistakes.** Asking for ICC(2,*) because a reviewer expects it.
ICC(2,*) treats pixels as a random sample from a pixel population, which
they aren’t. The numbers are usually close to ICC(3,\*) at high pixel
counts but the model is mis-specified - report both with
`variants = c("3_1", "3_k", "2_1", "2_k")` if asked.

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
#> * $results$loo
#> <rcicrely leave-one-out sensitivity>
#>   N producers:        30
#>   mean LOO r:         1.000
#>   SD:                 0.000
#>   median LOO r:       1.000
#>   MAD:                0.000
#>   flag rule:          r < mean - 2.50 * SD   (=> r < 1.000)
#>   flagged producers:  none
#> ---
#> * $results$icc
#> <rcicrely ICC>
#>   model:        ICC(3,*) / two-way mixed; pixels fixed, participants random
#>   N targets:    1024 pixels
#>   N raters:     30 participants
#> 
#>   ICC(3,1):      0.8788
#>   ICC(3,k):      0.9954
```

`$results$split_half`, `$results$loo`, `$results$icc` are the three
result objects. `plot(within_A)` lays out all three on one page.

## 6. Function reference: between-condition inference

### 6.1 `pixel_t_test()`

**What it does.** Vectorised Welch’s t per pixel between two condition
signal matrices.

**When to use it.** You want the raw t-map (e.g. for custom plotting,
threshold sensitivity analysis) without cluster-level inference.

``` r
t_vec <- pixel_t_test(sig_A, sig_B)
range(t_vec); length(t_vec)
#> [1] -25.20620  27.27447
#> [1] 1024
```

Welch’s t (unequal variances) is correct - 2IFC and Brief-RC conditions
can differ in N and variance.

### 6.2 `rel_cluster_test()`

**What it does.** Cluster-based permutation test with max-statistic FWER
control. Pixels above the threshold form clusters via 4-connectivity;
cluster mass = sum of t-values inside; null built by **stratified**
label permutation; max-mass approach controls family-wise error.

**When to use it.** “Where in the image do these two conditions differ?”
with proper multiple-comparison control across pixels.

``` r
clust <- rel_cluster_test(
  sig_A, sig_B,
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
#>   clusters found:    2 (2 significant)
#> 
#>  cluster_id direction  mass size p_value significant
#>           1       neg -9886  512       0        TRUE
#>           1       pos  9807  512       0        TRUE
```

``` r
plot(clust)
```

![Per-pixel Welch t-map; significant clusters
outlined.](tutorial_files/figure-html/plot-cluster-1.png)

Per-pixel Welch t-map; significant clusters outlined.

**Reading the result.** `$observed_t` is the raw t-map. `$clusters` is a
data.frame with `cluster_id`, `direction`, `mass`, `size`, `p_value`,
`significant`. `$pos_labels` / `$neg_labels` are integer matrices for
plotting overlays. `$null_distribution` exposes the per-permutation max
masses.

**Common mistakes.** Reading `cluster_threshold = 2.0` as “p \< 0.05” -
it is the t-cutoff for forming candidate clusters, not the cluster
significance threshold. Significance is the cluster’s mass relative to
the null. Trusting cluster p-values with `n_permutations < 1000`.

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
#>   n_boot:     500
#>   CI level:   95%
#>   Pearson r   = -0.996   [-0.994, -0.991]   SE = 0.001
#>   Euclidean   = 9.668   [9.293, 10.089]   SE = 0.207
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
#>   n_boot:     500
#>   CI level:   95%
#>   Pearson r   = -0.996   [-0.994, -0.991]   SE = 0.001
#>   Euclidean   = 9.668   [9.293, 10.089]   SE = 0.207
```

## 7. Brief-RC end-to-end

Schmitz, Rougier & Yzerbyt (2024) introduced Brief-RC 12: each trial
shows 12 noisy faces (6 oriented, 6 inverted, drawn from a shared noise
pool); participant picks one. The recorded format is one row per trial,
columns `participant_id`, `trial`, `stimulus` (chosen pool id),
`response` (`+1` if oriented, `-1` if inverted). `rcicrely`’s
implementation is native (no rcicr `_brief` functions, which don’t exist
in canonical rcicr v1.0.1).

``` r
# 1. Generate the noise pool (once per project) - heavy step.
rcicr::generateStimuli2IFC(
  base_face_files     = list(base = "data/base.jpg"),
  n_trials            = 600L,        # becomes the Brief-RC pool size
  img_size            = 256L,
  stimulus_path       = "data/stim",
  label               = "mystudy",
  ncores              = 1L,          # macOS-safe default
  return_as_dataframe = TRUE,        # capture the noise matrix!
  seed                = 1L
) |> saveRDS("data/noise_matrix.rds")  # cache for fast loading

# 2. Collect Brief-RC responses in the per-trial format.

# 3. Compute individual CIs.
res <- ci_from_responses_briefrc(
  responses       = my_brief_rc_responses,
  noise_matrix    = readRDS("data/noise_matrix.rds"),
  base_image_path = "data/base.jpg",
  method          = "briefrc12",
  scaling         = "matched"   # also rendered for plotting
)

# 4. Reliability assessment.
within  <- run_within(res$signal_matrix, seed = 1L)
between <- run_between(res_A$signal_matrix, res_B$signal_matrix,
                       seed = 1L)

# 5. Save rendered CI to PNG (visualization only - DO NOT use for stats)
if (requireNamespace("png", quietly = TRUE)) {
  for (j in seq_along(res$participants)) {
    img <- matrix(res$rendered_ci[, j], res$img_dims[1], res$img_dims[2])
    img <- pmin(pmax(img, 0), 1)
    png::writePNG(img, sprintf("plots/p%02d.png", j))
  }
}
```

**Brief-RC `infoVal`** is deferred (no canonical implementation in
either rcicr or rcicrely). It needs a reference distribution matched to
each producer’s trial count, not the pool size. The `rel_*` reliability
metrics do not depend on it. If you hand-roll an `infoVal` for Brief-RC
today, compute it on the **raw mask** (`res$signal_matrix`), not on
`res$rendered_ci` or any saved PNG (those carry the scaling step and
will give a wrong reference comparison).

## 8. Interpreting results and common pitfalls

**Sample size.** `rcicrely` aborts below 4 producers per condition and
warns below 30 (see Cone et al. 2021; who recommend N \>= 60 for stable
reliability assessment). At N = 10 the metrics are computable but their
CIs are wide.

**Split-half.** Think of `r_hh` ~ 0.3 as weak, 0.5 as moderate, 0.7+ as
strong producer-level consistency. `r_SB` is what to report when asked
about the full-sample CI.

**LOO flagging.** A flag is a starting point for investigation, not a
verdict. Use `flag_method = "mad"` if the SD rule is dominated by the
very outliers it tries to find.

**ICC(3,*) vs ICC(2,*).** Numerically close at high pixel counts but
ICC(3,\*) is the correctly-specified model. Don’t swap because of a
reviewer; report both with `variants = c("3_1", "3_k", "2_1", "2_k")` if
pressed.

**Cluster threshold != p-value.** `cluster_threshold = 2.0` is the
t-cutoff for forming clusters. Significance is decided by mass relative
to the null. The default is roughly two-tailed p \< 0.05 at large df,
but at small N it under-rejects - inspect the t histogram.

**Raw vs. rendered CIs.** This package’s reliability math is correct
**given** a raw signal matrix. Mode 1 (PNGs on disk) loads
`scaling(mask)` instead of `mask`; results will differ from what raw
responses would produce. The size of the difference depends on the
scaling that was applied at PNG-write time. Mode 2 is the safe path.

**`infoVal` correctness.**
[`rcicr::computeInfoVal2IFC()`](https://rdrr.io/pkg/rcicr/man/computeInfoVal2IFC.html)
extracts the raw `$ci` from its input CI-list internally, so the
canonical 2IFC path is safe regardless of `scaling`. The pitfall is in
hand-rolled `infoVal` implementations, for Brief-RC, or for any custom
code that consumes a bare matrix, which must be fed the raw mask.
Feeding rendered/scaled values silently distorts the reference
comparison and produces wrong significance verdicts.

## 9. References

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
- Schmitz, M., Rougier, M., & Yzerbyt, V. (2024). Introducing the brief
  reverse correlation: An improved tool to assess visual
  representations. *European Journal of Social Psychology*.
  <https://doi.org/10.1002/ejsp.3100>
- Shrout, P. E., & Fleiss, J. L. (1979). Intraclass correlations: Uses
  in assessing rater reliability. *Psychological Bulletin*, 86(2),
  420-428. <https://doi.org/10.1037/0033-2909.86.2.420>
