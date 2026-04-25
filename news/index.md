# Changelog

## rcicrely 0.2.2

### Bug fix: `cli` progress bar lookup in `rel_*()` permutation loops

[`rel_split_half()`](https://olivethree.github.io/rcicrely/reference/rel_split_half.md)
(and any other `rel_*()` that ticks a progress bar from inside a
`with_seed({ ... })` block) could abort mid-loop with
`Error in cli::cli_progress_update(id = id): Cannot find progress bar 'cli-XXXXX-NN'`.
The bar was created in the `rel_*()` frame but `progress_tick()`
resolved to a different parent frame (the captured `with_seed`
expression), so `cli`’s frame-keyed bar lookup missed.

Progress bars are now pinned to a dedicated package-private env
(`.rcicrely_progress_env`, `parent = baseenv()`) shared by
`progress_start()`, `progress_tick()`, and `progress_done()`. The env
outlives any single function frame and is the same object for all three
helpers, so the lookup is stable under `with_seed()`, nested calls,
knitr chunks, and non-interactive `R -e`. The env’s
[`baseenv()`](https://rdrr.io/r/base/environment.html) parent lets
`cli`’s `glue`-based template formatting resolve `::` and other base R
names (an earlier fix that pinned bars to the existing
`.rcicrely_session_state` env failed because that env has
`parent = emptyenv()` for warning-state isolation).

### Improvement: `read_noise_matrix()` accepts gzipped text

`read_noise_matrix("noise_matrix.txt.gz")` now works. Previously, files
with gzip magic bytes (`0x1f 0x8b`) were treated as `.rds` or `.RData`
only; gzipped whitespace-delimited matrices aborted with a misleading
“could not parse as binary R file” error. The function now falls through
to
[`data.table::fread()`](https://rdrr.io/pkg/data.table/man/fread.html)
(which handles gzip transparently) when neither
[`readRDS()`](https://rdrr.io/r/base/readRDS.html) nor
[`load()`](https://rdrr.io/r/base/load.html) succeeds on a gzip-magic,
non-RDX file. Useful when the noise matrix is large at full double
precision (~700 MB for 65536 x 600); gzip typically halves the on-disk
size with no API change required.

### Documentation: `infoval()` framing for Brief-RC

Clarifies how to read the per-producer infoVal z-score and what
trial-count matching does (and does not) say about Brief-RC.

- [`infoval()`](https://olivethree.github.io/rcicrely/reference/infoval.md)
  help and
  [`vignette("tutorial")`](https://olivethree.github.io/rcicrely/articles/tutorial.md)
  (Brief-RC end-to-end section) now state explicitly that infoVal is a
  **magnitude** statistic (Frobenius norm of the mask): it answers “is
  the mask larger than chance?” but not “is it pointing at the right
  pattern?”. Cross-paradigm comparisons (Brief-RC vs 2IFC absolute
  z-scores) need care because Brief-RC’s per-trial cognitive richness
  may produce more accurately localized masks without changing
  magnitude.
- Reframes the trial-count-matching rationale: Brief-RC producers are
  **not** perceptually exposed to fewer pool items than 2IFC producers
  (each Brief-RC trial integrates 12 noisy faces of context). What
  [`infoval()`](https://olivethree.github.io/rcicrely/reference/infoval.md)
  matches the reference to is the number of **recorded** mask
  contributions per producer, which is typically smaller than the pool
  size in Brief-RC and equals the pool size in 2IFC. This is a
  calibration choice, not a claim about cognitive exposure.
- No code changes. The simulator at
  `R/infoval.R::simulate_reference_norms()` mirrors `genMask()` (random
  `(stim, +/-1)` pairs, mean-by-stim collapse, divide by number of
  unique chosen stims, Frobenius norm); under random responding the
  chosen face is oriented vs inverted with 50/50 probability, so the
  reference null is correctly specified for both paradigms.

### Wording: “canonical” -\> “standard” / “original” / “upstream” / “typical”

Mechanical project-wide rephrasing in shipped docs (R/, vignettes/,
README.md, NEWS.md, tests/). No semantic change.

## rcicrely 0.1.1 (unreleased)

### Raw vs. rendered CIs - documentation overhaul + safety net

Reliability metrics operate on the **raw mask** (the participant’s
unscaled noise contribution). Reading PNGs from disk via
[`read_cis()`](https://olivethree.github.io/rcicrely/reference/read_cis.md)
silently loaded the *rendered* CI (`base + scaling(mask)`) and fed
`scaling(mask)` to downstream metrics, which can distort variance-based
statistics
([`rel_icc()`](https://olivethree.github.io/rcicrely/reference/rel_icc.md),
Euclidean dissimilarity,
[`pixel_t_test()`](https://olivethree.github.io/rcicrely/reference/pixel_t_test.md),
[`rel_cluster_test()`](https://olivethree.github.io/rcicrely/reference/rel_cluster_test.md)).
The standard 2IFC `infoVal` path
([`rcicr::computeInfoVal2IFC()`](https://rdrr.io/pkg/rcicr/man/computeInfoVal2IFC.html))
is unaffected because it extracts the raw `$ci` element from the rcicr
CI-list internally; hand-rolled `infoVal` implementations (Brief-RC,
custom code) are affected and need the raw mask. The mathematical core
of the package was always correct; this release sharpens the boundary so
users see the issue.

- [`read_cis()`](https://olivethree.github.io/rcicrely/reference/read_cis.md),
  [`extract_signal()`](https://olivethree.github.io/rcicrely/reference/extract_signal.md),
  [`load_signal_matrix()`](https://olivethree.github.io/rcicrely/reference/load_signal_matrix.md)
  now emit a once-per-session `cli` warning explaining the
  raw-vs-rendered distinction, which metrics are affected, and how to
  silence (`acknowledge_scaling = TRUE` argument or
  `options(rcicrely.silence_scaling_warning = TRUE)`).
- Every `rel_*()` and `run_*()` function now runs a cheap heuristic on
  entry (`looks_scaled()`) and emits a once-per-session warning if the
  input matrix’s dynamic range looks rendered rather than raw. Heuristic
  only; same silencing options.
- [`vignette("tutorial")`](https://olivethree.github.io/rcicrely/articles/tutorial.md)
  is restructured into a function-by-function user guide (~600 lines,
  was ~350). Adds a dedicated “Mental model: the signal matrix” chapter
  with a numerical demonstration of how scaling distorts variance-based
  metrics while leaving Pearson correlation intact, and an “Interpreting
  results” chapter that flags the `infoVal` pitfall.
- `README.md` adds a top-level “Important: raw vs. rendered CIs”
  callout.

### Brief-RC: visualization-only `$rendered_ci`

[`ci_from_responses_briefrc()`](https://olivethree.github.io/rcicrely/reference/ci_from_responses_briefrc.md)
gains a `scaling` argument (`"none"` (default) / `"matched"` /
`"constant"`) and a `scaling_constant` argument. The returned
`$signal_matrix` is **always** the raw mask, regardless of `scaling`.
When `scaling != "none"`, an additional `$rendered_ci` field carries
`base + scaling(mask)` per producer for visualisation - **not** for
downstream stats. Mirrors
[`ci_from_responses_2ifc()`](https://olivethree.github.io/rcicrely/reference/ci_from_responses_2ifc.md)’s
new `keep_rendered = TRUE` argument, which adds a `$rendered_ci` field
extracted from rcicr’s `$combined`.

### `rel_loo()` MAD-based outlier rule

- New `flag_method = c("sd", "mad")` argument (default `"sd"`, preserves
  previous behavior). `"mad"` flags producers whose LOO correlation
  falls below `median(r) - flag_threshold * mad(r)`, which is robust to
  the few atypical producers that often dominate the SD rule.
- `flag_threshold_sd` argument is renamed to `flag_threshold` (the
  multiplier applies to either SD or MAD). The old name still works as
  an alias so existing code does not break.
- Result object gains `$median_r`, `$mad_r`, `$flag_method`,
  `$flag_threshold` fields.

### Other

- New internal helpers: `looks_scaled()` heuristic
  (`R/scaling_diagnostics.R`); `warn_mode1_scaling()`,
  `warn_looks_scaled()`, `reset_session_warnings()` in `R/utils.R`.
- Roxygen overhaul: every exported function now has uniform
  `@description`, `@section Reading the result:`,
  `@section Common mistakes:` blocks for consistent help-page UX.

## rcicrely 0.1.0

Initial public release.

### Features

- **Within-condition reliability**
  - [`rel_split_half()`](https://olivethree.github.io/rcicrely/reference/rel_split_half.md) -
    permuted split-half with Spearman-Brown correction; handles odd N
    via per-permutation drop.
  - [`rel_loo()`](https://olivethree.github.io/rcicrely/reference/rel_loo.md) -
    leave-one-out sensitivity with 2.5-SD outlier flagging by default.
  - [`rel_icc()`](https://olivethree.github.io/rcicrely/reference/rel_icc.md) -
    ICC(3,1) and ICC(3,k) via direct mean-squares (two-way mixed model:
    pixels fixed, participants random); ICC(2,\*) available on request.
- **Between-condition inference**
  - [`pixel_t_test()`](https://olivethree.github.io/rcicrely/reference/pixel_t_test.md) -
    vectorised Welch’s t per pixel.
  - [`rel_cluster_test()`](https://olivethree.github.io/rcicrely/reference/rel_cluster_test.md) -
    cluster-based permutation test with 4-connectivity labelling and
    max-statistic FWER control; permutation is stratified (preserves
    N_A, N_B).
  - [`rel_dissimilarity()`](https://olivethree.github.io/rcicrely/reference/rel_dissimilarity.md) -
    bootstrap correlation + Euclidean distance with percentile CIs (no
    `boot` dependency).
- **Orchestrators**
  - [`run_within()`](https://olivethree.github.io/rcicrely/reference/run_within.md),
    [`run_between()`](https://olivethree.github.io/rcicrely/reference/run_between.md) -
    run every metric, wrap in an `rcicrely_report`.
- **Data loading**
  - [`read_cis()`](https://olivethree.github.io/rcicrely/reference/read_cis.md),
    [`extract_signal()`](https://olivethree.github.io/rcicrely/reference/extract_signal.md),
    [`load_signal_matrix()`](https://olivethree.github.io/rcicrely/reference/load_signal_matrix.md)
    for pre-computed CI images.
  - [`read_noise_matrix()`](https://olivethree.github.io/rcicrely/reference/read_noise_matrix.md)
    for rcicr `.Rdata` or plain-text noise matrices.
  - [`ci_from_responses_2ifc()`](https://olivethree.github.io/rcicrely/reference/ci_from_responses_2ifc.md) -
    wraps
    [`rcicr::batchGenerateCI2IFC()`](https://rdrr.io/pkg/rcicr/man/batchGenerateCI2IFC.html)
    with all known gotchas handled.
  - [`ci_from_responses_briefrc()`](https://olivethree.github.io/rcicrely/reference/ci_from_responses_briefrc.md) -
    native Brief-RC 12 implementation per Schmitz, Rougier & Yzerbyt
    (2024), matching `genMask()` exactly including duplicate-stimulus
    tie-breaking.
- **S3 methods** - `print`, `summary`, `plot` for every `rcicrely_*`
  result class plus `rcicrely_report`.

### Known limitations

- Brief-RC infoVal is deferred (needs a per-participant-trial-count
  reference distribution; out of scope for v0.1).
- Brief-RC 20 is unsupported;
  `ci_from_responses_briefrc(method = "briefrc20")` aborts with a clear
  message.
- BCa bootstrap CIs are future work; only percentile CIs in v0.1.
