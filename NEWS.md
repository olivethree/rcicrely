# rcicrely 0.1.1 (unreleased)

## Raw vs. rendered CIs - documentation overhaul + safety net

Reliability metrics operate on the **raw mask** (the participant's
unscaled noise contribution). Reading PNGs from disk via `read_cis()`
silently loaded the *rendered* CI (`base + scaling(mask)`) and fed
`scaling(mask)` to downstream metrics, which can distort
variance-based statistics (`rel_icc()`, Euclidean dissimilarity,
`pixel_t_test()`, `rel_cluster_test()`). The canonical 2IFC `infoVal`
path (`rcicr::computeInfoVal2IFC()`) is unaffected because it
extracts the raw `$ci` element from the rcicr CI-list internally;
hand-rolled `infoVal` implementations (Brief-RC, custom code) are
affected and need the raw mask. The mathematical core of the
package was always correct; this release sharpens the boundary so
users see the issue.

- `read_cis()`, `extract_signal()`, `load_signal_matrix()` now emit a
  once-per-session `cli` warning explaining the raw-vs-rendered
  distinction, which metrics are affected, and how to silence
  (`acknowledge_scaling = TRUE` argument or
  `options(rcicrely.silence_scaling_warning = TRUE)`).
- Every `rel_*()` and `run_*()` function now runs a cheap heuristic
  on entry (`looks_scaled()`) and emits a once-per-session warning
  if the input matrix's dynamic range looks rendered rather than raw.
  Heuristic only; same silencing options.
- `vignette("tutorial")` is restructured into a function-by-function
  user guide (~600 lines, was ~350). Adds a dedicated "Mental model:
  the signal matrix" chapter with a numerical demonstration of how
  scaling distorts variance-based metrics while leaving Pearson
  correlation intact, and an "Interpreting results" chapter that
  flags the `infoVal` pitfall.
- `README.md` adds a top-level "Important: raw vs. rendered CIs"
  callout.

## Brief-RC: visualization-only `$rendered_ci`

`ci_from_responses_briefrc()` gains a `scaling` argument
(`"none"` (default) / `"matched"` / `"constant"`) and a
`scaling_constant` argument. The returned `$signal_matrix` is
**always** the raw mask, regardless of `scaling`. When
`scaling != "none"`, an additional `$rendered_ci` field carries
`base + scaling(mask)` per producer for visualisation - **not** for
downstream stats. Mirrors `ci_from_responses_2ifc()`'s new
`keep_rendered = TRUE` argument, which adds a `$rendered_ci` field
extracted from rcicr's `$combined`.

## `rel_loo()` MAD-based outlier rule

- New `flag_method = c("sd", "mad")` argument (default `"sd"`,
  preserves previous behavior). `"mad"` flags producers whose LOO
  correlation falls below `median(r) - flag_threshold * mad(r)`,
  which is robust to the few atypical producers that often dominate
  the SD rule.
- `flag_threshold_sd` argument is renamed to `flag_threshold` (the
  multiplier applies to either SD or MAD). The old name still works
  as an alias so existing code does not break.
- Result object gains `$median_r`, `$mad_r`, `$flag_method`,
  `$flag_threshold` fields.

## Other

- New internal helpers: `looks_scaled()` heuristic
  (`R/scaling_diagnostics.R`); `warn_mode1_scaling()`,
  `warn_looks_scaled()`, `reset_session_warnings()` in `R/utils.R`.
- Roxygen overhaul: every exported function now has uniform
  `@description`, `@section Reading the result:`,
  `@section Common mistakes:` blocks for consistent help-page UX.

# rcicrely 0.1.0

Initial public release.

## Features

- **Within-condition reliability**
  - `rel_split_half()` - permuted split-half with Spearman-Brown
    correction; handles odd N via per-permutation drop.
  - `rel_loo()` - leave-one-out sensitivity with 2.5-SD outlier
    flagging by default.
  - `rel_icc()` - ICC(3,1) and ICC(3,k) via direct mean-squares
    (two-way mixed model: pixels fixed, participants random);
    ICC(2,*) available on request.
- **Between-condition inference**
  - `pixel_t_test()` - vectorised Welch's t per pixel.
  - `rel_cluster_test()` - cluster-based permutation test with
    4-connectivity labelling and max-statistic FWER control;
    permutation is stratified (preserves N_A, N_B).
  - `rel_dissimilarity()` - bootstrap correlation + Euclidean
    distance with percentile CIs (no `boot` dependency).
- **Orchestrators**
  - `run_within()`, `run_between()` - run every metric, wrap in an
    `rcicrely_report`.
- **Data loading**
  - `read_cis()`, `extract_signal()`, `load_signal_matrix()` for
    pre-computed CI images.
  - `read_noise_matrix()` for rcicr `.Rdata` or plain-text noise
    matrices.
  - `ci_from_responses_2ifc()` - wraps `rcicr::batchGenerateCI2IFC()`
    with all known gotchas handled.
  - `ci_from_responses_briefrc()` - native Brief-RC 12
    implementation per Schmitz, Rougier & Yzerbyt (2024), matching
    `genMask()` exactly including duplicate-stimulus tie-breaking.
- **S3 methods** - `print`, `summary`, `plot` for every
  `rcicrely_*` result class plus `rcicrely_report`.

## Known limitations

- Brief-RC infoVal is deferred (needs a per-participant-trial-count
  reference distribution; out of scope for v0.1).
- Brief-RC 20 is unsupported; `ci_from_responses_briefrc(method =
  "briefrc20")` aborts with a clear message.
- BCa bootstrap CIs are future work; only percentile CIs in v0.1.
