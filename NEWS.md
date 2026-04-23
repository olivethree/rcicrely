# rcicrely 0.1.0 (unreleased)

Initial public release.

## Features

- **Within-condition reliability**
  - `rel_split_half()` — permuted split-half with Spearman–Brown
    correction; handles odd N via per-permutation drop.
  - `rel_loo()` — leave-one-out sensitivity with 2.5-SD outlier
    flagging by default.
  - `rel_icc()` — ICC(3,1) and ICC(3,k) via direct mean-squares
    (two-way mixed model: pixels fixed, participants random);
    ICC(2,*) available on request.
- **Between-condition inference**
  - `pixel_t_test()` — vectorised Welch's t per pixel.
  - `rel_cluster_test()` — cluster-based permutation test with
    4-connectivity labelling and max-statistic FWER control;
    permutation is stratified (preserves N_A, N_B).
  - `rel_dissimilarity()` — bootstrap correlation + Euclidean
    distance with percentile CIs (no `boot` dependency).
- **Orchestrators**
  - `run_within()`, `run_between()` — run every metric, wrap in an
    `rcicrely_report`.
- **Data loading**
  - `read_cis()`, `extract_signal()`, `load_signal_matrix()` for
    pre-computed CI images.
  - `read_noise_matrix()` for rcicr `.Rdata` or plain-text noise
    matrices.
  - `ci_from_responses_2ifc()` — wraps `rcicr::batchGenerateCI2IFC()`
    with all known gotchas handled.
  - `ci_from_responses_briefrc()` — native Brief-RC 12
    implementation per Schmitz, Rougier & Yzerbyt (2024), matching
    `genMask()` exactly including duplicate-stimulus tie-breaking.
- **S3 methods** — `print`, `summary`, `plot` for every
  `rcicrely_*` result class plus `rcicrely_report`.

## Known limitations

- Brief-RC infoVal is deferred (needs a per-participant-trial-count
  reference distribution; out of scope for v0.1).
- Brief-RC 20 is unsupported; `ci_from_responses_briefrc(method =
  "briefrc20")` aborts with a clear message.
- BCa bootstrap CIs are future work; only percentile CIs in v0.1.
