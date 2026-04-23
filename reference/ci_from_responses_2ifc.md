# Compute individual 2IFC CIs from trial-level responses

Thin wrapper around
[`rcicr::batchGenerateCI2IFC()`](https://rdrr.io/pkg/rcicr/man/batchGenerateCI2IFC.html)
that handles the known rcicr-integration gotchas and exposes a uniform
return shape across the 2IFC and Brief-RC paths.

Use this when you have 2IFC trial-level responses and want the package
to compute individual CIs ready for the reliability metrics.

## Usage

``` r
ci_from_responses_2ifc(
  responses,
  rdata_path,
  baseimage = NULL,
  participant_col = "participant_id",
  stimulus_col = "stimulus",
  response_col = "response",
  scaling = "autoscale",
  keep_rendered = FALSE,
  targetpath = tempfile("rcicrely_2ifc_"),
  save_as_png = FALSE
)
```

## Arguments

- responses:

  A data frame with one row per trial. Must contain `participant_col`
  (producer id), the column named by `stimulus_col` (stimulus id,
  integer – index into the rcicr noise pool), and `response_col` with
  values in `{-1, +1}`.

- rdata_path:

  Path to the `.Rdata` file produced by
  [`rcicr::generateStimuli2IFC()`](https://rdrr.io/pkg/rcicr/man/generateStimuli2IFC.html).

- baseimage:

  Base image label used at stimulus-generation time
  (`baseimage = "base"` is the convention in the test harness; whatever
  you passed to
  `generateStimuli2IFC(base_face_files = list(<label> = path))`). If
  `NULL`, tries to read the single base stored in `rdata_path`.

- participant_col, stimulus_col, response_col:

  Column names in `responses`.

- scaling:

  rcicr scaling option; one of `"autoscale"`, `"independent"`,
  `"constant"`, `"none"`. Passed through to
  [`rcicr::batchGenerateCI2IFC()`](https://rdrr.io/pkg/rcicr/man/batchGenerateCI2IFC.html).
  Scaling only affects the `$combined` (rendered) image, not `$ci` – so
  the returned `$signal_matrix` is the raw mask regardless of this
  argument.

- keep_rendered:

  If `TRUE`, also extract rcicr's `$combined` image (base + scaled
  noise) per producer and stack into `$rendered_ci`. Default `FALSE`.
  Use only for visualisation; do not feed `$rendered_ci` to `rel_*`
  functions.

- targetpath:

  Where rcicr writes PNGs (defaults to an auto-cleaned tempdir so the
  working directory isn't polluted).

- save_as_png:

  Whether rcicr writes individual CI PNGs. Defaults to `FALSE` for speed
  in a pure reliability pipeline.

## Value

A list with

- `signal_matrix` – pixels x participants raw mask (always raw,
  regardless of `scaling`).

- `rendered_ci` – pixels x participants, present only when
  `keep_rendered = TRUE`. Visualisation only.

- `participants` – character vector of participant ids in the column
  order of `signal_matrix`.

- `img_dims` – integer `c(nrow, ncol)`.

- `scaling` – the `scaling` value passed to rcicr.

- `rcicr_result` – the raw rcicr list (one element per producer).

## Details

What the wrapper does for you:

- Attaches `foreach`, `tibble`, `dplyr` at runtime (`%dopar%` /
  `tribble` / `%>%` are not namespace-prefixed inside rcicr).

- Matches the `.Rdata` extension case-insensitively (rcicr writes
  lowercase on some filesystems).

- Extracts the per-participant CI **noise component** (`$ci`) from the
  rcicr result list and stacks it into a pixels x participants signal
  matrix – already base-subtracted, ready for `rel_*()`.

- When `keep_rendered = TRUE`, also extracts the rendered `$combined`
  image and stacks it as `$rendered_ci` for visualisation (display only,
  not for downstream stats).

`rcicr` must be installed (it's a Suggests dep; install with
`remotes::install_github("rdotsch/rcicr")` if needed).

## Reading the result

- `$signal_matrix` is the **raw mask** (rcicr's `$ci` per producer).
  This is what every `rel_*` function expects.

- `$rendered_ci`, present when `keep_rendered = TRUE`, is the
  `base + scaling(mask)` image rcicr would have written to PNG.
  **Visualisation only.**

- `$rcicr_result` is the raw return value of
  [`rcicr::batchGenerateCI2IFC()`](https://rdrr.io/pkg/rcicr/man/batchGenerateCI2IFC.html)
  in case you need anything else (e.g. the per-producer base,
  scaled-only ci).

## Common mistakes

- Passing `$rendered_ci` to `rel_*` functions – it carries the scaling
  step, which distorts variance-based metrics.

- Forgetting that the warning about Mode 1 also applies if you later
  read `$rendered_ci` back from PNG instead of using `$signal_matrix`
  directly.

## See also

[`ci_from_responses_briefrc()`](https://olivethree.github.io/rcicrely/reference/ci_from_responses_briefrc.md),
[`rcicr::batchGenerateCI2IFC()`](https://rdrr.io/pkg/rcicr/man/batchGenerateCI2IFC.html)

## Examples

``` r
if (FALSE) { # \dontrun{
res <- ci_from_responses_2ifc(
  responses  = my_responses,
  rdata_path = "data/rcicr_stimuli.Rdata",
  baseimage  = "base"
)
dim(res$signal_matrix)  # n_pixels x n_participants
} # }
```
