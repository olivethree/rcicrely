# Compute individual Brief-RC CIs from trial-level responses

Native implementation of Schmitz, Rougier & Yzerbyt (2024)'s Brief-RC
mask. **Does not** call any `rcicr::*_brief` function – those do not
exist in canonical rcicr v1.0.1. Only rcicr's noise-pattern pool (the
`stimuli` object inside an `.Rdata` from `generateStimuli2IFC()`) is
reused; the mask is computed in pure R.

Use this when you have Brief-RC 12 trial-level responses and want the
package to produce per-producer noise masks ready for the reliability
metrics.

## Usage

``` r
ci_from_responses_briefrc(
  responses,
  rdata_path = NULL,
  noise_matrix = NULL,
  base_image_path,
  participant_col = "participant_id",
  stimulus_col = "stimulus",
  response_col = "response",
  method = c("briefrc12", "briefrc20"),
  scaling = c("none", "matched", "constant"),
  scaling_constant = NULL
)
```

## Arguments

- responses:

  Data frame with one row per trial. Must contain the columns named by
  `participant_col`, `stimulus_col`, `response_col`. `response` values
  must be in `{-1, +1}`. Brief-RC 12 trial structure: `stimulus` =
  chosen pool id; `response = +1` if the oriented version was chosen,
  `-1` if the inverted version was chosen. Unchosen faces are not
  recorded.

- rdata_path, noise_matrix:

  Exactly one must be supplied. Provide `rdata_path` to read the noise
  matrix from an rcicr `.Rdata` (`stimuli` object), or pass a pre-loaded
  `noise_matrix` directly.

- base_image_path:

  Path to the base face image (PNG / JPEG). Used to validate the noise
  matrix shape and (when `scaling` is not `"none"`) to render the
  visualisation-only `$rendered_ci` field.

- participant_col, stimulus_col, response_col:

  Column names.

- method:

  One of `"briefrc12"`. `"briefrc20"` is reserved for a future release
  and currently aborts with a clear message.

- scaling:

  Visualisation-only scaling for the optional `$rendered_ci` field. One
  of `"none"` (default; no rendered_ci returned), `"matched"` (stretch
  mask to base range, then add to base – Schmitz's default), or
  `"constant"` (multiply mask by `scaling_constant`, then add to base).
  The mathematical `$signal_matrix` returned by this function is always
  the raw unscaled mask, regardless of this argument.

- scaling_constant:

  Numeric multiplier used when `scaling = "constant"`. Ignored
  otherwise.

## Value

A list with

- `signal_matrix` – pixels x participants raw mask (always raw, even
  when `scaling != "none"`).

- `rendered_ci` – pixels x participants, present only when
  `scaling != "none"`. Visualisation only.

- `participants` – character vector of participant ids.

- `img_dims` – integer `c(nrow, ncol)`.

- `scaling` – the scaling option that was used.

## Details

Formula (Schmitz's `genMask()` exactly):

    X    <- data.table(response, stim)
    X    <- X[, .(response = mean(response)), stim]  # collapse duplicates
    mask <- (noiseMatrix[, X$stim] %*% X$response) / length(X$response)

The `length(X$response)` denominator is the number of **unique** pool
ids chosen by that participant, not the raw trial count. If a
participant chooses the same stimulus on two trials with opposite
responses, those two cancel.

## Reading the result

- `$signal_matrix` is the **raw mask** per producer – the object every
  `rel_*` function expects. Always pass this (and only this) to
  reliability metrics or to any external infoVal computation.

- `$rendered_ci`, when present, is `base + scaling(mask)` per producer.
  **Visualization only.** Save to PNG, plot for inspection, but do not
  feed into reliability metrics – the scaling step distorts
  variance-based statistics. Also do not feed it to a hand-rolled
  Brief-RC `infoVal` (no canonical implementation exists in either rcicr
  or rcicrely): the reference distribution would compare against a
  different scale.

- `$participants` and `$img_dims` are convenience metadata.

## Common mistakes

- Passing the "expanded" 12-rows-per-trial response format. Brief-RC 12
  is one row per trial; `stimulus` is the chosen pool id and `response`
  is `+1` (oriented) or `-1` (inverted). Schmitz 2024 sec.3.1.2.

- Using `$rendered_ci` for downstream stats. It exists only because the
  user often wants to save a PNG.

- Asking for `method = "briefrc20"` – that's reserved for a future
  release and aborts.

## References

Schmitz, M., Rougier, M., & Yzerbyt, V. (2024). Introducing the brief
reverse correlation: an improved tool to assess visual representations.
*Behavior Research Methods*.

## See also

[`ci_from_responses_2ifc()`](https://olivethree.github.io/rcicrely/reference/ci_from_responses_2ifc.md),
[`run_within()`](https://olivethree.github.io/rcicrely/reference/run_within.md),
[`run_between()`](https://olivethree.github.io/rcicrely/reference/run_between.md)

## Examples

``` r
if (FALSE) { # \dontrun{
res <- ci_from_responses_briefrc(
  responses       = my_responses,
  rdata_path      = "data/rcicr_stimuli.Rdata",
  base_image_path = "data/base.jpg"
)
# signal_matrix goes into rel_*; rendered_ci is for plotting only.
rel_split_half(res$signal_matrix, seed = 1L)
} # }
```
