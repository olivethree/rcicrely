# Load the noise matrix from an rcicr `.Rdata`, an `.rds`, or a text file

Produces the pixels x n_trials numeric matrix of individual noise
patterns that
[`ci_from_responses_briefrc()`](https://olivethree.github.io/rcicrely/reference/ci_from_responses_briefrc.md)
needs. Use this when you want to load (or reconstruct) the noise pool
once and pass it to multiple Brief-RC CI computations without
re-loading.

Three input formats are supported:

1.  **rcicr `.Rdata`** (the output of
    [`rcicr::generateStimuli2IFC()`](https://rdrr.io/pkg/rcicr/man/generateStimuli2IFC.html)).
    The rdata does not store the noise matrix directly; it stores
    `stimuli_params` (per-trial parameter rows) and `p` (the sinusoid /
    gabor basis). This loader reconstructs each trial's noise pattern
    via
    [`rcicr::generateNoiseImage()`](https://rdrr.io/pkg/rcicr/man/generateNoiseImage.html)
    and stacks them. Requires `rcicr` to be installed.

2.  **R `.rds`** previously saved from
    `rcicr::generateStimuli2IFC(return_as_dataframe = TRUE)`. No
    reconstruction needed.

3.  **Plain text** (whitespace- or comma-delimited) matrix, one column
    per trial.

The loader dispatches on file magic bytes (not extension; rcicr writes
lowercase `.Rdata` on some filesystems, which trips case-sensitive
matching).

## Usage

``` r
read_noise_matrix(path, baseimage = NULL, stimuli_object = "stimuli")
```

## Arguments

- path:

  Path to the file.

- baseimage:

  For rcicr `.Rdata` inputs: which base label to reconstruct noise for.
  Defaults to the first key in `base_faces` if there is exactly one.

- stimuli_object:

  For `.Rdata` files that contain a pre-saved `stimuli` object (rare,
  rcicr does not save one by default), the name to look up. Defaults to
  `"stimuli"`.

## Value

A numeric pixels x n_trials matrix.

## Reading the result

Numeric pixels x n_trials matrix. Each column is one trial's noise
pattern (a sinusoidal or gabor combination of the basis `p`). Pool ids
referenced by Brief-RC `responses$stimulus` are column indices into this
matrix.

## Common mistakes

- Pointing at the rcicr `.Rdata` and assuming it contains a `stimuli`
  object, it doesn't. The loader reconstructs by calling
  [`rcicr::generateNoiseImage()`](https://rdrr.io/pkg/rcicr/man/generateNoiseImage.html)
  per row, which is correct but slow (~4s per 120 trials at 64x64).

- Caching the captured return of
  `rcicr::generateStimuli2IFC(return_as_dataframe = TRUE)` to `.rds`
  once and reusing it is much faster than reconstructing from the rdata
  each run.

- Multiple stale `.rdata` accumulating in a directory because rcicr
  timestamps its filenames; pick the most recent by mtime if you have
  several.

## See also

[`ci_from_responses_briefrc()`](https://olivethree.github.io/rcicrely/reference/ci_from_responses_briefrc.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# From an rcicr rdata:
nm <- read_noise_matrix("data/rcicr_stimuli.Rdata")

# From a cached rds you saved earlier:
nm <- read_noise_matrix("data/noise_matrix.rds")
} # }
```
