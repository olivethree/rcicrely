# Read a directory of CI images into a pixels x participants matrix

Reads every PNG / JPEG in `dir` (non-recursive), converts each to
grayscale via ITU-R BT.709 luminance weights
(`0.2126 R + 0.7152 G + 0.0722 B`), and returns a numeric matrix with
one column per file, sorted alphabetically by filename.

Use this when you have one CI image per producer already saved on disk
and you want a numeric matrix you can pass to
[`extract_signal()`](https://olivethree.github.io/rcicrely/reference/extract_signal.md)
or
[`load_signal_matrix()`](https://olivethree.github.io/rcicrely/reference/load_signal_matrix.md)
for the base subtraction step.

## Usage

``` r
read_cis(dir, pattern = "\\.(png|jpe?g)$", acknowledge_scaling = FALSE)
```

## Arguments

- dir:

  Directory containing the CI images.

- pattern:

  Regular expression used to select files. Defaults to any PNG / JPG /
  JPEG (case-insensitive).

- acknowledge_scaling:

  Set to `TRUE` to silence the once-per-session warning about Mode 1
  returning rendered (scaled) pixel data. Recommended only when you
  generated the PNGs yourself with
  `rcicr::generateCI2IFC(scaling = "none")` (or equivalent) and the PNG
  truly carries unscaled noise.

## Value

A numeric matrix with `nrow = n_pixels` and `ncol = n_participants`.
Carries an `img_dims` attribute.

## Details

All images must share the same dimensions; mixed sizes abort. Column
names are the filename without the extension.

## Reading the result

The returned matrix has `nrow = n_pixels` and `ncol = n_participants`.
Row-major pixel ordering follows `as.vector(img)` of the underlying
image matrix. The integer `c(nrow, ncol)` of the original image is
stamped on the result as an `img_dims` attribute, which downstream
functions read automatically.

## Common mistakes

- Forgetting to subtract the base. The matrix returned here contains the
  **rendered CI** (base + scaled noise). Pass it through
  [`extract_signal()`](https://olivethree.github.io/rcicrely/reference/extract_signal.md)
  before any reliability metric.

- Mixing CIs from different image sizes in one directory – this aborts
  with a clear error rather than silently misaligning.

- Trusting reliability numbers without first reading the **Note** below
  on raw vs. rendered CIs.

## Raw vs. rendered CIs

PNG pixel values are necessarily what was rendered to disk:
`base + scaling(mask)`. Subtracting the base in
[`extract_signal()`](https://olivethree.github.io/rcicrely/reference/extract_signal.md)
therefore yields `scaling(mask)`, not the raw mask. Variance-based
reliability metrics
([`rel_icc()`](https://olivethree.github.io/rcicrely/reference/rel_icc.md),
Euclidean half of
[`rel_dissimilarity()`](https://olivethree.github.io/rcicrely/reference/rel_dissimilarity.md),
[`pixel_t_test()`](https://olivethree.github.io/rcicrely/reference/pixel_t_test.md),
[`rel_cluster_test()`](https://olivethree.github.io/rcicrely/reference/rel_cluster_test.md))
are sensitive to the scaling transform; only correlation-based metrics
([`rel_split_half()`](https://olivethree.github.io/rcicrely/reference/rel_split_half.md),
[`rel_loo()`](https://olivethree.github.io/rcicrely/reference/rel_loo.md),
correlation half of
[`rel_dissimilarity()`](https://olivethree.github.io/rcicrely/reference/rel_dissimilarity.md))
survive a single uniform linear scaling unmodified – and even those
break under per-CI `"matched"`-style scaling. The canonical 2IFC
`infoVal` path
([`rcicr::computeInfoVal2IFC()`](https://rdrr.io/pkg/rcicr/man/computeInfoVal2IFC.html))
extracts the raw `$ci` internally from the rcicr CI-list and is **not**
affected; hand-rolled `infoVal` implementations (Brief-RC, custom code)
**are**, and should be fed the raw mask. Prefer
[`ci_from_responses_2ifc()`](https://olivethree.github.io/rcicrely/reference/ci_from_responses_2ifc.md)
/
[`ci_from_responses_briefrc()`](https://olivethree.github.io/rcicrely/reference/ci_from_responses_briefrc.md)
when raw responses are available; both return the raw mask. See
[`vignette("tutorial", package = "rcicrely")`](https://olivethree.github.io/rcicrely/articles/tutorial.md)
chapter 3 for the full discussion.

## See also

[`extract_signal()`](https://olivethree.github.io/rcicrely/reference/extract_signal.md),
[`load_signal_matrix()`](https://olivethree.github.io/rcicrely/reference/load_signal_matrix.md),
[`ci_from_responses_2ifc()`](https://olivethree.github.io/rcicrely/reference/ci_from_responses_2ifc.md),
[`ci_from_responses_briefrc()`](https://olivethree.github.io/rcicrely/reference/ci_from_responses_briefrc.md)

## Examples

``` r
if (FALSE) { # \dontrun{
cis <- read_cis("data/cis_condition_A/")
signal <- extract_signal(cis, base_image_path = "data/base.jpg")
} # }
```
