# Load an image-based face-region mask from a PNG or JPEG file

Reads a binary mask image from disk and returns a logical vector of
length `prod(img_dims)` in column-major order — the format the rest of
the package expects. White / light pixels (above `threshold`) become
`TRUE`, dark pixels become `FALSE`. Use this for image-based masks
created with
[webmorphR::mask_oval()](https://debruine.github.io/webmorphR/), painted
in GIMP, drawn with PowerPoint shapes, or any other tool that produces a
binary PNG/JPEG.

Companion to
[`face_mask()`](https://olivethree.github.io/rcicrely/reference/face_mask.md)
(which generates parametric masks programmatically). Either function
returns the same logical-vector shape and either is accepted by
[`infoval()`](https://olivethree.github.io/rcicrely/reference/infoval.md)'s
`mask` argument or by row-subsetting on a signal matrix.

## Usage

``` r
load_face_mask(path, threshold = 0.5, expected_dims = NULL)
```

## Arguments

- path:

  Path to a PNG or JPEG mask image. Required to be the same dimensions
  as the analysis image.

- threshold:

  Numeric in `[0, 1]`. Pixels with luminance above this threshold are
  `TRUE`. Default `0.5` (mid-grey).

- expected_dims:

  Optional integer `c(nrow, ncol)`. When set, the function aborts if the
  loaded image's dimensions do not match — useful for catching a
  wrong-resolution mask before it silently corrupts a downstream
  computation.

## Value

Logical vector of length `prod(img_dims)`, column-major.

## See also

[`face_mask()`](https://olivethree.github.io/rcicrely/reference/face_mask.md),
[`infoval()`](https://olivethree.github.io/rcicrely/reference/infoval.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Hand-painted mask: white = include, black = exclude
fm <- load_face_mask("masks/oval_256.png", expected_dims = c(256L, 256L))
iv <- infoval(signal_matrix, noise_matrix, trial_counts,
              mask = fm, iter = 1000L, seed = 1L)
} # }
```
