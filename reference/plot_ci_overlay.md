# Render a group CI as a translucent agreement-map overlay on a base image

Plots the producer-mean classification image (or the per-pixel *t*-map
from
[`agreement_map_test()`](https://olivethree.github.io/rcicrely/reference/agreement_map_test.md))
as a diverging-palette translucent layer composited over the base face
image. Useful as a publication-grade single figure that combines
anatomical context (the base) with the empirical signal (the CI /
agreement map). Optionally overlays significance contours from an
[`agreement_map_test()`](https://olivethree.github.io/rcicrely/reference/agreement_map_test.md)
result.

Pure base graphics — no `grid` / `ggplot2` dependency. Returns invisibly
the composed `nrow x ncol x 3` numeric raster so it can be saved to disk
via [`png::writePNG()`](https://rdrr.io/pkg/png/man/writePNG.html) /
[`jpeg::writeJPEG()`](https://rdrr.io/pkg/jpeg/man/writeJPEG.html) if
needed.

## Usage

``` r
plot_ci_overlay(
  signal_matrix,
  base_image,
  img_dims = NULL,
  test = NULL,
  mask = NULL,
  threshold = NULL,
  alpha_max = 0.7,
  palette = c("diverging"),
  contour_col = "black",
  contour_lwd = 1,
  main = NULL
)
```

## Arguments

- signal_matrix:

  Pixels x participants; the producer mean is computed and rendered as
  the heatmap layer. Alternatively, a numeric vector of length
  `prod(img_dims)` (e.g. `agreement_map_test(...)$observed_t`) which is
  rendered directly.

- base_image:

  Either a numeric matrix (`nrow x ncol`, grayscale, values in 0-1) or a
  path to a PNG/JPEG file. When it's a path, the image is read as
  grayscale (BT.709 luminance for RGB inputs) the same way as
  [`read_cis()`](https://olivethree.github.io/rcicrely/reference/read_cis.md).

- img_dims:

  Integer `c(nrow, ncol)`. If `NULL`, inferred from
  `attr(signal_matrix, "img_dims")` or, when `signal_matrix` is a
  vector, the dimensions of `base_image`.

- test:

  Optional
  [`agreement_map_test()`](https://olivethree.github.io/rcicrely/reference/agreement_map_test.md)
  result. When supplied, significance contours (around
  `significant_mask`) are drawn on top of the heatmap.

- mask:

  Optional logical vector restricting the visible overlay region; pixels
  outside the mask render as the bare base image.

- threshold:

  Optional numeric. Pixels with `|signal| < threshold` are rendered as
  fully transparent (only the base shows through). Default `NULL` (no
  threshold).

- alpha_max:

  Numeric in `[0, 1]`. Maximum opacity at the peak of the signal
  magnitude. Default 0.7.

- palette:

  Character. Currently only `"diverging"` is implemented
  (red-white-blue). Future-proofed by argument.

- contour_col, contour_lwd:

  Significance-contour colour and line width when `test` is supplied.
  Defaults: `"black"`, `1.0`.

- main:

  Optional plot title.

## Value

Invisibly the composed `nrow x ncol x 3` raster. The plot is drawn on
the active device as a side effect.

## See also

[`agreement_map_test()`](https://olivethree.github.io/rcicrely/reference/agreement_map_test.md),
[`face_mask()`](https://olivethree.github.io/rcicrely/reference/face_mask.md)
