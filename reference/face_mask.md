# Oval face-region mask for a square image

Returns a logical vector of length `prod(img_dims)` marking an
elliptical face region (or a face sub-region) centred on the image. Pass
to
[`infoval()`](https://olivethree.github.io/rcicrely/reference/infoval.md)
to restrict both observed and reference Frobenius norms to the masked
region; pass to any `rel_*()` function via row-subsetting
(`signal_matrix[mask, ]`) to compute reliability or discriminability on
a single anatomical region.

Five regions are supported:

- `"full"` (default): the full face oval (Schmitz, Rougier, & Yzerbyt
  2024 geometry).

- `"eyes"`: two small ellipses at the typical eye positions.

- `"nose"`: a narrow vertical ellipse along the midline.

- `"mouth"`: a wide-and-short ellipse below centre.

- `"upper_face"`, `"lower_face"`: the top and bottom halves of the full
  face oval, useful for crude top/bottom comparisons.

All region geometries are heuristic approximations matched to a typical
centred face on a square base image (e.g., 256x256 male face from the
KDEF). For non-default base images, use `centre`, `half_width`,
`half_height` to tune the full-face oval; the sub-region geometries
scale relative to that ellipse.

## Usage

``` r
face_mask(
  img_dims,
  region = c("full", "eyes", "nose", "mouth", "upper_face", "lower_face"),
  centre = c(0.5, 0.5),
  half_width = 0.35,
  half_height = 0.45
)
```

## Arguments

- img_dims:

  Integer `c(nrow, ncol)`.

- region:

  Character. One of `"full"` (default), `"eyes"`, `"nose"`, `"mouth"`,
  `"upper_face"`, `"lower_face"`.

- centre:

  Numeric `c(row, col)` in (0, 1) coordinates. Default `c(0.5, 0.5)`.
  Defines the centre of the full-face oval; sub- region positions are
  calibrated relative to this centre.

- half_width:

  Full-face ellipse horizontal half-axis, in (0, 1) fraction of image
  width. Default 0.35.

- half_height:

  Full-face ellipse vertical half-axis, in (0, 1) fraction of image
  height. Default 0.45.

## Value

Logical vector of length `prod(img_dims)`, column-major to match the
package's image vectorisation convention.

## References

Schmitz, M., Rougier, M., & Yzerbyt, V. (2024). Introducing the brief
reverse correlation: an improved tool to assess visual representations.
*European Journal of Social Psychology*.
[doi:10.1002/ejsp.3100](https://doi.org/10.1002/ejsp.3100)

## See also

[`infoval()`](https://olivethree.github.io/rcicrely/reference/infoval.md),
[`plot_agreement_map()`](https://olivethree.github.io/rcicrely/reference/plot_agreement_map.md)

## Examples

``` r
full  <- face_mask(c(128L, 128L))                   # ~0.49 of image
eyes  <- face_mask(c(128L, 128L), region = "eyes")  # ~0.02 of image
mouth <- face_mask(c(128L, 128L), region = "mouth") # ~0.04 of image
c(full = mean(full), eyes = mean(eyes), mouth = mean(mouth))
#>       full       eyes      mouth 
#> 0.48754883 0.01538086 0.01440430 
```
