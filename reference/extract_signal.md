# Subtract the base image from each column of a CI matrix

The pixel values returned by
[`read_cis()`](https://olivethree.github.io/rcicrely/reference/read_cis.md)
contain the shared base image plus a (typically scaled) noise
contribution. Reliability computations must operate on the **noise**
alone, otherwise the shared base inflates inter-participant
correlations. This function reads the base image, converts to grayscale
if needed, validates dimensions, and returns `cis - base` column-wise.

## Usage

``` r
extract_signal(cis, base_image_path, acknowledge_scaling = FALSE)
```

## Arguments

- cis:

  A numeric matrix from
  [`read_cis()`](https://olivethree.github.io/rcicrely/reference/read_cis.md)
  (pixels x participants).

- base_image_path:

  Path to the base face image. PNG or JPEG.

- acknowledge_scaling:

  Set to `TRUE` to silence the once-per-session warning. See
  [`read_cis()`](https://olivethree.github.io/rcicrely/reference/read_cis.md).

## Value

A numeric matrix the same shape as `cis`, base-subtracted.

## Reading the result

Numeric matrix the same shape as `cis`. Column names propagate. The
`img_dims` attribute is preserved (or filled in from the base image if
missing).

## Common mistakes

- Passing a base image with different dimensions, aborts with a clear
  pixel-count mismatch.

- Treating the result as the raw mask. The output is `scaling(mask)` if
  the input PNGs were rendered with any scaling. See **Raw vs. rendered
  CIs** below.

## Raw vs. rendered CIs

PNG pixel values are necessarily what was rendered to disk:
`base + scaling(mask)`. Subtracting the base in `extract_signal()`
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
survive a single uniform linear scaling unmodified, and even those break
under per-CI `"matched"`-style scaling. The standard 2IFC `infoVal` path
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

[`read_cis()`](https://olivethree.github.io/rcicrely/reference/read_cis.md),
[`load_signal_matrix()`](https://olivethree.github.io/rcicrely/reference/load_signal_matrix.md)
