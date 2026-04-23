# Read a directory of CI images and extract the base-subtracted signal

Convenience wrapper composing
[`read_cis()`](https://olivethree.github.io/rcicrely/reference/read_cis.md)
and
[`extract_signal()`](https://olivethree.github.io/rcicrely/reference/extract_signal.md).
Returns a pixels x participants signal matrix in one call. Use this when
you don't need to intervene between reading and extracting (e.g. masking
pixels, swapping the base).

## Usage

``` r
load_signal_matrix(
  dir,
  base_image_path,
  pattern = "\\.(png|jpe?g)$",
  acknowledge_scaling = FALSE
)
```

## Arguments

- dir:

  Directory of CI images (as for
  [`read_cis()`](https://olivethree.github.io/rcicrely/reference/read_cis.md)).

- base_image_path:

  Path to the base face image used at stimulus-generation time.

- pattern:

  Regex for file selection (see
  [`read_cis()`](https://olivethree.github.io/rcicrely/reference/read_cis.md)).

- acknowledge_scaling:

  Set to `TRUE` to silence the once-per-session warning. See
  [`read_cis()`](https://olivethree.github.io/rcicrely/reference/read_cis.md).

## Value

A numeric pixels x participants matrix, base-subtracted.

## Reading the result

Pixels x participants numeric matrix, base-subtracted. Carries an
`img_dims` attribute. Column names are the filename without the
extension.

## Common mistakes

- Pointing `base_image_path` at the wrong file (e.g. a generated CI PNG
  instead of the base face used at stimulus-generation time). The
  dimensions check will catch a wholly wrong image but not a
  wrong-base-of-the-right-size.

- Treating the result as raw mask. See **Note** below.

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
survive a single uniform linear scaling unmodified, and even those break
under per-CI `"matched"`-style scaling. The canonical 2IFC `infoVal`
path
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
[`extract_signal()`](https://olivethree.github.io/rcicrely/reference/extract_signal.md)

## Examples

``` r
if (FALSE) { # \dontrun{
signal <- load_signal_matrix(
  dir             = "data/cis_condition_A/",
  base_image_path = "data/base.jpg"
)
} # }
```
