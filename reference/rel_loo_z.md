# Z-scored leave-one-out influence (accessor)

Convenience accessor that returns a data frame of producer ids and their
z-scored LOO influence, ordered from most-influential (lowest, most
negative `z_score`) to least. Accepts either a signal matrix (runs
[`rel_loo()`](https://olivethree.github.io/rcicrely/reference/rel_loo.md)
under the hood) or an existing `rcicrely_loo` result object (cheap, no
recomputation).

Use this when you want the ordered ranking of producer influence without
the full `rcicrely_loo` list — e.g., for joining against producer-level
metadata, passing to `dplyr`/`ggplot2`, or a tidy Supplementary table.

## Usage

``` r
rel_loo_z(x, ...)
```

## Arguments

- x:

  Either a `pixels x participants` signal matrix or an object of class
  `rcicrely_loo` (as returned by
  [`rel_loo()`](https://olivethree.github.io/rcicrely/reference/rel_loo.md)).

- ...:

  Passed to
  [`rel_loo()`](https://olivethree.github.io/rcicrely/reference/rel_loo.md)
  when `x` is a signal matrix (e.g. `flag_threshold`, `flag_method`).
  Ignored when `x` is already a result object.

## Value

A data frame with columns `participant_id`, `correlation`, `z_score`,
`flag`, sorted by `z_score` ascending.

## See also

[`rel_loo()`](https://olivethree.github.io/rcicrely/reference/rel_loo.md)
