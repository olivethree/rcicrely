# Compute only the empirical-null distribution for split-half

Standalone access to the null-distribution simulation used by
`rel_split_half(null = ...)`. Useful when the user wants to precompute
and reuse a null across multiple observed analyses (e.g., when comparing
across conditions with the same producer count).

## Usage

``` r
rel_split_half_null(
  n_participants,
  n_pixels,
  null = c("permutation", "random_responders"),
  noise_matrix = NULL,
  n_permutations = 2000L,
  seed = NULL,
  pid = NULL
)
```

## Arguments

- n_participants:

  Integer. Number of producers per side of the split. Should match the
  observed analysis's `ncol(signal_matrix)`.

- n_pixels:

  Integer. Number of pixels in the signal vector (used by
  `null = "permutation"` to size the random columns).

- null:

  Either `"permutation"` or `"random_responders"`.

- noise_matrix:

  Required for `"random_responders"`. Pixels x pool-size numeric matrix.

- n_permutations:

  Integer. Default 2000.

- seed:

  Optional integer.

- pid:

  Optional internal progress-bar id (used when called from inside
  [`rel_split_half()`](https://olivethree.github.io/rcicrely/reference/rel_split_half.md));
  end users should leave at `NULL`.

## Value

Numeric vector of length `n_permutations` containing per-iteration
`r_hh` values under the chosen null.
