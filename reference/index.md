# Package index

## Loading data

Read CI images, build the signal matrix.

- [`load_signal_matrix()`](https://olivethree.github.io/rcicrely/reference/load_signal_matrix.md)
  : Read a directory of CI images and extract the base-subtracted signal

- [`read_cis()`](https://olivethree.github.io/rcicrely/reference/read_cis.md)
  : Read a directory of CI images into a pixels x participants matrix

- [`extract_signal()`](https://olivethree.github.io/rcicrely/reference/extract_signal.md)
  : Subtract the base image from each column of a CI matrix

- [`read_noise_matrix()`](https://olivethree.github.io/rcicrely/reference/read_noise_matrix.md)
  :

  Load the noise matrix from an rcicr `.Rdata`, an `.rds`, or a text
  file

## CI from raw responses

Compute individual CIs from trial-level response data.

- [`ci_from_responses_2ifc()`](https://olivethree.github.io/rcicrely/reference/ci_from_responses_2ifc.md)
  : Compute individual 2IFC CIs from trial-level responses
- [`ci_from_responses_briefrc()`](https://olivethree.github.io/rcicrely/reference/ci_from_responses_briefrc.md)
  : Compute individual Brief-RC CIs from trial-level responses

## Within-condition reliability

How consistent is a single condition’s CI?

- [`rel_split_half()`](https://olivethree.github.io/rcicrely/reference/rel_split_half.md)
  : Permuted split-half reliability with Spearman-Brown correction
- [`rel_loo()`](https://olivethree.github.io/rcicrely/reference/rel_loo.md)
  : Leave-one-out sensitivity
- [`rel_icc()`](https://olivethree.github.io/rcicrely/reference/rel_icc.md)
  : Intraclass correlation coefficients via direct mean squares
- [`run_within()`](https://olivethree.github.io/rcicrely/reference/run_within.md)
  : Run every within-condition reliability metric

## Between-condition inference

Is condition A’s CI distinguishable from condition B’s?

- [`pixel_t_test()`](https://olivethree.github.io/rcicrely/reference/pixel_t_test.md)
  : Vectorised pixel-wise Welch t-test
- [`rel_cluster_test()`](https://olivethree.github.io/rcicrely/reference/rel_cluster_test.md)
  : Cluster-based permutation test with max-statistic FWER control
- [`rel_dissimilarity()`](https://olivethree.github.io/rcicrely/reference/rel_dissimilarity.md)
  : Representational dissimilarity with bootstrap CIs
- [`run_between()`](https://olivethree.github.io/rcicrely/reference/run_between.md)
  : Run every between-condition discriminability metric
