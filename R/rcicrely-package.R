#' rcicrely: reliability assessment for reverse-correlation CIs
#'
#' Within-condition consistency (permuted split-half, leave-one-out,
#' ICC(3,\*)) and between-condition discriminability (Welch pixel t,
#' cluster-based permutation with max-statistic FWER control, bootstrap
#' dissimilarity) for reverse-correlation classification images.
#'
#' Works with both standard 2IFC (via the canonical `rcicr` package) and
#' Brief-RC 12 (Schmitz, Rougier & Yzerbyt, 2024, implemented natively
#' in this package). Operates directly on the pixel-level signal
#' produced by the original producers, so reliability does not depend
#' on a second-phase trait-rating study.
#'
#' See `vignette("tutorial", package = "rcicrely")` for a complete
#' end-to-end walkthrough.
#'
#' @keywords internal
"_PACKAGE"
