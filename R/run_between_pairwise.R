#' Run all pairwise between-condition comparisons across K conditions
#'
#' @description
#' Generalises [run_between()] from a 2-condition comparison to K
#' conditions: runs [rel_cluster_test()] and [rel_dissimilarity()] on
#' every K-choose-2 pair and applies a family-wise error correction
#' across pairs.
#'
#' @section FWER scope:
#' Cluster-level *p*-values within each pair are already
#' max-statistic FWER-controlled by [rel_cluster_test()]. The
#' Holm/Bonferroni adjustment in `run_between_pairwise()` controls
#' family-wise error across the K-choose-2 *pair comparisons* (a
#' second layer above the cluster test's internal control), not
#' over individual pixels or clusters.
#'
#' For each pair, the statistic carried into the across-pairs
#' adjustment is the **minimum cluster-level *p*-value within that
#' pair** (i.e., the most-significant cluster's *p*). Within-pair
#' cluster *p*-values are not re-adjusted: they remain the
#' max-statistic FWER-controlled values from `rel_cluster_test()`.
#' A pair with no clusters at all contributes `p_min = 1.0` so the
#' Holm ordering is well-defined.
#'
#' @param signal_matrices Named list of pixels x participants signal
#'   matrices, one per condition. Names become condition labels in
#'   the output. Must have at least 2 elements.
#' @param fwer One of `"holm"` (default), `"bonferroni"`, or
#'   `"none"`. The across-pairs family-wise correction applied to
#'   the per-pair minimum cluster *p*-value.
#' @param img_dims Integer `c(nrow, ncol)`. If `NULL`, inferred
#'   from the `img_dims` attribute on the first matrix.
#' @param paired Logical. When `TRUE`, all pairs use the paired
#'   variant. All matrices must have identical column names; rows
#'   correspond to the same producer across matrices.
#' @param method Cluster-test method. Forwarded to
#'   [rel_cluster_test()]. Default `"threshold"`.
#' @param n_permutations,n_boot,cluster_threshold,alpha,ci_level
#'   Forwarded to the per-pair `rel_cluster_test()` and
#'   `rel_dissimilarity()` calls.
#' @param dissim_null Forwarded to [rel_dissimilarity()] as `null`.
#'   Default `"none"` (skip the per-pair dissimilarity null to
#'   keep wall time bounded).
#' @param mask Optional logical vector forwarded to both
#'   per-pair calls.
#' @param seed Optional integer.
#' @param progress Show `cli` progress bars per per-pair call.
#' @param acknowledge_scaling Logical. Forwarded.
#' @return Object of class `rcicrely_pairwise_report` with:
#' * `$pairs`, data.frame with one row per K-choose-2 pair:
#'   columns `pair_id`, `cond_a`, `cond_b`, `n_clusters`,
#'   `p_min` (minimum within-pair cluster *p*-value), `p_adj_pair`
#'   (FWER-adjusted across pairs under `fwer`), `significant`
#'   (logical: `p_adj_pair < alpha`), `euclidean`,
#'   `euclidean_normalised`.
#' * `$results`, named list of per-pair `rel_cluster_test` and
#'   `rel_dissimilarity` results (keyed by `pair_id`).
#' * `$conditions`, the input names.
#' * `$fwer`, the correction method used.
#' * `$alpha`, the across-pairs alpha.
#' @seealso [run_between()], [rel_cluster_test()], [rel_dissimilarity()]
#' @export
run_between_pairwise <- function(signal_matrices,
                                 fwer                = c("holm",
                                                          "bonferroni",
                                                          "none"),
                                 img_dims            = NULL,
                                 paired              = FALSE,
                                 method              = c("threshold",
                                                          "tfce"),
                                 n_permutations      = 2000L,
                                 n_boot              = 2000L,
                                 cluster_threshold   = 2.0,
                                 alpha               = 0.05,
                                 ci_level            = 0.95,
                                 dissim_null         = c("none",
                                                          "permutation"),
                                 mask                = NULL,
                                 seed                = NULL,
                                 progress            = TRUE,
                                 acknowledge_scaling = FALSE) {
  fwer        <- match.arg(fwer)
  method      <- match.arg(method)
  dissim_null <- match.arg(dissim_null)

  if (!is.list(signal_matrices) || is.null(names(signal_matrices)) ||
        any(!nzchar(names(signal_matrices)))) {
    cli::cli_abort(
      "{.arg signal_matrices} must be a named list of signal matrices."
    )
  }
  if (length(signal_matrices) < 2L) {
    cli::cli_abort(c(
      "{.arg signal_matrices} must contain at least 2 conditions.",
      "*" = "Got {length(signal_matrices)}."
    ))
  }
  conditions <- names(signal_matrices)
  for (nm in conditions) {
    validate_signal_matrix(signal_matrices[[nm]], name = nm)
  }
  if (is.null(img_dims)) {
    img_dims <- attr(signal_matrices[[1L]], "img_dims")
  }

  pair_idx <- utils::combn(length(conditions), 2L, simplify = FALSE)
  pair_ids <- vapply(pair_idx, function(p) {
    sprintf("%s_vs_%s", conditions[p[1L]], conditions[p[2L]])
  }, character(1L))

  per_pair <- vector("list", length(pair_idx))
  names(per_pair) <- pair_ids

  pair_rows <- data.frame(
    pair_id              = pair_ids,
    cond_a               = vapply(pair_idx, function(p) conditions[p[1L]],
                                  character(1L)),
    cond_b               = vapply(pair_idx, function(p) conditions[p[2L]],
                                  character(1L)),
    n_clusters           = NA_integer_,
    p_min                = NA_real_,
    p_adj_pair           = NA_real_,
    significant          = NA,
    euclidean            = NA_real_,
    euclidean_normalised = NA_real_,
    stringsAsFactors     = FALSE
  )

  for (i in seq_along(pair_idx)) {
    p   <- pair_idx[[i]]
    a_m <- signal_matrices[[p[1L]]]
    b_m <- signal_matrices[[p[2L]]]
    if (isTRUE(paired)) {
      validate_paired_matrices(a_m, b_m,
                                name_a = conditions[p[1L]],
                                name_b = conditions[p[2L]])
    }

    ct <- rel_cluster_test(
      a_m, b_m,
      img_dims            = img_dims,
      method              = method,
      paired              = paired,
      n_permutations      = n_permutations,
      cluster_threshold   = cluster_threshold,
      alpha               = alpha,
      mask                = mask,
      seed                = seed,
      progress            = progress,
      acknowledge_scaling = acknowledge_scaling
    )
    dr <- rel_dissimilarity(
      a_m, b_m,
      paired              = paired,
      n_boot              = n_boot,
      ci_level            = ci_level,
      null                = dissim_null,
      mask                = mask,
      seed                = seed,
      progress            = progress,
      acknowledge_scaling = acknowledge_scaling
    )

    per_pair[[i]] <- list(cluster_test = ct, dissimilarity = dr)

    cl <- ct$clusters
    pair_rows$n_clusters[i] <- if (is.null(cl) || nrow(cl) == 0L) {
      0L
    } else nrow(cl)
    p_min <- if (is.null(cl) || nrow(cl) == 0L) {
      1.0
    } else min(cl$p_value, na.rm = TRUE)
    pair_rows$p_min[i] <- p_min
    pair_rows$euclidean[i]            <- dr$euclidean
    pair_rows$euclidean_normalised[i] <- dr$euclidean_normalised
  }

  pair_rows$p_adj_pair <- switch(
    fwer,
    holm       = stats::p.adjust(pair_rows$p_min, method = "holm"),
    bonferroni = stats::p.adjust(pair_rows$p_min, method = "bonferroni"),
    none       = pair_rows$p_min
  )
  pair_rows$significant <- pair_rows$p_adj_pair < alpha

  structure(
    list(
      pairs            = pair_rows,
      results          = per_pair,
      conditions       = conditions,
      fwer             = fwer,
      alpha            = alpha,
      paired           = isTRUE(paired),
      method           = method
    ),
    class            = c("rcicrely_pairwise_report", "rcicrely_report"),
    rcicrely_version = utils::packageVersion("rcicrely")
  )
}
