## S3 constructors for result classes. None are exported, users get
## these back from `rel_*()`, `run_within()`, `run_between()`.
## Print / summary / plot methods live in R/plot_methods.R.

#' Build an S3 result with shared metadata
#'
#' Centralised wrapper that attaches `rcicrely_version` to every
#' result object (§1.8). Future fields with the same lifetime
#' (e.g., a build date) belong here.
#'
#' @keywords internal
#' @noRd
new_rcicrely_object <- function(payload, subclass) {
  structure(
    payload,
    class            = c(subclass, "rcicrely_result"),
    rcicrely_version = utils::packageVersion("rcicrely")
  )
}

new_rcicrely_split_half <- function(r_hh, r_sb, ci_95, ci_95_sb,
                                    distribution, n_participants,
                                    n_permutations,
                                    null              = "none",
                                    null_distribution = NULL,
                                    r_hh_null_p95     = NA_real_,
                                    r_hh_excess       = NA_real_,
                                    r_sb_excess       = NA_real_) {
  new_rcicrely_object(
    list(
      r_hh              = r_hh,
      r_sb              = r_sb,
      ci_95             = ci_95,
      ci_95_sb          = ci_95_sb,
      distribution      = distribution,
      null              = null,
      null_distribution = null_distribution,
      r_hh_null_p95     = r_hh_null_p95,
      r_hh_excess       = r_hh_excess,
      r_sb_excess       = r_sb_excess,
      n_participants    = n_participants,
      n_permutations    = n_permutations
    ),
    subclass = "rcicrely_split_half"
  )
}

new_rcicrely_loo <- function(correlations, z_scores,
                             mean_r, sd_r,
                             median_r, mad_r,
                             threshold, flagged, summary_df,
                             flag_method, flag_threshold) {
  new_rcicrely_object(
    list(
      correlations   = correlations,
      z_scores       = z_scores,
      mean_r         = mean_r,
      sd_r           = sd_r,
      median_r       = median_r,
      mad_r          = mad_r,
      threshold      = threshold,
      flagged        = flagged,
      summary_df     = summary_df,
      flag_method    = flag_method,
      flag_threshold = flag_threshold
    ),
    subclass = "rcicrely_loo"
  )
}

new_rcicrely_icc <- function(icc_3_1, icc_3_k, icc_2_1, icc_2_k,
                             ms_rows, ms_cols, ms_error,
                             n_raters, n_targets, model, variants) {
  new_rcicrely_object(
    list(
      icc_3_1   = icc_3_1,
      icc_3_k   = icc_3_k,
      icc_2_1   = icc_2_1,
      icc_2_k   = icc_2_k,
      ms_rows   = ms_rows,
      ms_cols   = ms_cols,
      ms_error  = ms_error,
      n_raters  = n_raters,
      n_targets = n_targets,
      model     = model,
      variants  = variants
    ),
    subclass = "rcicrely_icc"
  )
}

new_rcicrely_cluster_test <- function(observed_t, clusters,
                                      null_distribution, img_dims,
                                      pos_labels, neg_labels,
                                      cluster_threshold, alpha,
                                      n_permutations,
                                      n_participants_a,
                                      n_participants_b,
                                      method = "threshold",
                                      paired = FALSE,
                                      tfce_map = NULL,
                                      tfce_pmap = NULL,
                                      tfce_significant_mask = NULL,
                                      tfce_H = NA_real_,
                                      tfce_E = NA_real_,
                                      tfce_n_steps = NA_integer_) {
  new_rcicrely_object(
    list(
      observed_t             = observed_t,
      method                 = method,
      paired                 = paired,
      clusters               = clusters,
      null_distribution      = null_distribution,
      img_dims               = img_dims,
      pos_labels             = pos_labels,
      neg_labels             = neg_labels,
      cluster_threshold      = cluster_threshold,
      alpha                  = alpha,
      n_permutations         = n_permutations,
      n_participants_a       = n_participants_a,
      n_participants_b       = n_participants_b,
      tfce_map               = tfce_map,
      tfce_pmap              = tfce_pmap,
      tfce_significant_mask  = tfce_significant_mask,
      tfce_H                 = tfce_H,
      tfce_E                 = tfce_E,
      tfce_n_steps           = tfce_n_steps
    ),
    subclass = "rcicrely_cluster_test"
  )
}

new_rcicrely_dissim <- function(correlation, euclidean,
                                euclidean_normalised,
                                boot_cor, boot_dist,
                                ci_cor, ci_dist,
                                boot_se_cor, boot_se_dist,
                                n_boot, ci_level,
                                n_pixels,
                                null              = "none",
                                null_distribution = NULL,
                                d_null_p95        = NA_real_,
                                d_z               = NA_real_,
                                d_ratio           = NA_real_,
                                paired            = FALSE) {
  new_rcicrely_object(
    list(
      correlation          = correlation,
      euclidean            = euclidean,
      euclidean_normalised = euclidean_normalised,
      boot_cor             = boot_cor,
      boot_dist            = boot_dist,
      ci_cor               = ci_cor,
      ci_dist              = ci_dist,
      boot_se_cor          = boot_se_cor,
      boot_se_dist         = boot_se_dist,
      n_boot               = n_boot,
      ci_level             = ci_level,
      n_pixels             = n_pixels,
      null                 = null,
      null_distribution    = null_distribution,
      d_null_p95           = d_null_p95,
      d_z                  = d_z,
      d_ratio              = d_ratio,
      paired               = paired
    ),
    subclass = "rcicrely_dissim"
  )
}

new_rcicrely_infoval <- function(infoval, norms, reference,
                                 ref_median, ref_mad, trial_counts,
                                 mask, iter, n_pool, seed) {
  new_rcicrely_object(
    list(
      infoval      = infoval,
      norms        = norms,
      reference    = reference,
      ref_median   = ref_median,
      ref_mad      = ref_mad,
      trial_counts = trial_counts,
      mask         = mask,
      iter         = iter,
      n_pool       = n_pool,
      seed         = seed
    ),
    subclass = "rcicrely_infoval"
  )
}

new_rcicrely_report <- function(results, method, img_dims = NULL) {
  structure(
    list(
      results  = results,
      method   = method,
      img_dims = img_dims
    ),
    class            = "rcicrely_report",
    rcicrely_version = utils::packageVersion("rcicrely")
  )
}
