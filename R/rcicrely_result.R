## S3 constructors for result classes. None are exported, users get
## these back from `rel_*()`, `run_within()`, `run_between()`.
## Print / summary / plot methods live in R/plot_methods.R.

new_rcicrely_split_half <- function(r_hh, r_sb, ci_95, ci_95_sb,
                                    distribution, n_participants,
                                    n_permutations) {
  structure(
    list(
      r_hh            = r_hh,
      r_sb            = r_sb,
      ci_95           = ci_95,
      ci_95_sb        = ci_95_sb,
      distribution    = distribution,
      n_participants  = n_participants,
      n_permutations  = n_permutations
    ),
    class = c("rcicrely_split_half", "rcicrely_result")
  )
}

new_rcicrely_loo <- function(correlations, mean_r, sd_r,
                             median_r, mad_r,
                             threshold, flagged, summary_df,
                             flag_method, flag_threshold) {
  structure(
    list(
      correlations   = correlations,
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
    class = c("rcicrely_loo", "rcicrely_result")
  )
}

new_rcicrely_icc <- function(icc_3_1, icc_3_k, icc_2_1, icc_2_k,
                             ms_rows, ms_cols, ms_error,
                             n_raters, n_targets, model, variants) {
  structure(
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
    class = c("rcicrely_icc", "rcicrely_result")
  )
}

new_rcicrely_cluster_test <- function(observed_t, clusters,
                                      null_distribution, img_dims,
                                      pos_labels, neg_labels,
                                      cluster_threshold, alpha,
                                      n_permutations,
                                      n_participants_a,
                                      n_participants_b) {
  structure(
    list(
      observed_t         = observed_t,
      clusters           = clusters,
      null_distribution  = null_distribution,
      img_dims           = img_dims,
      pos_labels         = pos_labels,
      neg_labels         = neg_labels,
      cluster_threshold  = cluster_threshold,
      alpha              = alpha,
      n_permutations     = n_permutations,
      n_participants_a   = n_participants_a,
      n_participants_b   = n_participants_b
    ),
    class = c("rcicrely_cluster_test", "rcicrely_result")
  )
}

new_rcicrely_dissim <- function(correlation, euclidean,
                                boot_cor, boot_dist,
                                ci_cor, ci_dist,
                                boot_se_cor, boot_se_dist,
                                n_boot, ci_level) {
  structure(
    list(
      correlation   = correlation,
      euclidean     = euclidean,
      boot_cor      = boot_cor,
      boot_dist     = boot_dist,
      ci_cor        = ci_cor,
      ci_dist       = ci_dist,
      boot_se_cor   = boot_se_cor,
      boot_se_dist  = boot_se_dist,
      n_boot        = n_boot,
      ci_level      = ci_level
    ),
    class = c("rcicrely_dissim", "rcicrely_result")
  )
}

new_rcicrely_report <- function(results, method, img_dims = NULL) {
  structure(
    list(
      results  = results,
      method   = method,
      img_dims = img_dims
    ),
    class = c("rcicrely_report")
  )
}
