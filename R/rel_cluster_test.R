#' Cluster-based permutation test with max-statistic FWER control
#'
#' @description
#' Pixel-level discriminability test between two condition signal
#' matrices. Use this to identify **where** in the image two
#' conditions' CIs differ, with family-wise error control across
#' pixels. Pair with [rel_dissimilarity()] for an overall magnitude. Implements the Maris & Oostenveld (2007) approach with
#' Nichols & Holmes (2002) max-statistic family-wise error control:
#'
#' 1. Compute observed pixel-wise Welch t.
#' 2. Threshold at `|t| > cluster_threshold` (separately for positive
#'    and negative tails); label connected components with
#'    **4-connectivity**.
#' 3. For each observed cluster, cluster mass = sum of t-values in the
#'    cluster (not pixel count).
#' 4. Build the null distribution via `n_permutations` **stratified**
#'    label permutations (each preserves `(N_A, N_B)`). For each, find
#'    clusters and record max positive mass and max negative mass.
#' 5. Observed cluster's p-value = fraction of null max-masses
#'    (matching sign) that equal or exceed the observed mass in
#'    absolute value.
#'
#' Stratification (sec.7.5): unstratified permutation lets the two group
#' sizes drift, which changes Welch degrees of freedom under the null
#' and biases the max-mass distribution. We preserve `(N_A, N_B)`
#' exactly.
#'
#' @param signal_matrix_a,signal_matrix_b Pixels x participants,
#'   base-subtracted.
#' @param img_dims Integer `c(nrow, ncol)`. Can be inferred if the
#'   signal matrices carry an `img_dims` attribute.
#' @param n_permutations Number of label permutations. 2000 default.
#' @param cluster_threshold Absolute t-value threshold for forming
#'   clusters. Default 2.0 (~ two-tailed p < 0.05 at large df).
#' @param alpha Significance level for flagging clusters. Default
#'   0.05.
#' @param seed Optional integer; RNG state restored on exit.
#' @param progress Show a `cli` progress bar.
#' @section Reading the result:
#' * `$observed_t`, per-pixel Welch t vector for the observed
#'   condition labelling.
#' * `$clusters`, data.frame with `cluster_id`, `direction`
#'   (`"pos"`/`"neg"`), `mass`, `size`, `p_value`, `significant`. One
#'   row per supra-threshold cluster, sorted by direction then mass.
#' * `$null_distribution`, list with `$pos` and `$neg` vectors of
#'   per-permutation max masses, useful for plotting the null.
#' * `$pos_labels`, `$neg_labels`, integer matrices the same shape
#'   as the image, where each non-zero value labels a cluster. Drives
#'   the contour overlay in the result's `plot()` method.
#' * `$cluster_threshold`, `$alpha`, `$n_permutations`,
#'   `$n_participants_a`, `$n_participants_b`, metadata.
#'
#' @section Common mistakes:
#' * Reading `cluster_threshold` (default 2.0) as a p-value cutoff. It
#'   is the t-value above which pixels join a cluster; significance
#'   is decided afterwards by comparing cluster mass to the null.
#' * Trusting cluster significance with `n_permutations < 1000`
#'   (tail probabilities are noisy at low iter counts).
#' * Permuting **pixels** instead of condition labels (the function
#'   does the right thing internally; this is just a warning if you
#'   reimplement it yourself).
#'
#' @section Reliability metrics expect raw masks:
#' Welch t and cluster mass are variance-based and sensitive to any
#' scaling. If `signal_matrix_a` / `_b` came from rendered PNGs,
#' results are distorted. See
#' `vignette("tutorial", package = "rcicrely")` chapter 3.
#'
#' @return Object of class `rcicrely_cluster_test`. Fields described
#'   in **Reading the result** above.
#' @seealso [rel_dissimilarity()], [run_between()]
#' @references
#' Maris, E., & Oostenveld, R. (2007). Nonparametric statistical
#' testing of EEG- and MEG-data. *Journal of Neuroscience Methods*,
#' 164(1), 177-190. \doi{10.1016/j.jneumeth.2007.03.024}
#'
#' Nichols, T. E., & Holmes, A. P. (2002). Nonparametric permutation
#' tests for functional neuroimaging: a primer with examples.
#' *Human Brain Mapping*, 15(1), 1-25. \doi{10.1002/hbm.1058}
#' @export
rel_cluster_test <- function(signal_matrix_a,
                             signal_matrix_b,
                             img_dims          = NULL,
                             n_permutations    = 2000L,
                             cluster_threshold = 2.0,
                             alpha             = 0.05,
                             seed              = NULL,
                             progress          = TRUE) {
  validate_two_signal_matrices(signal_matrix_a, signal_matrix_b)
  if (looks_scaled(signal_matrix_a) || looks_scaled(signal_matrix_b)) {
    warn_looks_scaled("signal_matrix_a / _b")
  }
  n_pix <- nrow(signal_matrix_a)
  n_a <- ncol(signal_matrix_a)
  n_b <- ncol(signal_matrix_b)
  n_total <- n_a + n_b

  if (is.null(img_dims)) {
    attr_a <- attr(signal_matrix_a, "img_dims")
    if (!is.null(attr_a)) {
      img_dims <- attr_a
    } else {
      side <- sqrt(n_pix)
      if (side != as.integer(side)) {
        cli::cli_abort(
          "Cannot infer {.arg img_dims}, pixel count {n_pix} is \\
           not a perfect square. Pass {.arg img_dims} explicitly."
        )
      }
      img_dims <- c(as.integer(side), as.integer(side))
    }
  }
  img_dims <- validate_img_dims(img_dims, n_pix)

  n_permutations <- as.integer(n_permutations)
  if (n_permutations < 100L) {
    cli::cli_warn(
      "{.arg n_permutations} = {n_permutations} is low for cluster \\
       p-values; consider >= 1000."
    )
  }

  # ---- observed stats ---------------------------------------------------

  observed_t <- pixel_t_test(signal_matrix_a, signal_matrix_b)
  obs <- find_clusters(observed_t, img_dims, cluster_threshold)

  # ---- build null via stratified label permutation ---------------------

  combined <- cbind(signal_matrix_a, signal_matrix_b)
  null_pos <- numeric(n_permutations)
  null_neg <- numeric(n_permutations)

  # row-variance helper that takes a precomputed rowMeans for speed
  row_var <- function(x, m, n_cols) {
    rowSums((x - m)^2) / (n_cols - 1L)
  }

  pid <- progress_start(n_permutations, "cluster permutation",
                        show = progress)
  on.exit(progress_done(pid), add = TRUE)

  with_seed(seed, {
    for (i in seq_len(n_permutations)) {
      perm <- sample.int(n_total)
      idx_a <- perm[seq_len(n_a)]
      idx_b <- perm[(n_a + 1L):n_total]
      a_mat <- combined[, idx_a, drop = FALSE]
      b_mat <- combined[, idx_b, drop = FALSE]

      ma <- rowMeans(a_mat)
      mb <- rowMeans(b_mat)
      va <- row_var(a_mat, ma, n_a)
      vb <- row_var(b_mat, mb, n_b)
      se <- sqrt(va / n_a + vb / n_b)
      t_perm <- (ma - mb) / se
      t_perm[!is.finite(t_perm)] <- 0

      cl <- find_clusters(t_perm, img_dims, cluster_threshold)
      null_pos[i] <- if (length(cl$pos_masses) > 0L)
        max(cl$pos_masses) else 0
      null_neg[i] <- if (length(cl$neg_masses) > 0L)
        min(cl$neg_masses) else 0

      progress_tick(pid)
    }
  })

  # ---- p-values ---------------------------------------------------------

  n_pos <- length(obs$pos_masses)
  n_neg <- length(obs$neg_masses)
  clusters_df <- data.frame(
    cluster_id  = c(seq_len(n_pos), seq_len(n_neg)),
    direction   = c(rep("pos", n_pos), rep("neg", n_neg)),
    mass        = c(obs$pos_masses, obs$neg_masses),
    size        = c(
      if (n_pos) tabulate(obs$pos_labels[obs$pos_labels != 0L]) else
        integer(0L),
      if (n_neg) tabulate(obs$neg_labels[obs$neg_labels != 0L]) else
        integer(0L)
    ),
    p_value     = c(
      if (n_pos)
        vapply(obs$pos_masses,
               function(m) mean(null_pos >= m),
               numeric(1L))
      else numeric(0L),
      if (n_neg)
        vapply(obs$neg_masses,
               function(m) mean(null_neg <= m),
               numeric(1L))
      else numeric(0L)
    ),
    stringsAsFactors = FALSE
  )
  clusters_df$significant <- clusters_df$p_value < alpha
  clusters_df <- clusters_df[
    order(clusters_df$direction, -abs(clusters_df$mass)), ,
    drop = FALSE
  ]
  rownames(clusters_df) <- NULL

  new_rcicrely_cluster_test(
    observed_t         = observed_t,
    clusters           = clusters_df,
    null_distribution  = list(pos = null_pos, neg = null_neg),
    img_dims           = img_dims,
    pos_labels         = obs$pos_labels,
    neg_labels         = obs$neg_labels,
    cluster_threshold  = cluster_threshold,
    alpha              = alpha,
    n_permutations     = n_permutations,
    n_participants_a   = n_a,
    n_participants_b   = n_b
  )
}
