#' Cluster-based permutation test with family-wise error control
#'
#' @description
#' Pixel-level discriminability test between two condition signal
#' matrices. Returns either (a) the classical threshold-based cluster
#' test (Maris & Oostenveld, 2007) with cluster mass as the statistic,
#' or (b) threshold-free cluster enhancement (TFCE; Smith & Nichols,
#' 2009), which integrates across thresholds and removes the arbitrary
#' cluster-forming cutoff. Both paths use the same stratified
#' permutation scheme and the maximum-statistic approach to family-
#' wise error control (Nichols & Holmes, 2002).
#'
#' Pair with [rel_dissimilarity()] when you also want an overall
#' magnitude summary.
#'
#' @details
#' Common scaffold (both methods):
#'
#' 1. Compute observed pixel-wise Welch t (see [pixel_t_test()]).
#' 2. Build the null via `n_permutations` **stratified** label
#'    permutations, each preserving `(N_A, N_B)` exactly.
#' 3. Max-statistic over the null controls FWER in the strong sense.
#'
#' Method-specific step:
#'
#' * `method = "threshold"` (default): threshold at
#'   `|t| > cluster_threshold`, label connected components with
#'   4-connectivity, cluster mass = sum of t-values within the
#'   cluster. Per observed cluster, p-value = fraction of null
#'   max-masses (matching sign) that equal or exceed the observed
#'   mass in absolute value.
#' * `method = "tfce"`: enhance the observed t-map into a TFCE map
#'   (integral over thresholds of `size^E * h^H dh`, pos and neg
#'   tails enhanced separately and recombined with sign preserved;
#'   H = 2.0, E = 0.5 default). Compute the same enhancement for each
#'   permuted t-map and record `max(|TFCE|)`. Per-pixel p-value =
#'   fraction of null max-TFCE values that equal or exceed the
#'   observed `|TFCE|` at that pixel. No free threshold parameter.
#'
#' @param signal_matrix_a,signal_matrix_b Pixels x participants,
#'   base-subtracted.
#' @param img_dims Integer `c(nrow, ncol)`. Can be inferred if the
#'   signal matrices carry an `img_dims` attribute.
#' @param method One of `"threshold"` (default) or `"tfce"`. The
#'   threshold-based method requires `cluster_threshold`; TFCE does
#'   not.
#' @param paired Logical. `FALSE` (default) for a between-subjects
#'   design: Welch t per pixel and stratified label permutation for
#'   the null. `TRUE` for a within-subjects design: paired t per
#'   pixel on `A - B` and random sign-flip of matched pairs for the
#'   null. When `paired = TRUE`, the two matrices must have the
#'   same number of columns and, if named, identical column names.
#' @param n_permutations Number of stratified label permutations
#'   (between-subjects) or sign-flip permutations (paired). Default
#'   2000.
#' @param cluster_threshold Absolute t-value threshold for forming
#'   clusters under `method = "threshold"`. Default 2.0. Ignored
#'   under `method = "tfce"`.
#' @param tfce_H,tfce_E TFCE height and extent exponents. Defaults
#'   match Smith & Nichols (2009): `H = 2.0`, `E = 0.5`. Ignored
#'   under `method = "threshold"`.
#' @param tfce_n_steps Number of thresholds in the TFCE integration
#'   grid. Default 100. Finer grids cost more per permutation.
#'   Ignored under `method = "threshold"`.
#' @param alpha Significance level. Default 0.05.
#' @param mask Optional logical vector of length `n_pixels`
#'   restricting cluster / TFCE inference to a region (e.g., from
#'   [face_mask()] or [load_face_mask()]). Implementation uses the
#'   **zero-out pattern**, not row-subsetting: pixels where
#'   `!mask` have their per-pixel t set to 0 (in both observed and
#'   permutation t-maps) before cluster identification. This way
#'   the 2D image structure is preserved for 4-connectivity logic,
#'   but masked-out pixels can never join a cluster
#'   (`|t| = 0 < cluster_threshold`) and contribute 0 to TFCE.
#'   Apply the same mask in any companion analyses for consistency.
#' @param seed Optional integer; RNG state restored on exit.
#' @param progress Show a `cli` progress bar.
#'
#' @section Reading the result (threshold method):
#' * `$observed_t`, per-pixel Welch t vector.
#' * `$clusters`, data.frame with `cluster_id`, `direction`
#'   (`"pos"`/`"neg"`), `mass`, `size`, `p_value`, `significant`.
#' * `$null_distribution$pos`, `$neg`, per-permutation max masses.
#' * `$pos_labels`, `$neg_labels`, integer matrices for plotting.
#' * `$cluster_threshold`, `$alpha`, `$n_permutations`,
#'   `$n_participants_a`, `$n_participants_b`, `$method = "threshold"`.
#'
#' @section Reading the result (TFCE method):
#' * `$observed_t`, per-pixel Welch t vector (before enhancement).
#' * `$tfce_map`, per-pixel signed TFCE values.
#' * `$tfce_pmap`, per-pixel p-values against the max-TFCE null.
#' * `$tfce_significant_mask`, logical vector flagging pixels with
#'   `p < alpha`.
#' * `$null_distribution$max_abs_tfce`, per-permutation max|TFCE|.
#' * `$tfce_H`, `$tfce_E`, `$tfce_n_steps`, `$alpha`,
#'   `$n_permutations`, `$n_participants_a`, `$n_participants_b`,
#'   `$method = "tfce"`.
#'
#' @section Common mistakes:
#' * Reading `cluster_threshold` as a p-value cutoff. It is the
#'   t-value above which pixels join a cluster; significance is
#'   decided afterwards by comparing cluster mass to the null.
#' * Trusting cluster significance with `n_permutations < 1000`.
#' * Permuting **pixels** instead of condition labels (the function
#'   does the right thing internally; this is a warning if you
#'   reimplement it yourself).
#'
#' @section Reliability metrics expect raw masks:
#' Welch t and cluster mass / TFCE are variance-based and sensitive
#' to any scaling. If `signal_matrix_a` / `_b` came from rendered
#' PNGs, results are distorted. See
#' `vignette("tutorial", package = "rcicrely")` chapter 3.
#'
#' @return Object of class `rcicrely_cluster_test`. Fields described
#'   in **Reading the result**.
#' @seealso [rel_dissimilarity()], [run_between()]
#' @references
#' Maris, E., & Oostenveld, R. (2007). Nonparametric statistical
#' testing of EEG- and MEG-data. *Journal of Neuroscience Methods*,
#' 164(1), 177-190. \doi{10.1016/j.jneumeth.2007.03.024}
#'
#' Nichols, T. E., & Holmes, A. P. (2002). Nonparametric permutation
#' tests for functional neuroimaging: a primer with examples.
#' *Human Brain Mapping*, 15(1), 1-25. \doi{10.1002/hbm.1058}
#'
#' Smith, S. M., & Nichols, T. E. (2009). Threshold-free cluster
#' enhancement: addressing problems of smoothing, threshold
#' dependence and localisation in cluster inference. *NeuroImage*,
#' 44(1), 83-98. \doi{10.1016/j.neuroimage.2008.03.061}
#' @export
rel_cluster_test <- function(signal_matrix_a,
                             signal_matrix_b,
                             img_dims          = NULL,
                             method            = c("threshold", "tfce"),
                             paired            = FALSE,
                             n_permutations    = 2000L,
                             cluster_threshold = 2.0,
                             tfce_H            = 2.0,
                             tfce_E            = 0.5,
                             tfce_n_steps      = 100L,
                             alpha             = 0.05,
                             mask              = NULL,
                             seed              = NULL,
                             progress          = TRUE) {
  method <- match.arg(method)
  validate_two_signal_matrices(signal_matrix_a, signal_matrix_b)
  if (isTRUE(paired)) {
    validate_paired_matrices(signal_matrix_a, signal_matrix_b)
  }
  if (looks_scaled(signal_matrix_a) || looks_scaled(signal_matrix_b)) {
    warn_looks_scaled("signal_matrix_a / _b")
  }
  n_pix   <- nrow(signal_matrix_a)
  n_a     <- ncol(signal_matrix_a)
  n_b     <- ncol(signal_matrix_b)
  n_total <- n_a + n_b

  img_dims <- resolve_img_dims(img_dims, signal_matrix_a, n_pix)
  n_permutations <- as.integer(n_permutations)
  if (n_permutations < 100L) {
    cli::cli_warn(
      "{.arg n_permutations} = {n_permutations} is low for cluster \\
       p-values; consider >= 1000."
    )
  }

  # `mask` for cluster_test uses a zero-out (rather than drop) pattern:
  # we keep the full 2D image so 4-connectivity / TFCE work, but zero
  # out per-pixel t-values where !mask. Masked-out pixels then never
  # join a cluster (|t|=0 < threshold) and never contribute to TFCE
  # (0^H = 0 in the integrand).
  if (!is.null(mask)) {
    if (!is.logical(mask) || length(mask) != n_pix) {
      cli::cli_abort(c(
        "{.arg mask} must be a logical vector of length {.code nrow}.",
        "*" = "{.arg mask}: {length(mask)}",
        "*" = "{.code nrow(signal_matrix_a)}: {n_pix}"
      ))
    }
    if (sum(mask) < 4L) {
      cli::cli_abort(
        "{.arg mask} selects too few pixels ({sum(mask)})."
      )
    }
  }

  # ---- observed stats --------------------------------------------------
  observed_t <- pixel_t_test(signal_matrix_a, signal_matrix_b,
                             paired = paired)
  if (!is.null(mask)) observed_t[!mask] <- 0

  if (method == "threshold") {
    obs <- find_clusters(observed_t, img_dims, cluster_threshold)
  } else {
    obs_tfce <- tfce_enhance(observed_t, img_dims,
                             H = tfce_H, E = tfce_E,
                             n_steps = tfce_n_steps)
  }

  # ---- permutation null ------------------------------------------------
  # For paired designs, the null schema is sign-flip on per-producer
  # difference rows. For between-subjects, stratified label permutation
  # over the concatenated matrix.
  row_var <- function(x, m, n_cols) {
    rowSums((x - m)^2) / (n_cols - 1L)
  }

  if (isTRUE(paired)) {
    diff_mat <- signal_matrix_a - signal_matrix_b  # pixels x N
    n_prod   <- ncol(diff_mat)
    perm_tmap <- function() {
      signs <- sample(c(-1L, 1L), n_prod, replace = TRUE)
      # flip sign on each producer's difference contribution
      d_flipped <- sweep(diff_mat, 2L, signs, `*`)
      mean_d <- rowMeans(d_flipped)
      var_d  <- rowSums((d_flipped - mean_d)^2) / (n_prod - 1L)
      se     <- sqrt(var_d / n_prod)
      t_perm <- mean_d / se
      t_perm[!is.finite(t_perm)] <- 0
      t_perm
    }
  } else {
    combined <- cbind(signal_matrix_a, signal_matrix_b)
    perm_tmap <- function() {
      perm  <- sample.int(n_total)
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
      t_perm
    }
  }

  pid <- progress_start(n_permutations,
                        if (method == "tfce") "cluster permutation (TFCE)"
                        else "cluster permutation",
                        show = progress)
  on.exit(progress_done(pid), add = TRUE)

  if (method == "threshold") {
    null_pos <- numeric(n_permutations)
    null_neg <- numeric(n_permutations)
    with_seed(seed, {
      for (i in seq_len(n_permutations)) {
        t_perm <- perm_tmap()
        if (!is.null(mask)) t_perm[!mask] <- 0
        cl <- find_clusters(t_perm, img_dims, cluster_threshold)
        null_pos[i] <- if (length(cl$pos_masses) > 0L)
          max(cl$pos_masses) else 0
        null_neg[i] <- if (length(cl$neg_masses) > 0L)
          min(cl$neg_masses) else 0
        progress_tick(pid)
      }
    })
  } else {
    null_max_abs_tfce <- numeric(n_permutations)
    with_seed(seed, {
      for (i in seq_len(n_permutations)) {
        t_perm <- perm_tmap()
        if (!is.null(mask)) t_perm[!mask] <- 0
        tfce_perm <- tfce_enhance(t_perm, img_dims,
                                  H = tfce_H, E = tfce_E,
                                  n_steps = tfce_n_steps)
        null_max_abs_tfce[i] <- max(abs(tfce_perm))
        progress_tick(pid)
      }
    })
  }

  # ---- build result ---------------------------------------------------
  if (method == "threshold") {
    clusters_df <- assemble_clusters_df(obs, null_pos, null_neg, alpha)
    return(new_rcicrely_cluster_test(
      observed_t         = observed_t,
      method             = "threshold",
      clusters           = clusters_df,
      null_distribution  = list(pos = null_pos, neg = null_neg),
      img_dims           = img_dims,
      pos_labels         = obs$pos_labels,
      neg_labels         = obs$neg_labels,
      cluster_threshold  = cluster_threshold,
      alpha              = alpha,
      n_permutations     = n_permutations,
      n_participants_a   = n_a,
      n_participants_b   = n_b,
      paired             = isTRUE(paired),
      tfce_map           = NULL,
      tfce_pmap          = NULL,
      tfce_significant_mask = NULL,
      tfce_H             = NA_real_,
      tfce_E             = NA_real_,
      tfce_n_steps       = NA_integer_
    ))
  }

  # tfce branch: per-pixel p-values against max-|TFCE| null
  abs_tfce <- abs(obs_tfce)
  p_map <- vapply(
    abs_tfce,
    function(v) (sum(null_max_abs_tfce >= v) + 1L) /
                (length(null_max_abs_tfce) + 1L),
    numeric(1L)
  )
  sig_mask <- p_map < alpha

  new_rcicrely_cluster_test(
    observed_t            = observed_t,
    method                = "tfce",
    clusters              = data.frame(
      cluster_id  = integer(0L),
      direction   = character(0L),
      mass        = numeric(0L),
      size        = integer(0L),
      p_value     = numeric(0L),
      significant = logical(0L),
      stringsAsFactors = FALSE
    ),
    null_distribution     = list(max_abs_tfce = null_max_abs_tfce),
    img_dims              = img_dims,
    pos_labels            = matrix(0L, img_dims[1L], img_dims[2L]),
    neg_labels            = matrix(0L, img_dims[1L], img_dims[2L]),
    cluster_threshold     = NA_real_,
    alpha                 = alpha,
    n_permutations        = n_permutations,
    n_participants_a      = n_a,
    n_participants_b      = n_b,
    paired                = isTRUE(paired),
    tfce_map              = obs_tfce,
    tfce_pmap             = p_map,
    tfce_significant_mask = sig_mask,
    tfce_H                = tfce_H,
    tfce_E                = tfce_E,
    tfce_n_steps          = as.integer(tfce_n_steps)
  )
}

#' Resolve img_dims from arg / attr / sqrt(n_pix)
#'
#' @keywords internal
#' @noRd
resolve_img_dims <- function(img_dims, signal_matrix, n_pix) {
  if (is.null(img_dims)) {
    attr_dims <- attr(signal_matrix, "img_dims")
    if (!is.null(attr_dims)) {
      img_dims <- attr_dims
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
  validate_img_dims(img_dims, n_pix)
}

#' Build the clusters data frame from find_clusters() output and
#' per-permutation max-mass vectors.
#'
#' @keywords internal
#' @noRd
assemble_clusters_df <- function(obs, null_pos, null_neg, alpha) {
  n_pos <- length(obs$pos_masses)
  n_neg <- length(obs$neg_masses)
  clusters_df <- data.frame(
    cluster_id = c(seq_len(n_pos), seq_len(n_neg)),
    direction  = c(rep("pos", n_pos), rep("neg", n_neg)),
    mass       = c(obs$pos_masses, obs$neg_masses),
    size       = c(
      if (n_pos) tabulate(obs$pos_labels[obs$pos_labels != 0L])
      else integer(0L),
      if (n_neg) tabulate(obs$neg_labels[obs$neg_labels != 0L])
      else integer(0L)
    ),
    p_value    = c(
      if (n_pos)
        vapply(obs$pos_masses,
               function(m) mean(null_pos >= m), numeric(1L))
      else numeric(0L),
      if (n_neg)
        vapply(obs$neg_masses,
               function(m) mean(null_neg <= m), numeric(1L))
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
  clusters_df
}
