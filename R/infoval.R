#' Per-producer informational value with trial-count-matched reference
#'
#' @description
#' Computes a z-scored informational value (infoVal) for each
#' producer's classification image, using a reference distribution
#' matched to that producer's trial count. Handles both 2IFC and
#' Brief-RC paradigms with a single function: the difference is
#' entirely in what the user passes as `noise_matrix`, not in how
#' the statistic is computed.
#'
#' @details
#' The observed statistic is the Frobenius norm of producer j's
#' classification image mask, optionally restricted to a logical
#' `mask` (e.g. an oval face region via [face_mask()]):
#'
#' ```
#' norm_j = sqrt(sum( signal_matrix[mask, j]^2 ))
#' ```
#'
#' The reference distribution for producer j is built by simulating
#' `iter` random masks at the same trial count `trial_counts[j]`,
#' mirroring the `genMask()` construction exactly (Schmitz, Rougier,
#' & Yzerbyt, 2024):
#'
#' 1. Sample `trial_counts[j]` stimulus ids uniformly from
#'    `1:n_pool` with replacement.
#' 2. Sample `trial_counts[j]` responses uniformly from `{-1, +1}`.
#' 3. Collapse `response` by `stim` via `mean()` (this is the
#'    duplicate-stim rule used by Brief-RC).
#' 4. Build the mask: `noise_matrix[, unique_stims] %*%
#'    mean_response / n_unique_stims`.
#' 5. Apply `mask` (if supplied) and compute the Frobenius norm.
#'
#' The per-producer z-score is
#' `(norm_j - median(ref)) / mad(ref)`, where the reference
#' distribution is the one matched to producer j's trial count.
#' Producers sharing a trial count share a reference; this makes the
#' simulation efficient when several producers have the same number
#' of trials.
#'
#' Choice of statistic: the **Frobenius norm** is the correct form
#' per Schmitz, Rougier, & Yzerbyt (2019, comment) and the 2020
#' erratum. The original `rcicr::computeInfoVal2IFC()` already uses it
#' (`norm(matrix(target_ci[["ci"]]), "f")`); this function matches
#' that convention for both 2IFC and Brief-RC data.
#'
#' Trial-count matching: the original
#' `rcicr::generateReferenceDistribution2IFC()` builds the reference
#' using the full pool size (`n_trials = ncol(noise_matrix)`). For
#' 2IFC this is appropriate because every producer typically responds
#' to every pool item, so per-producer trial count equals the pool
#' size. For Brief-RC the producer's *cognitive exposure* per trial is
#' richer than 2IFC (12 noisy faces vs 2), but only one chosen
#' `(stim, +/-1)` is recorded per trial, so the number of mask
#' contributions is `n_trials` and is typically smaller than
#' `n_pool`. A pool-size reference is therefore a mismatch for
#' Brief-RC and biases infoVal downward. `infoval()` closes this gap
#' by keying the reference on the actual per-producer recorded trial
#' count.
#'
#' What infoVal measures (and what it does not). The Frobenius norm
#' is a **magnitude** statistic, not a target-alignment one: it
#' answers "is the mask larger than chance?" but not "is it pointing
#' at the right pattern?". Two consequences worth flagging when
#' interpreting per-producer z-scores:
#'
#' * **Cross-paradigm comparisons need care.** Brief-RC and 2IFC are
#'   placed on the same z-scale by `infoval()` (per-producer trial-
#'   count-matched reference, identical mask construction), but the
#'   *cognitive* processes generating those masks differ. A producer
#'   who benefits from Brief-RC's richer per-trial context might
#'   produce a more accurately localized but not necessarily larger
#'   mask, and the magnitude metric will not reward that. Treat
#'   absolute Brief-RC vs 2IFC z-score differences with caution.
#' * **Stability and discriminability are different questions.**
#'   `rel_split_half()` / `rel_icc()` ask whether the signal is
#'   stable; `rel_cluster_test()` / `rel_dissimilarity()` ask whether
#'   conditions are separable. These are the right complements to a
#'   magnitude-based infoVal.
#'
#' @param signal_matrix Pixels x participants numeric matrix of raw
#'   masks (as returned by [ci_from_responses_2ifc()] or
#'   [ci_from_responses_briefrc()]).
#' @param noise_matrix Pixels x pool-size numeric matrix of noise
#'   patterns (for 2IFC: the `stimuli` object returned by
#'   `rcicr::generateStimuli2IFC(..., return_as_dataframe = TRUE)`;
#'   for Brief-RC: the same pool). Row count must match
#'   `signal_matrix`.
#' @param trial_counts Named integer vector of trial counts per
#'   producer. Names must match `colnames(signal_matrix)`.
#' @param iter Reference-distribution Monte Carlo size. Default
#'   10000 (matches rcicr convention).
#' @param mask Optional logical vector of length
#'   `nrow(signal_matrix)`. When supplied, both observed and
#'   reference norms are computed on the masked pixel subset.
#' @param cache_path Optional path to an `.rds` file. When set and
#'   the file exists, the reference distributions are loaded from
#'   it (keyed on trial count); when set but the file does not exist,
#'   the computed distributions are written there after running. This
#'   accelerates repeated runs on the same `(noise_matrix, iter, mask,
#'   seed)` configuration. The cache stores only reference norms,
#'   never observed data.
#' @param seed Optional integer; RNG state restored on exit.
#' @param progress Show a `cli` progress bar.
#' @return Object of class `rcicrely_infoval` with fields:
#' * `$infoval`, per-producer z-scores (named numeric vector).
#' * `$norms`, per-producer observed Frobenius norms.
#' * `$reference`, named list of reference norm vectors, keyed by
#'   trial count.
#' * `$ref_median`, `$ref_mad`, named numeric vectors (centre / spread)
#'   per trial count.
#' * `$trial_counts`, `$mask`, `$iter`, `$n_pool`, `$seed`, metadata.
#' @seealso [face_mask()], [ci_from_responses_2ifc()],
#'   [ci_from_responses_briefrc()].
#' @references
#' Brinkman, L., Goffin, S., van de Schoot, R., van Haren, N. E. M.,
#' Dotsch, R., & Aarts, H. (2019). Quantifying the informational
#' value of classification images. *Behavior Research Methods*,
#' 51(5), 2059-2073. \doi{10.3758/s13428-019-01232-2}
#'
#' Schmitz, M., Rougier, M., & Yzerbyt, V. (2019). Comment on
#' "Quantifying the informational value of classification images": A
#' miscomputation of the infoVal metric. *Behavior Research Methods*.
#' \doi{10.3758/s13428-019-01295-1}
#'
#' Schmitz, M., Rougier, M., Yzerbyt, V., Brinkman, L., & Dotsch, R.
#' (2020). Erratum to: Comment on "Quantifying the informational value
#' of classification images": Miscomputation of infoVal metric was a
#' minor issue and is now corrected. *Behavior Research Methods*, 52,
#' 1800-1801. \doi{10.3758/s13428-020-01367-7}
#'
#' Schmitz, M., Rougier, M., & Yzerbyt, V. (2024). Introducing the
#' brief reverse correlation: an improved tool to assess visual
#' representations. *European Journal of Social Psychology*.
#' \doi{10.1002/ejsp.3100}
#' @export
infoval <- function(signal_matrix,
                    noise_matrix,
                    trial_counts,
                    iter        = 10000L,
                    mask        = NULL,
                    cache_path  = NULL,
                    seed        = NULL,
                    progress    = TRUE) {
  if (!is.matrix(signal_matrix) || !is.numeric(signal_matrix)) {
    cli::cli_abort("{.arg signal_matrix} must be a numeric matrix.")
  }
  if (!is.matrix(noise_matrix) || !is.numeric(noise_matrix)) {
    cli::cli_abort("{.arg noise_matrix} must be a numeric matrix.")
  }
  if (nrow(signal_matrix) != nrow(noise_matrix)) {
    cli::cli_abort(c(
      "Row counts of {.arg signal_matrix} and {.arg noise_matrix} \\
       must match.",
      "*" = "signal: {nrow(signal_matrix)} pixels",
      "*" = "noise: {nrow(noise_matrix)} pixels"
    ))
  }
  n_pix  <- nrow(signal_matrix)
  n_pool <- ncol(noise_matrix)

  if (is.null(colnames(signal_matrix))) {
    colnames(signal_matrix) <- sprintf("p%03d", seq_len(ncol(signal_matrix)))
  }
  producers <- colnames(signal_matrix)

  if (is.null(names(trial_counts))) {
    cli::cli_abort(
      "{.arg trial_counts} must be a named integer vector; names \\
       must match {.code colnames(signal_matrix)}."
    )
  }
  missing_ids <- setdiff(producers, names(trial_counts))
  if (length(missing_ids) > 0L) {
    cli::cli_abort(c(
      "Missing trial counts for {length(missing_ids)} producer{?s}:",
      "*" = "{.val {missing_ids}}"
    ))
  }
  trial_counts <- as.integer(trial_counts[producers])

  if (!is.null(mask)) {
    if (!is.logical(mask) || length(mask) != n_pix) {
      cli::cli_abort(
        "{.arg mask} must be a logical vector of length {n_pix}."
      )
    }
  }

  iter <- as.integer(iter)
  if (iter < 100L) {
    cli::cli_warn(
      "{.arg iter} = {iter} is very low; MAD will be unstable."
    )
  }

  # ---- observed norms -------------------------------------------------
  norms <- vapply(
    seq_len(ncol(signal_matrix)),
    function(j) {
      x <- signal_matrix[, j]
      if (!is.null(mask)) x <- x[mask]
      sqrt(sum(x * x))
    },
    numeric(1L)
  )
  names(norms) <- producers

  # ---- reference distributions, keyed on unique trial counts ---------
  unique_n_trials <- sort(unique(trial_counts))

  cache <- NULL
  cache_hit <- FALSE
  if (!is.null(cache_path) && file.exists(cache_path)) {
    cache <- readRDS(cache_path)
    if (is.list(cache) &&
          !is.null(cache$reference) &&
          !is.null(cache$iter) &&
          cache$iter == iter &&
          setequal(names(cache$reference), as.character(unique_n_trials)) &&
          identical(cache$n_pool, n_pool) &&
          identical(cache$mask_sig, mask_signature(mask))) {
      cache_hit <- TRUE
    }
  }

  if (cache_hit) {
    reference <- cache$reference[as.character(unique_n_trials)]
  } else {
    reference <- stats::setNames(
      vector("list", length(unique_n_trials)),
      as.character(unique_n_trials)
    )
    total_sims <- iter * length(unique_n_trials)
    pid <- progress_start(total_sims, "infoval reference", show = progress)
    on.exit(progress_done(pid), add = TRUE)
    with_seed(seed, {
      for (n_t in unique_n_trials) {
        reference[[as.character(n_t)]] <-
          simulate_reference_norms(noise_matrix, n_t, iter,
                                   mask, pid)
      }
    })
    if (!is.null(cache_path)) {
      saveRDS(list(
        reference = reference,
        iter      = iter,
        n_pool    = n_pool,
        mask_sig  = mask_signature(mask)
      ), cache_path)
    }
  }

  ref_median <- vapply(reference, stats::median, numeric(1L))
  ref_mad    <- vapply(reference, stats::mad,    numeric(1L))

  # ---- per-producer z-scores -----------------------------------------
  iv <- vapply(
    seq_along(norms),
    function(j) {
      key <- as.character(trial_counts[j])
      (norms[j] - ref_median[key]) / ref_mad[key]
    },
    numeric(1L)
  )
  names(iv) <- producers

  new_rcicrely_infoval(
    infoval      = iv,
    norms        = norms,
    reference    = reference,
    ref_median   = ref_median,
    ref_mad      = ref_mad,
    trial_counts = stats::setNames(trial_counts, producers),
    mask         = mask,
    iter         = iter,
    n_pool       = n_pool,
    seed         = seed
  )
}

#' Simulate `iter` random Frobenius norms for a given trial count
#'
#' Matches Schmitz's `genMask()` construction (random sign per stim,
#' mean-by-stim collapse, divide by number of unique chosen stims,
#' Frobenius norm of the resulting mask), but samples stim ids
#' **without replacement** when the trial count fits in the pool —
#' the typical case for both 2IFC (every producer sees every pool
#' item once) and standard Brief-RC (every producer sees a random
#' subset of distinct pool items). Sampling with replacement was the
#' v0.2.0 default and produced systematically inflated reference
#' norms for 2IFC, biasing per-producer infoVal z-scores
#' downward (Study 1 replication revealed this on real data;
#' the fix landed in v0.2.1). Replacement is still used when
#' `n_trials > n_pool`, where it is the only feasible sampling
#' scheme.
#'
#' @keywords internal
#' @noRd
simulate_reference_norms <- function(noise_matrix, n_trials, iter,
                                     mask = NULL, pid = NULL) {
  n_pool <- ncol(noise_matrix)
  use_replace <- n_trials > n_pool
  norms  <- numeric(iter)
  for (k in seq_len(iter)) {
    stim <- sample.int(n_pool, n_trials, replace = use_replace)
    resp <- sample(c(-1, 1), n_trials, replace = TRUE)
    # Schmitz's genMask duplicate-stim rule: mean(response) by stim,
    # then mask = noise[, unique_stim] %*% mean_resp / n_unique_stim.
    uniq <- sort(unique(stim))
    if (length(uniq) < length(stim)) {
      idx <- match(stim, uniq)
      sums <- tabulate(idx, length(uniq))
      wts  <- as.vector(tapply(resp, idx, sum)) / sums
    } else {
      # fast path: no duplicates
      ord  <- order(stim)
      uniq <- stim[ord]
      wts  <- resp[ord]
    }
    mask_vec <- noise_matrix[, uniq, drop = FALSE] %*% wts / length(uniq)
    if (!is.null(mask)) mask_vec <- mask_vec[mask]
    norms[k] <- sqrt(sum(mask_vec * mask_vec))
    progress_tick(pid)
  }
  norms
}

#' Cheap signature of a mask vector for cache identity checking
#'
#' @keywords internal
#' @noRd
mask_signature <- function(mask) {
  if (is.null(mask)) return("none")
  # sum + length is cheap and collides only on masks with identical
  # count. Combined with iter + n_pool it's enough for a cache check.
  sprintf("%d/%d", sum(mask), length(mask))
}
