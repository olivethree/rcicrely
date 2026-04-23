## S3 print / summary / plot methods for every rcicrely_* class.
##
## These are registered via roxygen @exportS3Method tags so R CMD
## check is happy even though none of the generics are imported from
## another package  --  they all come from `base`.

# ---- rcicrely_split_half ---------------------------------------------------

#' @export
print.rcicrely_split_half <- function(x, ...) {
  cat("<rcicrely split-half reliability>\n")
  cat(sprintf("  N producers:          %d\n", x$n_participants))
  cat(sprintf("  n_permutations:       %d\n", x$n_permutations))
  cat(sprintf("  mean per-split r:     %.3f  [%.3f, %.3f]\n",
              x$r_hh, x$ci_95[1], x$ci_95[2]))
  cat(sprintf("  Spearman-Brown r_SB:  %.3f  [%.3f, %.3f]\n",
              x$r_sb, x$ci_95_sb[1], x$ci_95_sb[2]))
  invisible(x)
}

#' @export
summary.rcicrely_split_half <- function(object, ...) {
  print(object, ...)
  cat("  (95% CI via percentile method on permutation distribution.)\n")
  invisible(object)
}

#' @export
plot.rcicrely_split_half <- function(x, ...,
                                     main = "Split-half reliability") {
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  graphics::hist(
    x$distribution, breaks = 40,
    col = "grey80", border = "white",
    main = main,
    xlab = "per-split Pearson r"
  )
  graphics::abline(v = x$r_hh,   col = "red",  lwd = 2)
  graphics::abline(v = x$ci_95,  col = "red",  lwd = 1, lty = 2)
  graphics::legend(
    "topleft",
    legend = c(
      sprintf("r_hh = %.3f", x$r_hh),
      sprintf("95%% CI = [%.3f, %.3f]", x$ci_95[1], x$ci_95[2])
    ),
    bty = "n", text.col = "red"
  )
  invisible(x)
}

# ---- rcicrely_loo ----------------------------------------------------------

#' @export
print.rcicrely_loo <- function(x, ...) {
  cat("<rcicrely leave-one-out sensitivity>\n")
  cat(sprintf("  N producers:        %d\n", length(x$correlations)))
  cat(sprintf("  mean LOO r:         %.3f\n", x$mean_r))
  cat(sprintf("  SD:                 %.3f\n", x$sd_r))
  cat(sprintf("  flag threshold:     r < %.3f  (mean - %.2f SD)\n",
              x$threshold,
              (x$mean_r - x$threshold) / x$sd_r))
  if (length(x$flagged) == 0L) {
    cat("  flagged producers:  none\n")
  } else {
    cat(sprintf(
      "  flagged producers:  %s  (%d)\n",
      paste(x$flagged, collapse = ", "),
      length(x$flagged)
    ))
  }
  invisible(x)
}

#' @export
summary.rcicrely_loo <- function(object, ...) {
  print(object, ...)
  cat("\n  Per-participant LOO correlation (sorted):\n")
  srt <- sort(object$correlations)
  for (nm in names(srt)) {
    flag <- if (nm %in% object$flagged) "  *FLAG*" else ""
    cat(sprintf("    %-20s %.4f%s\n", nm, srt[nm], flag))
  }
  invisible(object)
}

#' @export
plot.rcicrely_loo <- function(x, ..., main = "Leave-one-out sensitivity") {
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  srt <- sort(x$correlations)
  graphics::par(mar = c(8, 4, 4, 2) + 0.1)
  cols <- ifelse(names(srt) %in% x$flagged, "red", "grey40")
  graphics::barplot(
    srt, las = 2, col = cols, border = NA,
    main = main, ylab = "Correlation with full-sample group CI",
    ylim = c(min(srt, x$threshold) - 0.01, 1)
  )
  graphics::abline(h = x$mean_r,    col = "blue", lty = 2)
  graphics::abline(h = x$threshold, col = "red",  lty = 2)
  graphics::legend(
    "bottomright", bty = "n",
    legend = c(sprintf("mean = %.3f", x$mean_r),
               sprintf("threshold = %.3f", x$threshold)),
    text.col = c("blue", "red")
  )
  invisible(x)
}

# ---- rcicrely_icc ----------------------------------------------------------

#' @export
print.rcicrely_icc <- function(x, ...) {
  cat("<rcicrely ICC>\n")
  cat(sprintf("  model:        %s\n", x$model))
  cat(sprintf("  N targets:    %d pixels\n",       x$n_targets))
  cat(sprintf("  N raters:     %d participants\n", x$n_raters))
  cat("\n")
  pad <- function(label, value) {
    sprintf("  %-14s %.4f\n", label, value)
  }
  if ("3_1" %in% x$variants) cat(pad("ICC(3,1):",  x$icc_3_1))
  if ("3_k" %in% x$variants) cat(pad("ICC(3,k):",  x$icc_3_k))
  if ("2_1" %in% x$variants) cat(pad("ICC(2,1):",  x$icc_2_1))
  if ("2_k" %in% x$variants) cat(pad("ICC(2,k):",  x$icc_2_k))
  invisible(x)
}

#' @export
summary.rcicrely_icc <- function(object, ...) {
  print(object, ...)
  cat(sprintf("\n  MS_rows (targets)  = %.5g\n", object$ms_rows))
  cat(sprintf("  MS_cols (raters)   = %.5g\n",    object$ms_cols))
  cat(sprintf("  MS_error           = %.5g\n",    object$ms_error))
  invisible(object)
}

#' @export
plot.rcicrely_icc <- function(x, ..., main = "ICC") {
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  vals <- c(
    "ICC(3,1)" = x$icc_3_1,
    "ICC(3,k)" = x$icc_3_k,
    "ICC(2,1)" = x$icc_2_1,
    "ICC(2,k)" = x$icc_2_k
  )
  vals <- vals[!is.na(vals)]
  graphics::par(mar = c(5, 5, 4, 2) + 0.1)
  graphics::barplot(
    vals, col = "steelblue", border = NA,
    main = main, ylab = "ICC value",
    ylim = c(min(0, min(vals)), max(1, max(vals)))
  )
  graphics::abline(h = 0, col = "grey50")
  invisible(x)
}

# ---- rcicrely_cluster_test -------------------------------------------------

#' @export
print.rcicrely_cluster_test <- function(x, ...) {
  cat("<rcicrely cluster-based permutation test>\n")
  cat(sprintf(
    "  N_A = %d, N_B = %d\n",
    x$n_participants_a, x$n_participants_b
  ))
  cat(sprintf("  image dims:        %d x %d\n",
              x$img_dims[1], x$img_dims[2]))
  cat(sprintf("  cluster threshold: |t| > %.2f\n", x$cluster_threshold))
  cat(sprintf("  n_permutations:    %d\n",  x$n_permutations))
  cat(sprintf("  alpha:             %.2f\n", x$alpha))
  n_sig <- sum(x$clusters$significant)
  cat(sprintf("  clusters found:    %d (%d significant)\n",
              nrow(x$clusters), n_sig))
  if (nrow(x$clusters) > 0L) {
    cat("\n")
    print(x$clusters, row.names = FALSE, digits = 3)
  }
  invisible(x)
}

#' @export
summary.rcicrely_cluster_test <- function(object, ...) {
  print(object, ...)
  invisible(object)
}

#' @export
plot.rcicrely_cluster_test <- function(x, ...,
                                       main = "Cluster t-map") {
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  tmap <- matrix(x$observed_t, x$img_dims[1L], x$img_dims[2L])
  # symmetric colour scale
  rng <- max(abs(x$observed_t), na.rm = TRUE)
  graphics::image(
    t(tmap[nrow(tmap):1L, ]),
    col = grDevices::hcl.colors(256L, "RdBu", rev = TRUE),
    zlim = c(-rng, rng),
    main = main, axes = FALSE,
    useRaster = TRUE
  )
  # overlay significant clusters as contours
  sig_pos_ids <- x$clusters$cluster_id[
    x$clusters$direction == "pos" & x$clusters$significant
  ]
  sig_neg_ids <- x$clusters$cluster_id[
    x$clusters$direction == "neg" & x$clusters$significant
  ]
  add_contour <- function(labels, ids, col) {
    if (length(ids) == 0L) return(invisible())
    mask <- matrix(labels %in% ids, x$img_dims[1L], x$img_dims[2L])
    graphics::contour(
      t(mask[nrow(mask):1L, ]),
      levels = 0.5, add = TRUE, drawlabels = FALSE,
      col = col, lwd = 2
    )
  }
  add_contour(x$pos_labels, sig_pos_ids, "black")
  add_contour(x$neg_labels, sig_neg_ids, "black")
  invisible(x)
}

# ---- rcicrely_dissim -------------------------------------------------------

#' @export
print.rcicrely_dissim <- function(x, ...) {
  cat("<rcicrely representational dissimilarity>\n")
  cat(sprintf("  n_boot:     %d\n", x$n_boot))
  cat(sprintf("  CI level:   %.0f%%\n", x$ci_level * 100))
  cat(sprintf(
    "  Pearson r   = %.3f   [%.3f, %.3f]   SE = %.3f\n",
    x$correlation, x$ci_cor[1], x$ci_cor[2], x$boot_se_cor
  ))
  cat(sprintf(
    "  Euclidean   = %.3f   [%.3f, %.3f]   SE = %.3f\n",
    x$euclidean, x$ci_dist[1], x$ci_dist[2], x$boot_se_dist
  ))
  invisible(x)
}

#' @export
summary.rcicrely_dissim <- function(object, ...) {
  print(object, ...)
  invisible(object)
}

#' @export
plot.rcicrely_dissim <- function(x, ..., main = "Bootstrap distributions") {
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  graphics::par(mfrow = c(1L, 2L))

  graphics::hist(x$boot_cor, breaks = 40, col = "grey80",
                 border = "white",
                 main = "Pearson r", xlab = "r")
  graphics::abline(v = x$correlation, col = "red",  lwd = 2)
  graphics::abline(v = x$ci_cor,      col = "red",  lty = 2)

  graphics::hist(x$boot_dist, breaks = 40, col = "grey80",
                 border = "white",
                 main = "Euclidean distance", xlab = "distance")
  graphics::abline(v = x$euclidean, col = "red",  lwd = 2)
  graphics::abline(v = x$ci_dist,   col = "red",  lty = 2)

  graphics::mtext(main, outer = TRUE, line = -1.5, cex = 1.1)
  invisible(x)
}

# ---- rcicrely_report -------------------------------------------------------

#' @export
print.rcicrely_report <- function(x, ...) {
  cat(sprintf("<rcicrely report: %s>\n", x$method))
  for (nm in names(x$results)) {
    cat("---\n")
    cat(sprintf("* $results$%s\n", nm))
    print(x$results[[nm]])
  }
  invisible(x)
}

#' @export
summary.rcicrely_report <- function(object, ...) {
  cat(sprintf("<rcicrely report: %s>\n", object$method))
  for (nm in names(object$results)) {
    cat(sprintf("\n-- $results$%s --\n", nm))
    summary(object$results[[nm]], ...)
  }
  invisible(object)
}

#' @export
plot.rcicrely_report <- function(x, ...) {
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  n <- length(x$results)
  # one row per result; use layout() for clean stacking
  graphics::par(mfrow = c(ceiling(n / 2L), 2L))
  for (nm in names(x$results)) {
    plot(x$results[[nm]], main = nm)
  }
  invisible(x)
}
