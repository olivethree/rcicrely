## S3 print / summary / plot methods for every rcicrely_* class.
##
## These are registered via roxygen @exportS3Method tags so R CMD
## check is happy even though none of the generics are imported from
## another package, they all come from `base`.
##
## Plot methods aim for publication-grade base R aesthetics: minimal
## axis chrome, soft gridlines, viridis-friendly palettes, shaded
## confidence regions where applicable. The helpers below set these
## defaults consistently across all `plot.rcicrely_*()` methods.

#' Set publication-style par() defaults. Use inside plot methods.
#' @keywords internal
#' @noRd
set_pub_par <- function(mar = c(4, 4, 2.5, 1) + 0.1) {
  graphics::par(
    mar     = mar,
    mgp     = c(2.0, 0.5, 0),
    tcl     = -0.3,
    las     = 1,
    cex.lab = 0.9,
    cex.axis= 0.8,
    cex.main= 1.0,
    font.main = 1,
    bty     = "n",
    family  = ""
  )
}

#' Soft horizontal grid behind a plot. Call after `plot()` setup but
#' before drawing the data.
#' @keywords internal
#' @noRd
add_soft_grid <- function(side = "h", at = NULL,
                          col = "grey92", lwd = 0.7) {
  if (is.null(at)) at <- pretty(graphics::par("usr")[if (side == "h") 3:4 else 1:2])
  if (side == "h") {
    graphics::abline(h = at, col = col, lwd = lwd)
  } else {
    graphics::abline(v = at, col = col, lwd = lwd)
  }
}

#' Light L-frame axes (left + bottom only).
#' @keywords internal
#' @noRd
add_clean_axes <- function(x_at = NULL, y_at = NULL,
                           x_labels = TRUE, y_labels = TRUE,
                           xlab = NULL, ylab = NULL) {
  graphics::axis(1, at = x_at, labels = x_labels,
                 col = "grey50", col.axis = "grey20", lwd = 0.7)
  graphics::axis(2, at = y_at, labels = y_labels,
                 col = "grey50", col.axis = "grey20", lwd = 0.7)
  if (!is.null(xlab)) graphics::title(xlab = xlab, col.lab = "grey20")
  if (!is.null(ylab)) graphics::title(ylab = ylab, col.lab = "grey20")
}

#' Subtle filled CI band on a 1D x-axis.
#' @keywords internal
#' @noRd
ci_shade <- function(xlo, xhi, col = "#377EB8", alpha = 0.15) {
  usr <- graphics::par("usr")
  graphics::rect(xlo, usr[3], xhi, usr[4],
                 col = grDevices::adjustcolor(col, alpha.f = alpha),
                 border = NA)
}

# ---- rcicrely_split_half ---------------------------------------------------

#' @export
print.rcicrely_split_half <- function(x, ...) {
  warn_known_regression(x)
  cat("<rcicrely split-half reliability>\n")
  cat(sprintf("  N producers:          %d\n", x$n_participants))
  cat(sprintf("  n_permutations:       %d\n", x$n_permutations))
  cat(sprintf("  mean per-split r:     %.3f  [%.3f, %.3f]\n",
              x$r_hh, x$ci_95[1], x$ci_95[2]))
  cat(sprintf("  Spearman-Brown r_SB:  %.3f  [%.3f, %.3f]\n",
              x$r_sb, x$ci_95_sb[1], x$ci_95_sb[2]))
  null <- if (is.null(x$null)) "none" else x$null
  if (!identical(null, "none")) {
    cat(sprintf("  null:                 %s\n", null))
    cat(sprintf("  null p95 (r_hh):      %.3f\n", x$r_hh_null_p95))
    cat(sprintf("  excess r_hh:          %+.3f\n", x$r_hh_excess))
    cat(sprintf("  excess r_SB:          %+.3f\n", x$r_sb_excess))
  }
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
  set_pub_par()

  h <- graphics::hist(x$distribution, breaks = 40, plot = FALSE)
  graphics::plot.new()
  graphics::plot.window(
    xlim = range(h$breaks),
    ylim = c(0, max(h$counts) * 1.08)
  )
  add_soft_grid("h")

  # Shaded 95% CI band, then the histogram bars on top
  ci_shade(x$ci_95[1], x$ci_95[2], col = "#377EB8", alpha = 0.18)
  graphics::rect(
    h$breaks[-length(h$breaks)], 0,
    h$breaks[-1L], h$counts,
    col    = "#9ECAE1",
    border = "white", lwd = 0.5
  )

  # Reference lines: r_hh (mean), r_sb (Spearman-Brown projected)
  graphics::abline(v = x$r_hh, col = "#08519C", lwd = 2)
  graphics::abline(v = x$r_sb, col = "#D62728", lwd = 2, lty = 2)

  add_clean_axes(xlab = "per-split Pearson r", ylab = "frequency")
  graphics::title(main = main, col.main = "grey10")
  graphics::mtext(
    sprintf("N = %d producers,  %d permutations",
            x$n_participants, x$n_permutations),
    side = 3, line = 0.2, cex = 0.78, col = "grey40"
  )

  graphics::legend(
    "topleft",
    legend = c(
      sprintf("r_hh = %.3f", x$r_hh),
      sprintf("r_SB  = %.3f  (Spearman-Brown)", x$r_sb),
      sprintf("95%% CI on r = [%.3f, %.3f]",
              x$ci_95[1], x$ci_95[2])
    ),
    bty = "n", cex = 0.8,
    text.col = c("#08519C", "#D62728", "grey30")
  )
  invisible(x)
}

# ---- rcicrely_loo ----------------------------------------------------------

#' @export
print.rcicrely_loo <- function(x, ...) {
  warn_known_regression(x)
  cat("<rcicrely leave-one-out influence screening>\n")
  cat(sprintf("  N producers:        %d\n", length(x$correlations)))
  cat(sprintf("  flag rule:          %s (threshold = %.2f)\n",
              x$flag_method, x$flag_threshold))

  # Lead with the informative quantity: z-scored ordering.
  cat("\n  z-scored influence (most influential first):\n")
  ordered <- order(x$z_scores)
  n_show  <- min(length(ordered), 5L)
  for (i in ordered[seq_len(n_show)]) {
    pid  <- names(x$z_scores)[i]
    flag <- if (pid %in% x$flagged) "  *FLAG*" else ""
    cat(sprintf(
      "    %-20s  z = %+6.2f  (r_loo = %.4f)%s\n",
      pid, x$z_scores[i], x$correlations[i], flag
    ))
  }
  if (length(ordered) > n_show) {
    cat(sprintf("    ... (%d more)\n", length(ordered) - n_show))
  }

  if (length(x$flagged) == 0L) {
    cat("\n  flagged producers:  none\n")
  } else {
    cat(sprintf(
      "\n  flagged producers:  %s  (%d)\n",
      paste(x$flagged, collapse = ", "),
      length(x$flagged)
    ))
  }

  cat("\n  Note: r_loo values are near 1 by construction because the\n")
  cat("  full-sample and leave-one-out means share (N-1)/N of their\n")
  cat("  data. Use z_scores (relative ordering), not r_loo levels,\n")
  cat("  to interpret influence. This is a diagnostic, not a\n")
  cat("  reliability statistic.\n")
  invisible(x)
}


#' @export
summary.rcicrely_loo <- function(object, ...) {
  print(object, ...)
  cat("\n  Full per-participant table (sorted by z_score):\n")
  for (i in seq_len(nrow(object$summary_df))) {
    row  <- object$summary_df[i, ]
    flag <- if (row$flag) "  *FLAG*" else ""
    cat(sprintf(
      "    %-20s  z = %+6.2f  (r_loo = %.4f)%s\n",
      row$participant_id, row$z_score, row$correlation, flag
    ))
  }
  invisible(object)
}

#' @export
plot.rcicrely_loo <- function(x, ...,
                              main = "Leave-one-out influence (z-scored)") {
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  set_pub_par(mar = c(6.5, 4, 2.8, 1) + 0.1)

  srt   <- sort(x$z_scores)
  flagd <- names(srt) %in% x$flagged
  ylim  <- c(min(srt, -x$flag_threshold) - 0.4,
             max(srt, 1)                  + 0.4)

  bar_x <- graphics::barplot(
    srt, las = 2, col = NA, border = NA,
    axes = FALSE, ylim = ylim,
    main = "", names.arg = rep("", length(srt))
  )
  add_soft_grid("h")

  cols <- ifelse(flagd, "#D62728", "#5B9BD5")
  graphics::rect(
    bar_x - 0.4, 0, bar_x + 0.4, srt,
    col = cols, border = NA
  )
  graphics::abline(h = 0,                  col = "grey50", lwd = 1)
  graphics::abline(h = -x$flag_threshold,  col = "#D62728",
                   lwd = 1, lty = 2)

  add_clean_axes(y_at = pretty(ylim),
                 ylab = "z-score (LOO influence)")
  graphics::title(main = main, col.main = "grey10")

  # Producer ids on x-axis, vertical
  graphics::axis(
    1, at = bar_x, labels = names(srt),
    las = 2, cex.axis = 0.65, col = "grey50",
    col.axis = "grey20", lwd = 0.7
  )

  graphics::legend(
    "topleft", bty = "n", cex = 0.78,
    legend = c(
      sprintf("flag threshold (z = -%.2f)", x$flag_threshold),
      sprintf("flagged producers: %d/%d",
              length(x$flagged), length(x$z_scores))
    ),
    text.col = c("#D62728", "grey20")
  )
  invisible(x)
}

# ---- rcicrely_icc ----------------------------------------------------------

#' @export
print.rcicrely_icc <- function(x, ...) {
  warn_known_regression(x)
  cat("<rcicrely ICC>\n")
  cat(sprintf("  model:        %s\n", x$model))
  cat(sprintf("  N targets:    %d pixels\n",       x$n_targets))
  cat(sprintf("  N raters:     %d participants\n", x$n_raters))
  pad <- function(label, value) {
    sprintf("  %-14s %.4f\n", label, value)
  }
  # Lead with ICC(3,1): single-producer reliability, resolution-comparable.
  if ("3_1" %in% x$variants) {
    cat("\n  Primary: ICC(3,1) (single-producer reliability)\n")
    cat(pad("ICC(3,1):", x$icc_3_1))
  }
  if ("3_k" %in% x$variants) {
    cat("\n  Secondary: ICC(3,k) (group-mean reliability)\n")
    cat(pad("ICC(3,k):", x$icc_3_k))
    if (isTRUE(x$n_targets > 50000L)) {
      cat("  Note: ICC(3,k) approaches 1 at large pixel counts (resolution\n")
      cat("  asymptote). For cross-resolution comparisons, report ICC(3,1).\n")
    }
  }
  shown_2 <- FALSE
  if ("2_1" %in% x$variants) {
    if (!shown_2) {
      cat("\n  Two-way-random variants (for reviewer comparability only):\n")
      shown_2 <- TRUE
    }
    cat(pad("ICC(2,1):", x$icc_2_1))
  }
  if ("2_k" %in% x$variants) {
    if (!shown_2) {
      cat("\n  Two-way-random variants (for reviewer comparability only):\n")
      shown_2 <- TRUE
    }
    cat(pad("ICC(2,k):", x$icc_2_k))
  }
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
plot.rcicrely_icc <- function(x, ..., main = "Intraclass correlation") {
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  set_pub_par(mar = c(4, 5.5, 2.8, 1) + 0.1)

  vals <- c(
    "ICC(3,1)" = x$icc_3_1,
    "ICC(3,k)" = x$icc_3_k,
    "ICC(2,1)" = x$icc_2_1,
    "ICC(2,k)" = x$icc_2_k
  )
  vals <- vals[!is.na(vals)]
  vals <- rev(vals) # so first variant ends up at top of horizontal bar

  # Horizontal layout reads better with two-line labels
  xlim <- c(min(0, min(vals)) - 0.05, max(1, max(vals)) + 0.05)
  bar_y <- graphics::barplot(
    vals, horiz = TRUE, las = 1,
    col = NA, border = NA, axes = FALSE,
    xlim = xlim, names.arg = rep("", length(vals))
  )
  add_soft_grid("v")

  pal <- c("ICC(3,1)" = "#08519C", "ICC(3,k)" = "#3182BD",
           "ICC(2,1)" = "#D9D9D9", "ICC(2,k)" = "#BDBDBD")
  cols <- pal[names(vals)]

  graphics::rect(
    pmin(0, vals), bar_y - 0.4,
    pmax(0, vals), bar_y + 0.4,
    col = cols, border = NA
  )
  graphics::abline(v = 0, col = "grey60", lwd = 0.8)

  add_clean_axes(x_at = pretty(xlim),
                 xlab = "ICC value")
  graphics::axis(
    2, at = bar_y, labels = names(vals), las = 1,
    cex.axis = 0.85, col = "grey50", col.axis = "grey20",
    lwd = 0.7
  )
  graphics::title(main = main, col.main = "grey10")
  graphics::mtext(
    sprintf("%d targets x %d raters", x$n_targets, x$n_raters),
    side = 3, line = 0.2, cex = 0.78, col = "grey40"
  )

  # Value annotations at the end of each bar
  graphics::text(
    vals, bar_y, labels = sprintf("%.3f", vals),
    pos = ifelse(vals >= 0, 4, 2), cex = 0.78,
    col = "grey20", offset = 0.3
  )
  invisible(x)
}

# ---- rcicrely_cluster_test -------------------------------------------------

#' @export
print.rcicrely_cluster_test <- function(x, ...) {
  warn_known_regression(x)
  method <- if (is.null(x$method)) "threshold" else x$method
  if (method == "tfce") {
    cat("<rcicrely cluster-based permutation test (TFCE)>\n")
    cat(sprintf(
      "  N_A = %d, N_B = %d\n",
      x$n_participants_a, x$n_participants_b
    ))
    cat(sprintf("  image dims:        %d x %d\n",
                x$img_dims[1], x$img_dims[2]))
    cat(sprintf("  TFCE parameters:   H = %.2f, E = %.2f, n_steps = %d\n",
                x$tfce_H, x$tfce_E, x$tfce_n_steps))
    cat(sprintf("  n_permutations:    %d\n",  x$n_permutations))
    cat(sprintf("  alpha:             %.2f\n", x$alpha))
    n_sig <- sum(x$tfce_significant_mask)
    cat(sprintf("  pixels below alpha: %d / %d\n",
                n_sig, length(x$tfce_pmap)))
    if (n_sig > 0L) {
      cat(sprintf("  min p-value:       %.4f\n", min(x$tfce_pmap)))
    }
    return(invisible(x))
  }
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
                                       main = NULL) {
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  method <- if (is.null(x$method)) "threshold" else x$method
  set_pub_par(mar = c(1, 1, 3, 6) + 0.1)

  pal <- grDevices::hcl.colors(256L, "RdBu", rev = TRUE)

  if (method == "tfce") {
    if (is.null(main)) main <- "TFCE map (signed, FWE-corrected)"
    tfce_mat <- matrix(x$tfce_map, x$img_dims[1L], x$img_dims[2L])
    rng <- max(abs(x$tfce_map), na.rm = TRUE)
    if (!is.finite(rng) || rng == 0) rng <- 1
    zlim <- c(-rng, rng)

    graphics::image(
      seq_len(x$img_dims[2L]), seq_len(x$img_dims[1L]),
      t(tfce_mat[nrow(tfce_mat):1L, ]),
      col       = pal,
      zlim      = zlim,
      main      = main, axes = FALSE,
      xlab      = "", ylab = "",
      asp       = x$img_dims[1L] / x$img_dims[2L],
      useRaster = TRUE
    )
    graphics::box(col = "grey80", lwd = 0.5)

    sig_mat <- matrix(x$tfce_significant_mask,
                      x$img_dims[1L], x$img_dims[2L])
    if (any(sig_mat)) {
      graphics::contour(
        seq_len(x$img_dims[2L]), seq_len(x$img_dims[1L]),
        t(sig_mat[nrow(sig_mat):1L, ]),
        levels = 0.5, add = TRUE, drawlabels = FALSE,
        col = "grey10", lwd = 1.5
      )
    }
    add_colour_bar(zlim, pal, label = "TFCE value (signed)")
    graphics::mtext(
      sprintf("alpha = %.2f, %d permutations, H = %.1f, E = %.1f",
              x$alpha, x$n_permutations, x$tfce_H, x$tfce_E),
      side = 3, line = 0.2, cex = 0.78, col = "grey40"
    )
    return(invisible(x))
  }

  if (is.null(main)) main <- "Cluster-based permutation t-map"
  tmap <- matrix(x$observed_t, x$img_dims[1L], x$img_dims[2L])
  rng  <- max(abs(x$observed_t), na.rm = TRUE)
  zlim <- c(-rng, rng)

  graphics::image(
    seq_len(x$img_dims[2L]), seq_len(x$img_dims[1L]),
    t(tmap[nrow(tmap):1L, ]),
    col       = pal,
    zlim      = zlim,
    main      = main, axes = FALSE,
    xlab      = "", ylab = "",
    asp       = x$img_dims[1L] / x$img_dims[2L],
    useRaster = TRUE
  )
  graphics::box(col = "grey80", lwd = 0.5)

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
      seq_len(x$img_dims[2L]), seq_len(x$img_dims[1L]),
      t(mask[nrow(mask):1L, ]),
      levels = 0.5, add = TRUE, drawlabels = FALSE,
      col = col, lwd = 1.5
    )
  }
  add_contour(x$pos_labels, sig_pos_ids, "grey10")
  add_contour(x$neg_labels, sig_neg_ids, "grey10")

  add_colour_bar(zlim, pal, label = "Welch t")
  n_sig <- sum(x$clusters$significant)
  graphics::mtext(
    sprintf("|t| threshold = %.1f, alpha = %.2f, %d perms, %d/%d clusters significant",
            x$cluster_threshold, x$alpha, x$n_permutations,
            n_sig, nrow(x$clusters)),
    side = 3, line = 0.2, cex = 0.78, col = "grey40"
  )
  invisible(x)
}

# ---- rcicrely_dissim -------------------------------------------------------

#' @export
print.rcicrely_dissim <- function(x, ...) {
  warn_known_regression(x)
  cat("<rcicrely representational dissimilarity>\n")
  cat(sprintf("  n_boot:               %d\n", x$n_boot))
  cat(sprintf("  CI level:             %.0f%%\n", x$ci_level * 100))
  cat(sprintf("  n_pixels:             %d\n",
              if (is.null(x$n_pixels)) NA_integer_ else x$n_pixels))
  cat("\n  Primary: Euclidean distance between group-mean CIs\n")
  cat(sprintf(
    "    Euclidean             = %.3f   [%.3f, %.3f]   SE = %.3f\n",
    x$euclidean, x$ci_dist[1], x$ci_dist[2], x$boot_se_dist
  ))
  if (!is.null(x$euclidean_normalised)) {
    cat(sprintf(
      "    Euclidean / sqrt(n)   = %.4f   (resolution-normalised)\n",
      x$euclidean_normalised
    ))
  }
  null <- if (is.null(x$null)) "none" else x$null
  if (!identical(null, "none")) {
    cat(sprintf(
      "\n  Empirical null (%s%s):\n",
      null,
      if (isTRUE(x$paired)) "; sign-flip" else "; condition-label"
    ))
    cat(sprintf("    null p95 (Euclidean)  = %.3f\n", x$d_null_p95))
    cat(sprintf("    z-equivalent          = %+.2f\n", x$d_z))
    cat(sprintf("    observed / null med   = %.2fx\n", x$d_ratio))
  }
  cat("\n  [Deprecated - will be removed in v0.3]\n")
  cat(sprintf(
    "    Pearson r             = %.3f   [%.3f, %.3f]   SE = %.3f\n",
    x$correlation, x$ci_cor[1], x$ci_cor[2], x$boot_se_cor
  ))
  cat("    Note: correlation between two base-subtracted CIs has a\n")
  cat("    positive baseline from shared image-domain structure and is\n")
  cat("    not a clean similarity score. Prefer Euclidean distance.\n")
  invisible(x)
}

#' @export
summary.rcicrely_dissim <- function(object, ...) {
  print(object, ...)
  invisible(object)
}

#' @export
plot.rcicrely_dissim <- function(x, ...,
                                 main = "Bootstrap distributions") {
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  graphics::par(mfrow = c(1L, 2L), oma = c(0, 0, 2, 0))

  hist_ci <- function(values, observed, ci, xlab, panel_title,
                      col_band) {
    set_pub_par()
    h <- graphics::hist(values, breaks = 40, plot = FALSE)
    graphics::plot.new()
    graphics::plot.window(
      xlim = range(h$breaks),
      ylim = c(0, max(h$counts) * 1.1)
    )
    add_soft_grid("h")
    ci_shade(ci[1], ci[2], col = col_band, alpha = 0.20)
    graphics::rect(
      h$breaks[-length(h$breaks)], 0,
      h$breaks[-1L], h$counts,
      col    = grDevices::adjustcolor(col_band, 0.55),
      border = "white", lwd = 0.5
    )
    graphics::abline(v = observed, col = col_band, lwd = 2)
    add_clean_axes(xlab = xlab, ylab = "frequency")
    graphics::title(main = panel_title, col.main = "grey10")
  }

  hist_ci(x$boot_dist, x$euclidean, x$ci_dist,
          xlab = "Euclidean distance",
          panel_title = "Euclidean (primary)",
          col_band = "#08519C")
  hist_ci(x$boot_cor, x$correlation, x$ci_cor,
          xlab = "Pearson r",
          panel_title = "Pearson r (deprecated)",
          col_band = "#969696")

  graphics::mtext(
    sprintf("%s   |   n_boot = %d, %.0f%% CI",
            main, x$n_boot, x$ci_level * 100),
    outer = TRUE, line = 0.4, cex = 0.95, col = "grey20"
  )
  invisible(x)
}

#' Compare Euclidean dissimilarity intervals across multiple
#' contrasts on a single plot
#'
#' @description
#' Side-by-side comparison of `rel_dissimilarity()` results: one
#' point per contrast, with its 95% percentile bootstrap CI as a
#' horizontal range. Useful for paper figures showing whether two
#' contrasts (e.g., Trust vs Friendly, Dominant vs Competent) have
#' overlapping or non-overlapping CIs without forcing the reader to
#' read four numbers from a table.
#'
#' @param ... Named `rcicrely_dissim` objects (use `name = obj`
#'   syntax). Names become the y-axis labels.
#' @param metric `"euclidean"` (default) or `"euclidean_normalised"`.
#' @param main Plot title.
#' @return Invisibly a data frame of (label, observed, ci_low,
#'   ci_high) for further use.
#' @seealso [rel_dissimilarity()]
#' @export
#' @examples
#' \dontrun{
#' d_AB <- rel_dissimilarity(sig_A, sig_B, n_boot = 500L, seed = 1L)
#' d_AC <- rel_dissimilarity(sig_A, sig_C, n_boot = 500L, seed = 1L)
#' compare_dissimilarity("A vs B" = d_AB, "A vs C" = d_AC)
#' }
compare_dissimilarity <- function(...,
                                  metric = c("euclidean",
                                             "euclidean_normalised"),
                                  main   = "Between-condition Euclidean distance") {
  metric <- match.arg(metric)
  contrasts <- list(...)
  if (length(contrasts) == 0L) {
    cli::cli_abort("Pass at least one named {.cls rcicrely_dissim} object.")
  }
  if (is.null(names(contrasts)) || any(names(contrasts) == "")) {
    cli::cli_abort("All arguments must be named (e.g. `\"A vs B\" = d_AB`).")
  }
  if (!all(vapply(contrasts, inherits, logical(1L),
                  what = "rcicrely_dissim"))) {
    cli::cli_abort("All arguments must be {.cls rcicrely_dissim} objects.")
  }

  # Pull observed value, CI bounds, and bootstrap distribution per contrast
  if (metric == "euclidean") {
    obs <- vapply(contrasts, function(d) d$euclidean, numeric(1L))
    ci  <- t(vapply(contrasts, function(d) d$ci_dist, numeric(2L)))
    boots <- lapply(contrasts, function(d) d$boot_dist)
    xlab  <- "Euclidean distance"
  } else {
    obs <- vapply(contrasts, function(d) d$euclidean_normalised, numeric(1L))
    n   <- vapply(contrasts, function(d) d$n_pixels, numeric(1L))
    ci  <- t(vapply(contrasts, function(d) {
      d$ci_dist / sqrt(d$n_pixels)
    }, numeric(2L)))
    boots <- lapply(contrasts, function(d) d$boot_dist / sqrt(d$n_pixels))
    xlab  <- "Euclidean / sqrt(n_pixels)"
  }
  labels <- names(contrasts)
  k <- length(labels)

  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  set_pub_par(mar = c(4, max(8, max(nchar(labels)) * 0.55), 2.8, 1) + 0.1)

  xlim <- range(c(0, ci, obs))
  xlim[2] <- xlim[2] * 1.05
  graphics::plot.new()
  graphics::plot.window(xlim = xlim, ylim = c(0.5, k + 0.5))
  add_soft_grid("v")

  cols <- grDevices::hcl.colors(max(k, 2L), "Dark 3")[seq_len(k)]
  ys <- seq_len(k)

  # Density "violin" lookalike: rescaled bootstrap density per row
  for (i in seq_len(k)) {
    dens <- stats::density(boots[[i]])
    sc   <- 0.35 / max(dens$y)
    graphics::polygon(
      c(dens$x, rev(dens$x)),
      c(ys[i] + dens$y * sc, rev(ys[i] - dens$y * sc)),
      col    = grDevices::adjustcolor(cols[i], 0.18),
      border = NA
    )
  }

  # CI bars + observed points
  graphics::segments(
    ci[, 1], ys, ci[, 2], ys,
    col = cols, lwd = 4, lend = "butt"
  )
  graphics::points(obs, ys, pch = 19, col = cols, cex = 1.4)
  graphics::points(obs, ys, pch = 21, col = "white", bg = cols,
                   cex = 1.4, lwd = 0.8)

  add_clean_axes(x_at = pretty(xlim), xlab = xlab)
  graphics::axis(2, at = ys, labels = labels, las = 1,
                 cex.axis = 0.85, col = "grey50",
                 col.axis = "grey20", lwd = 0.7)
  graphics::title(main = main, col.main = "grey10")
  graphics::mtext(
    sprintf("Points = observed; bars = 95%% percentile CI; shaded = bootstrap density"),
    side = 3, line = 0.2, cex = 0.75, col = "grey40"
  )

  invisible(data.frame(
    label    = labels,
    observed = obs,
    ci_low   = ci[, 1],
    ci_high  = ci[, 2],
    stringsAsFactors = FALSE
  ))
}

# ---- rcicrely_report -------------------------------------------------------

#' @export
print.rcicrely_pairwise_report <- function(x, ...) {
  warn_known_regression(x)
  cat(sprintf(
    "<rcicrely pairwise report: %d conditions, %d pairs, fwer = %s>\n",
    length(x$conditions), nrow(x$pairs), x$fwer
  ))
  cat(sprintf("  alpha (across pairs): %.3f\n", x$alpha))
  cat(sprintf("  cluster method:       %s\n", x$method))
  if (isTRUE(x$paired)) cat("  design:               paired (within-subjects)\n")
  cat("\n")
  pr <- x$pairs
  pr$pair <- sprintf("%s vs %s", pr$cond_a, pr$cond_b)
  pr$p_min       <- formatC(pr$p_min, format = "f", digits = 4)
  pr$p_adj_pair  <- formatC(pr$p_adj_pair, format = "f", digits = 4)
  pr$euclidean   <- formatC(pr$euclidean, format = "f", digits = 3)
  pr$euclidean_normalised <-
    formatC(pr$euclidean_normalised, format = "f", digits = 4)
  pr$sig <- ifelse(pr$significant, "*", "")
  print(pr[, c("pair", "n_clusters", "p_min", "p_adj_pair", "sig",
               "euclidean", "euclidean_normalised")],
        row.names = FALSE)
  cat(sprintf(
    "\n  p_adj_pair: minimum within-pair cluster p, %s-adjusted across %d pairs.\n",
    x$fwer, nrow(x$pairs)
  ))
  invisible(x)
}

#' @export
print.rcicrely_report <- function(x, ...) {
  warn_known_regression(x)
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

# ---- rcicrely_infoval ------------------------------------------------------

#' @export
print.rcicrely_infoval <- function(x, ...) {
  warn_known_regression(x)
  cat("<rcicrely informational value>\n")
  cat(sprintf("  producers:          %d\n", length(x$infoval)))
  cat(sprintf("  unique trial counts: %s\n",
              paste(sort(unique(x$trial_counts)), collapse = ", ")))
  cat(sprintf("  noise-matrix pool:  %d\n", x$n_pool))
  cat(sprintf("  iter (per n_trials): %d\n", x$iter))
  cat(sprintf("  mask:               %s\n",
              if (is.null(x$mask)) "none"
              else sprintf("logical vector, %d pixels inside",
                           sum(x$mask))))
  n_sig <- sum(x$infoval > 1.96)
  cat(sprintf("  producers above z = 1.96: %d / %d\n",
              n_sig, length(x$infoval)))
  cat("\n  per-producer z (top 5 and bottom 5):\n")
  ord <- order(-x$infoval)
  show_ix <- if (length(ord) > 10L) c(ord[1:5], ord[(length(ord) - 4L):length(ord)])
             else ord
  for (i in show_ix) {
    cat(sprintf("    %-20s  z = %+6.2f  (norm = %.4f, n_trials = %d)\n",
                names(x$infoval)[i], x$infoval[i],
                x$norms[i], x$trial_counts[i]))
  }
  invisible(x)
}

#' @export
summary.rcicrely_infoval <- function(object, ...) {
  print(object, ...)
  cat("\n  reference median / MAD by trial count:\n")
  for (nm in names(object$ref_median)) {
    cat(sprintf("    n_trials = %s:   median = %.4f, MAD = %.4f\n",
                nm, object$ref_median[nm], object$ref_mad[nm]))
  }
  invisible(object)
}

#' @export
plot.rcicrely_infoval <- function(x, ...,
                                  main = "Informational value (per-producer z)") {
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  set_pub_par(mar = c(6.5, 4, 2.8, 1) + 0.1)

  srt  <- sort(x$infoval)
  ylim <- c(min(srt, -1) - 0.5, max(srt, 2.5) + 0.5)

  bar_x <- graphics::barplot(
    srt, las = 2, col = NA, border = NA, axes = FALSE,
    ylim = ylim, main = "", names.arg = rep("", length(srt))
  )
  add_soft_grid("h")

  cols <- ifelse(
    srt >  1.96, "#08519C",
    ifelse(srt >  0, "#9ECAE1",
           ifelse(srt > -1.96, "#FDBB84", "#D62728"))
  )
  graphics::rect(
    bar_x - 0.4, 0, bar_x + 0.4, srt,
    col = cols, border = NA
  )
  graphics::abline(h = 0,    col = "grey50", lwd = 1)
  graphics::abline(h = 1.96, col = "#D62728", lwd = 1, lty = 2)
  graphics::abline(h = -1.96, col = "#D62728", lwd = 1, lty = 2)

  add_clean_axes(y_at = pretty(ylim),
                 ylab = "infoVal z-score")
  graphics::axis(
    1, at = bar_x, labels = names(srt),
    las = 2, cex.axis = 0.65, col = "grey50",
    col.axis = "grey20", lwd = 0.7
  )
  graphics::title(main = main, col.main = "grey10")

  n_above <- sum(x$infoval > 1.96)
  graphics::mtext(
    sprintf("median z = %+.2f, %d/%d producers above z = 1.96",
            stats::median(x$infoval), n_above,
            length(x$infoval)),
    side = 3, line = 0.2, cex = 0.78, col = "grey40"
  )
  graphics::legend(
    "topleft", bty = "n", cex = 0.78,
    legend = c("z > 1.96 (clear signal)",
               "z > 0",
               "z < 0 (low-rank / idiosyncratic)",
               "z < -1.96"),
    fill = c("#08519C", "#9ECAE1", "#FDBB84", "#D62728"),
    border = NA
  )
  invisible(x)
}
