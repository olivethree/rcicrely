# 04_generate_oliveira_2019_examples.R
#
# Pre-renders two figures from Oliveira, Garcia-Marques, Dotsch, &
# Garcia-Marques (2019), Study 1, embedded by `vignettes/tutorial.Rmd`:
#
#   1. cluster_dominant_vs_competent.png      (full image)
#   2. cluster_trust_vs_dominant_eyes_mouth.png (region-restricted)
#
# Heavy: noise reconstruction + cluster permutations on 65,536 pixels
# x 1,000 permutations. Run once and commit the resulting PNGs under
# `inst/extdata/oliveira_2019_examples/`.
#
# Requires `temp/Study 1/` (developer-only). Not part of the
# user-facing build; data-raw/ is .Rbuildignore-d.

library(rcicrely)

study_dir <- file.path("temp", "Study 1")
paths <- list(
  responses = file.path(study_dir, "data",        "study1data.csv"),
  rdata     = file.path(study_dir, "R scripts",
                        "rcic_seed_1_time_fev_05_2015_03_17.Rdata")
)
stopifnot(all(file.exists(unlist(paths))))

out_dir <- file.path("inst", "extdata", "oliveira_2019_examples")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# -------- Load -----------------------------------------------------
responses <- read.csv2(paths$responses, stringsAsFactors = FALSE)
rdata_env <- new.env()
load(paths$rdata, envir = rdata_env)

# -------- Reconstruct noise matrix from legacy rcicr 0.3.0 format --
build_noise_matrix_legacy <- function(env, base = "male") {
  s        <- env$s
  params   <- env$stimuli_params[[base]]
  n_pix    <- prod(dim(s$sinusoids)[1:2])
  n_layers <- dim(s$sinusoids)[3]
  n_trials <- nrow(params)
  sin_flat <- matrix(s$sinusoids, n_pix, n_layers)
  idx_flat <- matrix(s$sinIdx,    n_pix, n_layers)
  noise <- matrix(0, n_pix, n_trials)
  for (i in seq_len(n_trials)) {
    noise[, i] <- rowSums(sin_flat * params[i, idx_flat])
  }
  attr(noise, "img_dims") <- c(env$img_size, env$img_size)
  noise
}

message("Reconstructing 65,536 x 300 noise matrix ...")
noise_matrix <- build_noise_matrix_legacy(rdata_env, "male")

# Base face: 256 x 256 grayscale matrix in [0, 1]
base_face <- rdata_env$base_faces[["male"]]
stopifnot(identical(dim(base_face), c(256L, 256L)))

# -------- Build per-trait signal matrices --------------------------
build_signal_matrix <- function(noise_matrix, responses_subset) {
  subjects <- unique(responses_subset$subject)
  out <- matrix(NA_real_, nrow = nrow(noise_matrix), ncol = length(subjects))
  colnames(out) <- as.character(subjects)
  for (j in seq_along(subjects)) {
    sub <- responses_subset[responses_subset$subject == subjects[j], , drop = FALSE]
    sub <- sub[order(sub$trial), ]
    out[, j] <- (noise_matrix[, sub$stimulus] %*% sub$response) / nrow(sub)
  }
  attr(out, "img_dims") <- attr(noise_matrix, "img_dims")
  attr(out, "source")   <- "raw"
  out
}

message("Building per-trait signal matrices ...")
sm_dominant  <- build_signal_matrix(noise_matrix,
                                    responses[responses$trait == "DOMINANT",  ])
sm_competent <- build_signal_matrix(noise_matrix,
                                    responses[responses$trait == "COMPETENT", ])
sm_trust     <- build_signal_matrix(noise_matrix,
                                    responses[responses$trait == "TRUST",     ])

img_dims <- c(256L, 256L)

# -------- Overlay helper: t-map (signed) + base face + sig contours
#
# `plot_ci_overlay()` ships in the package but its `test` argument
# accepts only `rcicrely_agreement_map_test` results. Cluster-test
# results carry significance via `pos_labels` / `neg_labels` plus
# `$clusters$significant`. This helper renders the same composite
# (translucent diverging heatmap over the grayscale base) and then
# draws contours around the significant clusters.
plot_cluster_overlay <- function(clust, base, main = NULL,
                                 alpha_max = 0.7) {
  img_dims <- as.integer(clust$img_dims)
  t_vec    <- clust$observed_t

  rng <- max(abs(t_vec), na.rm = TRUE)
  if (!is.finite(rng) || rng == 0) rng <- 1
  norm_vec <- pmax(-1, pmin(1, t_vec / rng))
  norm_mat <- matrix(norm_vec, img_dims[1], img_dims[2])
  alpha    <- abs(norm_mat) * alpha_max

  pos  <- norm_mat > 0
  neg  <- norm_mat < 0
  fg_r <- ifelse(pos, 0.85, ifelse(neg, 0.10, 0))
  fg_g <- ifelse(pos, 0.10, ifelse(neg, 0.20, 0))
  fg_b <- ifelse(pos, 0.10, ifelse(neg, 0.85, 0))

  composed <- array(0, dim = c(img_dims[1], img_dims[2], 3L))
  composed[, , 1] <- (1 - alpha) * base + alpha * fg_r
  composed[, , 2] <- (1 - alpha) * base + alpha * fg_g
  composed[, , 3] <- (1 - alpha) * base + alpha * fg_b
  composed[composed < 0] <- 0
  composed[composed > 1] <- 1

  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  graphics::par(mar = c(0.5, 0.5, if (is.null(main)) 0.5 else 2.8, 0.5))
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, img_dims[2]),
                        ylim = c(0, img_dims[1]),
                        asp = 1, xaxs = "i", yaxs = "i")
  graphics::rasterImage(composed, 0, 0, img_dims[2], img_dims[1],
                        interpolate = FALSE)

  # Significant-cluster contours
  sig_mat <- matrix(0L, img_dims[1], img_dims[2])
  if (NROW(clust$clusters) > 0L) {
    for (i in seq_len(NROW(clust$clusters))) {
      if (!isTRUE(clust$clusters$significant[i])) next
      cid <- clust$clusters$cluster_id[i]
      dir <- clust$clusters$direction[i]
      lbl <- if (dir == "pos") clust$pos_labels else clust$neg_labels
      sig_mat[lbl == cid] <- 1L
    }
  }
  if (any(sig_mat > 0)) {
    graphics::contour(
      x = seq_len(img_dims[2]) - 0.5,
      y = seq_len(img_dims[1]) - 0.5,
      z = t(sig_mat[nrow(sig_mat):1L, , drop = FALSE]),
      levels = 0.5, drawlabels = FALSE, add = TRUE,
      col = "black", lwd = 1.4
    )
  }
  if (!is.null(main)) {
    graphics::title(main = main, line = 1, cex.main = 1.0,
                    font.main = 1)
  }
  invisible(composed)
}

# -------- Figure A: Dominant vs Competent (full image) -------------
message("Cluster test: Dominant vs Competent (full image) ...")
clust_DC <- rel_cluster_test(
  sm_dominant, sm_competent,
  img_dims          = img_dims,
  cluster_threshold = 2.0,
  n_permutations    = 1000L,
  alpha             = 0.05,
  seed              = 1L,
  progress          = FALSE
)

png_A <- file.path(out_dir, "cluster_dominant_vs_competent.png")
grDevices::png(png_A, width = 720, height = 720, res = 96)
plot_cluster_overlay(
  clust_DC, base_face,
  main = "Dominant vs Competent (Oliveira et al., 2019, Study 1)"
)
grDevices::dev.off()
message("  -> ", png_A)

# -------- Figure B: Trust vs Dominant, region-restricted -----------
zero_fill <- function(sm, keep) {
  out <- sm
  out[!keep, ] <- 0
  attr(out, "img_dims") <- attr(sm, "img_dims")
  attr(out, "source")   <- attr(sm, "source")
  out
}

m_eyes  <- face_mask(img_dims, region = "eyes")
m_mouth <- face_mask(img_dims, region = "mouth")

message("Cluster test: Trust vs Dominant, eyes ...")
clust_eyes <- rel_cluster_test(
  zero_fill(sm_trust,    m_eyes),
  zero_fill(sm_dominant, m_eyes),
  img_dims          = img_dims,
  cluster_threshold = 2.0,
  n_permutations    = 1000L,
  alpha             = 0.05,
  seed              = 2L,
  progress          = FALSE
)

message("Cluster test: Trust vs Dominant, mouth ...")
clust_mouth <- rel_cluster_test(
  zero_fill(sm_trust,    m_mouth),
  zero_fill(sm_dominant, m_mouth),
  img_dims          = img_dims,
  cluster_threshold = 2.0,
  n_permutations    = 1000L,
  alpha             = 0.05,
  seed              = 3L,
  progress          = FALSE
)

# `plot.rcicrely_cluster_test()` does `par(no.readonly = TRUE)` on
# entry and restores on exit; that clobbers the `mfg` panel pointer
# under `mfrow`, so two consecutive plots land in the same panel.
# Workaround: render each panel to its own PNG, then composite via
# rasterImage on a fresh device.
tmp_eyes  <- tempfile(fileext = ".png")
tmp_mouth <- tempfile(fileext = ".png")

grDevices::png(tmp_eyes, width = 720, height = 720, res = 96)
plot_cluster_overlay(clust_eyes, base_face,
                     main = "Trust vs Dominant: eyes")
grDevices::dev.off()

grDevices::png(tmp_mouth, width = 720, height = 720, res = 96)
plot_cluster_overlay(clust_mouth, base_face,
                     main = "Trust vs Dominant: mouth")
grDevices::dev.off()

img_eyes  <- png::readPNG(tmp_eyes)
img_mouth <- png::readPNG(tmp_mouth)

png_B <- file.path(out_dir, "cluster_trust_vs_dominant_eyes_mouth.png")
grDevices::png(png_B, width = 1440, height = 720, res = 96)
op <- graphics::par(mfrow = c(1, 2), mar = c(0, 0, 0, 0))
graphics::plot.new()
graphics::rasterImage(img_eyes,  0, 0, 1, 1)
graphics::plot.new()
graphics::rasterImage(img_mouth, 0, 0, 1, 1)
graphics::par(op)
grDevices::dev.off()
unlink(c(tmp_eyes, tmp_mouth))
message("  -> ", png_B)

message("\nDone. Output files:")
print(list.files(out_dir, full.names = TRUE))
