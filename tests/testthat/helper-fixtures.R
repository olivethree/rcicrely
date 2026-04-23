## Shared fixtures for rcicrely tests. Kept small (32 x 32 = 1024 pixels,
## 10–12 producers) so every test file runs in well under a second.

make_sig <- function(n_pix = 1024L, n_p = 12L,
                     signal_strength = 0.2, seed = 1L) {
  set.seed(seed)
  # half the pixels get the signal; the other half doesn't
  out <- matrix(rnorm(n_pix * n_p, sd = 0.05), n_pix, n_p)
  signal_mask <- c(rep(1, n_pix / 2L), rep(0, n_pix / 2L))
  out <- out + signal_strength * outer(signal_mask,
                                       stats::runif(n_p, 0.7, 1.3))
  attr(out, "img_dims") <- c(as.integer(sqrt(n_pix)),
                             as.integer(sqrt(n_pix)))
  out
}

make_sig_pair <- function(n_pix = 1024L, n_p = 12L, seed = 1L) {
  # A emphasises the top half, B emphasises the bottom half
  side <- as.integer(sqrt(n_pix))
  top_mask <- c(rep(1, n_pix / 2L), rep(0, n_pix / 2L))
  bot_mask <- rev(top_mask)

  set.seed(seed)
  noise_a <- matrix(rnorm(n_pix * n_p, sd = 0.05), n_pix, n_p)
  a <- noise_a + 0.2 * outer(top_mask, stats::runif(n_p, 0.7, 1.3))
  attr(a, "img_dims") <- c(side, side)

  noise_b <- matrix(rnorm(n_pix * n_p, sd = 0.05), n_pix, n_p)
  b <- noise_b + 0.2 * outer(bot_mask, stats::runif(n_p, 0.7, 1.3))
  attr(b, "img_dims") <- c(side, side)

  list(a = a, b = b, img_dims = c(side, side))
}

write_tiny_png <- function(path, data) {
  if (!requireNamespace("png", quietly = TRUE)) {
    testthat::skip("png not installed")
  }
  png::writePNG(data, path)
  path
}
