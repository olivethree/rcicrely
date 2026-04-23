## Internal helpers for cluster-based permutation testing.
## Pure R, no compiled code. Connected-component labelling via
## breadth-first search with 4-connectivity (CLAUDE.md sec.7.5  --  the more
## conservative choice than 8-connectivity).

#' Label connected components in a 2D boolean mask, 4-connectivity
#'
#' Returns an integer matrix the same shape as `mask`. 0 = background
#' (cell was `FALSE`); 1, 2, ... = component ids in the order they were
#' first encountered during a row-major scan.
#'
#' Classic BFS. Queue is a preallocated integer vector with head/tail
#' pointers so we don't pay for `c()` growth.
#'
#' @param mask Logical `nrow x ncol` matrix.
#' @return Integer matrix, same shape.
#' @keywords internal
#' @noRd
label_components_4conn <- function(mask) {
  stopifnot(is.matrix(mask), is.logical(mask))
  nr <- nrow(mask)
  nc <- ncol(mask)
  labels <- matrix(0L, nrow = nr, ncol = nc)

  # preallocate queue large enough to hold every masked pixel
  max_queue <- sum(mask)
  if (max_queue == 0L) return(labels)

  qrow <- integer(max_queue)
  qcol <- integer(max_queue)

  next_label <- 0L
  for (c0 in seq_len(nc)) {
    for (r0 in seq_len(nr)) {
      if (!mask[r0, c0] || labels[r0, c0] != 0L) next
      next_label <- next_label + 1L
      # seed
      head <- 1L
      tail <- 1L
      qrow[1L] <- r0
      qcol[1L] <- c0
      labels[r0, c0] <- next_label
      while (head <= tail) {
        r <- qrow[head]
        cc <- qcol[head]
        head <- head + 1L
        # 4 neighbours
        for (k in seq_len(4L)) {
          nr_n <- r + switch(k, -1L, 1L, 0L, 0L)
          nc_n <- cc + switch(k, 0L, 0L, -1L, 1L)
          if (nr_n < 1L || nr_n > nr || nc_n < 1L || nc_n > nc) next
          if (!mask[nr_n, nc_n] || labels[nr_n, nc_n] != 0L) next
          labels[nr_n, nc_n] <- next_label
          tail <- tail + 1L
          qrow[tail] <- nr_n
          qcol[tail] <- nc_n
        }
      }
    }
  }
  labels
}


#' Given a t-map and a threshold, find clusters (positive and negative)
#' and compute their masses.
#'
#' @param t_vec Numeric vector of t values.
#' @param img_dims Integer `c(nrow, ncol)`.
#' @param threshold Cluster-forming threshold (positive scalar).
#' @return List with:
#'   * `pos_labels`, `neg_labels`  --  integer matrices (`img_dims`).
#'   * `pos_masses`, `neg_masses`  --  numeric vectors of per-cluster
#'     masses. Positive masses are sums of t-values inside positive
#'     clusters (so always > threshold x cluster_size). Negative
#'     masses are sums (so always < -threshold x cluster_size).
#' @keywords internal
#' @noRd
find_clusters <- function(t_vec, img_dims, threshold) {
  tmat <- matrix(t_vec, nrow = img_dims[1L], ncol = img_dims[2L])
  pos_mask <- tmat > threshold
  neg_mask <- tmat < -threshold

  pos_labels <- label_components_4conn(pos_mask)
  neg_labels <- label_components_4conn(neg_mask)

  pos_ids <- seq_len(max(pos_labels))
  pos_masses <- vapply(
    pos_ids,
    function(i) sum(tmat[pos_labels == i]),
    numeric(1L)
  )
  if (length(pos_ids) == 0L) pos_masses <- numeric(0L)

  neg_ids <- seq_len(max(neg_labels))
  neg_masses <- vapply(
    neg_ids,
    function(i) sum(tmat[neg_labels == i]),
    numeric(1L)
  )
  if (length(neg_ids) == 0L) neg_masses <- numeric(0L)

  list(
    pos_labels = pos_labels,
    neg_labels = neg_labels,
    pos_masses = pos_masses,
    neg_masses = neg_masses
  )
}
