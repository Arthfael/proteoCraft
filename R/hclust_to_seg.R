#' hclust_to_seg
#'
#' @description 
#' Great function built with the help of chatGPT to extract segments and labels from an hclust efficiently, for the purpose of plotting them in a heatmap.
#' 
#' @param hc A hclust object.
#'
#' @export

hclust_to_seg <- function(hc) {
  stopifnot(inherits(hc, "hclust"))
  #
  n <- length(hc$order)
  merge <- hc$merge
  height <- hc$height
  n_merge <- nrow(merge)
  n_nodes <- n + n_merge
  #
  # Allocate
  x <- numeric(n_nodes)
  y <- numeric(n_nodes)
  #
  # Leaves
  x[hc$order] <- seq_len(n)
  y[1L:n] <- 0
  #
  # Precompute absolute indices
  left_idx  <- ifelse(merge[, 1L] < 0, -merge[, 1L], n + merge[, 1L])
  right_idx <- ifelse(merge[, 2L] < 0, -merge[, 2L], n + merge[, 2L])
  #
  # Compute internal node positions (iterative dependency)
  i_seq <- seq_len(n_merge)
  y[n + i_seq] <- hc$height
  for (i in i_seq) {
    x[n + i] <- (x[left_idx[i]] + x[right_idx[i]]) / 2
    y[n + i] <- height[i]
  }
  #
  # ---- build segments ----
  # For each merge, we add 3 segments: left vertical, right vertical, top horizontal
  parent_idx <- n + i_seq
  segs <- data.frame(x = numeric(3 * n_merge),
                     y = numeric(3 * n_merge),
                     xend = numeric(3 * n_merge),
                     yend = numeric(3 * n_merge))
  # vertical left
  segs$x[3L*i_seq - 2L] <- x[left_idx]
  segs$y[3L*i_seq - 2L] <- y[left_idx]
  segs$xend[3L*i_seq - 2L] <- x[left_idx]
  segs$yend[3L*i_seq - 2L] <- y[parent_idx]
  # vertical right
  segs$x[3L*i_seq - 1L] <- x[right_idx]
  segs$y[3L*i_seq - 1L] <- y[right_idx]
  segs$xend[3L*i_seq - 1L] <- x[right_idx]
  segs$yend[3L*i_seq - 1L] <- y[parent_idx]
  # horizontal top
  segs$x[3L*i_seq]    <- x[left_idx]
  segs$y[3L*i_seq]    <- y[parent_idx]
  segs$xend[3L*i_seq] <- x[right_idx]
  segs$yend[3L*i_seq] <- y[parent_idx]
  # Labels
  labs <- data.frame(x = x[hc$order],
                     y = y[hc$order],
                     label = hc$labels[hc$order])
  #
  list(segments = segs,
       labels = labs)
}
