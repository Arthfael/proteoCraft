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
  y[1:n] <- 0
  #
  # Precompute absolute indices
  left_idx  <- ifelse(merge[, 1] < 0, -merge[, 1], n + merge[, 1])
  right_idx <- ifelse(merge[, 2] < 0, -merge[, 2], n + merge[, 2])
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
  segs$x[3*i_seq - 2] <- x[left_idx]
  segs$y[3*i_seq - 2] <- y[left_idx]
  segs$xend[3*i_seq - 2] <- x[left_idx]
  segs$yend[3*i_seq - 2] <- y[parent_idx]
  # vertical right
  segs$x[3*i_seq - 1] <- x[right_idx]
  segs$y[3*i_seq - 1] <- y[right_idx]
  segs$xend[3*i_seq - 1] <- x[right_idx]
  segs$yend[3*i_seq - 1] <- y[parent_idx]
  # horizontal top
  segs$x[3*i_seq]    <- x[left_idx]
  segs$y[3*i_seq]    <- y[parent_idx]
  segs$xend[3*i_seq] <- x[right_idx]
  segs$yend[3*i_seq] <- y[parent_idx]
  # Labels
  labs <- data.frame(x = x[hc$order],
                     y = y[hc$order],
                     label = hc$labels[hc$order])
  #
  list(segments = segs,
       labels = labs)
}
