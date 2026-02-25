#' .ggSave_worker
#' 
#' Worker function within Volcano.plot() for saving evaluated ggplots.
#'
#' @param x A list with 4 slots:\cr
#' - Plot: the evaluated ggplot to save\cr
#' - Path: the directory in which to save the plot\cr
#' - Ttl: name of the destination file (without extension)\cr
#' - Ext: extension of the destination file

.ggSave_worker <- function(x) {
  # Even though the plot is evaluated, we must load ggrepel explicitly on the cluster,
  # otherwise the plot is saved without the ggrepel labels!
  require(ggrepel)
  suppressMessages({
    ggplot2::ggsave(paste0(x$Path, "/", x$Ttl, ".", x$Ext), x$Plot,
                    dpi = 300, width = 10, height = 10, units = "in")
  })
}