#' plotEval
#' 
#' @description 
#' A simple wrapup for evaluation a ggplot, allowing to save it in a parallel manner.
#' By default this now outputs a gtable, but may also output a full ggplot (if ggplotify is TRUE).
#' 
#' @param plot The ggplot to evaluate. 
#' @param ggplotify Alternate (classic) method using ggplotify, producing a more complete (larger) evaluated ggplot.
#' 
#' @returns
#' Evaluated ggplot.
#' 
#' @export

plotEval <- function(plot,
                     ggplotify = FALSE) {
  if (ggplotify) {
    ggplotify::as.ggplot(ggplotify::as.grob(plot))
  } else {
    ggplot2::ggplot_gtable(ggplot2::ggplot_build(plot))
  }
}
