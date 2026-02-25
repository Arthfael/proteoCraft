#' plotEval
#' 
#' @description 
#' A simple wrapup for evaluation a ggplot, allowing to save it in a parallel manner.
#' 
#' @param plot ggplot to evaluate.
#' 
#' @returns
#' Evaluated ggplot.
#' 
#' @export

plotEval <- function(plot) {
  ggplotify::as.ggplot(ggplotify::as.grob(plot))
}
