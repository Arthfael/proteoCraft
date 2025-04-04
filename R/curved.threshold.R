#' curved.threshold
#'
#' @description
#' Placeholder!!! Should be replaced eventually by a real curved threshold (as in SAM).
#' A function to create a curved threshold for significance in volcano plots.
#' Can be used to calculate the local threshold value for a given x. 
#' As input to ggplot2's stat_function, can also be used to draw said threshold.
#' 
#' @param x The value on the horizontal axis (usually a log2 ratio).
#' @param s0 Tuned SAM s0 fudging parameter
#' @param fdr Acceptable FDR values.
#' 
#' @examples
#' dat$Theshold <- thresh(seq(-xm, xm, 0.1))
#' Add the threshold graph to the plot:
#' plot <- plot + stat_function(fun = thresh, color = "red") + xlim(-xm, xm)
#' 

curved.threshold <- function(x,
                             s0,
                             fdr) {
  return(x)
}
