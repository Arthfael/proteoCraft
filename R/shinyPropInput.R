#' shinyPropInput
#'
#' @description 
#' Convenience wrapper for shiny::numericInput for use in DT data table columns.
#' This is specifically to query users for a proportion (between 0 and 1 by steps of 0.001).
#' 
#' @param startVal Starting values of each proportion.
#' @param width The width of the input, e.g. "400px", or "100%"; see htmltools::validateCssUnit().
#' @param root Root of the created input names, default = "Proportion". Input names generated will be of the form "root___i", where i is the 1:n row number.
#' 
#' @examples
#' # Use within a shiny server as below:
#' 
#' server <- function(input, output, session) {
#' 
#'   Table$Proportion <- shinyPropInput(Table$Proportion)
#' 
#' }
#' 
#' @export

shinyPropInput <- function(startVal,
                           width = "300px",
                           root = "Proportion") {
  l <- length(startVal)
  stopifnot(l > 0)
  # Pre-process values
  startVal <- suppressWarnings(as.numeric(startVal))
  startVal[which(is.na(startVal))] <- 1
  startVal[which(startVal > 1)] <- 1
  startVal[which(startVal < 0)] <- 0
  rg <- seq_len(l)
  ids <- paste0(root, "___", as.character(rg))
  inputs <- sapply(rg, function(i) {
    as.character(shiny::numericInput(ids[i],
                                     "",
                                     startVal[i],
                                     0,
                                     1,
                                     0.001,
                                     width))
  })
  inputs
}
