#' shinyTextInput
#'
#' @description 
#' Convenience wrapper for shiny::textInput for use in DT data table columns.
#' 
#' @param startVal Starting values of each text input.
#' @param root Root of the created input names. Names generated will be of the form "root___i", where i is the 1:n row number.
#' @param width The width of the input, e.g. "400px", or "100%"; see htmltools::validateCssUnit().
#' 
#' @examples
#' # Use within a shiny server as below:
#' 
#' server <- function(input, output, session) {
#' 
#'   Table$Free_input <- shinyTextInput("Enter something...",
#'                                      "FreeInput")
#' 
#' }
#' 
#' @export

shinyTextInput <- function(startVal,
                           root,
                           width) {
  l <- length(startVal)
  stopifnot(l > 0)
  rg <- seq_len(l)
  ids <- paste0(root, "___", as.character(rg))
  inputs <- sapply(rg, function(i) {
    as.character(shiny::textInput(ids[i],
                                  "",
                                  startVal[i],
                                  width))
  })
  inputs
}
