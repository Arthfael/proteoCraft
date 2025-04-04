#' shinyCheckInput
#'
#' @description 
#' Convenience wrapper for shiny::checkboxInput for use in DT data table columns.
#' 
#' @param startVal Starting values of each checkbox.
#' @param root Root of the created input names. Names generated will be of the form "root___i", where i is the 1:n row number.
#' 
#' @examples
#' # Use within a shiny server as below:
#' 
#' server <- function(input, output, session) {
#' 
#' ...
#' 
#'   Table$CheckBox <- shinyCheckInput(Table$UseThisSample)
#' 
#' ...
#' 
#' }
#' 
#' @export

shinyCheckInput <- function(startVal,
                            root = "check") {
  l <- length(startVal)
  stopifnot(l > 0)
  rg <- seq_len(l)
  ids <- paste0(root, "___", as.character(rg))
  inputs <- sapply(rg, function(i) {
    as.character(shiny::checkboxInput(ids[i],
                                      NULL,
                                      as.logical(startVal[i])))
  })
  inputs
}
