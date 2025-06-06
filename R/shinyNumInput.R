#' shinyNumInput
#'
#' @description 
#' A function to create numeric input columns for use in DT data table columns.
#' 
#' @param startVal Starting values of each numeric value in the column.
#' @param min Minimum allowed value
#' @param max Maximum allowed value
#' @param step Interval to use when stepping between min and max
#' @param naRep What to replace NA values with? Default = NA; set e.g. to 1 to use 1 as default replacement for incorrect initial inputs.
#' @param width The width of the input, e.g. '400px', or '100%'; see htmltools::validateCssUnit().
#' @param root Root of the created input names. Names generated will be of the form "root___i", where i is the 1:n row number.
#' 
#' @examples
#' # Use within a shiny server as below:
#' 
#' server <- function(input, output, session) {
#' 
#'   Table$Counts <- shinyNumInput(Table$Counts, 0, Inf, 1, root = "Counts")
#'   Table$Intensity <- shinyNumInput(Table$Counts, root = "Intensity")
#' 
#' }
#' 
#' @export

shinyNumInput <- function(startVal,
                          min = NA,
                          max = NA,
                          step = NA,
                          naRep = NA,
                          width = "300px",
                          root) {
  l <- length(startVal)
  stopifnot(l > 0)
  startVal <- suppressWarnings(as.numeric(startVal))
  naRep <- suppressWarnings(as.numeric(naRep))
  if (!is.na(naRep)) { startVal[which(is.na(startVal))] <- naRep }
  rg <- seq_len(l)
  ids <- paste0(root, "___", as.character(rg))
  inputs <- sapply(rg, function(i) {
    as.character(shiny::numericInput(ids[i],
                                     "",
                                     startVal[i],
                                     0,
                                     width = width))
  })
  inputs
}
