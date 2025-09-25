#' shinyFDInput
#'
#' @description 
#' A function to create column fill-down buttons for use in DT data table columns.
#' Based on shiny::actionButton.
#' 
#' @param root Root of the created input names. Names generated will be of the form "root___i___FD", where i is the 1:n row number.
#' @param l Number of rows of the table.
#' @param useIcon Ignored, only kept because of old code!
#' @param width The width of the input, e.g. "400px", or "100\%"; see htmltools::validateCssUnit().\cr Default = "15px"
#' 
#' @examples
#' # Use within a shiny server as below:
#' 
#' lev <- unique(df$Choice)
#' width <- paste0(30*max(c(nchar(as.character(lev)), 2)), "px")
#' nr <- nrow(df)
#' server <- function(input, output, session) {
#' 
#'   Table <- df
#'   Table$Choice <- shinySelectInput(Table$Choice,
#'                                    "Choice",
#'                                    lev,
#'                                    width)
#'   Table$Choice___FD <- shinyFDInput("Choice",
#'                                     nr,
#'                                     TRUE)
#' 
#' }
#' 
#' @export

shinyFDInput <- function(root,
                         l,
                         useIcon = FALSE,
                         width = "15px") {
  stopifnot(l > 0)
  rg <- seq_len(l)
  ids <- paste0(root, "___", as.character(rg), "___FD")
  if ((length(width) != 1)||(is.na(width))||(is.null(width))) {
    width <- paste0(15, "px")
  }
  inputs <- vapply(rg, function(i) {
    as.character(shinyWidgets::actionBttn(ids[i],
                                          "",
                                          shiny::icon("arrow-down"),
                                          size = "xs",
                                          color = "primary",
                                          style = "pill",
                                          width = width))
  }, "a")
  inputs
}
