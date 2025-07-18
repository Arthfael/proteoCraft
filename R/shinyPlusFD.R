#' shinyPlusFD
#'
#' @description 
#' A function to create a "+" button column fill-down button for use in DT data table columns.
#' Could be used for many different things depending on how you write the shiny server, here we use it for filling down with increment.
#' Based on shiny::actionButton.
#' 
#' @param root Root of the created input names. Names generated will be of the form "root___i___INCR", where i is the 1:n row number.
#' @param l Number of rows of the table.
#' @param buttonLabel Which text to display on the button. Default = "+"
#' @param stem Input suffix. Default = "INCR"
#' @param width The width of the input, e.g. "400px", or "100%"; see htmltools::validateCssUnit().\cr Default = "15px"
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
#'   Table$Choice___INCR <- shinyPlusFD("Choice",
#'                                      nr)
#' 
#' }
#' 
#' @export

shinyPlusFD <- function(root,
                        l,
                        buttonLabel = "+",
                        stem = "INCR",
                        width = "15px") {
  stopifnot(l > 0)
  rg <- seq_len(l)
  ids <- paste0(root, "___", as.character(rg), "___", stem)
  if ((length(width) != 1)||(is.na(width))||(is.null(width))) {
    width <- paste0((nchar(buttonLabel) + 1)*10, "px")
  }
  inputs <- vapply(rg, function(i) {
    if (!nchar(buttonLabel)) {
      rs <- shinyWidgets::actionBttn(ids[i],
                                     "",
                                     shiny::icon("arrow-down"),
                                     size = "xs",
                                     color = "primary",
                                     style = "pill",
                                     width = width)
    } else {
      rs <- shinyWidgets::actionBttn(ids[i],
                                     buttonLabel,
                                     size = "xs",
                                     color = "primary",
                                     style = "pill",
                                     width = width)
    }
    return(as.character(rs))
  }, "a")
  inputs
}
