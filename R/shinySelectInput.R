#' shinySelectInput
#'
#' @description 
#' Convenience wrapper for shiny::selectInput for use in DT data table columns.
#' 
#' @param startVal Starting values of each text input.
#' @param root Root of the created input names. Names generated will be of the form "root___i", where i is the 1:n row number.
#' @param levels Allowed values.
#' @param width The width of the input, e.g. "400px", or "100\%".\cr Default = "600px"
#' @param NAs If TRUE (default), NA selections are allowed (and NA is included as another additional allowed level). If FALSE, they are automatically replaced with the first authorized level at initialization.
#' 
#' @examples
#' # Use within a shiny server as below:
#' 
#' lev <- unique(df$Choice)
#' width <- paste0(30*max(c(nchar(as.character(lev)), 2)), "px")
#' server <- function(input, output, session) {
#' 
#'   Table <- df
#'   Table$Choice <- shinySelectInput(Table$Choice,
#'                                    "Choice",
#'                                    lev,
#'                                    width)
#' 
#' }
#' 
#' @export

shinySelectInput <- function(startVal,
                             root,
                             levels,
                             width = "600px",
                             NAs = TRUE) {
  l <- length(startVal)
  stopifnot(l > 0)
  levels <- unique(levels)
  stopifnot(length(levels) > 0)
  if (NAs) {
    startVal[which(!startVal %in% levels)] <- NA
    levels <- unique(c(levels, NA))
  } else {
    levels <- levels[which(!is.na(levels))]
    stopifnot(length(levels) > 0)
    startVal[which(!startVal %in% levels)] <- levels[1]
  }
  rg <- seq_len(l)
  ids <- paste0(root, "___", as.character(rg))
  if ((length(width) != 1)||(is.na(width))||(is.null(width))) {
    width <- paste0(600, "px")
  }
  inputs <- vapply(rg, function(i) {
    as.character(shiny::selectInput(ids[i],
                                    "",
                                    levels,
                                    startVal[i],
                                    FALSE,
                                    FALSE,
                                    width = width))
  }, "a")
  inputs
}
