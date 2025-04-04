#' chooseFrom
#'
#' @description
#' A wrapper function to allow easy selection from a range of options.
#' Completely useless really... what was I thinking.
#' It does in a more complicated way what the wrapped function did...
#' Candidate for deletion...
#' 
#' @param title The caption displayed at the top of the choices.
#' @param choices Available items to choose from.
#' @param default Integer, the choice by default (default = 1).
#' @param multiple Are multiple choices allowed? Default = FALSE
#'
#' @examples
#' temp <- chooseFrom("What do you want?", c("Peace", "War"), 2) # (Please choose 1!)
#' 
#' @export

chooseFrom <- function(title,
                       choices,
                       default = 1,
                       multiple = FALSE) {
  warning("I am a useless function, please phase me out...")
  default <- as.integer(default)
  stopifnot(is.integer(default),
            default %in% 1:length(choices))
  msg <- paste0(title, "\n",  paste(paste0(" - ", 1:length(choices), ": ", choices), collapse = "\n"))
  if (multiple) {
    svDialogs::dlgList(choices, default, TRUE, title)
  } else {
    x <- suppressWarnings(as.integer(svDialogs::dlg_input(msg, default)$res))
    while ((is.na(x))||(!x %in% 1:length(choices))) {
      x <- suppressWarnings(as.integer(svDialogs::dlg_input(msg, default)$res))
    }
  }
  return(choices[x])
}
