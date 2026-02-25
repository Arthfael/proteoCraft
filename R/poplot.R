#' poplot
#'
#' @description 
#' A wrapup for printing a ggplot in a popup window (or not: can also just do in local window). Should work on all devices (some strange people like Macs.)
#' 
#' @param plot ggplot to print.
#' @param height Height of the popup. Default = 10.
#' @param width Width of the popup. Default = 15.
#' @param new.window TRUE by default, set to FALSE to not plot in a popup window.
#' @param max.gr.dev 62 by default, the maximum number of graphical devices which the system will allow to be open simultaneously. I could not find the information anywhere but checked on 3 different PCs and the result was always 62.
#'
#' @export

poplot <- function(plot,
                   height = 10,
                   width = 10,
                   new.window = TRUE,
                   max.gr.dev = 62) {
  if (max.gr.dev < 2) {
    warning("Ignoring \"max.gr.dev\" as the value isn't valid!")
    max.gr.dev <- 62
  }
  if (new.window) {
    dl <- dev.list()
    if (length(dl) > max.gr.dev-2) {
      warning(paste0("Closing ", max(c(2, max.gr.dev))-1, " graphical windows to avoid exceeding the maximum allowed number!"))
      graphics.off()
    }
    if (.Platform$OS.type == "windows") {
      grDevices::windows(width = width, height = height)
    } else {
      if (.Platform$OS.type != "unix") { warning("What system is this? Defaulting to x11()") }
      grDevices::x11(width = width, height = height)
    }
  }
  print(plot)
}
