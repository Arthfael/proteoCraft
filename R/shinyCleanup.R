#' shinyCleanup
#' 
#' @description 
#' A function to clean up after a shiny app finishes running, so that the next one is not affected by leftover from the previous.
#' I do not have time for designing this type of function, so must admit that this is 100% from chatGPT.
#' 
#' @param max_iter Max number of iterations, default = 100.
#' 
#' @returns
#' This function does not return anything.
#' 
#' @export

shinyCleanup <- function(max_iter = 100L) {
  # Ensure any running app exits
  #try(shiny::stopApp(), silent = TRUE) # Nope, this may be breaking the flow of code execution...
  # Drain async/reactive event loop
  # idle_rounds <- 0L
  # while ((idle_rounds < 5L)&&(max_iter)) {
  #   ran <- later::run_now(0L) # ... or this could be what's breaking the flow...
  #   idle_rounds <- if (!ran) { idle_rounds + 1L } else { 0L }
  #   Sys.sleep(0.01)
  #   max_iter <- max_iter - 1L
  # }
  # final GC to drop observer closures
  Sys.sleep(0.05)
  gc(full = TRUE)
  invisible(NULL)
}
