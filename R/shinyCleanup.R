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

shinyCleanup <- function(max_iter = 100) {
  # Ensure any running app exits
  try(shiny::stopApp(), silent = TRUE)
  # Drain async/reactive event loop
  idle_rounds <- 0
  while ((idle_rounds < 5)&&(max_iter)) {
    ran <- later::run_now(0)
    idle_rounds <- if (!ran) { idle_rounds + 1 } else { 0 }
    Sys.sleep(0.01)
    max_iter <- max_iter - 1
  }
  # final GC to drop observer closures
  gc(full = TRUE)
  invisible(NULL)
}
