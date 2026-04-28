#' .bind_worker
#' 
#' Worker function to bind variables to the environment of a function which is to be run in parallel within a larger function.
#'
#' @param fun Function to bind.
#' @param data Named list of variables used by the function; this should include:\cr - the function's argument,\cr - any other package function called by fun's own code using fun() instead of pkg::fun() syntax (usually this would exclusively include other proteoCraft functions such as is.all.good, listMelt etc...)
#' @param envir Environment in which to bind.
#' 
#' @export

.bind_worker <- function(fun,
                         data = list(),
                         parent = baseenv()) {
  e <- new.env(parent = parent)
  list2env(data, e)
  environment(fun) <- e
  return(fun)
}
