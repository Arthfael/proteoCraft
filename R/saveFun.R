#' saveFun
#'
#' @description
#' Function to save objects quickly.
#' The backend used is qs2::qs_savem.
#' 
#' @param x The object to save.
#' @param file The path/file to save to.
#' 
#' @examples
#' a <- "Hello world!"
#' proteoCraft::saveFun(a, "test.qs2")
#' rm(a)
#' proteoCraft::loadFun("test.qs2")
#' print(a)
#' 
#' @export

saveFun <- function(x, file) {
  #if (.Platform$OS.type == "windows") {
    if (!exists(deparse(substitute(x)), envir = .GlobalEnv)) { error("Object doesn't exist!") }
    tmp <- paste0("qs2::qs_savem(", deparse(substitute(x)),
                  ", file = '", file, "', nthreads = max(c(parallel::detectCores()-1, 1)))")
    #cat(tmp)
    eval(parse(text = tmp), envir = .GlobalEnv)
  # }
  # if (.Platform$OS.type == "unix") {
  #   fastSave::save.lbzip2(x, file = file)
  # }
}
