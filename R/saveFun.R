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
#' saveFun(a, "test.qs2")
#' rm(a)
#' loadFun("test.qs2")
#' print(a)
#' 
#' @export

saveFun <- function(x, file) {
  sz <- as.numeric(object.size(x))
  N <- parallel::detectCores()
  max_threads <- ceiling(N/2)
  nThreads <- if (sz < 1e8) { 1L } else {
    if (sz >= 1e9) { max_threads } else { round(1L+(max_threads-1L)*(log10(sz)-8L)) }
  }
  #if (.Platform$OS.type == "windows") {
    if (!exists(deparse(substitute(x)), envir = .GlobalEnv)) { error("Object doesn't exist!") }
    tmp <- paste0("qs2::qs_savem(", deparse(substitute(x)),
                  ", file = '", file, "', nthreads = ", nThreads, ")")
    # Note: keep qs_savem() here! Even though we are saving only one object, if we replace qs_savem with qs_save above,
    # reloading doesn't assign properly the variable to an object with the original name!
    #cat(tmp)
    eval(parse(text = tmp), envir = .GlobalEnv)
  # }
  # if (.Platform$OS.type == "unix") {
  #   fastSave::save.lbzip2(x, file = file)
  # }
}
