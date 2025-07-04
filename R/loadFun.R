#' loadFun
#'
#' @description
#' Function to reload object(s) saved using saveFun() or saveImgFun().
#' The backend used is qs2::qs_savem under Windows and fastSave::save.lbzip2 under Linux.
#' 
#' @param file The path/file to save to.
#' @param tryClassic In past versions, we used the classic R save/save.image/load functions. If loading a backup fails, should we try to re-load it using these functions? Default = FALSE.
#' 
#' @examples
#' a <- "Hello world!"
#' proteoCraft::saveFun(a, "test.qs2")
#' rm(a)
#' proteoCraft::loadFun("test.qs2")
#' print(a)
#' 
#' @export

loadFun <- function(file, tryClassic = FALSE) {
  tryCatch({
    if (.Platform$OS.type == "windows") {
      tst <- try(qs2::qs_readm(file, env = globalenv(), nthreads = max(c(parallel::detectCores()-1, 1))), silent = TRUE)
    }
    if (.Platform$OS.type == "unix") {
      tst <- try(fastSave::load.lbzip2(file, envir = globalenv()), silent = TRUE)
    }
    if (("try-error" %in% class(tst))&&(tryClassic)) { load(file, envir = globalenv()) } # This is for historic reasons
  }, warning = "There was at least one error when reloading the backup!")
  library(proteoCraft) # Reload package
}
