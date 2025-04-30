require(proteoCraft)
saveFun %<o% save
saveImgFun %<o% save.image
loadFun %<o% load
if (.Platform$OS.type == "windows") {
  if (!require("qs2", quietly = TRUE)) { install.packages("qs2") }
  if (require(qs2)) {
    if (exists("cran_req")) { cran_req <- unique(c(cran_req, "qs")) }
    library(qs2)
    saveFun %<o% function(x, file) {
      if (!exists(deparse(substitute(x)), envir = .GlobalEnv)) { error("Object doesn't exist!") }
      tmp <- paste0("qs2::qs_savem(", deparse(substitute(x)),
                    ", file = '", file, "', nthreads = max(c(parallel::detectCores()-1, 1)))")
      #cat(tmp)
      eval(parse(text = tmp), envir = .GlobalEnv)
    }
    saveImgFun <- function(file) { # This one adapted from https://github.com/qsbase/qs2/issues/new?template=Blank+issue
      obj <- base::ls(envir = .GlobalEnv)
      if (exists(".obj", envir = .GlobalEnv)) {
        obj <- unique(c(".obj", obj))
        obj <- obj[which(sapply(obj, exists, envir = .GlobalEnv))]
      }
      obj <- grep("^[A-Za-z\\.][A-Za-z\\.0-9_]*$", obj, value = TRUE)
      do.call(qs2::qs_savem,
              c(lapply(obj, as.symbol),
                file = file,
                nthreads = parallel::detectCores()-1)
      )
    }
    loadFun %<o% function(file) {
      tst <- try(qs2::qs_readm(file, env = globalenv(), nthreads = max(c(parallel::detectCores()-1, 1))), silent = TRUE)
      if ("try-error" %in% class(tst)) { load(file, envir = globalenv()) }
      require(proteoCraft) # This is because we have to remove the 
    }
  }
}
if (.Platform$OS.type == "unix") {
  # Unfortunately, this does not work on Windows :(
  #One can dream...
  if (!require("fastSave", quietly = TRUE)) { devtools::install_github("barkasn/fastSave") }
  if (require(fastSave)) {
    library(fastSave)
    saveFun %<o% save.lbzip2
    saveImgFun %<o% save.image.lbzip2
    loadFun %<o% function(file) {
      tst <- try(load.lbzip2(file, envir = globalenv()), silent = TRUE)
      if ("try-error" %in% class(tst)) { load(file, envir = globalenv()) }
    }
  }
}
