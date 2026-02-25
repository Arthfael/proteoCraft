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
    saveImgFun <- function(file) {
      obj <- base::ls(envir = .GlobalEnv)
      if (exists(".obj", envir = .GlobalEnv)) {
        obj <- unique(c(".obj", obj))
        obj <- obj[which(sapply(obj, exists, envir = .GlobalEnv))]
      }
      obj <- grep("^[A-Za-z\\.][A-Za-z\\.0-9_]*$", obj, value = TRUE)
      # Very important:
      # Some object classes are connections and not meant to be reloaded!
      # If you get the `Error in xml_ns.xml_document(x) : external pointer is not valid` error when reloading,
      # then you are presumably saving objects which should not be.
      # Add their class to the vector below of excluded classes which should not be saved.
      # 
      obj <- obj[which(vapply(obj, function(x) {
        sum(c("cluster", # e.g. parClust
              "connection", # any connection
              "rdocx") # created by package officer
            %in% class(get(x)))
      }, 1) == 0)]
      #
      # Alternate way to using do.call:
      # cmd <- paste0("qs2::qs_savem(",
      #               paste(obj, collapse = ", "),
      #               ", file = \"", file, "\", nthreads = max(c(parallel::detectCores()-1, 1)))")
      # eval(parse(text = cmd))
      do.call(qs2::qs_savem,
              c(lapply(obj, as.symbol),
                file = file,
                nthreads = max(c(parallel::detectCores()-1, 1)))
      )
    }
    loadFun %<o% function(file) {
      tryCatch({
        tst <- try(qs2::qs_readm(file, env = globalenv(), nthreads = max(c(parallel::detectCores()-1, 1))), silent = TRUE)
        if ("try-error" %in% class(tst)) { load(file, envir = globalenv()) }
        library(proteoCraft) # Reload package
      }, warning = "There was at least one error when reloading the backup!")
    }
  }
}
# if (.Platform$OS.type == "unix") {
#   # Unfortunately, this does not work on Windows :(
#   #One can dream...
#   if (!require("fastSave", quietly = TRUE)) { devtools::install_github("barkasn/fastSave") }
#   if (require(fastSave)) {
#     library(fastSave)
#     saveFun %<o% save.lbzip2
#     saveImgFun %<o% save.image.lbzip2
#     loadFun %<o% function(file) {
#       tryCatch({
#         tst <- try(load.lbzip2(file, envir = globalenv()), silent = TRUE)
#         if ("try-error" %in% class(tst)) { load(file, envir = globalenv()) }
#       }, warning = "There was at least one error when reloading the backup!")
#     }
#   }
# }
