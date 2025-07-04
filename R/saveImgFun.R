#' saveImgFun
#'
#' @description
#' Function to save all objects in the environment quickly.
#' The backend used is qs2::qs_savem under Windows and fastSave::save.lbzip2 under Linux.
#' 
#' @param file The path/file to save to.
#' 
#' @examples
#' proteoCraft::saveImgFun("Backup.qs2")
#' proteoCraft::loadFun("Backup.qs2")
#' 
#' @export

saveImgFun <- function(file) {
  if (.Platform$OS.type == "windows") {
    obj <- base::ls(envir = .GlobalEnv)
    if (exists(".obj", envir = .GlobalEnv)) {
      obj <- unique(c(".obj", obj))
      obj <- obj[which(sapply(obj, exists, envir = .GlobalEnv))]
    }
    obj <- grep("^[A-Za-z\\.][A-Za-z\\.0-9_]*$", obj, value = TRUE)
    # Very important:
    # If you get the `Error in xml_ns.xml_document(x) : external pointer is not valid` error when reloading,
    # then you are presumably saving objects which should not be: some object classes are connections and not meant to be reloaded!
    # New package dependencies may add new object classes with connections which will cause new issues!
    # Check saving some only of the objects in obj until you identify which one is causing the issue,
    # then add its class to the vector below of excluded classes which should not be saved.
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
  if (.Platform$OS.type == "unix") {
    fastSave::save.image.lbzip2(file = file)
  }
}
