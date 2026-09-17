#' saveImgFun
#'
#' @description
#' Function to save all objects in the environment quickly.
#' The backend used is qs2::qs_savem under Windows and fastSave::save.lbzip2 under Linux.
#' 
#' @param file The path/file to save to.
#' @param N.clust A limit on the number of vCPUs to use. If left as NULL (default), uses the number of available clusters - 1, to a minimum of 1. This is probably not a good idea for large sessions however!
#' @param N.reserved Default = 1. Number of reserved vCPUs the function is not to use. Note that for obvious reasons it will always use at least one.
#' 
#' @examples
#' saveImgFun("Backup.qs2")
#' loadFun("Backup.qs2")
#' 
#' @export

saveImgFun <- function(file,
                       N.clust,
                       N.reserved = 1L) {
  TESTING <- FALSE
  #
  #DefArg(saveImgFun); TESTING <- TRUE
  misFun <- if (TESTING) {
    # Note:
    # This is not a perfect alternative to missing but will work in most cases, unless x matches a function imported by a package 
    \(x) { return(!exists(deparse(substitute(x)))) }
  } else { missing }
  #
  dc <- parallel::detectCores()
  if (misFun(N.reserved)) { N.reserved <- 1L }
  nMax <- max(c(dc - N.reserved, 1L))
  if (misFun(N.clust)) { N.clust <- nMax } else {
    if (N.clust > nMax) {
      warning("More cores specified than allowed, I will ignore the specified number! You should always leave at least one free for other processes, see the \"N.reserved\" argument.")
      N.clust <- nMax
    }
  }
  # if (.Platform$OS.type == "windows") {
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
  obj <- obj[which(vapply(obj, \(x) {
    inherits(get(x),
                 c("cluster", # e.g. parClust
                 "connection", # any connection
                 "rdocx") # created by package officer
    )
  }, 1L) == 0L)]
  #
  # Alternate way to using do.call:
  # cmd <- paste0("qs2::qs_savem(",
  #               paste(obj, collapse = ", "),
  #               ", file = \"", file, "\", nthreads = N.clust)")
  # eval(parse(text = cmd))
  do.call(qs2::qs_savem,
          c(lapply(obj, as.symbol),
            file = file,
            nthreads = N.clust)
  )
  # }
  # if (.Platform$OS.type == "unix") {
  #   fastSave::save.image.lbzip2(file = file)
  # }
}
