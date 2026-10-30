targDir <- "B:/home/aRmel_package/proteoCraft/inst/extdata/Sources"
fls <- list.files(targDir, "\\.[rR]$", full.names = TRUE)

flsCode <- setNames(lapply(fls, readLines), fls)
libTst <- vapply(fls, \(fl) { #fl <- fls[1L]
  sum(grepl("(library|require)\\(", flsCode[[fl]])) > 0L
}, TRUE)
w <- which(libTst)
cat(paste(paste0(" - ", fls[w], "\n"), collapse = ""))
