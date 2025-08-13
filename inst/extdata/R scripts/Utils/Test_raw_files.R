wd <- rstudioapi::selectDirectory(path = "C:/")
setwd(wd)
tst <- grep("ScanHeadsman-1\\.2", list.dirs("C:/", recursive = FALSE), value = TRUE)
ScanHdsMnTst <- length(tst) > 0
if (ScanHdsMnTst) { ScanHdsMnLoc <- normalizePath(tst[1], winslash = "/") } else {
  warning("No valid version of ScanHeadsman-1.2 could be located in the system!")
}

# Test with ScanHeadsman
raws <- grep("\\.raw$", list.files(), value = TRUE)
for (raw in raws) { #raw <- raws[1]
  cmd <- paste0("\"", ScanHdsMnLoc, "/ScanHeadsman.exe\" \"", raw, "\" -n -u -m=1 -t=55")
  #cat(paste0(cmd, "\n"))
  system(cmd)
}

# Test with rawr
if (!suppressWarnings(require(rawrr))) {
  devtools::install_github("cpanse/rawrr")
  #install.packages('http://fgcz-ms.uzh.ch/~cpanse/rawrr_0.2.1.tar.gz', repo = NULL)
  installRawFileReaderDLLs(sourceUrl = .thermofisherlsmsUrl())
  rawrr::installRawrrExe()
}
require(rawrr)
require(parallel)
dc <- detectCores()
N.clust <- max(c(dc-1, 1))
cl <- makeCluster(N.clust, type = "SOCK")
exports <- list("raws", "wd")
clusterExport(cl, exports, envir = environment())
clusterCall(cl, function() library(rawrr))
tic <- parLapply(cl, raws, function(raw) {
  x2 <- try(readChromatogram(paste0(wd, "/", raw), type = "tic"), silent = TRUE)
  if (!"try-error" %in% class(x2)) {
    x2 <- data.frame("Raw file" = raw,
                     "Retention time" = as.numeric(x2$times),
                     Intensity = as.numeric(x2$intensities),
                     check.names = FALSE)
  }
  return(x2)
})
stopCluster(cl)
setNames(tic, raws)
tst <- setNames(sapply(tic, function(x) { !"try-error" %in% class(x) }), raws)
print(tst)
