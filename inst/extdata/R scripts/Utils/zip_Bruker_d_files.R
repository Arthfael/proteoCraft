# Mass zip Bruker .d folders to avoid having to do it manually, which just kills me...
require(rstudioapi)
require(proteoCraft)
require(openxlsx2)

homePath %<o% paste0(normalizePath(Sys.getenv("HOME"), winslash = "/"), "/R/proteoCraft")
locFl <- paste0(homePath, "/Default_locations.xlsx")
locs <- openxlsx2::read_xlsx(locFl)
dflt <- archDir <- locs$Path[match("Archive folder", locs$Folder)]

if ((exists("targDir"))&&(dir.exists(targDir))) { dflt <- targDir }
targDir <- rstudioapi::selectDirectory(path = dflt)
#
allDirs <- list.dirs(targDir, full.names = TRUE, recursive = TRUE)
dDirs <- grep("\\.d$", allDirs, value = TRUE)
dDirs <- svDialogs::dlg_list(dDirs, dDirs, TRUE, "Selected Bruker .d directories to zip...")$res
if (length(dDirs)) {
  # We don't want any .d directory nested inside another - should never happen, but you never know.
  # If one is ever a sub-directory of another, we exclude the former.
  nc_dDir <- setNames(nchar(dDirs), dDirs)
  tst <- vapply(dDirs, function(x) {
    rest <- dDirs[which(dDirs != x)]
    w <- which(nc_dDir[rest] < nc_dDir[x])
    if (!length(w)) { return(FALSE) }
    rest <- rest[w]
    tmp <- vapply(rest, function(y) {
      substr(x, 1, nc_dDir[y]) == y
    }, TRUE)
    return(sum(tmp) > 0)
  }, TRUE)
  dDirs <- dDirs[which(!tst)]
  #
  RPath <- as.data.frame(library()$results)
  RPath <- normalizePath(RPath$LibPath[match("proteoCraft", RPath$Package)], winslash = "/")
  libPath <- paste0(RPath, "/proteoCraft")
  parSrc <- paste0(libPath, "/extdata/R scripts/Sources/make_check_Cluster.R")
  #rstudioapi::documentOpen(parSrc)
  source(parSrc, local = FALSE)
  #
  parLapply(parClust, dDirs, function(dr) { #dr <- dDirs[1]
    parDr <- dirname(dr)
    setwd(parDr)
    nm <- basename(dr)
    zip(paste0(dr, ".zip"),
        nm)
    return()
  })
}
