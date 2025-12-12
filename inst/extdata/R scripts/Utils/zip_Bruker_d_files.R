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
#
RPath <- as.data.frame(library()$results)
RPath <- normalizePath(RPath$LibPath[match("proteoCraft", RPath$Package)], winslash = "/")
libPath <- paste0(RPath, "/proteoCraft")
if (!exists("N.clust")) { N.clust <- max(c(round(parallel::detectCores()*0.95)-1, 1)) }
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
