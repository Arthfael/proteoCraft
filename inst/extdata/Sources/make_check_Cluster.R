# Let's check this cluster without doing anything with it
#parallel::stopCluster(parClust)
if (!exists("N.clust")) { N.clust <- max(c(round(parallel::detectCores()*0.95)-1L, 1L)) }
N.clust %<o% N.clust
# Check cluster
tst1 <- exists("parClust")
if (tst1) {
  tst2 <- clustCheck(parClust)
  if (!tst2) {
    try({
      parallel::stopCluster(parClust)
      rm(parClust)
    }, silent = TRUE)
    tst1 <- FALSE
  }
}
if (!tst1) {
  parClust %<o% parallel::makeCluster(N.clust, type = "SOCK")
}
#invisible(parallel::clusterCall(parClust, \() { NULL }))

if (!exists("validLogicPar")) {
  # Boolean functions to check parameter values
  Src <- paste0(libPath, "/extdata/Sources/parBooleans.R")
  #rstudioapi::documentOpen(Src)
  source(Src)
}

# Optional code for loading packages onto the cluster - set loadpack to TRUE to run it
# (can take a while to run)
if (!validLogicPar("loadpack")) { loadpack <- FALSE }
if (loadpack) {
  pck <- c()
  if (exists("cran_req")) { pck <- base::union(pck, cran_req) }
  if (exists("bioc_req")) { pck <- base::union(pck, bioc_req) }
  if (length(pck)) {
    parallel::clusterExport(parClust, "pck", envir = environment())
    invisible(parallel::clusterCall(parClust, \() { lapply(pck, requireNamespace) }))
  }
}

# Internal, for when testing: re-run script to load local versions of proteoCraft functions
if (!exists("useLocalFun")) { useLocalFun %<o% FALSE }
if (useLocalFun && exists("locPkgPath") && dir.exists(locPkgPath)) {
  pkgload::load_all(locPkgPath, quiet = TRUE)
  parallel::clusterExport(parClust, "locPkgPath", envir = environment())
  isDev <- unlist(parallel::clusterCall(parClust, \() { pkgload::is_dev_package("proteoCraft") }))
  if (!all(isDev)) {
    invisible(parallel::clusterCall(parClust, \() {
      pkgload::load_all(locPkgPath, quiet = TRUE)
      return(NULL)
    }))
  }
} else {
  invisible(parallel::clusterCall(parClust, \() { requireNamespace("proteoCraft") })) # otherwise this is handled by parSrc
}
