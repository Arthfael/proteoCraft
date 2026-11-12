# Let's check this cluster without doing anything with it
#parallel::stopCluster(parClust)
if (!exists("N.clust")) { N.clust <- max(c(round(parallel::detectCores()*0.95)-1L, 1L)) }
N.clust %<o% N.clust
# Check cluster
clusTst <- exists("parClust")
if (clusTst) {
  tst2 <- clustCheck(parClust)
  if (!tst2) {
    try({
      parallel::stopCluster(parClust)
      rm(parClust)
    }, silent = TRUE)
    clusTst <- FALSE
  }
}
if (!clusTst) {
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
  pkgs <- c()
  if (exists("cran_req")) { pkgs <- base::union(pkgs, cran_req) }
  if (exists("bioc_req")) { pkgs <- base::union(pkgs, bioc_req) }
  if (length(pkgs)) {
    parallel::clusterExport(parClust, "pkgs", envir = environment())
    invisible(parallel::clusterCall(parClust, \() { lapply(pkgs, requireNamespace) }))
  }
}

# Internal, for when testing: re-run script to load local versions of proteoCraft functions
if (!exists("useLocalFun")) { useLocalFun %<o% FALSE }
if (useLocalFun && exists("locPkgPath") && dir.exists(locPkgPath)) {
  pkgload::load_all(locPkgPath, quiet = TRUE)
  parallel::clusterExport(parClust, "locPkgPath", envir = environment())
  # Below: check commented because, even if this is already a dev version, we MUST check that we have the latest versions!
  #isDev <- unlist(parallel::clusterCall(parClust, \() { pkgload::is_dev_package("proteoCraft") }))
  #if (!all(isDev)) {
    invisible(parallel::clusterCall(parClust, \() {
      pkgload::load_all(locPkgPath, quiet = TRUE)
      return(NULL)
    }))
  #}
} else {
  invisible(parallel::clusterCall(parClust, \() { require(proteoCraft) }))
}
