# Let's check this cluster without doing anything with it
#stopCluster(parClust)
require(parallel)
if (!exists("N.clust")) { N.clust <- max(c(round(parallel::detectCores()*0.95)-1L, 1L)) }
N.clust %<o% N.clust
a <- 1
tst <- try(parallel::clusterExport(parClust, "a", envir = environment()), silent = TRUE)
if (inherits(tst, "try-error")) {
  if (exists("parClust")) { try(parallel::stopCluster(parClust), silent = TRUE) }
  parClust %<o% parallel::makeCluster(N.clust, type = "SOCK")
}
#invisible(parallel::clusterCall(parClust, \() { NULL }))
#
# Optional code for loading packages onto the cluster - set loadpack to TRUE to run it
# (can take a while to run)
if (!validLogicPar("loadpack")) { loadpack <- FALSE }
if (loadpack) {
  pck <- c()
  if (exists("cran_req")) { pck <- union(pck, cran_req) }
  if (exists("bioc_req")) { pck <- union(pck, bioc_req) }
  if (length(pck)) {
    parallel::clusterExport(parClust, "pck", envir = environment())
    invisible(parallel::clusterCall(parClust, \() { lapply(pck, \(x) { library(x, character.only = TRUE) }) }))
  }
}
