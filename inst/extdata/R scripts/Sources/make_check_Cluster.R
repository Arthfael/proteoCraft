# Let's check this cluster without doing anything with it
require(parallel)
if (!exists("N.clust")) { N.clust <- max(c(round(parallel::detectCores()*0.95)-1, 1)) }
a <- 1
tst <- try(parallel::clusterExport(parClust, "a", envir = environment()), silent = TRUE)
if ("try-error" %in% class(tst)) {
  if (exists("parClust")) { try(parallel::stopCluster(parClust), silent = TRUE) }
  parClust %<o% parallel::makeCluster(N.clust, type = "SOCK")
}
