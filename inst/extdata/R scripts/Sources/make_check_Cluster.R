# Let's check this cluster without doing anything with it
if (!exists("N.clust")) { N.clust <- max(c(round(parallel::detectCores()*0.95)-1, 1)) }
a <- 1
tst <- try(clusterExport(parClust, "a", envir = environment()), silent = TRUE)
if ("try-error" %in% class(tst)) {
  if (exists("parClust")) { try(stopCluster(parClust), silent = TRUE) }
  parClust %<o% makeCluster(N.clust, type = "SOCK")
}
