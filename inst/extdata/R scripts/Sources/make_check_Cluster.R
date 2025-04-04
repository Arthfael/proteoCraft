# Let's check this cluster without doing anytthing with it 
a <- 1
tst <- try(clusterExport(parClust, "a", envir = environment()), silent = TRUE)
if ("try-error" %in% class(tst)) {
  if (exists("parClust")) { try(stopCluster(parClust), silent = TRUE) }
  parClust %<o% makeCluster(N.clust, type = "SOCK")
}
