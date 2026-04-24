# Regularly called cleanup source:
# - Cleans up the environment
# - Updates ScriptPath
# - Saves a backup
# - Optionally destroys/recreates the cluster
rm(list = ls()[which(!ls() %in% .obj)])
clustTst <- try({
  if ((validLogicPar("stopClust"))&&(stopClust)) {
    stopCluster(parClust)
  } else {
    invisible(clusterCall(parClust, \(x) { rm(list = ls());gc() }))
  }
}, silent = TRUE)
Script <- readLines(ScriptPath)
gc()
saveImgFun(BckUpFl)
#loadFun(BckUpFl)
if (!inherits(clustTst, "try-error")) { source(parSrc, local = FALSE) }
