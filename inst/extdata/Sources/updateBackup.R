# Regularly called cleanup source:
# - Cleans up the environment
# - Updates ScriptPath
# - Saves a backup
# - Checks the cluster and fixes it if it existed
inEnv <- ls()
.obj <- base::intersect(.obj, inEnv) # New step to avoid keeping in memory objects which have been removed but used to be remanent, (e.g. because I changed my mind):
rm(list = base::setdiff(inEnv, .obj))
Script <- readr::read_lines(ScriptPath)
gc()
saveImgFun(BckUpFl)
#loadFun(BckUpFl)
if ("parClust" %in% .obj) { source(parSrc) }
