# Re-load a project from a renv lock file
#
# This is meant for analysis reproducibility, to restore an older R environment. Before running this script:
#  - Check that you have the correct version of R selected in RStudio (comparing against the old rsession info output).
#  - Create the folder you will want to work from.
# Then run this script and select the folder you just created.
# If this was not already a project folder, a project will be initiated and the R session will restart.
# Browse to the lock file you want to reload and select it, the old renv will be reloaded.

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
# This script will restart the R session, which will interrupt the sequence of commands. #
#                                       |        |                                       #
#                                       |        |                                       #
#                                       |        |                                       #
#                                      \          /                                      #
#                                       \        /                                       #
#                                        \      /                                        #
#                                         \    /                                         #
#                                          \  /                                          #
#                                           \/                                           #
#                                  Run it line by line!                                  #
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #

require(renv)
wd <- rstudioapi::selectDirectory("Select project directory")
#if (!dir.exists(wd)) { dir.create(wd, recursive = TRUE) }
setwd(wd)
#proteoCraft::openwd()
nuEnv <- (!file.exists("renv.lock"))
if (nuEnv) { renv::init() }
lockFl <- rstudioapi::selectFile("Select renv lock file to restore", filter = "renv lock files (*.lock)")
docsPath <- normalizePath(Sys.getenv("HOME"), winslash = "/")
if (grepl("^~", lockFl)) { lockFl <- gsub("^~", docsPath, lockFl) }
renv::restore(lockfile = lockFl, prompt = FALSE)
renv::status()
renv::project()
require(proteoCraft)
renv %<o% TRUE
