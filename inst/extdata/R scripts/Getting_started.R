# Useful commands

# To clean the previous version cached by pak
if (!require(pak)) {
  install.packages("pak")
}
require(pak)
#a <- pak::cache_summary()
#cat(a$cachepath, "\n")
#pak::cache_delete("proteoCraft")

# To install the latest version
try({
  unloadNamespace("proteoCraft")
  remove.packages("proteoCraft")
}, silent = TRUE)
pak::pkg_install("Arthfael/proteoCraft", upgrade = TRUE, ask = FALSE)
# Alternative way to install:
#devtools::install_github("Arthfael/proteoCraft", upgrade = TRUE)

# Load the package
library(proteoCraft)

# It is always a good idea to configure the package before a fresh installation,
# or even after any update: this will move a copy of the latest analysis workflows to your temporary analysis folder
proteoCraft::Configure()
#proteoCraft::load_Bckp()

homePath %<o% paste0(normalizePath(Sys.getenv("HOME"), winslash = "/"), "/R/proteoCraft")

# Run workflow of interest:
# =========================
# 
# Choose below the workflow you want to run
#
# Full proteomics with replicates
wrkflwSrc <- paste0(homePath, "/Regulation analysis - detailed script.R")
source(wrkflwSrc)
#
# Full proteomics without replicates
wrkflwSrc <- paste0(homePath, "/No replicates analysis - detailed script.R")
source(wrkflwSrc)
#
# Modified peptides only, with replicates
wrkflwSrc <- paste0(homePath, "/Regulation analysis - detailed script_pepOnly.R")
source(wrkflwSrc)
#
# Histones analysis
wrkflwSrc <- paste0(homePath, "/Histones_analysis_-_DiaNN_FragPipe_Skyline_or_alphaDIA_input.R")
source(wrkflwSrc)
