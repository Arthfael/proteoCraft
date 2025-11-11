# This script is meant as a top level script allowing users to...
#  - ... install and configure the package,
#  - ... choose and run the workflow of interest.

if (!require(pak)) {
  install.packages("pak")
}
require(pak)
# To clean the previous version cached by pak
#a <- pak::cache_summary()
#cat(a$cachepath, "\n")
#pak::cache_delete("proteoCraft")

# To install the latest version
tmp <- installed.packages()
if ("proteoCraft" %in% tmp[, 1]) {
  unloadNamespace("proteoCraft")
  remove.packages("proteoCraft")
}
pak::pkg_install("Arthfael/proteoCraft", upgrade = TRUE, ask = FALSE)
# Alternative way to install:
#devtools::install_github("Arthfael/proteoCraft", upgrade = TRUE)

# Load the package
library(proteoCraft)

# It is always a good idea to configure the package before a fresh installation,
# or even after any update: this will move a copy of the latest analysis workflows to your temporary analysis folder
proteoCraft::Configure()
#proteoCraft::Configure(TRUE) # To update ontologies, good to do once in a while but slower...
#proteoCraft::load_Bckp()

homePath %<o% paste0(normalizePath(Sys.getenv("HOME"), winslash = "/"), "/R/proteoCraft")

# Run workflow of interest:
# =========================
# 
# Choose below the analysis workflow you want to run:
require(svDialogs)
WrkFlwsTbl <- data.frame(Name = c("Protein Groups - with replicates",
                                  "Protein Groups - no replicates",
                                  "Modified peptidoforms - with replicates",
                                  "Histone tails"),
                         Script = c("Regulation analysis - detailed script",
                                    "No replicates analysis - detailed script",
                                    "Regulation analysis - detailed script_pepOnly",
                                    "Histones_analysis_-_DiaNN_FragPipe_Skyline_or_alphaDIA_input"))
myWrkFlw %<o% dlg_list(WrkFlwsTbl$Name, WrkFlwsTbl$Name[1], title = "Choose analysis workflow")$res
wrkflwSrc %<o% paste0(homePath, "/", WrkFlwsTbl$Script[match(myWrkFlw, WrkFlwsTbl$Name)], ".R")
#rstudioapi::documentOpen(wrkflwSrc) # Inspect current workflow script
#rstudioapi::documentOpen(Src) # In case a specific sub-source throws an error
source(wrkflwSrc)


