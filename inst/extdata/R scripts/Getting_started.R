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
#rstudioapi::documentOpen(wrkflwSrc)
source(wrkflwSrc)


