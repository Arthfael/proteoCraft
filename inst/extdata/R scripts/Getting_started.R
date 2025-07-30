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
