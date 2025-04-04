# Generate run info text from raw file
library(svDialogs)

# You need to have ScanHeadsman.exe
# Specify the install directory below:
ScanHdsMnLoc <- "...User_Home/Downloads/ScanHeadsman-1.2.20200730"


filt <- matrix(c("Thermo raw file", "*.raw"), ncol = 2)
if (exists("dr")) { dfltdir <- dr } else { dfltdir <- "...Delivery_Folder" }
fl <- choose.files(paste0(dfltdir, "/.raw"), "Select a single raw file", multi = FALSE, filt, 1)
fl <- normalizePath(fl, winslash = "/")
dr <- dirname(fl)
setwd(dr)
flnm <- basename(fl)
#openwd()
#cmd <- paste0("\"", ScanHdsMnLoc, "/ScanHeadsman.exe\" -i=\"", fl, "\" -n -m=1 -u -t=55")
#cmd <- paste0("\"", ScanHdsMnLoc, "/ScanHeadsman.exe\" -i=\"", fl, "\" -n -u -m=1 -d=None -t=55")
cmd <- paste0("\"", ScanHdsMnLoc, "/ScanHeadsman.exe\" \"", fl, "\" -n -u -m=1 -t=55")
cat(paste0(cmd, "\n"))
system(cmd)
#system(paste0("open \"", methfl, "\""))
