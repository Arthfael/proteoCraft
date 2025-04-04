# Script to find all local raw files for the purpose of identifying which folders need cleaning up

TargDir <- "D:/Data/"
setwd(TargDir)
fls <- list.files(recursive = TRUE)
raws <- grep("\\.raw$", fls, value = TRUE, ignore.case = TRUE)
dirs <- unique(dirname(raws))
report <- paste0("The following directories contain raw files which are candidates for deletion:",
                 paste0("\n - ", dirs, collapse = ""), "\n")
cat(report)
write(report, "D:/Detected raw files.txt")
system("open \"D:/Detected raw files.txt\"")
