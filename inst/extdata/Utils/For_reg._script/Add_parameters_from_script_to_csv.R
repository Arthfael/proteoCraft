# Detect parameters from latest Regulation analysis script and add them to an existing parameters file
scriptFl <- choose.files("Regulation analysis - detailed script.R", multi = FALSE)
paramFl <- choose.files("Parameters.csv", multi = FALSE)
script <- readLines(scriptFl)
script <- grep("#", script, value = TRUE, invert = TRUE)
param <- read.csv(paramFl, header = FALSE)

pat <- "Param\\$[A-Za-z0-9\\._]+"
New.param <- unique(grep(pat, unlist(strsplit(grep(pat, script, value = TRUE), "[^\\$A-Za-z0-9\\._]+")), script, value = TRUE))
New.param <- gsub("^Param\\$", "", New.param)
tst <- grep("Param\\$", New.param, value = TRUE)
stopifnot(length(tst) == 0)

New.param <- New.param[which(!New.param %in% param$V1)]
New.param
if (length(New.param)) {
  param1 <- param[grep("[D,d]eprecated", param$V2, invert = TRUE),]
  param2 <- param[grep("[D,d]eprecated", param$V2),]
  New.param <- data.frame(V1 = New.param, V2 = "", V3 = "")
  param <- rbind(param1, New.param, param2)
  write.table(param, gsub("\\.csv$", "_Completed.csv", paramFl), sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)
}