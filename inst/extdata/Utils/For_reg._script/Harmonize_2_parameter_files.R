options(stringsAsFactors = FALSE)
library(proteoCraft)
library(svDialogs)

rm(list =ls())

# Get local work directory:
wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
Sys.sleep(1) # Required because sometimes R will attempt to set the word directory too quickly and fail.
setwd(wd)

# Compare two parameter files
filt <- matrix(c("Parameters csv file for regulation analysis script", "Parameters.csv"), nrow = 1,
               dimnames = list("proteoCraft", NULL))
# Choose first parameters file
fl1 <- choose.files("D:\\groups_temp\\.", "Select 1st Parameters.csv file", filters = filt)
wd1 <- dirname(fl1)
setwd(wd1)
# Choose second parameters file
fl2 <- choose.files(paste0(gsub("/", "\\\\", wd1), "Select 2nd Parameters.csv file", "\\."), filters = filt)
wd2 <- dirname(fl2)
if (wd2 == wd1) {
  warning("The second file is the same as the first, choose another one!!!")
  while(wd2 == wd1) {
    fl2 <- choose.files(paste0(gsub("/", "\\\\", wd1), "Select 2nd Parameters.csv file", "\\."), filters = filt)
    wd2 <- dirname(fl2)
  }
}

fl1c <- read.delim(fl1, header = FALSE, sep = ",")
fl1p <- Param.load(fl1, WD_detect = FALSE)
fl2c <- read.delim(fl2, header = FALSE, sep = ",")
fl2p <- Param.load(fl2, WD_detect = FALSE)

params <- data.frame(Name = unique(c(fl1c[,1], fl2c[,1])))
params$Descr <- lapply(params$Name, function(x) { c(fl1c[match(x, fl1c[,1]),3], fl2c[match(x, fl2c[,1]),3]) })
params$Descr_f1 <- sapply(params$Descr, function(x) { !is.na(x[[1]]) })
params$Descr_f2 <- sapply(params$Descr, function(x) { !is.na(x[[2]]) })
params$Descr_conflict <- (params$Descr_f1 + params$Descr_f2 == 2)&sapply(params$Descr, function(x) { !identical(x[[1]], x[[2]]) })
params$Descr_consensus <- ""
Wh <- which(params$Descr_f1 + params$Descr_f2 == 1)
params$Descr_consensus[Wh] <- sapply(params$Descr[Wh], function(x) { unique(x[which(!is.na(x))]) })
Wh <- which((params$Descr_consensus == "")&(!params$Descr_conflict))
params$Descr_consensus[Wh] <- sapply(params$Descr[Wh], unique)
Wh <- which(params$Descr_conflict)
if (length(Wh)) {
  for (wh in Wh) { #wh <- Wh[1]
    print(paste0("Conflict found for parameter \"", params$Name[wh], "\""))
    descr_opt <- params$Descr[[wh]]
    whv2u <- dlgList(descr_opt)$res
    params$Descr_consensus[wh] <- whv2u
  }
} else { print("Yay! No conflicts were found!") }

params2 <- params1 <- params0 <- data.frame(Name = params$Name, Value = "", Description = params$Descr_consensus)
wh1 <- which(params1$Name %in% fl1c[,1])
wh2 <- which(params2$Name %in% fl2c[,1])
params1$Value[wh1] <- fl1c[match(params1$Name[wh1], fl1c[,1]) ,2]
params2$Value[wh2] <- fl2c[match(params2$Name[wh2], fl2c[,1]) ,2]
write.table(params1, paste0(wd1, "/Parameters_unif.csv"), row.names = FALSE, col.names = FALSE, quote = TRUE, sep = ",", qmethod = "double")
write.table(params2, paste0(wd2, "/Parameters_unif.csv"), row.names = FALSE, col.names = FALSE, quote = TRUE, sep = ",", qmethod = "double")
tst1 <- Param.load(paste0(wd1, "/Parameters_unif.csv"), WD_detect = FALSE)
k1 <- colnames(tst1)[which(colnames(tst1) %in% colnames(fl1p))]
w1 <- which(sapply(k1, function(k) { tst1[[k]] != fl1p[[k]] }))
tst2 <- Param.load(paste0(wd2, "/Parameters_unif.csv"), WD_detect = FALSE)
k2 <- colnames(tst2)[which(colnames(tst2) %in% colnames(fl2p))]
w2 <- which(sapply(k2, function(k) { tst2[[k]] != fl2p[[k]] }))
if (length(w1)||length(w2)) { stop("Corruption! Investigate!") } else {
  print("Both unified parameter files reload correctly, no corruption detected.")
}

write.nu.ref <- dlg_message("Should I write a new default based on the consensus?", "yesno")$res
write.nu.ref <- c(FALSE, TRUE)[match(write.nu.ref, c("no", "yes"))]
if (write.nu.ref) {
  test <- lapply(params0$Name, function(x) {
    m <- match(x, params0$Name)
    return(unique(params1$Value[m], params2$Value[m]))
  })
  w <- which(sapply(test, length) == 1)
  params0$Value[w] <- as.character(test[w])
  write.table(params0, "...Destination/Parameters_unif.csv", row.names = FALSE, quote = TRUE, sep = ",", qmethod = "double")
}
