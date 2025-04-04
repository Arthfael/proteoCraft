# A script to combine 2 DiaNN output folders into one,
# for the purpose of combining the results into a single data analysis using the scripts from this package.
# Useful if, for instance, you searched different files from the same experiment with different parameters but want to combine the outputs.
# This could be easily rewritten to allow for N different folders, rather than exactly 2.

require(rstudioapi)
require(data.table)
try(setDTthreads(threads = parallel::detectCores()-1), silent = TRUE)


dflt <- "C:"
if ((exists("dir1"))&&(dir.exists(dir1))) { dflt <- dir1 }
fl2lg <- fl1lg <- normalizePath(choose.files(paste0(dflt, "/*.log.txt"),
                                             "Select 1st DiaNN log file", FALSE), winslash = "/")
dir1 <- dirname(fl1lg)
while (fl2lg ==  fl1lg) {
  dflt <- dir1
  if ((exists("dir2"))&&(nchar(paste0(dir2, "")))&&(dir.exists(dir2))) { dflt <- dir2 }
  dflt <- gsub("/[^/]+$", "", dflt)
  fl2lg <- normalizePath(choose.files(paste0(dflt, "/*.log.txt"),
                                      "Select 2nd DiaNN log file", FALSE), winslash = "/")
  dir2 <- dirname(fl2lg)
}
dstDir <- NA
while ((is.na(dstDir))||(dstDir == dir1)||(dstDir == dir2)) { dstDir <- selectDirectory("Select destination directory", path = gsub("/[^/]+$", "/.*", dir2)) }
setwd(dstDir)

# Parse and edit calls to create a fake .log.txt
dirs <- logFls <- logs <- calls <- list()
for (i in 1:2) { #i <- 1
  nm <- paste0("DiaNN", i)
  dirs[[nm]] <- get(paste0("dir", i))
  logFls[[nm]] <- get(paste0("fl", i, "lg"))
  tmp <- readLines(logFls[[nm]])
  logs[[nm]] <- tmp
  tmp <- grep("^ *diann\\.exe ", tmp, value = TRUE)[1]
  tmp <- unlist(strsplit(tmp, " +--"))
  tmp <- tmp[2:length(tmp)]
  tmp <- strsplit(tmp, " +")
  l <- sapply(tmp, length)
  tmp[which(l == 1)] <- lapply(tmp[which(l == 1)], function(x) { c(x, "") })
  tmp <- data.frame(Arg = sapply(tmp, function(x) { x[[1]] }),
                    Val = sapply(tmp, function(x) { x[[2]] }))
  calls[[nm]] <- tmp
}
files <- unique(unlist(lapply(paste0("DiaNN", 1:2), function(nm) {
  tmp <- calls[[nm]]
  tmp$Val[which(tmp$Arg == "f")]
})))
fastas <- unique(unlist(lapply(paste0("DiaNN", 1:2), function(nm) {
  tmp <- calls[[nm]]
  tmp$Val[which(tmp$Arg == "fasta")]
})))
fastas <- fastas[which(nchar(fastas) > 0)]
libs <- unique(unlist(lapply(paste0("DiaNN", 1:2), function(nm) {
  tmp <- calls[[nm]]
  tmp$Val[which(tmp$Arg == "lib")]
})))
libs <- libs[which(nchar(libs) > 0)]
varmods <- unique(unlist(lapply(paste0("DiaNN", 1:2), function(nm) {
  tmp <- calls[[nm]]
  tmp$Val[which(tmp$Arg == "var-mod")]
})))
reports <- setNames(lapply(paste0("DiaNN", 1:2), function(nm) {
  tmp <- calls[[nm]]
  tmp <- tmp$Val[which(tmp$Arg == "out")]
  gsub("\\\\", "/", tmp)
}), paste0("DiaNN", 1:2))
FDRs <- as.numeric(unique(unlist(sapply(paste0("DiaNN", 1:2), function(nm) {
  tmp <- calls[[nm]]
  tmp$Val[which(tmp$Arg == "qvalue")]
}))))
if (length(FDRs) > 1) { warning("Do we really want to combine datasets at different FDRs?!") }

# Create fake DiaNN-like call to put into a fake log.txt file
call <- calls[[1]]
call$Val[which(call$Arg == "out")] <- gsub("/", "\\\\",paste0(dstDir, "/report.tsv"))
call$Val[which(call$Arg == "call")] <- ""
wf <- which(call$Arg == "f")
wbf <- 1:(min(wf)-1)
waf <- 1:nrow(call)
waf <- waf[which(!waf %in% c(wf, wbf))]
call <- rbind(call[wbf, ],
              data.frame(Arg = "f", Val = files),
              call[waf,])
wf <- which(call$Arg == "fasta")
wbf <- 1:(min(wf)-1)
waf <- 1:nrow(call)
waf <- waf[which(!waf %in% c(wf, wbf))]
call <- rbind(call[wbf, ],
              data.frame(Arg = "fasta", Val = fastas),
              call[waf,])
wf <- which(call$Arg == "var-mod")
wbf <- 1:(min(wf)-1)
waf <- 1:nrow(call)
waf <- waf[which(!waf %in% c(wf, wbf))]
call <- rbind(call[wbf, ],
              data.frame(Arg = "var-mod", Val = varmods),
              call[waf,])
Call <- paste0("diann.exe ", paste(paste0(" --", apply(call, 1, function(x) {
  x <- unlist(x)
  x <- x[which(x != "")]
  if (length(x) > 1) { x <- paste(x, collapse = " ") }
  return(x)
})), collapse = ""))
Call <- gsub(" --lib --", " --", Call)
Call <- gsub(" --fasta --", " --", Call)
tmp <- as.data.frame(t(sapply(strsplit(gsub("^var-mod UniMod:", "", varmods), ","), unlist)))
tmp <- apply(tmp, 1, function(x) {
  paste0("Modification ", x[[1]], " with mass delta ", x[[2]], " at ", x[[3]], " will be considered as variable")
})
Call <- c(Call, tmp)
Call <- c(Call, paste0("Output will be filtered at ", gsub("^qvalue +", "", min(FDRs)), " FDR"))
write(Call, paste0(dstDir, "/report.log.txt")) # Fake DiaNN log with just the minimum to allow the conversion function to work
#system(paste0("open \"", dstDir, "/report.log.txt"))

# Process reports and created single merged one
PSMs <- list()
for (i in 1:2) { #i <- 1 #i <- 2
  nm <- paste0("DiaNN", i)
  tmp <- reports[[nm]]
  tmp2 <- gsub("/[^/]+$", "", tmp)
  fixPths <- dirs[[nm]] != tmp2
  if (fixPths) {
    warning(paste0("Input ", i, ": paths have been modified or files moved since DiaNN was run!"))
    tmp <- paste0(dirs[[nm]], "/", gsub(".*/", "", tmp))
  }
  tmp2 <- as.data.frame(fread(tmp))
  #fls <- data.frame(Original = unique(tmp2$File.Name))
  #fls$Exist <- file.exists(fls$Original)
  PSMs[[nm]] <- list(Table = tmp2, PSMs_file = tmp)
}
lapply(PSMs, function(x) { unique(x$Table$File.Name) })
PSMsTbls <- lapply(PSMs, function(x) { x$Table })
allPSMs <- plyr::rbind.fill(PSMsTbls)
unique(allPSMs$File.Name)
tst <- aggregate(1:nrow(allPSMs), list(allPSMs$File.Name), length)
tst$Group.1 <- gsub(".*\\\\", "", tst$Group.1)
colnames(tst) <- c("File", "PSMs")
View(tst)
data.table::fwrite(allPSMs, paste0(dstDir, "/report.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write(c(paste0("1st folder = ", dir1),
        paste0("2nd folder = ", dir2),
        ""), paste0(dstDir, "/Original_DiaNN_folders.tsv"))
#tst <- read.delim(paste0(dstDir, "/report.tsv")); identical(allPSMs, tst) # Check that re-import produces the same file!
msg <- paste0("Done!\nNow you can treat folder:\n\t---> ", normalizePath(dstDir), "\nas a (minimal) DiaNN output folder for the purpose of running this package's data analysis scripts!\n")
cat(msg)
