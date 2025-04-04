require(parallel)
require(XML)
dSrcs <- c("...MS_File_Archive",
           "...Projects_Doc_Folder/Mass_Spec/Acquired data_v2",
           "D:/Data/Projects")

dSrcs <- dSrcs[which(dir.exists(dSrcs))]

# Create parallel processing cluster
N.clust <- detectCores()-1
a <- 1
tst <- try(clusterExport(parClust, "a", envir = environment()), silent = TRUE)
if ("try-error" %in% class(tst)) {
  try(stopCluster(parClust), silent = TRUE)
  parClust <- makeCluster(N.clust, type = "SOCK")
}

dFls <- parLapply(parClust, dSrcs, function(dr) { 
  x <- list.dirs(dr)
  grep("_[0-9]+\\.d$", x, value = TRUE)
})
names(dFls) <- dSrcs

dFls2 <- data.frame(File = unlist(dFls))
dFls2$"Bruker run" <- as.integer(gsub(".*_|\\.d$", "", dFls2$File))
dFls2$execLog <- parLapply(parClust, dFls2$File, function(dFl) { #dFl <- dDrs[1]
  sbdr <- grep("\\.m$", list.dirs(dFl), invert = TRUE, value = TRUE)
  sbdr <- sbdr[which(sbdr != dFl)]
  grep("execution-log", list.files(sbdr, full.names = TRUE), value = TRUE)
})
dFls2$execLog_ok <- sapply(dFls2$execLog, length) == 1
w <- which(dFls2$execLog_ok)
dFls2$startPress <- NA
dFls2$startPress[w] <- parSapply(parClust, dFls2$execLog[w], function(lgFl) { #lgFl <- dFls2$execLog[w[1]]
  lg <- readLines(unlist(lgFl))
  writeClipboard(lg)
  ... now I need to figure out which pressure to extract
  
  
})

# Get pressure profile
# Unfortunately not doable easily
# Step 1 - convert to mzML
deer <- list(ParsDir = "C:/Program Files/ThermoRawFileParser",
             RawConvertDir = "C:/Program Files/RawConverter_x64")
MSConvertInst <- ("C:/Program Files/ProteoWizard"%in% list.dirs("C:/Program Files", recursive = FALSE))
if (MSConvertInst) {
  MSConvertDirs <- grep("/ProteoWizard/ProteoWizard [^/]+", list.dirs("C:/Program Files/ProteoWizard", recursive = FALSE), value = TRUE)
  if (!length(MSConvertDirs)) { MSConvertInst <- FALSE } else {
    
    
    
    MSConvertVers <- as.data.frame(t(sapply(strsplit(gsub(".*/ProteoWizard ", "", MSConvertDirs), "\\."), unlist)))
    MSConvertVers <- MSConvertVers[order(MSConvertVers$V1, MSConvertVers$V2, MSConvertVers$V3, MSConvertVers$V4, decreasing = TRUE),]
    MSConvertVers <- MSConvertVers[1,]
    MSConvertDir <- paste0("C:/Program Files/ProteoWizard/ProteoWizard ", paste(MSConvertVers, collapse = "."))
    deer$MSConvertDir <- MSConvertDir
  }
}
write(fls, file = paste0(wd, "/tmp_MS_files.txt"))
precRecal <- FALSE
cmd <- paste0("\"", deer$MSConvertDir, "/msconvert.exe\" -f \"", wd, "/tmp_MS_files.txt\" -o \"",
              wd, "\" --mzML --64 -z --filter \"peakPicking true 1-\"")
#cat(cmd)
system(cmd)
unlink("tmp_MS_files.txt")
# Step 2 - extract pressure profile from mzMLs
mzMLs <- gsub("\\.raw$", ".mzML", fls, ignore.case = TRUE)
#tst <- grabMzmlData(mzMLs[1], "everything")
#View(tst$metadata)
pressScript <- paste0("C:/Users/", Sys.getenv("USERNAME"), "/Downloads/mzML_pressure_to_csv.py")
if (!file.exists(pressScript)) {
  download.file("https://gist.githubusercontent.com/caetera/0921b33f0c6201a538436906cc965fff/raw/d1af134fe228ce6a23be3e5bc3c49d20a8447ab2/mzML_pressure_to_csv.py",
                pressScript)
}
cmd <- paste0("python \"", pressScript, "\" ", paste0("\"", mzMLs, "\"", collapse = " "))
#cat(cmd)
system(cmd)
# Step 3 - read and plot
pressFls <- gsub("\\.mzML$", "_Pump_Pressure_1.csv", mzMLs)
press <- parLapply(parClust, pressFls, function(fl) {
  x <- read.csv(fl)
  fl <- gsub("_Pump_Pressure_1\\.csv$", ".raw", fl)
  colnames(x) <- c("Retention time", "Pressure")
  x$"Raw file" <- fl
  x$"Raw file name" <- basename(fl)
  x$"Raw file path" <- dirname(fl)
  return(x)
})
press <- plyr::rbind.fill(press)
press$`Raw file` <- factor(press$`Raw file`, levels = fls)
press$`Raw file name` <- factor(press$`Raw file name`, levels = basename(fls))
offset <- max(press$Pressure)/3
nc <- nchar(floor(offset))
offset <- ceiling(offset/10^(nc-1))*10^(nc-1)
press$Press <- press$Pressure + offset*(match(press$`Raw file`, fls)-1)
plot <- ggplot(press, aes(x = `Retention time`, y = Press, colour = `Raw file name`)) +
  geom_line(aes(group = `Raw file name`)) +
  theme_bw() #+ facet_grid(`Raw file name`~.)
proteoCraft::poplot(plot, 12, 22)
