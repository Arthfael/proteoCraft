if (!require("BiocManager", quietly = TRUE)) { install.packages("BiocManager") }
if (!require("mzID", quietly = TRUE)) { BiocManager::install("mzID") }
require(mzID)
require(rstudioapi)

wd <- selectDirectory(path = "...Search_Folder")
setwd(wd)

fasta <- choose.files(paste0(wd, "/*.fasta"), multi = FALSE)
mzMLs <- grep("\\.mzML$", list.files(), value = TRUE)
MSGFplus <- "C:\\MSGFPlus_v20230112\\MSGFPlus.jar"
for (mzML in mzMLs) { #mzML <- mzMLs[2]
  cmd <- paste0('java -Xmx3500M -jar "', MSGFplus, '" -s "',
                normalizePath(wd), "\\\\", mzML, '" -d "', fasta,'" -inst 1 -t 20ppm -ti -1,2 -ntt 2 -tda 1 -o "',
                normalizePath(wd), "\\\\", gsub("\\.mzML$", ".mzid", mzML),
                '" -conf "', gsub("/", "\\\\", wd), '" -maxMissedCleavages 2')
  cat(cmd)
  write(cmd, paste0(wd, "/MSGF command line.bat"))
  tst <- system(paste0('"', wd, '/MSGF command line.bat"'))
}

mzIDs <- lapply(mzMLs, function(mzML) { mzID(gsub("\\.mzML$", ".mzid", mzML)) })
names(mzIDs) <- mzMLs
Evs <- setNames(lapply(names(mzIDs), function(x) { #x <- names(mzIDs)[1]
  Ev <- evidence(mzIDs[[x]])
  Ev$File <- x
  Ev$Sequence <- gsub("\\[?[+-][0-9]+", "", gsub("^Pep_", "", Ev$peptide_ref))
  return(Ev)
}), mzMLs)
#View(Evs$`VLX_S1-DDA.mzML`)
PSMs <- setNames(lapply(names(mzIDs), function(x) { #x <- names(mzIDs)[1]
  psms <- id(mzIDs[[x]])
  psms$File <- x
  psms$Sequence <- gsub("\\[?[+-][0-9]+", "", gsub("^Pep_", "", psms$peptide_ref))
  return(psms)
}), mzMLs)

# Plot
MeltEvs <- plyr::rbind.fill(Evs)
MeltPSMs <- plyr::rbind.fill(PSMs)
MeltPSMs$isdecoy <- MeltEvs$isdecoy[match(MeltPSMs$peptide_ref, MeltEvs$peptide_ref)]
require(ggplot2)
ttl <- "-log10(E-value)"
plot <- ggplot(MeltPSMs) +
  geom_density(stat = "density", aes(-log10(`ms-gf:evalue`), colour = isdecoy)) +
  facet_wrap(~File+passthreshold) + theme_bw() + ggtitle(ttl)
proteoCraft::poplot(plot)
ttl <- "PSM-level P-value"
plot <- ggplot(MeltPSMs) +
  geom_density(stat = "density", aes(`ms-gf:qvalue`, colour = isdecoy)) +
  facet_wrap(~File+passthreshold) + theme_bw() + ggtitle(ttl)
proteoCraft::poplot(plot)
ttl <- "Peptide-level P-value"
plot <- ggplot(MeltPSMs) +
  geom_density(stat = "density", aes(`ms-gf:pepqvalue`, colour = isdecoy)) +
  facet_wrap(~File+passthreshold) + theme_bw() + ggtitle(ttl)
proteoCraft::poplot(plot)
ttl <- "Score"
plot <- ggplot(MeltPSMs) +
  geom_density(stat = "density", aes(`ms-gf:rawscore`, colour = isdecoy)) +
  facet_wrap(~File+passthreshold) + theme_bw() + ggtitle(ttl)
proteoCraft::poplot(plot)
# Conclusion: I need to filter, cannot take at face value the results!!!
Thr <- 0.05
FiltMeltPSMs <- MeltPSMs[which(MeltPSMs$`ms-gf:pepqvalue` <= Thr),]
FiltMeltEvs <- MeltEvs[which(MeltEvs$peptide_ref %in% FiltMeltPSMs$peptide_ref),]
ZPSMs <- FiltMeltPSMs[grep("Z", FiltMeltPSMs$Sequence),]
ZEvs <- FiltMeltEvs[grep("Z", FiltMeltEvs$Sequence),]
