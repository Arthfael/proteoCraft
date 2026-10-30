require(proteoCraft)


wd <- rstudioapi::selectDirectory(path = "...Search_Folder")
setwd(wd)

#fasta <- choose.files(paste0(wd, "/*.fasta"), multi = FALSE)
fasta <- rstudioapi::selectFile("Select fasta file",
                                path = paste0(defltdir, "/*"),
                                filter = "fasta file (*.fasta|*.fas|*.fa|*.faa|*.fasta.fas|*.txt)")

mzMLs <- grep("\\.mzML$", list.files(), value = TRUE)
MSGFplus <- "C:\\MSGFPlus_v20230112\\MSGFPlus.jar"
for (mzML in mzMLs) { #mzML <- mzMLs[2]
  cmd <- paste0('java -Xmx3500M -jar "', MSGFplus, '" -s "',
                normalizePath(wd), "\\\\", mzML, '" -d "', fasta,'" -inst 1 -t 20ppm -ti -1,2 -ntt 2 -tda 1 -o "',
                normalizePath(wd), "\\\\", sub("\\.mzML$", ".mzid", mzML),
                '" -conf "', gsub("/", "\\\\", wd), '" -maxMissedCleavages 2')
  cat(cmd)
  write(cmd, paste0(wd, "/MSGF command line.bat"))
  tst <- system(paste0('"', wd, '/MSGF command line.bat"'))
}

mz_IDs <- lapply(mzMLs, \(mzML) { mzID::mzID(sub("\\.mzML$", ".mzid", mzML)) })
names(mz_IDs) <- mzMLs
Evs <- setNames(lapply(names(mz_IDs), \(x) { #x <- names(mz_IDs)[1L]
  Ev <- mzID::evidence(mz_IDs[[x]])
  Ev$File <- x
  Ev$Sequence <- gsub("\\[?[+-][0-9]+", "", sub("^Pep_", "", Ev$peptide_ref))
  return(Ev)
}), mzMLs)
#View(Evs$`VLX_S1-DDA.mzML`)
PSMs <- setNames(lapply(names(mz_IDs), \(x) { #x <- names(mz_IDs)[1L]
  psms <- mzID::id(mz_IDs[[x]])
  psms$File <- x
  psms$Sequence <- gsub("\\[?[+-][0-9]+", "", sub("^Pep_", "", psms$peptide_ref))
  return(psms)
}), mzMLs)

# Plot
MeltEvs <- plyr::rbind.fill(Evs)
MeltPSMs <- plyr::rbind.fill(PSMs)
MeltPSMs$isdecoy <- MeltEvs$isdecoy[match(MeltPSMs$peptide_ref, MeltEvs$peptide_ref)]
ttl <- "-log10(E-value)"
plot <- ggplot2::ggplot(MeltPSMs) +
  ggplot2::geom_density(stat = "density",
                        ggplot2::aes(-log10(`ms-gf:evalue`), colour = isdecoy)) +
  ggplot2::facet_wrap(~File+passthreshold) +
  ggplot2::theme_bw() + ggplot2::ggtitle(ttl)
poplot(plot)
ttl <- "PSM-level P-value"
plot <- ggplot2::ggplot(MeltPSMs) +
  ggplot2::geom_density(stat = "density",
                        ggplot2::aes(`ms-gf:qvalue`, colour = isdecoy)) +
  ggplot2::facet_wrap(~File + passthreshold) +
  ggplot2::theme_bw() + ggplot2::ggtitle(ttl)
poplot(plot)
ttl <- "Peptide-level P-value"
plot <- ggplot2::ggplot(MeltPSMs) +
  ggplot2::geom_density(stat = "density",
                        ggplot2::aes(`ms-gf:pepqvalue`, colour = isdecoy)) +
  ggplot2::facet_wrap(~File + passthreshold) +
  ggplot2::theme_bw() + ggplot2::ggtitle(ttl)
poplot(plot)
ttl <- "Score"
plot <- ggplot2::ggplot(MeltPSMs) +
  ggplot2::geom_density(stat = "density",
                        ggplot2::aes(`ms-gf:rawscore`, colour = isdecoy)) +
  ggplot2::facet_wrap(~File + passthreshold) +
  ggplot2::theme_bw() + ggplot2::ggtitle(ttl)
poplot(plot)
# Conclusion: I need to filter, cannot take at face value the results!!!
Thr <- 0.05
FiltMeltPSMs <- MeltPSMs[which(MeltPSMs$`ms-gf:pepqvalue` <= Thr),]
FiltMeltEvs <- MeltEvs[which(MeltEvs$peptide_ref %in% FiltMeltPSMs$peptide_ref),]
ZPSMs <- FiltMeltPSMs[grep("Z", FiltMeltPSMs$Sequence),]
ZEvs <- FiltMeltEvs[grep("Z", FiltMeltEvs$Sequence),]
