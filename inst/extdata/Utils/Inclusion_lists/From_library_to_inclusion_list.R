# A script to predict retention time from iRT and write an inclusion list
#
# This script takes the following inputs:
# - DIA-NN tsv library
# - DIA-NN report from real identifications with the desired chromatographic gradient
# - Fasta file for filtering peptides in the list with proteins of interest
#
require(reshape2)
require(ggplot2)
require(proteoCraft)

wd <- rstudioapi::selectDirectory(path = "...Search_Folder")
setwd(wd)
InclDir <- paste0(wd, "/Inclusion lists")
if (!dir.exists(InclDir)) { dir.create(InclDir, recursive = TRUE) }

filt <- matrix(c("DIA-NN tsv library", "*.tsv"), ncol = 2)
LibFl <- choose.files(paste0(wd, "/*.tsv"), "Select a DIA-NN predicted library file (tsv format only!):", multi = FALSE, filt, 1)
#LibFl <- "predictedlib_tsv_CRF.tsv"

Library <- read.delim(LibFl)
#View(Library[1:100,])

PepLib <- aggregate(Library$Tr_recalibrated, list(Library$ModifiedPeptide), function(x) { list(unique(x)) })
if (max(sapply(PepLib$x, length)) > 1) {
  warning("Some peptides have multiple predicted retention times!")
  PepLib$x <- sapply(PepLib$x, mean)
} else { PepLib$x <- unlist(PepLib$x) }
colnames(PepLib) <- c("Modified sequence", "iRT")
tmp <- aggregate(Library$PrecursorCharge, list(Library$ModifiedPeptide), function(x) { list(unique(x)) })
PepLib$Charges <- tmp$x[match(PepLib$`Modified sequence`, tmp$Group.1)]
PepLib$Proteins <- Library$ProteinGroup[match(PepLib$`Modified sequence`, Library$ModifiedPeptide)]

filt <- matrix(c("DIA-NN report tsv file", "*.tsv"), ncol = 2)
ObsFl <- choose.files(paste0(wd, "/*.tsv"), "Select a DIA-NN report file containing peptide identifications from a sample with the desired gradient:", multi = FALSE, filt, 1)

Obs <- read.delim(ObsFl, check.names = FALSE)
w <- which(Obs$Modified.Sequence %in% PepLib$`Modified sequence`)

#if (!require("xcms", character.only = TRUE)) { BiocManager::install("xcms", update = FALSE) }
#require("xcms", character.only = TRUE)

Obs$iRT <- PepLib$iRT[match(Obs$Modified.Sequence, PepLib$`Modified sequence`)]
Obs <- Obs[which(!is.na(Obs$iRT)),]

# View overlapping data and fit model
fit1 <- lm(RT ~ poly(iRT, 3), Obs) # Alternatively, it may be possible to use a LOESS fit
xrange <- range(Obs$iRT)
xseq <- seq(from = xrange[1], to = xrange[2], length = 80)
pred <- predict(fit1, newdata = data.frame(iRT = xseq), se = TRUE)
y <- pred$fit
ci <- pred$se.fit*qt(0.95/2 + 0.5, pred$df)
ymin <- y - ci
ymax <- y + ci
lm.DF <- data.frame(iRT = xseq, RT = y, RT_min = ymin, RT_max = ymax, SE = pred$se.fit)
ttl1 <- "iRT vs real RT"
plot1 <- ggplot(Obs, aes(x = iRT, y = RT)) + geom_point(size = 0.1, alpha = 0.1) +
  geom_smooth(data = lm.DF, aes(x = iRT, y = RT, ymin = RT_min, ymax = RT_min), stat="identity") +
  theme_bw() + ggtitle(ttl1, subtitle = "(3rd degree polynomial fit)")
poplot(plot1)
ggsave(paste0(InclDir, "/", ttl1, ".jpeg"), plot1, dpi = 300)

# Predict real RT for library peptides
tmp <- predict(fit1, newdata = PepLib[, "iRT", drop = FALSE], se = TRUE)
PepLib$"predicted RT" <- tmp$fit
w <- which(PepLib$`Modified sequence` %in% Obs$Modified.Sequence)
PepLib$Error <- NA
PepLib$Error[w] <- PepLib$`predicted RT`[w] - PepLib$iRT[w]
plot1a <- ggplot(PepLib[w,], aes(x = iRT, y = `predicted RT`), color = "red") + geom_point(size = 0.1, alpha = 0.1) +
  theme_bw() + ggtitle(ttl1, subtitle = "(3rd degree polynomial fit)")
#poplot(plot1a)
Obs$Error <- Obs$RT - PepLib$`predicted RT`[match(Obs$Modified.Sequence, PepLib$`Modified sequence`)]
ttl2 <- "Estimation of prediction error"
plot2 <- ggplot(Obs) + geom_density(stat = "density", aes(x = Error)) +
  theme_bw() + ggtitle(ttl2)
poplot(plot2)
ggsave(paste0(InclDir, "/", ttl2, ".jpeg"), plot2, dpi = 300)

# Filter library for protein groups
filt <- matrix(c("FASTA file", "*.fasta"), ncol = 2)
FiltFastaFl <- choose.files(paste0(wd, "/*.fasta"), "Select a fasta file of proteins for which you want to write the inclusion list", multi = FALSE, filt, 1)
FiltFastaFl <- normalizePath(FiltFastaFl, winslash = "/")
FiltFasta <- Format.DB(FiltFastaFl)
tst <- which(sapply(strsplit(PepLib$Proteins, ";"), function(x) { sum(x %in% FiltFasta$`Protein ID`) }) > 0)
if (!length(tst)) {
  warning("Not a single match to proteins of interest found! No inclusion list will be written!")
} else {
  FiltPepLib <- PepLib[,]
  FiltPepLib <- FiltPepLib[which(FiltPepLib$`predicted RT` > 0),]
  #tmp <- melt(setNames(FiltPepLib$Charges, 1:nrow(FiltPepLib)))
  #colnames(tmp) <- c("Charge", "row")
  tmp <- data.frame(Charge = sapply(FiltPepLib$Charges, min), row = 1:nrow(FiltPepLib))
  tmp$row <- as.integer(tmp$row)
  tmp <- tmp[which(tmp$Charge >= 2),]
  tmp <- tmp[which(tmp$Charge <= 3),]
  kol <- c("Modified sequence", "iRT", "predicted RT", "Proteins")
  tmp[, kol] <- FiltPepLib[tmp$row,  kol]
  tmp$row <- NULL
  FiltPepLib <- tmp
  w <- which(Library$ModifiedPeptide %in% FiltPepLib$`Modified sequence`)
  tmp <- aggregate(Library$PrecursorMz[w], list(Library$ModifiedPeptide[w], Library$PrecursorCharge[w]), function(x) { list(unique(x)) })
  if (max(sapply(tmp$x, length)) > 1) {
    warning("Some peptides have multiple predicted retention times!")
    tmp$x <- sapply(tmp$x, mean)
  } else { tmp$x <- unlist(tmp$x) }
  FiltPepLib$"M/Z" <- apply(FiltPepLib[, c("Modified sequence", "Charge")], 1, function(x) {
    tmp$x[which((tmp$Group.1 == x[[1]])&(tmp$Group.2 == x[[2]]))]
  })
  FiltPepLib$Sequence <- gsub("\\([^\\)]+\\)", "", FiltPepLib$`Modified sequence`)
  nc <- nchar(FiltPepLib$Sequence)
  FiltPepLib <- FiltPepLib[which(nc >= 7),]
  FiltPepLib <- FiltPepLib[which(nc <= 25),]
  FiltPepLib2 <- FiltPepLib[which(FiltPepLib$`Modified sequence` == FiltPepLib$Sequence),]
  
  # Write inclusion list
  inclLst <- data.frame("Mass [m/z]" = FiltPepLib2$`M/Z`,
                        "Formula [M]" = "",
                        "Species" = "",
                        "CS [z]" = FiltPepLib2$Charge,
                        "Polarity" = "Positive",
                        "Start [min]" = FiltPepLib2$`predicted RT`-2.5,
                        "End [min]" = FiltPepLib2$`predicted RT`+2.5,
                        "NCE" = 28,
                        "Comment" = FiltPepLib2$`Modified sequence`,
                        check.names = FALSE)
  write.csv(inclLst, paste0(InclDir, "/Incl._list_", gsub(":", "-", Sys.time()), ".csv"), row.names = FALSE, na = "")
}

