options(stringsAsFactors = FALSE)
setwd("...Search_Folder/Test_PRM/combined/txt")
require(reshape2)
require(magrittr)
require(proteoCraft)
require(ggplot2)

db <- Format.DB("D:/Fasta_databases/Homo_sapiens_-_uniprot-proteome_UP000005640 - 20210413.fasta")

RTWindow <- 2 #min; this one is the length around previously observed RT where we will try to look for the precursor
RTWin2 <- 5 #sec; this one is the scale at which we look for overlaps - to optimize the number of precursors to include
PrecPerProt <- 2 # Aiming for this number of precursors per protein of interest; these will be included even if there end up being more overlapping than allowed by MaxPrec
MaxPrec <- 10 # Maximum precursors per unit of time for the second pass (expanding list once PrecPerProt has been attained)
MQ.load(pep = FALSE)
Zmin <- 2
Zmax <- 3

MQ.load(pep = FALSE)

raw <- unique(ev$`Raw file`)
aggregate(ev$`Modified sequence`, list(ev$`Raw file`), function(x) { length(unique(x)) })
ev <- ev[grep("^general_Full_dd", ev$`Raw file`),]
ev <- ev[which(ev$`Potential contaminant` == ""),]
ev <- ev[which(ev$Reverse == ""),]
ev <- ev[which(ev$Intensity > 0),]
ev <- ev[which((ev$Charge >= Zmin)&(ev$Charge <= Zmax)),]
pep <- set_colnames(aggregate(ev$Intensity, list(ev$`Modified sequence`), sum), c("Modified sequence", "Intensity"))
pep <- pep[order(pep$Intensity, decreasing = TRUE),]
pep$id <- 1:nrow(pep) # Since peptides are sorted by intensity, sorting by IDs is the same as sorting by intensity, which is convenient later down
pep[, c("Mass [m/z]", "CS [z]")] <- ev[match(pep$`Modified sequence`, ev$`Modified sequence`), c("m/z", "Charge")]
temp <- aggregate(ev$id, list(ev$`Modified sequence`), paste, collapse = ";")
pep$"Evidence IDs" <- temp$x
temp <- aggregate(ev$`Retention time`, list(ev$`Modified sequence`), mean)
pep$"Retention time" <- temp$x
pep$"Start [min]" <- pep$"Retention time"-RTWindow/2
pep$"End [min]" <- pep$"Retention time"+RTWindow/2
pep$Proteins <- ev$Proteins[match(pep$`Modified sequence`, ev$`Modified sequence`)]
# Filter pep to include only unique peptides
pep <- pep[grep(";", pep$Proteins, invert = TRUE),]
ev <- ev[which(ev$`Modified sequence` %in% pep$`Modified sequence`),]

temp <- melt(setNames(strsplit(PG$`Evidence IDs`, ";"), PG$id))
temp <- temp[which(temp$value %in% ev$id),]
PG <- PG[which(PG$id %in% temp$L1),]
PG$`Evidence IDs` <- lapply(strsplit(PG$`Evidence IDs`, ";"), function(x) { x[which(x %in% ev$id)] })
PG$`Peptide IDs` <- lapply(PG$`Evidence IDs`, function(x) { #x <- unlist(PG$`Evidence IDs`[1])
  sq <- ev$`Modified sequence`[match(x, ev$id)]
  x <- sort(pep$id[match(sq, pep$`Modified sequence`)], decreasing = FALSE)
  return(setNames(x, 1:length(x)))
})
temp <- melt(setNames(PG$"Peptide IDs", PG$id))
temp <- aggregate(temp$L1, list(temp$value), list)
pep$"Protein group IDs" <- temp$x[match(pep$id, temp$Group.1)]

pep$Include <- FALSE
pep$Step <- 0

protList <- sapply(strsplit(PG$`Protein IDs`, ";")[1:100], function(x) { x[[1]] }) #This could be any other list
protList2 <- db$`Full ID`[match(protList, db$`Protein ID`)]
#writeClipboard(protList2) #(for filtering in Skyline)
gr <- grsep2(protList, PG$`Protein IDs`)
# process:
# - up to PrecPerProt peptides for proteins of interest
# - add peptides from those proteins until running out of peptides/space
# - add peptides from other proteins until running out of peptides/space
#
### Step 1: get up to PrecPerProt precursors per target protein
temp <- unique(unlist(sapply(PG$`Peptide IDs`[gr], function(x) { x[1:min(c(length(x), PrecPerProt))] })))
pep$Include[match(temp, pep$id)] <- TRUE
w <- which(pep$Include)
tst <- data.frame(Time = (1:(ceiling(max(pep$`Retention time`))*60/RTWin2))*RTWin2/60)
tst$Targets <- sapply(tst$Time, function(x) { sum((pep$`Start [min]`[w] <= x)&(pep$`End [min]`[w] >= x)) })
tst$Targets_Step1 <- tst$Targets
pep$Step[w] <- 1
### Step 2: add precursors from the same proteins so long that they do not cause us to exceed the max number of precursors per second
temp <- unique(unlist(PG$`Peptide IDs`[gr]))
temp <- sort(temp[which(!pep$Include[match(temp, pep$id)])], decreasing = FALSE)
tst$Targets_Step2 <- 0
for (tmp in temp) { #tmp <- temp[1]
  mtch <- match(tmp, pep$id)
  w1 <- which((tst$Time >= pep$`Start [min]`[mtch])&(tst$Time <= pep$`End [min]`[mtch]))
  tmp2 <- tst$Targets[w1]
  if (max(tmp2) < MaxPrec) {
    tst$Targets_Step2[w1] <- tst$Targets_Step2[w1]+1
    tst$Targets[w1] <- apply(tst[w1, paste0("Targets_Step", 1:2)], 1, sum)
    pep$Include[mtch] <- TRUE
    pep$Step[mtch] <- 2
  }
}
### Step 3: add precursors from other proteins so long that they do not cause us to exceed the max number of precursors per second
temp <- sort(pep$id[which(!pep$Include)], decreasing = FALSE)
tst$Targets_Step3 <- 0
for (tmp in temp) { #tmp <- temp[1]
  mtch <- match(tmp, pep$id)
  w1 <- which((tst$Time >= pep$`Start [min]`[mtch])&(tst$Time <= pep$`End [min]`[mtch]))
  tmp2 <- tst$Targets[w1]
  if (max(tmp2) < MaxPrec) {
    tst$Targets_Step3[w1] <- tst$Targets_Step3[w1]+1
    tst$Targets[w1] <- apply(tst[w1, paste0("Targets_Step", 1:3)], 1, sum)
    pep$Include[mtch] <- TRUE
    pep$Step[mtch] <- 3
  }
}
###
tst <- tst[which(tst$Targets > 0),]
tst2 <- melt(tst[, c("Time", paste0("Targets_Step", 1:3))], id.vars = "Time")
tst2$Step <- factor(as.numeric(gsub("^Targets_Step", "", tst2$variable)), levels = 3:1)
plot <- ggplot(tst2) + geom_bar(stat = "identity", aes(x = Time, y = value, fill = Step)) + theme_bw()
poplot(plot)

InclList <- pep[which(pep$Include),]
InclList <- InclList[order(InclList$Step, decreasing = FALSE),]
InclList$Polarity <- "Positive"
InclList$Polarity <- "Positive"
InclList$"Formula [M]" <- NA
InclList$"Species" <- NA
InclList$NCE <- 27
InclList$Comment <- InclList$`Modified sequence`
InclList <- InclList[, c("Mass [m/z]", "Formula [M]", "Species", "CS [z]", "Polarity", "Start [min]", "End [min]", "NCE", "Comment")]
write.csv(InclList, paste0("Incl. list ", gsub(":", "-", Sys.time()), ".csv"), row.names = FALSE, na = "")
openwd()

