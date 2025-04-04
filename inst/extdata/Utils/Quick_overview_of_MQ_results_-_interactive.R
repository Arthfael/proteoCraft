# Quick overview of MQ results
options(stringsAsFactors = FALSE)

library(gtools)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(parallel)
library(rstudioapi)
library(magrittr)
library(svDialogs)
library(plotly)
library(htmlwidgets)
library(tibble)
library(openxlsx)
library(proteoCraft)

wd <- selectDirectory(path = "...Search_Folder")
Sys.sleep(1) # Required because sometimes R will attempt to set the word directory too quickly and fail.
setwd(wd)

#load("Backup.RData")
AnalysisParam <- data.frame(Folder = wd, check.names = FALSE)

# Type of experiment
msg <- "Select a type of experiment:\n - 1: Discovery -> no comparisons \n - 2: Regulation -> ratio analysis (up and down)\n - 3: Pull-down (incl. BioID) -> ratio analysis (choice between up only or up and down)"
AnalysisType <- as.numeric(dlg_input(msg)$res)
while(!AnalysisType %in% 1:3) { AnalysisType <- as.numeric(dlg_input(msg)$res) }
AnalysisType <- c("Discovery", "Regulation", "Pull-down")[AnalysisType]
AnalysisParam$"Type" <- AnalysisType
MakeRatios <- FALSE
if (AnalysisType %in% c("Regulation", "Pull-down")) {
  filt <- matrix(c("Samples map", "*.csv"), ncol = 2)
  msg <- "Select csv file with at least 2 columns:\n - \"MQ.exp\": experiment names used for the MaxQuant search\n - \"New name\": new sample names\nOther optional columns (required for ratio analysis):\n - \"Reference\": is this a reference sample for ratios calculations?\n - \"Ratios group\": which ratios group is this sample part of?"
  SamplesMap_file <- choose.files(paste0(wd, "/SamplesMap.csv"), msg,
                                  multi = FALSE, filter = filt)
  SamplesMap <- read.delim(SamplesMap_file, check.names = FALSE)
  if (ncol(SamplesMap) == 1) { SamplesMap <- read.delim(SamplesMap_file, check.names = FALSE, sep = ",") }
  if (ncol(SamplesMap) == 1) { SamplesMap <- read.csv(SamplesMap_file, check.names = FALSE) }
  if (ncol(SamplesMap) == 1) { SamplesMap <- read.csv(SamplesMap_file, check.names = FALSE, sep = "\t") }
  if (nrow(SamplesMap) == 1) {
    warning("Only one sample, skipping ratios analysis!")
  } else {
    w <- which(!c("Ratios group", "Reference") %in% colnames(SamplesMap))
    if (length(w)) {
      warning(paste0("Samples map is missing column", c("", "s")[length(w)], " ", paste(c("\"Ratios group\"", "\"Reference\"")[w], collapse = " and "), "!"))
    } else { MakeRatios <- TRUE }
  }
  AnalysisParam$"Ratios analysis" <- MakeRatios
  if (MakeRatios) {
    # In absence of replicates we can only work with arbitrary thresholds:
    msg <- "Enter fold change threshold value (log2, absolute)"
    RatiosThresh <- as.numeric(dlg_input(msg, 1)$res)
    while((!is.numeric(RatiosThresh))||(RatiosThresh < 0)||(is.na(RatiosThresh))) {
      RatiosThresh <- as.numeric(dlg_input(msg, 1)$res)
    }
    AnalysisParam$"Ratios analysis - threshold" <- RatiosThresh
    if (AnalysisType == "Pull-down") {
      RatiosThresh_2sided <- c(TRUE, FALSE)[match(dlg_message("Are you interested in down-regulated protein groups as well?", "yesno")$res, c("yes", "no"))]
    } else { RatiosThresh_2sided <- TRUE}
    AnalysisParam$"Ratios analysis - threshold is two-sided" <- RatiosThresh_2sided
  }
}
if (AnalysisType == "Pull-down") {
  IsBioID <- c(TRUE, FALSE)[match(dlg_message("Is this a BioID (and variants) experiment?", "yesno")$res, c("yes", "no"))]
  if (IsBioID) { AnalysisParam$"Type - advanced" <- "BioID" }
} else { IsBioID <- FALSE }
# Labelling
LabelType <- c("LFQ", "Isobaric")[as.numeric(dlgInput("Are samples Label-free (1) or Isobarically-labelled (2)?", 1)$res)]
if (LabelType == "Isobaric") {
  Labels <- gsub("^Reporter\\.intensity\\.", "", grep("^Reporter\\.intensity\\.", colnames(ev), value = TRUE))
}
AnalysisParam$"Label type" <- LabelType
# Number of threads
N.clust <- as.numeric(dlg_input("How many vCPUs am I allowed to use for this analysis?", 5)$res)
if (N.clust > detectCores()-1) {
  warning("You cannot specify more vCPUs than available, reducing number to max -1!")
  N.clust <- min(N.clust, detectCores()-1)
}
AnalysisParam$"N. of threads" <- N.clust
# Custom protein groups
custPGs <- c(TRUE, FALSE)[match(dlg_message("Do you want to provide a table of a priori known protein groups?", "yesno")$res, c("yes", "no"))]
if (custPGs) {
  filt <- matrix(c("Custom protein groups file", "*.csv"), ncol = 2)
  custPGs_file <- normalizePath(choose.files("Custom_PGs.csv", filters = filt), winslash = "/")
  custPGs <- read.delim(custPGs_file, check.names = FALSE)
  if (colnames(custPGs)[1] != "Leading protein IDs") { custPGs <- read.delim(custPGs_file, check.names = FALSE, sep = ",") }
  if (colnames(custPGs)[1] != "Leading protein IDs") { custPGs <- read.csv(custPGs_file, check.names = FALSE) }
  if (colnames(custPGs)[1] != "Leading protein IDs") { custPGs <- read.csv(custPGs_file, check.names = FALSE, sep = "\t") }
  if (colnames(custPGs)[1] != "Leading protein IDs") {
    warning("I could not make sense of that file! Skipping.")
    custPGs <- NA
  }
} else { custPGs <- NA }
#
# Min. number of peptidoforms for discovery/quantitation
NPep <- as.numeric(dlg_input("How many peptidoforms should a protein have at least to be quantified and reported?",
                             c(2, 2, 1)[match(AnalysisType, c("Discovery", "Regulation", "Pull-down"))])$res)
AnalysisParam$"N. of peptidoforms for quantitation" <- NPep
# Classes of peptides eligible for quantitation
pepclasses <- c("Unique peptide IDs", "Razor peptide IDs", "Peptide IDs")
msg <- paste0("Select peptides eligible for quantitation:\n", paste(" - ", 1:3, ": ", pepclasses, collapse = "\n"))
Pep4Quant <- as.numeric(dlg_input(msg, 2)$res)
while (!Pep4Quant %in% 1:3) { Pep4Quant <- as.numeric(dlg_input(msg, 2)$res) }
Pep4Quant <- pepclasses[Pep4Quant]
AnalysisParam$"Peptides eligible for quantitation" <- Pep4Quant
# Normalize?
NormalizePG <- c(TRUE, FALSE)[match(dlg_message("Do you want to normalize protein groups-level quantitative values?", "yesno")$res, c("yes", "no"))]
AnalysisParam$"Normalisation - Protein Groups" <- NormalizePG
#
# Proteome ruler
msg <- "Do you want to calculate Proteome Ruler estimates of copy number per cell? (y/n)"
def <- c("y", "y", "n")[match(AnalysisType, c("Discovery", "Regulation", "Pull-down"))]
protrul <- dlg_input(msg, def)$res
while (!protrul %in% c("y", "n")) { protrul <- dlg_input(msg, def)$res }
protrul <- c(TRUE, FALSE)[match(protrul, c("y", "n"))]
AnalysisParam$"Proteome ruler calculated" <- protrul
# Use match-between-runs peptides for coverage?
removeMBR <- c(TRUE, FALSE)[match(dlg_message("Should we remove Match-Between-Runs identifications from the coverage analysis?", "yesno")$res, c("yes", "no"))]
AnalysisParam$"Proteins list: remove match-between-runs" <- removeMBR
# Fasta databases
SSH_on <- FALSE
fastas <- readLines("parameters.txt")
fastas <- grep("\\.fasta$", fastas, value = TRUE, ignore.case = TRUE)
fastas <- unlist(strsplit(gsub("\\\\", "/", gsub("^Fasta file\t", "", fastas)), ";"))
fastas <- data.frame(Full = fastas, Name = basename(fastas), Dir = dirname(fastas))
fastas$Loc <- c("Local", "Cluster")[sapply(fastas$Dir, function(x) { grepl("^/nfs/", x) }) + 1 ]
fastas$Exists <- file.exists(fastas$Full)
fastas$ExistsHere <- file.exists(fastas$Name)
w1 <- which((fastas$ExistsHere)&(fastas$Loc == "Cluster"))
if (length(w1)) {
  require(tools)
  if (!SSH_on) {
    require(ssh)
    print("SSH session required")
    sshsess <- ssh_connect("bjoern22")
    while (class(sshsess) != "ssh_session") { sshsess <- ssh_connect("bjoern22") }
    #print(sshsess)
    SSH_on <- TRUE
  }
  tst1 <- md5sum(fastas$Name[w1])
  tst2 <- sapply(fastas$Full[w1], function(x) { #x <- fastas$Full[w1]
    gsub(" .*", "", capture.output(ssh_exec_wait(sshsess, paste0("md5sum ", x)))[1])
  })
  w1y <- w1[which(tst1 == tst2)]
  w1n <- w1[which(tst1 != tst2)]
  fastas$Exists[w1n] <- FALSE
  fastas$ExistsHere[w1y] <- fastas$Exists[w1y] <- TRUE
  fastas$Full[w1y] <- paste0(wd, "/", fastas$Name[w1y])
  fastas$Dir[w1y] <- wd
  fastas$Loc[w1y] <- "Local"
}
w2 <- which((!fastas$Exists)&(fastas$Loc == "Local"))
if (length(w2)) {
  tst <- (length(w2) == 1)+1
  msg <- paste0(c("Several", "One")[tst], " Fasta file", c("s", "")[tst], " used to search the data could not be located:\n",
                paste(paste0(" - \"", fastas$Full[w2], "\""), collapse = "\n"))
  stop(msg)
}
w3 <- which((!fastas$ExistsHere)&(fastas$Exists))
if (length(w3)) {
  file.copy(fastas$Full[w3], to = wd)
  fastas$ExistsHere[w3] <- TRUE
  fastas$Full[w3] <- paste0(wd, "/", fastas$Name[w3])
  fastas$Dir[w3] <- wd
}
w4 <- which((!fastas$ExistsHere)&(fastas$Loc == "Cluster"))
if (length(w4)) {
  if (!SSH_on) {
    require(ssh)
    print("I need to download locally some of the Fastas used for the search from the cluster, please open an SSH session.")
    sshsess <- ssh_connect("bjoern22")
    while (class(sshsess) != "ssh_session") { sshsess <- ssh_connect("bjoern22") }
    #print(sshsess)
    SSH_on <- TRUE
  }
  scp_download(sshsess, fastas$Full[w4], wd)
  fastas$ExistsHere[w4] <- TRUE
  fastas$Full[w4] <- paste0(wd, "/", fastas$Name[w4])
  fastas$Dir[w4] <- wd
}
fasta_types <- c("UNIPROTKB", "ENSEMBL", "REFSEQPROTEIN", "NCBI", "TAIR")
fastas$Data <- lapply(fastas$Full, readLines)
if (nrow(fastas) == 1) {
  mes <- "Of which type is the Fasta database?\n - 1: UniprotKB\n - 2: ENSEMBL\n - 3: REFSEQPROTEIN\n - 4: NCBI\n - 5: TAIR"
} else {
  mes <- paste0("Of which type is FASTA database \"", fastas$Full[1], "\"?\n - 1: UniprotKB\n - 2: ENSEMBL\n - 3: REFSEQPROTEIN\n - 4: NCBI\n - 5: TAIR")
}
fasta_type <- fasta_types[as.numeric(dlg_input(mes, 1)$res)]
db <- Format.DB(unlist(fastas$Data[[1]]), in.env = TRUE, mode = fasta_type)
if (nrow(fastas) > 1) {
  for (i in 2:nrow(fastas)) {
    mes <- paste0("Of which type is FASTA database \"", fastas$Full[i], "\"?\n - 1: UniprotKB\n - 2: ENSEMBL\n - 3: REFSEQPROTEIN\n - 4: NCBI\n - 5: TAIR")
    fasta_type <- fasta_types[as.numeric(dlg_input(mes, 1)$res)]
    tmp <- Format.DB(unlist(fastas$Data[[i]]), in.env = TRUE, mode = fasta_type)
    c1 <- colnames(tmp)[which(!colnames(tmp) %in% colnames(db))]
    c2 <- colnames(db)[which(!colnames(db) %in% colnames(tmp))]
    if (length(c1) > 0) { db[,c1] <- NA }
    if (length(c2) > 0) { tmp[,c2] <- NA }
    db <- rbind(db, tmp)
  }
}
for (i in 1:nrow(fastas)) { if (!file.exists(paste0(wd, "/", fastas$Name[i]))) { file.copy(fastas$Full[i], wd) } }
AnalysisParam$Fastas <- list(fastas$Full)
# Proteins of interest
protlist <- c(TRUE, FALSE)[match(dlg_message("Do you want to provide a list of proteins of interest?", "yesno")$res, c("yes", "no"))]
IDs.list <- list()
prot.names <- c()
protlistsource <- c()
if (protlist) {
  msg <- paste(c("Under which form will the proteins of interest be provided?", " - 1: fasta file(s)",
                 " - 2: Protein names list (txt file, with comma, semicolon or line separated individual names)",
                 " - 3: Manually enter protein names", " - 4: Manually enter protein accessions"), collapse = "\n")
  protlistsource <- as.numeric(dlg_input(msg, 1)$res)
  while (!protlistsource %in% c(1:4)) { protlistsource <- as.numeric(dlg_input(msg, 1)$res) }
  if (protlistsource == 1) {
    filt <- matrix(c("Fasta file", "*.fasta"), ncol = 2)
    tmp <- choose.files(default = c("*.fasta", "Proteins of interest.fasta")[file.exists("Proteins of interest.fasta")+1],
                        filters = filt)
    tmp <- setNames(lapply(tmp, Format.DB), paste0("Fasta #", 1:length(tmp)))
    tmp <- suppressMessages(melt(tmp))
    tmp$Sequence <- gsub("\\*$", "", tmp$Sequence)
    prot.names <- db$`Common Name`[unique(unlist(sapply(tmp$Sequence, function(x) { which(db$Sequence == x) })))]
  }
  if (protlistsource == 2) {
    filt <- matrix(c("Text file", "*.txt"), ncol = 2)
    tmp <- choose.files("Proteins of interest.txt", filters = filt, multi = FALSE)
    tmp <- unlist(strsplit(suppressWarnings(readLines(tmp)), " ?[,;] ?"))
  }
  if (protlistsource %in% 3:4) {
    tmp <- unlist(strsplit(dlg_input("Enter values (comma or semicolon separated)", filters = filt, multi = FALSE)$res, " ?[,;] ?"))
  }
  if (protlistsource %in% 2:3) { prot.names <- tmp[which(tmp %in% c(db$Name, db$`Common Name`))] }
  if (protlistsource == 4) { prot.names <- db$`Common Name`[match(tmp, db$"Protein ID")] }
  prot.names <- prot.names[which(!is.na(prot.names))]
  prot.names <- prot.names[which(nchar(prot.names) > 0)]
  w <- which(!prot.names %in% db$`Common Name`)
  if (length(w)) {
    prot.names[w] <- db$`Common Name`[which(db$Name %in% prot.names[w])]
  }
  prot.names <- unique(prot.names)
  IDs.list <- setNames(lapply(prot.names, function(x) { db$"Protein ID"[which(db$`Common Name` == x)] }), prot.names)
}
if (!is.na(list(Custom = custPGs))) {
  temp <- unlist(strsplit(custPGs$"Leading protein IDs", ";"))
  temp2 <- db$`Common Name`[match(temp, db$`Protein ID`)]
  temp2 <- setNames(lapply(unique(temp2), function(x) { temp[which(temp2 == x)] }), unique(temp2))
  w1 <- which(names(temp2) %in% names(IDs.list))
  w2 <- which(!names(temp2) %in% names(IDs.list))
  IDs.list[names(temp2)[w1]] <- temp2[w1]
  if (length(w2)) {
    IDs.list[names(temp2)[w2]] <- lapply(names(temp2)[w2], function(x) { unique(unlist(c(IDs.list[[x]], temp2[[x]]))) })
  }
  prot.names <- unique(c(prot.names, names(temp2)))
  protlistsource <- c(protlistsource, "Custom protein groups table")
}
protlist <- length(prot.names) > 0
AnalysisParam$"Proteins list: provided?" <- protlist
if (protlist) {
  temp <- db[which(db$`Protein ID` %in% unlist(IDs.list)),]
  temp2 <- rep("", nrow(temp)*3)
  temp2[(1:nrow(temp))*3-2] <- temp$Header
  temp2[(1:nrow(temp))*3-1] <- temp$Sequence
  write(temp2, "Proteins of interest.fasta")
  prot.IDs <- unlist(IDs.list)
  AnalysisParam$"Proteins list: protein" <- list(IDs.list)
  AnalysisParam$"Proteins list: names" <- list(prot.names)
  AnalysisParam$"Proteins list: source" <- protlistsource
}
# Annotate
setwd(wd)
Annotate <- c(TRUE, FALSE)[match(dlg_message("Do you want to annotate the database and protein groups tables?\n(This is required for GO terms enrichment analysis)", "yesno")$res, c("yes", "no"))]
if (Annotate) {
  if ("Parsed_annotations.RData" %in% list.files()) { load("Parsed_annotations.RData") } else {
    Parsed_annotations <- lapply(gsub("\\.fasta$", ".txt", fastas$Full), function(x) {
      if (file.exists(x)) { return(Format.DB_txt(readLines(x))) } else { return(NA) }
    })
    w <- which(lapply(Parsed_annotations, class) == "data.frame")
    if (!length(w)) { warning("No annotations files provided, skipping annotations!") } else {
      Parsed_annotations <- Parsed_annotations[w]
      temp1 <- Parsed_annotations[[1]]
      if (length(Parsed_annotations) > 1) {
        for (i in 2:length(Parsed_annotations)) {
          temp2 <- Parsed_annotations[[i]]
          kol1 <- colnames(temp1)
          kol2 <- colnames(temp2)
          kol1 <- kol1[which(!kol1 %in% kol2)]
          kol2 <- kol2[which(!kol2 %in% kol1)]
          temp1[,kol2] <- NA
          temp2[,kol1] <- NA
          temp1 <- rbind(temp1, temp2)
        }
      }
      Parsed_annotations <- temp1; rm(temp1)
      save(Parsed_annotations, file = "Parsed_annotations.RData")
    }
  }
}
Annotate <- exists("Parsed_annotations") # Update condition
AnalysisParam$Annotations <- Annotate
# GO terms enrichment analysis
if (Annotate) {
  GO_analysis <- c(TRUE, FALSE)[match(dlg_message("Do you want to perform GO terms enrichment analysis?", "yesno")$res, c("yes", "no"))]
  if (GO_analysis) {
    # Get pre-registered email for connection to DAVID web access
    DavidEmail <- dlg_input("Enter a valid email for DAVID web access\n(must be pre-registered)")$res
  }
} else { GO_analysis <- FALSE }
AnalysisParam$"GO terms enrichment analysis" <- GO_analysis
if (GO_analysis) { AnalysisParam$"GO terms enrichment analysis: DAVID email" <- DavidEmail }

# Load data
MQ.load(pep = FALSE, prot = FALSE)
ev$"Protein group IDs" <- NULL
ev$"Peptide ID" <- NULL
# Convert modified sequence column to old format + make PTMs table
temp <- cor_mod_seq(ev)
Modifs <- temp$PTMs
ev <- temp$Peptides

# Peptide modifications eligible for protein groups quantitation
Mod4Quant <- Modifs$Mark[grep("phospho|biot", Modifs$`Full name`, ignore.case = TRUE, invert = TRUE)] # This will have to be made interactive eventually
msg <- paste0("Which modifications are eligible for protein group quantitation?\n(enter comma-separated numbers, or 0 for all displayed)\n",
              paste(paste0(" - ", 1:length(Mod4Quant), ": ", Mod4Quant), collapse = "\n"))
temp <- suppressWarnings(as.numeric(unlist(strsplit(dlg_input(msg, 0)$res, ","))))
temp <- temp[which(!is.na(temp))]
while (!length(temp)) {
  temp <- suppressWarnings(as.numeric(unlist(strsplit(dlg_input(msg, 0)$res, ","))))
  temp <- temp[which(!is.na(temp))]
}
if ((length(temp) == 1)&&(temp == 0)) { temp <- 1:length(Mod4Quant) }
Mod4Quant <- Mod4Quant[temp]
AnalysisParam$"PTMs eligible for quantitation" <- list(Mod4Quant)

# Experiments
if ("Experiment" %in% colnames(ev)) {
  if (!exists("SamplesMap")) {
    filt <- matrix(c("Samples map", "*.csv"), ncol = 2)
    SamplesMap_file <- choose.files(paste0(wd, "/SamplesMap.csv"), "Select 2 column file, with two columns:\n - \"MQ.exp\": names used for the MaxQuant search\n - \"New name\": new sample names",
                                    multi = FALSE, filter = filt)
    while (!length(SamplesMap_file)) {
      SamplesMap_file <- choose.files(paste0(wd, "/SamplesMap.csv"), "Select 2 column file, with two columns:\n - \"MQ.exp\": names used for the MaxQuant search\n - \"New name\": new sample names",
                                      multi = FALSE, filter = filt)
    }
    SamplesMap <- read.delim(SamplesMap_file, check.names = FALSE)
    if (ncol(SamplesMap) == 1) { SamplesMap <- read.delim(SamplesMap_file, check.names = FALSE, sep = ",") }
    if (ncol(SamplesMap) == 1) { SamplesMap <- read.csv(SamplesMap_file, check.names = FALSE) }
    if (ncol(SamplesMap) == 1) { SamplesMap <- read.csv(SamplesMap_file, check.names = FALSE, sep = "\t") }
  }
  # Optional: rename samples:
  RenameSamples <- FALSE
  if ("New name" %in% colnames(SamplesMap)) {
    RenameSamples <- c(TRUE, FALSE)[match(dlg_message("Do you want to rename samples?", "yesno")$res, c("yes", "no"))]
    AnalysisParam$"Rename samples?" <- RenameSamples
    if (RenameSamples) {
      SamplesMap$Experiment <- SamplesMap$"New name"
      ev$Experiment <- SamplesMap$Experiment[match(ev$Experiment, SamplesMap$MQ.exp)]
    } else { SamplesMap$Experiment <- SamplesMap$MQ.exp }
  }
  SamplesMap <- SamplesMap[which(!is.na(SamplesMap$Reference)),]
  Exp <- SamplesMap$Experiment
  w1 <- which(ev$Experiment %in% Exp)
  w2 <- which(!ev$Experiment %in% Exp)
  if (length(w2)) {
    warning(paste0("Removing ", length(w2), " evidences from undefined experiments (check \"Reference\" column of Samples map)..."))
    ev <- ev[w1,]
  }
} else {
  Exp <- "Exp1"
  ev$Experiment <- Exp
  SamplesMap <- data.frame(Experiment = Exp, check.names = FALSE, "Ratios group" = 1, Reference = TRUE)
}
Venn_Obs <- Venn_Ratios <- FALSE
if (length(Exp) > 1) {
  # Venn diagrams
  Venn_Obs <- c(TRUE, FALSE)[match(dlg_message("Since there is more than 1 sample, we could create Venn diagrams of observed protein groups per sample. Should we?", "yesno")$res, c("yes", "no"))]
  AnalysisParam$"Venn diagrams: observed" <- Venn_Obs
  if (MakeRatios) {
    Venn_Ratios <- c(TRUE, FALSE)[match(dlg_message("Should we create Venn diagrams of regulated protein groups per sample?", "yesno")$res, c("yes", "no"))]
    AnalysisParam$"Venn diagrams: regulated (up)" <- Venn_Ratios
    AnalysisParam$"Venn diagrams: regulated (down)" <- Venn_Ratios & RatiosThresh_2sided
  }
}

# Filter ev
w <- which((is.na(ev$Reverse))|(ev$Reverse != "+"))
l <- length(which(!1:nrow(ev) %in% w))
if (l) {
  warning(paste0("Removing ", l, " reverse peptide evidences!"))
  revEv <- ev[-w,]
  ev <- ev[w, ]
}
w <- which((is.na(ev$"Potential contaminant"))|(ev$"Potential contaminant" != "+"))
l <- length(which(!1:nrow(ev) %in% w))
if (l) {
  warning(paste0("Removing ", l, " contaminant peptide evidences!"))
  contEv <- ev[-w, ]
  ev <- ev[w, ]
}
w <- which((is.all.good(ev$Intensity, 2))&(ev$Intensity > 0))
l <- length(which(!1:nrow(ev) %in% w))
if (l) {
  warning(paste0("Removing ", l, " peptide evidences with invalid or null intensity values!"))
  nullEv <- ev[-w, ]
  ev <- ev[w, ]
}
Exp <- Exp[which(Exp %in% ev$Experiment)] # Update experiments

# Create peptidoforms file
pep <- set_colnames(aggregate(ev$id, list(ev$"Modified sequence"), function(x) {paste(sort(x), collapse = ";")}),
                    c("Modified sequence", "Evidence IDs"))
# Create matches between ev and pep
ev_to_pep <- match(pep$"Modified sequence", ev$"Modified sequence") # Only to be used when all matches would return the same value, e.g. sequence
pep_to_ev <- match(ev$"Modified sequence", pep$"Modified sequence")# In the reverse direction, works for all
# Proteins column
pep[, c("Sequence", "Proteins")] <- ev[ev_to_pep,  c("Sequence", "Proteins")]
# PEP: for a peptide, the lowest of the PEP of individual matching evidences
temp <- aggregate(ev$PEP, list(ev$"Modified sequence"), function(x) {
  x <- x[which(!is.na(x))]
  if (length(x)) { return(min(x)) } else { return(NA) }
})
pep$PEP <-  temp$x[match(pep$"Modified sequence", temp$Group.1)]
# IDs
pep$id <- c(1:nrow(pep)) # Note: these are new IDs, so will not match the peptide IDs in other tables, e.g. protein groups or evidence
ev$"Mod. peptide ID" <- NULL
ev$"Peptide ID" <- pep$id[pep_to_ev]
# Amino acid counts
AA <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
for (aa in AA) {
  ev[[paste0(aa, " Count")]] <- nchar(ev$Sequence) - nchar(gsub(aa, "", ev$Sequence))
  pep[[paste0(aa, " Count")]] <- ev[ev_to_pep, paste0(aa, " Count")]
}
# Length
ev$Length <- nchar(ev$Sequence)
pep$Length <-  ev$Length[ev_to_pep]
# Intensity
temp <- aggregate(ev$Intensity, list(ev$"Modified sequence"), sum, na.rm = TRUE)
pep$Intensity <-  temp$x[match(pep$"Modified sequence", temp$Group.1)]
for (e in Exp) { #e <- Exp[1]
  temp <- ev[which(ev$Experiment == e),]
  if (LabelType == "Isobaric") { int.col <- paste0("Reporter intensity ", Labels) }
  if (LabelType == "LFQ") { int.col <- "Intensity" }
  temp <-  set_colnames(aggregate(temp[[int.col]], list(temp$"Modified sequence"), function(x) { sum(is.all.good(x)) }),
                        c("Modified sequence", int.col))
  pep[, paste0(int.col, " - ", e)] <- 0
  w <- which(pep$"Modified sequence" %in% temp$"Modified sequence")
  pep[w, paste0(int.col, " - ", e)] <- temp[match(pep$"Modified sequence"[w], temp$"Modified sequence"), int.col]
}

raw.files <- sort(unique(ev$"Raw file"))

## Peptides level:
# Intensity distribution:
temp <- pep[,c("Modified sequence", paste0(int.col, " - ", Exp))]
temp <- melt(temp, id.vars = "Modified sequence")
temp$variable <- gsub(topattern(paste0(int.col, " - ")), "", temp$variable)
temp$variable <- factor(temp$variable, levels = Exp)
temp$value <- log10(temp$value)
temp <- temp[which(is.all.good(temp$value, 2)),]
ttl <- "Density plot - Peptides level"
plot <- ggplot(temp) + geom_density(stat = "density", aes(x = value, colour = variable)) +
  facet_wrap(~variable, strip.position = "right") + ggtitle(ttl) +
  theme_bw() +
  xlab("log10(Peptides Intensity)") +
  ylab("Density")
poplot(plot)
dir <- "Workflow control"
if (!dir.exists(dir)) { dir.create(dir) }
setwd(paste0(wd, "/", dir))
ggsave(paste0(ttl, ".jpeg"), plot, dpi = 300, width = 10, height = 10, units = "in")
ggsave(paste0(ttl, ".pdf"), plot)
setwd(wd)
# Correlation:
if (length(Exp) > 1) {
  temp <- pep[,c("Modified sequence", paste0(int.col, " - ", Exp))]
  for (e in Exp) {
    temp[[paste0("log10(Intensity) - ", e)]] <- log10(temp[[paste0(int.col, " - ", e)]])
  }
  comb <- as.data.frame(combinations(length(Exp), 2, 1:length(Exp)))
  for (i in 1:2) { comb[,i] <- Exp[comb[,i]] } # Required to keep the original order of Exp
  temp2 <- temp[, c("Modified sequence", paste0("log10(Intensity) - ", comb[1,]))]
  temp2$X <- comb[1,1]
  temp2$Y <- comb[1,2]
  temp2$Comparison <- paste0(comb[1,1], " (X) vs ", comb[1,2], " (Y)")
  colnames(temp2)[2:3] <- c("log10 (X intensity)", "log10 (Y intensity)")
  if (length(Exp) > 2) {
    for (i in 2:nrow(comb)) {
      temp3 <- temp[, c("Modified sequence", paste0("log10(Intensity) - ", comb[i,]))]
      temp3$X <- comb[i,1]
      temp3$Y <- comb[i,2]
      temp3$Comparison <- paste0(comb[i,1], " (X) vs ", comb[i,2], " (Y)")
      colnames(temp3)[2:3] <- c("log10 (X intensity)", "log10 (Y intensity)")
      temp2 <- rbind(temp2, temp3)
    }
  }
  test <- apply(temp2[,c("log10 (X intensity)", "log10 (Y intensity)")], 1, function(x) { length(is.all.good(x)) }) == 2
  temp2 <- temp2[which(test),]
  temp2$X <- factor(temp2$X, levels = Exp)
  temp2$Y <- factor(temp2$Y, levels = Exp)
  temp3 <- as.data.frame(t(sapply(unique(temp2$Comparison), function(x) { #x <- unique(temp2$Comparison)[1]
    x1 <- temp2[which(temp2$Comparison == unlist(x)), c("log10 (X intensity)", "log10 (Y intensity)")]
    x1 <- x1$"log10 (Y intensity)"-x1$"log10 (X intensity)"
    return(setNames(c(x, paste0("Median = ", round(median(x1), 3)), paste0("S.D. = ", round(sd(x1), 3))),
                    c("Comparison", "Median", "SD")))
  })))
  temp3$X <- gsub(" \\(X\\).+", "", temp3$Comparison)
  temp3$Y <- gsub(".+\\(X\\) vs ", "", gsub(" \\(Y\\)$", "", temp3$Comparison))
  temp3$X <- factor(temp3$X, levels = Exp)
  temp3$Y <- factor(temp3$Y, levels = Exp)
  temp3$R <- sapply(strsplit(gsub(" \\([XY]\\)", "", temp3$Comparison), " vs "), function(x) {
    return(paste0("R = ", round(cor(pep[[paste0(int.col, " - ", x[[1]])]], pep[[paste0(int.col, " - ", x[[2]])]]), 3)))
  })
  x_min <- min(temp2$`log10 (X intensity)`)
  y_min <- min(temp2$`log10 (Y intensity)`)
  y_max <- max(temp2$`log10 (Y intensity)`)
  ttl <- "Correlation plot - Peptides level"
  plot <- ggplot(temp2) + geom_point(aes(x = `log10 (X intensity)`, y = `log10 (Y intensity)`), alpha = 0.25, size = 0.01) +
    geom_abline(intercept = 0, slope = 1, colour = "red") +
    geom_text(data = temp3, x = x_min, y = y_max - 0.01*(y_max - y_min), aes(label = R), hjust = 0, cex = 2.5) +
    geom_text(data = temp3, x = x_min, y = y_max - 0.06*(y_max - y_min), aes(label = Median), hjust = 0, cex = 2.5) +
    geom_text(data = temp3, x = x_min, y = y_max - 0.11*(y_max - y_min), aes(label = SD), hjust = 0, cex = 2.5) +
    facet_grid(Y~X) + coord_fixed(1) +
    theme_bw() + theme(strip.text.y = element_text(angle = 0, vjust = 0.5, hjust = 0),
                       strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0)) + ggtitle(ttl)
  poplot(plot)
  if (!dir.exists(dir)) { dir.create(dir) }
  setwd(paste0(wd, "/", dir))
  ggsave(paste0(ttl, ".jpeg"), plot, dpi = 300, width = 10, height = 10, units = "in")
  ggsave(paste0(ttl, ".pdf"), plot)
  setwd(wd)
}

save.image("Backup.RData")
#load("Backup.RData")

## Protein groups level
# Assemble protein groups
PG_assembly <- PG_assemble(pep, DB = db, Ev = ev, N.clust = N.clust, Custom_PGs = custPGs)
save(PG_assembly, file = "PG_assembly.RData")
#load("PG_assembly.RData")
PG <- PG_assembly$Protein.groups
pep <- PG_assembly$Peptides
db <- PG_assembly$Database
if ("Evidences" %in% names(PG_assembly)) { ev <- PG_assembly$Evidences }
w <- which(c("Organism_Full", "Organism") %in% colnames(db))
tstorg <- length(w) > 0
if (tstorg) {
  kol <- c("Organism_Full", "Organism")[w[1]]
  test <- sapply(strsplit(PG$`Protein IDs`, ";"), function(x) {
    paste(sort(unique(db[match(x, db$`Protein ID`), kol])), collapse = ";")
  })
  orgkol <- c("Organism", "Organism(s)")[(sum(grepl(";", test))>0)+1]
  PG[[orgkol]] <- test
}

# Some more columns
tmp <- strsplit(PG$"Leading protein IDs", ";")
for (i in c("No Isoforms", "Names", "Genes")) {
  if (i == "No Isoforms") { j <- i } else { j <- gsub("s$", "", i) }
  if (!j %in% colnames(db)) {
    j <- paste0(gsub("s$", "", i), c(" ID", " IDs", ""))
    w <- which(j %in% colnames(db))
    if (length(w)) { j <- j[w[1]] } else {
      warning(paste0("No near matching column name found for \"", i, "\" in the protein data base table."))
    }
  }
  if (length(j) == 1) {
    PG[[i]] <- sapply(tmp, function(x) {
      x <- unlist(x)
      m1 <- match(x, db$"Protein ID")
      if ("Full ID" %in% colnames(db)) {
        m1 <- data.frame(m1 = m1, m2 = match(x, db$"Full ID"))
        m1 <- apply(m1, 1, function(y) {
          y <- unique(y[which(!is.na(y))])
          if (!length(y)) { y <- "" }
          return(y)
        })
        m1 <- unlist(m1)
      }
      x <- db[m1, j]
      x[which(x %in% c("", " ", "NA", NA))] <- ""
      if (!length(x)) { x <- "" }
      if (i == "Genes") { x <- unique(x) }
      x <- x[which(x != "")]
      return(paste(x, collapse = ";"))
    })
  }
}

# Here, if this is a BioID type experiment, we also want to mark protein groups which have Biotin peptides:
if (IsBioID) {
  wbiot <- grep("biot", Modifs$"Full name", ignore.case = TRUE)
  l <- length(wbiot)
  if (length(wbiot)) {
    if (l == 1) { tmp <- Modifs$"Full name"[wbiot] } else { tmp <- paste0(paste(Modifs$"Full name"[wbiot[1:(l-1)]], collapse = "\", \""), "\" and \"", Modifs$"Full name"[wbiot[l]]) }
    warning(paste0("Modifications \"", tmp, "\" were detected as biotinylations, check that this is correct!"))
    AnalysisParam$"Type - advanced (mod)" <- list(Modifs$"Full name"[wbiot])
  } else {
    warning("I could not identify any biotinylated PTMs in the modifications table, did you include them in the search?")
    IsBioID <- FALSE
    AnalysisParam$"Type - advanced" <- "BioID... - except I could not find biotinylations among the PTMs!"
  }
}
if (IsBioID) {
  g <- grep(topattern(Modifs$Mark[wbiot], start = FALSE), pep$"Modified sequence")
  if (length(g)) {
    wpg <- unique(unlist(strsplit(pep$"Protein group IDs"[g], ";")))
    wpg <- which(PG$id %in% wpg)
    PG[["Biot. peptide IDs"]] <- ""
    temp <- setNames(strsplit(PG$"Peptide IDs", ";"), PG$id)
    temp <- melt(temp)
    temp <- temp[which(temp$value %in% pep$id[g]),]
    temp <- aggregate(temp$value, list(temp$L1), function(x) { paste(sort(x), collapse = ";") })
    PG[wpg, "Biot. peptide IDs"] <- temp$x[match(PG$id[wpg], temp$Group.1)]
    PG[["Biot. peptides count"]] <- sapply(strsplit(PG[["Biot. peptide IDs"]], ";"), length)
    PG[["Biot. peptides [%]"]] <- round(100*PG[["Biot. peptides count"]]/PG$"Peptides count", 1)
  } else {
    warning("I could not find any biotinylated peptides!")
    IsBioID <- FALSE
    AnalysisParam$"Type - advanced" <- "BioID... - except that no biotinylated peptides were detected!"
  }
}

# Number of evidences and peptides per sample:
temp_PG <- data.frame(id = PG$id, Accession1 = sapply(strsplit(PG$"Leading protein IDs", ";"), function(x) { unlist(x)[1] }))
temp_PG$Pep <- strsplit(PG$"Peptide IDs", ";")
temp_PG$Pep <- sapply(temp_PG$Pep, function(x) { pep$Sequence[match(as.numeric(x), pep$id)] })
temp_PG$Seq <- db$Sequence[match(temp_PG$Accession1, db$"Protein ID")]
PG$"Sequence coverage [%]" <- apply(temp_PG[, c("Seq", "Pep")], 1, function(x) {
  round(100*c(Coverage(x[[1]], x[[2]])), 1)
})
for (exp in Exp) { #exp <- Exp[1]
  kol <- paste0(sapply(c("Peptide", "Evidence"), function(x) { paste0(x, c("s count", " IDs")) }), " - ", exp)
  kolk <- grep("count", kol, value = TRUE)
  koli <- grep("IDs", kol, value = TRUE)
  kole <- grep("^Evidence", kol, value = TRUE)
  kolp <- grep("^Peptide", kol, value = TRUE)
  PG[, kolk] <- 0
  PG[, koli] <- ""
  PG[[paste0("Sequence coverage [%] - ", exp)]] <- 0
  w <- which(ev$Experiment == exp)
  if (length(w)) {
    e <- ev[w, , drop = FALSE]
    temp1 <- melt(setNames(strsplit(e$"Protein group IDs", ";"), e$"Peptide ID"))
    temp1 <- do.call(data.frame, aggregate(temp1$L1, list(temp1$value), function(x) {
      x <- unique(x)
      return(c(Count = length(x), List = list(x)))
    }))
    temp1$x.Count <- unlist(temp1$x.Count)
    temp1$Pasted <- sapply(temp1$x.List, function(x) { paste(sort(as.numeric(unlist(x))), collapse = ";") })
    temp1$Pepseq <- lapply(temp1$x.List, function(x) { unique(pep$Sequence[match(unlist(x), pep$id)]) }) # Unique here because different peptidoforms can have the same sequence
    temp2 <- melt(setNames(strsplit(e$"Protein group IDs", ";"), e$id))
    temp2 <- do.call(data.frame, aggregate(temp2$L1, list(temp2$value), function(x) {
      c(Count = length(x), IDs = paste(sort(as.numeric(x)), collapse = ";"))
    }))
    temp2$x.Count <- as.integer(temp2$x.Count)
    w <- which(PG$id %in% temp1$Group.1)
    m <- match(PG$id[w], temp1$Group.1)
    PG[w, kolp] <- temp1[m, c("x.Count", "Pasted")]
    temp_PG$Pep <- NA
    temp_PG$Pep[w] <- temp1$Pepseq[m]
    PG[w, paste0("Sequence coverage [%] - ", exp)] <- apply(temp_PG[w, c("Seq", "Pep")], 1, function(x) {
      round(100*c(Coverage(x[[1]], x[[2]])), 1)
    })
    w <- which(PG$id %in% temp2$Group.1)
    m <- match(PG$id[w], temp2$Group.1)
    PG[w, kole] <- temp2[m, c("x.Count", "x.IDs")]
    if (IsBioID) {
      kol <- paste0(sapply(paste0("Biot. ", c("peptide", "evidence")), function(x) { paste0(x, c("s count", " IDs")) }), " - ", exp)
      kolk <- grep("count", kol, value = TRUE)
      koli <- grep("IDs", kol, value = TRUE)
      kole <- grep("^Biot\\. evidence", kol, value = TRUE)
      kolp <- grep("^Biot\\. peptide", kol, value = TRUE)
      PG[, kolk] <- 0
      PG[, koli] <- ""
      g <- grep(topattern(Modifs$Mark[wbiot], start = FALSE), e$"Modified sequence")
      if (length(g)) {
        e <- e[g, , drop = FALSE]
        temp1 <- melt(setNames(strsplit(e$"Protein group IDs", ";"), e$"Peptide ID"))
        temp1 <- do.call(data.frame, aggregate(temp1$L1, list(temp1$value), function(x) {
          x <- unique(x)
          return(c(Count = length(x), IDs = paste(sort(as.numeric(x)), collapse = ";")))
        }))
        temp1$x.Count <- as.integer(temp1$x.Count)
        temp1 <- do.call(data.frame, temp1)
        temp2 <- melt(setNames(strsplit(e$"Protein group IDs", ";"), e$id))
        temp2 <- do.call(data.frame, aggregate(temp2$L1, list(temp2$value), function(x) {
          c(Count = length(x), IDs = paste(sort(as.numeric(x)), collapse = ";"))
        }))
        temp2$x.Count <- as.integer(temp2$x.Count)
        w <- which(PG$id %in% temp1$Group.1)
        m <- match(PG$id[w], temp1$Group.1)
        PG[w, kolp] <- temp1[m, c("x.Count", "x.IDs")]
        w <- which(PG$id %in% temp2$Group.1)
        m <- match(PG$id[w], temp2$Group.1)
        PG[w, kole] <- temp2[m, c("x.Count", "x.IDs")]
      }
    }
  }
}
rm(temp_PG, temp1, temp2)

# Simple LFQ quant:
temp <- TopN(3, PG, Pep4Quant, pep, "id",
             Pep.Intens.Nms = grep(topattern(paste0(int.col, " - ")), colnames(pep), value = TRUE),
             log.Pep.Intens = FALSE,
             Mods = Mod4Quant,
             Min.Pep.Nb = NPep, corr = "global", Out.Norm = FALSE)
PG.int.col <- paste0("Top3 log10(", int.col, ") - ")
colnames(temp) <- gsub(paste0("^Top3 log10 - ", int.col, " - "), PG.int.col, colnames(temp))
PG[, colnames(temp)] <- temp

# Normalize (Levenberg-Marquardt)
if (length(Exp) > 1) {
  if (NormalizePG) {
    g <- grep(topattern(PG.int.col), colnames(PG), value = TRUE)
    temp <- PG[, c("id", g)]
    m <- apply(temp[,g], 2, function(x) { median(is.all.good(x)) })
    M <- median(is.all.good(unlist(temp[,g])))
    temp[,g] <- sweep(temp[,g], 2, m, "-") + M
    temp <- AdvNorm.IL(temp[, c("id", g)], "id", exprs.col = g, exprs.log = TRUE)
    PG[, paste0("Norm. ", g)] <- temp[, paste0("AdvNorm.", g)]
    PG.int.col <- paste0("Norm. ", PG.int.col)
  }
}
# Average expression:
PG[[paste0("Mean ", gsub(" - $", "", PG.int.col))]] <- apply(PG[, grep(topattern(PG.int.col), colnames(PG), value = TRUE)],
                                                             1, mean, na.rm = TRUE)

# Intensity distribution:
dir <- "Workflow control"
temp <- PG[, c("Common Name (short)", paste0(PG.int.col, Exp))]
temp <- melt(temp, id.vars = "Common Name (short)")
temp$variable <- gsub(topattern(PG.int.col), "", temp$variable)
temp <- temp[which(is.all.good(temp$value, 2)),]
ttl <- "LFQ density plot - Protein groups level"
plot <- ggplot(temp) + geom_histogram(aes(x = value, fill = variable), bins = 100) +
  theme_bw() + theme(strip.text.y = element_text(angle = 0, vjust = 0.5, hjust = 0)) +
  ggtitle(ttl) + facet_grid(variable~.) +
  xlab("log10(Protein Groups LFQ)")
poplot(plot)
setwd(paste0(wd, "/", dir))
ggsave(paste0(ttl, ".jpeg"), plot, dpi = 300, width = 10, height = 10, units = "in")
ggsave(paste0(ttl, ".pdf"), plot)
setwd(wd)
# Correlation:
dir <- "Workflow control"
if (length(Exp) > 1) {
  temp <- PG[,c("Common Name (short)", paste0(PG.int.col, Exp))]
  comb <- combinations(length(Exp), 2, Exp) 
  temp2 <- temp[, c("Common Name (short)", paste0(PG.int.col, comb[1,]))]
  temp2$X <- comb[1,1]
  temp2$Y <- comb[1,2]
  temp2$Comparison <- paste0(comb[1,1], " (X) vs ", comb[1,2], " (Y)")
  colnames(temp2)[2:3] <- c("log10 (X LFQ)", "log10 (Y LFQ)")
  if (length(Exp) > 2) {
    for (i in 2:nrow(comb)) {
      temp3 <- temp[, c("Common Name (short)", paste0(PG.int.col, comb[i,]))]
      temp3$X <- comb[i,1]
      temp3$Y <- comb[i,2]
      temp3$Comparison <- paste0(comb[i,1], " (X) vs ", comb[i,2], " (Y)")
      colnames(temp3)[2:3] <- c("log10 (X LFQ)", "log10 (Y LFQ)")
      temp2 <- rbind(temp2, temp3)
    }
  }
  test <- apply(temp2[,c("log10 (X LFQ)", "log10 (Y LFQ)")], 1, function(x) { length(is.all.good(x)) }) == 2
  temp2 <- temp2[which(test),] 
  temp3 <- as.data.frame(t(sapply(unique(temp2$Comparison), function(x) { #x <- unique(temp2$Comparison)[1]
    x1 <- temp2[which(temp2$Comparison == unlist(x)), c("log10 (X LFQ)", "log10 (Y LFQ)")]
    x1 <- x1$"log10 (Y LFQ)"-x1$"log10 (X LFQ)"
    return(setNames(c(x, paste0("Median = ", round(median(x1), 3)), paste0("S.D. = ", round(sd(x1), 3))), c("Comparison", "Median", "SD")))
  })))
  temp3$R <- sapply(strsplit(gsub(" \\([XY]\\)", "", temp3$Comparison), " vs "), function(x) {
    return(paste0("R = ", round(cor(pep[[paste0(int.col, " - ", x[[1]])]], pep[[paste0(int.col, " - ", x[[2]])]]), 3)))
  })
  temp2$Comparison <- gsub(" vs ", "\nvs\n", temp2$Comparison)
  temp3$Comparison <- gsub(" vs ", "\nvs\n", temp3$Comparison)
  x_min <- min(temp2$`log10 (X LFQ)`)
  y_min <- min(temp2$`log10 (Y LFQ)`)
  y_max <- max(temp2$`log10 (Y LFQ)`)
  ttl <- "Correlation plot - Protein groups level"
  plot <- ggplot(temp2) + geom_point(aes(x = `log10 (X LFQ)`, y = `log10 (Y LFQ)`), alpha = 0.25, size = 0.01) +
    geom_abline(intercept = 0, slope = 1, colour = "red") +
    geom_text(data = temp3, x = x_min, y = y_max - 0.01*(y_max - y_min), aes(label = R), hjust = 0, cex = 2.5) +
    geom_text(data = temp3, x = x_min, y = y_max - 0.06*(y_max - y_min), aes(label = Median), hjust = 0, cex = 2.5) +
    geom_text(data = temp3, x = x_min, y = y_max - 0.11*(y_max - y_min), aes(label = SD), hjust = 0, cex = 2.5) +
    facet_grid(Y~X) + ggtitle(ttl) + coord_fixed(1) +
    theme_bw() + theme(strip.text.y = element_text(angle = 0, vjust = 0.5, hjust = 0),
                       strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0))
  poplot(plot, 12, 20)
  setwd(paste0(wd, "/", dir))
  ggsave(paste0(ttl, ".jpeg"), plot, dpi = 300, width = 10, height = 10, units = "in")
  ggsave(paste0(ttl, ".pdf"), plot)
  setwd(wd)
}
g <- grep(topattern(PG.int.col), colnames(PG), value = TRUE)
test <- setNames(sapply(g, function(x) { length(is.all.good(PG[[x]])) }), gsub(topattern(PG.int.col), "", g))
print(test)

# Calculate ratios
if (MakeRatios) {
  PG.rat.col <- "log2(Ratio) - "
  rat.grps <- unique(SamplesMap$`Ratios group`)
  rat.grps <- rat.grps[which(!is.na(rat.grps))]
  temp <- data.frame(ID = PG$id, Accession = sapply(strsplit(PG$"Leading protein IDs", ";"), function(x) { x[[1]] }))
  temp$Sequence <- db$Sequence[match(temp$Accession, db$"Protein ID")]
  for (i in rat.grps) { #i <- rat.grps[1]
    m <- SamplesMap[which(SamplesMap$`Ratios group` == i),]
    if (sum(c(TRUE, FALSE) %in% m$Reference) < 2) { stop("The reference column should include TRUE and FALSE values!") } else {
      ref <- apply(pep[, paste0(int.col, " - ", m$Experiment[which(m$Reference)]), drop = FALSE], 1, function(x) {
        log_ratio_av(log10(x))
      }) # This makes it a log reference => use "-"
      w <- which(!m$Reference)
      for (x in w) { #x <- w[1]
        pep[[paste0(PG.rat.col, m$Experiment[x])]] <- (log10(pep[[paste0(int.col, " - ", m$Experiment[x])]]) - ref)/log10(2)
        temprat <- melt(setNames(strsplit(pep$"Protein group IDs", ";"), pep[[paste0(PG.rat.col, m$Experiment[x])]]))
        temprat$L1 <- as.numeric(temprat$L1)
        temprat <- aggregate(temprat$L1, list(temprat$value), median, na.rm = TRUE)
        PG[[paste0(PG.rat.col, m$Experiment[x])]] <- NA
        w2 <- which(PG$id %in% temprat$Group.1)
        PG[w2, paste0(PG.rat.col, m$Experiment[x])] <- temprat$x[match(PG$id[w2], temprat$Group.1)]
        if (RatiosThresh_2sided) {
          PG[[paste0("Regulated - ", m$Experiment[x])]] <- ""
          wup <- which(PG[[paste0(PG.rat.col, m$Experiment[x])]] >= RatiosThresh)
          wdwn <- which(PG[[paste0(PG.rat.col, m$Experiment[x])]] <= -RatiosThresh)
          PG[wup, paste0("Regulated - ", m$Experiment[x])] <- "up"
          PG[wdwn, paste0("Regulated - ", m$Experiment[x])] <- "down"
        } else {
          PG[[paste0("Enriched - ", m$Experiment[x])]] <- PG[[paste0(PG.rat.col, m$Experiment[x])]] >= RatiosThresh
        }
        e <- ev[which(ev$Experiment == m$Experiment[x]),]
        temp2 <- melt(setNames(strsplit(e$"Protein group IDs", ";"), e$Sequence))
        temp2 <- aggregate(temp2$L1, list(temp2$value), list)
        temp$Peptides <- list(NA)
        temp$Coverage <- 0
        w2 <- which(temp$ID %in% temp2$Group.1)
        temp$Peptides[w2] <- temp2$x[match(temp$ID[w2], temp2$Group.1)]
        temp$Coverage[w2] <- apply(temp[w2, c("ID", "Sequence", "Peptides")], 1, function(x) {
          round(100*Coverage(setNames(x[[2]], x[[1]]), x[[3]]), 1)
        })
        PG[[paste0("Sequence coverage [%] - ", m$Experiment[x])]] <- temp$Coverage
        if (IsBioID) {
          g <- grep(topattern(Modifs$Mark[wbiot], start = FALSE), e$"Modified sequence")
          PG[[paste0("Biot. peptides count - ", m$Experiment[x])]] <- 0
          if (length(g)) {
            e <- e[g,]
            temp2 <- melt(setNames(strsplit(e$"Protein group IDs", ";"), e$"Modified sequence"))
            temp2 <- aggregate(temp2$L1, list(temp2$value), function(x) { length(unique(x)) })
            temp$Count <- 0
            w2 <- which(temp$ID %in% temp2$Group.1)
            temp$Count[w2] <- temp2$x[match(temp$ID[w2], temp2$Group.1)]
            PG[[paste0("Biot. peptides count - ", m$Experiment[x])]] <- temp$Count
          }
        }
      }
    }
  }
  #
  # Fold change filters:
  test <- grep(topattern(PG.rat.col), colnames(PG), value =  TRUE)
  test <- gsub(topattern(PG.rat.col), "", test)
  FC_filt <- setNames(lapply(test, function(x) {
    x <- PG[[paste0(PG.rat.col, x)]]
    if (RatiosThresh_2sided) { x <- abs(x) }
    which(x >= RatiosThresh)
  }), test)
  g <- grep(topattern(PG.rat.col), colnames(PG), value = TRUE)
  test <- setNames(sapply(g, function(x) { length(is.all.good(PG[[x]])) }), gsub(topattern(PG.rat.col), "", g))
  print(test)
  #
  # Ratios distribution:
  dir <- "Workflow control"
  if (MakeRatios) {
    temp <- PG[, c("Common Name (short)", grep(topattern(PG.rat.col), colnames(PG), value = TRUE))]
    temp <- melt(temp, id.vars = "Common Name (short)")
    temp$variable <- gsub(topattern(PG.rat.col), "", temp$variable)
    temp <- temp[which(is.all.good(temp$value, 2)),]
    ttl <- "Ratios density plot - Protein groups level"
    plot <- ggplot(temp) + geom_histogram(aes(x = value, fill = variable), bins = 100) +
      geom_vline(xintercept = RatiosThresh, colour = "red") +
      theme_bw() + theme(strip.text.y = element_text(angle = 0, vjust = 0.5, hjust = 0)) +
      ggtitle(ttl) + facet_grid(variable~.) +
      xlab("log2(Ratio)")
    if (RatiosThresh_2sided) { plot <- plot + geom_vline(xintercept = -RatiosThresh, colour = "red") }
    poplot(plot)
    setwd(paste0(wd, "/", dir))
    ggsave(paste0(ttl, ".jpeg"), plot, dpi = 300, width = 10, height = 10, units = "in")
    ggsave(paste0(ttl, ".pdf"), plot)
    setwd(wd)
  }
}

# Proteome ruler
if (protrul) {
  temp <- Prot.Ruler(PG, db, Expr.roots = c(PG.int.col, paste0("Mean ", PG.int.col)))
  PG <- temp$Protein.groups
  db <- temp$Database
}

save.image("Backup.RData")
#load("Backup.RData")

# Some summaries of the data:
Exp_summary <- MQ.summary(wd = wd, ev = ev, pg = PG, mods = setNames(Modifs$Mark, Modifs$"Full name"),
                          raw.files.order = raw.files, sc = min(c(60, round(length(raw.files)/length(Exp)))),
                          save = c("jpeg", "pdf"))
write.csv(Exp_summary, file = "Summary.csv", row.names = FALSE)
AA_biases <- AA_bias(ev, db)
write.csv(AA_biases, file = "Amino Acid composition biases.csv", row.names = FALSE)
View(AA_biases)

save.image("Backup.RData")
#load("Backup.RData")

# Coverage map of proteins of interest
nCharLim <- 50
if (protlist) {
  for (prnm in 1:length(prot.names)) { #prnm <- 1
    p <- prot.names[prnm]
    IDs <- IDs.list[[p]]
    if (length(IDs) == 0) {
      p <- toupper(dlgInput("Sorry, I could not find this protein name in the database, do you want to provide an alternate name?", "")$res)
      if (p != "") {
        prot.names[prnm] <- p
        IDs <- db$"Protein ID"[which(db$`Common Name` == p)]
        if (length(IDs) == 0) { warning("Really sorry, but I really cannot make sense of this protein name, skipping...") }
      } else { warning("Ok, fine, we'll skip then.") }
    }
    if (length(IDs)) {
      TMP <- ev
      if (removeMBR) { TMP <- TMP[grep("-MATCH$", TMP$Type, invert = TRUE),] }
      test <- setNames(strsplit(TMP$Proteins, ";"), TMP$id)
      test <- melt(test)
      w <- which(test$value %in% IDs)
      if (!length(w)) {
        if (length(IDs) > 1) {
          mess <- paste0("(s) ", paste(IDs[1:(length(IDs)-1)], collapse = ", "), " and ", IDs[length(IDs)])
        } else { mess <- paste0(" ", IDs) }
        warning(paste0("No peptides identified for protein accession", mess, "."))
      } else {
        w <- as.numeric(unique(test$L1[w]))
        TMP <- TMP[match(w, TMP$id),]
        for (id in IDs) { #id <- IDs[1]
          nm <- gsub(" - $", "", paste0(id, "_", db$"Common Name"[match(id, db$"Protein ID")]))
          if (nchar(nm) > nCharLim) { nm <- paste0(substr(nm, 1, nCharLim-3), "...") }
          P <- setNames(db$Sequence[match(id, db$"Protein ID")], nm)
          P2 <- paste0("_", P, "_")
          s <- data.frame(Seq = unique(TMP$Sequence))
          s$Matches <- sapply(s$Seq, function(x) {
            x <- unlist(strsplit(P2, x))
            return(nchar(x[1:(length(x)-1)]) + 1)
          })
          s2 <- data.frame("Modified sequence" = unique(TMP$"Modified sequence"), check.names = FALSE)
          s2$Sequence <- TMP$Sequence[match(s2$"Modified sequence", TMP$"Modified sequence")]
          s2$Matches <- s$Matches[match(s2$Sequence, s$Seq)]
          tempev <- setNames(lapply(Exp, function(e) { #e <- Exp[1]
            wh <- which(TMP$Experiment == e); if (length(wh)) {
              e <- TMP[wh,]
              wh <- which(is.all.good(log10(e$Intensity), 2)); if (length(wh)) {
                res <- set_colnames(aggregate(e$Intensity[wh], list(e$"Modified sequence"[wh]), sum),
                                    c("Modified sequence", "Intensity"))
                e <- set_colnames(aggregate(e$PEP[wh], list(e$"Modified sequence"[wh]), function(x) {
                  x <- is.all.good(x)
                  if (!length(x)) { x <- NA } else { x <- max(x) } # Conservative: taking the worst estimate of PEP
                  return(x)
                }), c("Modified sequence", "PEP"))
                res$PEP <- e$PEP; rm(e)
                res$"log10(Intensity)" <- log10(res$Intensity)
                res$Intensity <- NULL
                res$Matches <- s2$Matches[match(res$"Modified sequence", s2$"Modified sequence")]
              } else { res <- NA }
            } else { res <- NA }
            return(res)
          }), Exp)
          for (e in Exp) { #e <- Exp[1]
            tmp <- tempev[[e]]
            if (class(tmp) == "data.frame") {
              setwd(wd); suppressWarnings(dir.create("Coverage plots"))
              setwd(paste0(wd, "/Coverage plots"))
              ttl <- paste0("Coverage - ", nm, " - ", e, " (intensity)")
              Coverage(P, tmp$"Modified sequence", Mode = "Align2", title = ttl, save = c("jpeg", "pdf"),
                       intensities = tmp$`log10(Intensity)`)
              ttl <- paste0("Coverage - ", nm, " - ", e, " (-log10(PEP))")
              Coverage(P, tmp$"Modified sequence", Mode = "Align2", title = ttl, save = c("jpeg", "pdf"),
                       intensities = -log10(tmp$PEP), colscale = 8, na = "red")
              setwd(wd)
            }
          }
          # Correlation and ratio plots:
          setwd(wd); suppressWarnings(dir.create("Correlation plots")); suppressWarnings(dir.create("Ratio plots"))
          wh <- which(!sapply(names(tempev), function(x) { (class(tempev[x]) != " data.frame")&(is.na(tempev[x])) }))
          if (length(wh) > 1) {
            temp1 <- melt(tempev[wh])
            temp1$variable <- as.character(temp1$variable)
            temp1$Dummy <- apply(temp1[,c("Modified sequence", "L1")], 1, paste, collapse = "---")
            temp11 <- temp1[which(temp1$variable == "log10(Intensity)"),]
            temp12 <- temp1[which(temp1$variable == "Matches"),]
            temp11$Matches <- temp12$value[match(temp11$Dummy, temp12$Dummy)]
            temp1 <- temp11
            rm(temp11, temp12)
            comb <- combinations(length(Exp), 2, Exp)
            kount <- 0
            for (j in 1:nrow(comb)) { #j <- 1
              s1 <- temp1[which(temp1$L1 == comb[j,1]),]
              s2 <- temp1[which(temp1$L1 == comb[j,2]),]
              if (sum(c(nrow(s1) > 0, nrow(s2) > 0)) == 2) {
                wtst <- which(s1$"Modified sequence" %in% s2$"Modified sequence")
                if (length(wtst)) {
                  kount <- kount + 1
                  temp3 <- data.frame("Modified sequence" = s1$"Modified sequence"[which(s1$"Modified sequence" %in% s2$"Modified sequence")], check.names = FALSE)
                  temp3[,c("Matches" , "log10 (X LFQ)")] <- s1[match(temp3$"Modified sequence", s1$"Modified sequence"), c("Matches", "value")]
                  temp3$"log10 (Y LFQ)" <- s2$value[match(temp3$"Modified sequence", s2$"Modified sequence")]
                  temp3$Comparison <- paste0(comb[j,1], " (X) vs ", comb[j,2], " (Y)")
                  if (kount == 1) { temp2 <- temp3 } else { temp2 <- rbind(temp2, temp3) }
                }
              }
            }
            test <- apply(temp2[,c("log10 (X LFQ)", "log10 (Y LFQ)")], 1, function(x) { length(is.all.good(x)) }) == 2
            temp2 <- temp2[which(test),]
            temp3 <- as.data.frame(t(sapply(unique(temp2$Comparison), function(x) { #x <- unique(temp2$Comparison)[1]
              x1 <- temp2[which(temp2$Comparison == unlist(x)), c("log10 (X LFQ)", "log10 (Y LFQ)")]
              x1 <- x1$"log10 (Y LFQ)" - x1$"log10 (X LFQ)"
              return(setNames(c(x, paste0("Median = ", round(median(x1), 3)), paste0("S.D. = ", round(sd(x1), 3))), c("Comparison", "Median", "SD")))
            })))
            temp3$R <- sapply(strsplit(gsub(" \\([XY]\\)", "", temp3$Comparison), " vs "), function(x) {
              return(paste0("R = ", round(cor(pep[[paste0(int.col, " - ", x[[1]])]], pep[[paste0(int.col, " - ", x[[2]])]]), 3)))
            })
            temp2[, c("X", "Y")] <- Isapply(strsplit(temp2$Comparison, " vs "), unlist)
            temp3[, c("X", "Y")] <- Isapply(strsplit(temp3$Comparison, " vs "), unlist)
            temp2 <- temp2[order(temp2$Comparison),]
            temp3 <- temp3[order(temp3$Comparison),]
            temp2$Sequence <- gsub("\\([^\\)]+\\)|_", "", temp2$"Modified sequence")
            temp2$"C-terminal extent" <- temp2$Matches + nchar(temp2$Sequence) - 1
            temp2$"log10(Y/X)" <- temp2$`log10 (Y LFQ)` - temp2$`log10 (X LFQ)`
            x_min <- min(temp2$`log10 (X LFQ)`)
            y_min <- min(temp2$`log10 (Y LFQ)`)
            y_max <- max(temp2$`log10 (Y LFQ)`)
            ttl <- paste0("Correlation plot - ", nm)
            plot <- ggplot(temp2) + geom_point(aes(x = `log10 (X LFQ)`, y = `log10 (Y LFQ)`, colour = `C-terminal extent`),
                                               alpha = 0.25, size = 0.01) +
              geom_abline(intercept = 0, slope = 1, colour = "red") + coord_fixed(1) +
              scale_colour_gradient(low = "green", high = "red") +
              geom_text(data = temp3, x = x_min, y = y_max - 0.01*(y_max - y_min), aes(label = R), hjust = 0, cex = 2.5) +
              geom_text(data = temp3, x = x_min, y = y_max - 0.06*(y_max - y_min), aes(label = Median), hjust = 0, cex = 2.5) +
              geom_text(data = temp3, x = x_min, y = y_max - 0.11*(y_max - y_min), aes(label = SD), hjust = 0, cex = 2.5) +
              facet_grid(Y~X) +
              theme_bw() + theme(strip.text.y = element_text(angle = 0, vjust = 0.5, hjust = 0),
                                 strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0)) +
              ggtitle(ttl)
            poplot(plot)
            setwd(paste0(wd, "/Correlation plots"))
            ggsave(paste0(gsub(":|\\*|\\?|<|>|\\||/", "-", ttl), ".jpeg"), plot, dpi = 300, width = 10, height = 10, units = "in")
            ggsave(paste0(gsub(":|\\*|\\?|<|>|\\||/", "-", ttl), ".pdf"), plot)
            setwd(wd)
            temp2$X <- gsub("\\(X\\)", "(a)", temp2$X)
            temp2$Y <- gsub("\\(Y\\)", "(b)", temp2$Y)
            sk <- (max(temp2$`C-terminal extent`) - min(temp2$`C-terminal extent`))/(max(temp2$`log10(Y/X)`) - min(temp2$`log10(Y/X)`))
            ttl <- paste0("Ratio plot - ", nm)
            plot <- ggplot(temp2) + geom_segment(aes(x = Matches, xend = `C-terminal extent`, y = `log10(Y/X)`, yend = `log10(Y/X)`, colour = Matches),
                                                 size = 1) +
              coord_fixed(sk/3) +
              scale_colour_gradient(low = "green", high = "red") +
              scale_x_continuous(breaks = 50*(1:floor(nchar(P)/50))) +
              facet_grid(Y~X) + ylab("log10(b/a)") +
              theme_bw() + theme(strip.text.y = element_text(angle = 0, vjust = 0.5, hjust = 0),
                                 strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0)) +
              ggtitle(ttl)
            poplot(plot, 12, 20)
            setwd(paste0(wd, "/Ratio plots"))
            ggsave(paste0(gsub(":|\\*|\\?|<|>|\\||/", "-", ttl), ".jpeg"), plot, dpi = 300, width = 10, height = 10, units = "in")
            ggsave(paste0(gsub(":|\\*|\\?|<|>|\\||/", "-", ttl), ".pdf"), plot)
            setwd(wd)
          }
        }
      }
    }
  }
}

# LFQ plots - proteins sorted by abundance:
LFQsortDir <- paste0(wd, "/Protein Groups sorted by LFQ")
suppressWarnings(dir.create(LFQsortDir))
setwd(LFQsortDir)
tstorg2 <- c()
myColors <- setNames("black", "-")
if (tstorg) {
  tstorg2 <- unique(PG[[orgkol]])
  if (length(tstorg2) > 1) {
    myColors <- setNames(colorRampPalette(c("blue", "green"))(length(tstorg2)), tstorg2)
  }
}
for (e in Exp) { #e <- Exp[1]
  temp <- PG[, c(#"Leading protein IDs",
    "Protein IDs", "Label", "Common Name (short)", paste0(PG.int.col, e), "id")]
  colnames(temp)[which(colnames(temp) == paste0(PG.int.col, e))] <- "Top3 log10 Intensity"
  temp <- temp[which(is.all.good(temp$"Top3 log10 Intensity", 2)),]
  if (nrow(temp)) {
    temp <- temp[order(temp$"Top3 log10 Intensity", decreasing = TRUE),]
    temp$"Protein Group" <- temp$Label; temp$Label <- NULL
    test <- aggregate(temp$"Protein Group", list(temp$"Protein Group"), length)
    w <- which(test$x > 1)
    if (length(w)) {
      test <- test[w,]
      for (i in test$Group.1) {
        w <- which(temp$"Protein Group" == i)
        temp$"Protein Group"[w[2:length(w)]] <- paste0(temp$"Protein Group"[w[2:length(w)]], "_", 2:length(w))
      }
    }
    temp$"In list" <- "-"
    if (length(tstorg2) > 1) { temp$Category <- PG[match(temp$id, PG$id), orgkol] } else {
      temp$Category <- "-"
    }
    temp$"Protein Group" <- factor(temp$"Protein Group", levels = temp$"Protein Group")
    temp2 <- setNames(strsplit(temp$"Protein IDs", ";"), temp$id)
    temp2 <- melt(temp2)
    if (protlist) {
      w <- which(temp2$value %in% db$"Protein ID"[which(db$`Common Name` %in% prot.names)])
      m <- match(temp2$L1[w], temp$id)
      temp$"In list"[m] <- "+"
      temp$Category[m] <- c("+", "In list")[(length(tstorg2) > 1)+1]
    }
    myColors[[c("+", "In list")[(length(tstorg2) > 1)+1]]] <- "red"
    myAlpha <- setNames(c(0.25, 1), c("-", "+"))
    colScale <- scale_colour_manual(name = "Category", values = myColors)
    fillScale <- scale_fill_manual(name = "Category", values = myColors)
    alphaScale <- scale_alpha_manual(name = "In list", values = myAlpha)
    intmin <- floor(min(temp$"Top3 log10 Intensity"))
    intmax <- ceiling(max(temp$"Top3 log10 Intensity"))
    intscale <- intmax-intmin
    #xmax <- max(c(max(c(0, which(temp$Category == "+")))+round(nrow(temp)/15), nrow(temp)))
    xmax <- nrow(temp)*18/15
    ttl <- paste0("Protein Groups sorted by LFQ - ", e)
    w1 <- which(temp$"In list" == "+")
    w2 <- which(temp$"In list" != "+")
    plot <- ggplot(temp)
    if (length(tstorg2) > 1) {
      plot <- plot + geom_bar(stat = "identity", aes(x = `Protein Group`, y = `Top3 log10 Intensity`, fill = Category))
    } else {
      plot <- plot + geom_bar(stat = "identity", aes(x = `Protein Group`, y = `Top3 log10 Intensity`, alpha = `In list`,
                                                     fill = Category)) + alphaScale
    }
    plot <- plot +
      annotate("text", nrow(temp)/2, intmax*1.2+intscale*0.04, label = paste0(nrow(temp), " Protein Groups"), hjust = 0.5, ) +
      colScale + fillScale + theme_bw() + ggtitle(ttl) + ylab("log10(Top3 Intensity)") +
      coord_cartesian(xlim = c(1, xmax), ylim = c(intmin, intmax*1.2+intscale*0.05)) +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"), axis.text.x = element_blank(), axis.ticks = element_blank(),
            plot.margin = margin(r = 100)) + guides(alpha = FALSE)
    plot_ly <- ggplotly(plot, tooltip = "Protein Group")
    setwd(LFQsortDir) # For some reason, unless I do this the default selfcontained = TRUE argument gets ignored and
    # a folder with external resources is created for each html plot!
    saveWidget(partial_bundle(plot_ly), paste0(ttl, ".html"))
    system(paste0("open \"", LFQsortDir, "/", ttl, ".html"))
    plot <- plot +
      geom_text(data = temp[w2,], aes(x = `Protein Group`, y = `Top3 log10 Intensity` + intmax*0.01, colour = Category,
                                      label = `Protein Group`, alpha = `In list`), angle = 45, hjust = 0, cex = 3.5) +
      geom_text(data = temp[w1,], aes(x = `Protein Group`, y = `Top3 log10 Intensity` + intmax*0.01, colour = Category,
                                      label = `Protein Group`), angle = 45, hjust = 0, cex = 4, fontface = "bold" )
    poplot(plot, 12, 22)
    setwd(LFQsortDir)
    ggsave(filename = paste0(ttl, ".jpeg"), plot, dpi = 300, width = 10, height = 10, units = "in")
    ggsave(filename = paste0(ttl, ".pdf"), plot)
    setwd(wd)
  }
}
setwd(wd)

save.image("Backup.RData")
#load("Backup.RData")

# Annotate
setwd(wd)
if (Annotate) {
  kol <- colnames(Parsed_annotations)
  annot.col <- kol[which(!kol %in% c("Accession", "id", "ID", "Names", "Sequence", "MW (Da)"))]
  annot.col2 <- annot.col[which(!annot.col %in% colnames(db))]
  if (length(annot.col2)) {
    db[, annot.col2] <- NA
    w <- which(db$"Protein ID" %in% Parsed_annotations$Accession)
    db[w, annot.col2] <- Parsed_annotations[match(db$"Protein ID"[w], Parsed_annotations$Accession), annot.col2]
    write.csv(db, "Parsed, annotated search db.csv", row.names = FALSE)
    #db <- read.csv("Parsed, annotated search db.csv", row.names = FALSE)
  }
  temp <- melt(setNames(strsplit(PG$"Leading protein IDs", ";"), PG$id))
  temp[, annot.col] <- db[match(temp$value, db$"Protein ID"), annot.col]
  stopifnot(length(which(temp$value != temp$Accession)) == 0)
  for (i in annot.col) { temp[[i]] <- strsplit(as.character(temp[[i]]), ";") }
  temp <- aggregate(temp[, annot.col], list(temp$L1), function(x) { list(unique(unlist(x))) })
  for (i in annot.col) { temp[[i]] <- sapply(temp[[i]], function(x) { paste(sort(unique(unlist(x))), collapse = ";") }) }
  PG[, annot.col] <- temp[match(PG$id, temp$Group.1), annot.col]
}

## Perform GO terms enrichment analysis
if (GO_analysis) {
  PG$"First protein" <- sapply(strsplit(PG$`Leading protein IDs`, ";"), function(x) { unlist(x)[1] })
  filt <- setNames(lapply(Exp, function(e) { which(PG[[paste0(PG.int.col, e)]] > 0) }), Exp)
  for (n in names(filt)) { #n <- names(filt)[1]
    setwd(wd)
    temp <- GO_enrich_DAVID(Prot = PG[filt[[n]],], Prot_ID_col = "First protein", mode = "dataset",
                            Ref = PG$"First protein",
                            Prot_FC_root = paste0(PG.int.col, n), title.root = paste0("Bubble_plot - ", n, "_"),
                            save = c("jpeg", "pdf"), return = TRUE,
                            GO_FDR = c(0.1, 0.2, 0.3), True_Zscore = FALSE, subfolder = "GO enrichment analysis",
                            subfolderpertype = FALSE, plotly = TRUE, plotly_local = TRUE,
                            DAVID_email = DavidEmail)
  }
  if (MakeRatios) {
    db_obs <- db[which(db$"Protein ID" %in% unlist(strsplit(PG$"Leading protein IDs", ";"))),]
    setwd(wd)
    temp <- GO_enrich_DAVID(Prot = PG, Prot_ID_col = "Protein IDs", mode = "regulated",
                            filters = FC_filt, Prot_FC_root = PG.rat.col, title.root = "Bubble_plot_",
                            save = c("jpeg", "pdf"), return = TRUE,
                            GO_FDR = c(0.1, 0.2, 0.3), True_Zscore = TRUE, subfolder = "GO enrichment analysis",
                            subfolderpertype = FALSE, plotly = TRUE, plotly_local = TRUE,
                            DAVID_email = DavidEmail)
  }
}

# Venn diagrams
if (Venn_Obs) {
  test <- setNames(lapply(Exp, function(e) {
    x <- PG[[paste0(PG.int.col, e)]]
    which((!is.na(x))&(x > 0))
  }), Exp)
  w <- which(sapply(test, length) > 0)
  if (length(w) > 1) {
    require(VennDiagram)
    if (length(w) > 5) {
      msg <- paste0("Too many samples, select at least 2 and up to 5 to include in the Venn diagram (commas-separated):\n",
                    paste(paste0(" - ", 1:length(w), ": ", Exp[w]), collapse = "\n"))
      w1 <- suppressWarnings(as.numeric(unlist(strsplit(dlg_input(msg)$res, " *, *"))))
      while (sum(is.na(w1))||((length(w1) > 5)&&(length(w1) > 1))) {
        w1 <- suppressWarnings(as.numeric(unlist(strsplit(dlg_input(msg)$res, " *, *"))))
      }
      w <- w[w1]
    }
    print("Creating LFC Venn diagrams:")
    VennExp <- Exp[w]
    test <- test[VennExp]
    AnalysisParam$"Venn diagram (LFC) - samples" <- list(VennExp)
    setwd(wd); suppressWarnings(dir.create("Venn diagrams"))
    ttl <- "LFQ Venn diagram - global"
    venn.diagram(x = test, filename = paste0("Venn diagrams/", ttl, ".tiff"), col = "transparent", title = ttl,
                 fill = rainbow(n = length(test)),  alpha = 1/length(test),
                 cex = 1.5, fontface = "bold", main = gsub(" - ", "\n", ttl), margin = 0.1)  
    system(paste0("open \"Venn diagrams/", ttl, ".tiff", "\""))
  } else { warning("Skipping LFC Venn diagrams: not enough valid samples!") }
}
if (Venn_Ratios) {
  tst <- sapply(FC_filt, length)
  if (0 %in% tst) {
    w <- which(tst == 0)
    l <- length(w)
    tmp <- names(FC_filt)[which(tst == 0)]
    if (l > 1) { tmp <- paste0(paste(tmp[1:(l-1)], collapse = ", "), " and ", tmp[l]) }
    warning(paste0("No filtered proteins for samples ", tmp, ", they will be skipped!"))
  }
  fc_filt <- FC_filt[which(tst > 0)]
  if (length(fc_filt) > 1) {
    require(VennDiagram)
    if (length(fc_filt) > 5) {
      msg <- paste0("Too many samples, select at least 2 and up to 5 to include in the Venn diagram (commas-separated):\n",
                    paste(paste0(" - ", 1:length(fc_filt), ": ", names(fc_filt)), collapse = "\n"))
      w1 <- suppressWarnings(as.numeric(unlist(strsplit(dlg_input(msg)$res, " *, *"))))
      while (sum(is.na(w))||((length(w) > 5)&&(length(w) > 1))) {
        w1 <- suppressWarnings(as.numeric(unlist(strsplit(dlg_input(msg)$res, " *, *"))))
      }
      fc_filt <- fc_filt[w1]
    }
    print("Creating ratios Venn diagrams:")
    VennExp <- names(fc_filt)
    fc_filt <- fc_filt[VennExp]
    AnalysisParam$"Venn diagram (ratios) - samples" <- list(VennExp)
    setwd(wd); suppressWarnings(dir.create("Venn diagrams"))
    ttl <- "Ratios Venn diagram - global"
    venn.diagram(x = fc_filt, filename = paste0("Venn diagrams/", ttl, ".tiff"), col = "transparent", title = ttl,
                 fill = rainbow(n = length(fc_filt)),  alpha = 1/length(fc_filt),
                 cex = 1.5, fontface = "bold", main = gsub(" - ", "\n", ttl), margin = 0.1)  
    system(paste0("open \"Venn diagrams/", ttl, ".tiff", "\""))
  } else { warning("Skipping ratios Venn diagrams: not enough valid samples!") }
}

# PCA plots
## Process data
if (length(Exp) > 2) {
  setwd(wd); suppressWarnings(dir.create("PCA plots"))
  temp <- PG[, grep(topattern(PG.int.col), colnames(PG), value = TRUE)]
  colnames(temp) <- gsub(topattern(PG.int.col), "", colnames(temp))
  ## Impute data
  ## Here we have no way to decide between MAR/MCAR/MNAR,
  ## so we will instead replace every missing value with a random value drawn from a gaussian distribution of reduced m and sd
  m <- median(is.all.good(unlist(temp)))
  sd <- sd(is.all.good(unlist(temp)))
  for (i in colnames(temp)) {
    w <- which(!is.all.good(temp[[i]], 2))
    temp2 <- rnorm(length(w), m-3, sd/2)
    temp[w, i] <- temp2
  }
  nrm <- PG[[paste0("Mean ", gsub(" - $", "", PG.int.col))]]
  w <- which(is.all.good(nrm, 2))
  temp <- sweep(temp[w,], 1, nrm[w], "-")
  rownames(temp) <- PG$"Leading protein IDs"[w]
  ## 1/ Samples level:
  pc <- prcomp(t(temp), scale. = TRUE)
  scores <- as.data.frame(pc$x)
  pv <- round(100*(pc$sdev)^2 / sum(pc$sdev^2), 0)
  pv <- pv[which(pv > 0)]
  pv <- paste0("Components: ", paste(sapply(1:length(pv), function(x) {
    paste0("PC", x, ": ", pv[x], "%")
  }), collapse = ", "))
  scores$Sample <- rownames(scores)
  rownames(scores) <- NULL
  tst <- apply(scores[, c("PC1", "PC2")], 2, function(x) { length(unique(x)) })
  if (min(tst) > 1) {
    ttl <- "PCA plot - Samples-level"
    plot <- ggplot(scores) + geom_point(aes(x = PC1, y = PC2, colour = Sample)) +
      coord_fixed() + ggtitle(ttl, subtitle = pv) +
      geom_text_repel(aes(x = PC1, y = PC2, label = Sample, colour = Sample), size = 2.5) + 
      theme_bw() + theme(legend.position = "none")
    poplot(plot)
    if ("PC3" %in% colnames(scores)) {
      plot_lyPCA <- plot_ly(scores, x = ~PC1, y = ~PC2, z = ~PC3, color = ~Sample, text = ~Sample, type = "scatter3d",
                            mode = "markers")
      plot_lyPCA <- add_trace(plot_lyPCA, scores, x = ~PC1, y = ~PC2, z = ~PC3, type = "scatter3d", mode = "text",
                              showlegend = FALSE)
      plot_lyPCA <- layout(plot_lyPCA, title = ttl)
      saveWidget(partial_bundle(plot_lyPCA), paste0(wd, "/PCA plots/", ttl, ".html"))
      system(paste0("open \"", wd, "/PCA plots/", ttl, ".html"))
    } else { poplot(plot, width = 18) }
    ggsave(paste0(wd, "/PCA plots/", ttl, ".jpeg"), plot, dpi = 150)
    ggsave(paste0(wd, "/PCA plots/", ttl, ".pdf"), plot)
    ## 2/ Protein groups level:
    pc <- prcomp(temp, scale. = TRUE)
    scores <- as.data.frame(pc$x)
    pv <- round(100*(pc$sdev)^2 / sum(pc$sdev^2), 0)
    pv <- pv[which(pv > 0)]
    pv <- paste0("Components: ", paste(sapply(1:length(pv), function(x) {
      paste0("PC", x, ": ", pv[x], "%")
    }), collapse = ", "))
    scores$"Leading protein IDs" <- rownames(scores)
    rownames(scores) <- NULL
    scores[, c("Protein IDs", "Protein group")] <- PG[match(scores$"Leading protein IDs", PG$"Leading protein IDs"),
                                                      c("Protein IDs", "Label")]
    scores$Alpha <- (scores$PC1^2 + scores$PC2^2)
    scores$Direction <- apply(temp[w,], 1, function(x) {
      wh <- which(is.all.good(10^x, 2))
      return(weighted.mean(c(1:length(colnames(temp)))[wh], 10^x[wh]))
    })
    scores$Class <- ""
    breaks <- 1:length(Exp)
    labels <- Exp
    if (protlist) {
      g1 <- grsep(prot.IDs, x = scores$"Protein IDs")
      if (length(g1)) {
        g2 <- grsep(prot.IDs, x = scores$"Protein IDs", invert = TRUE)
        scores2 <- scores[g1,]
        scores <- scores[g2,]
        scores2$Alpha <- 1
        scores2$Class <- length(Exp)+1
      }
    }
    ttl <- "PCA plot - Protein groups (PG-level)"
    plot <- ggplot(scores) + geom_point(aes(x = PC1, y = PC2, alpha = Alpha, colour = Direction, text = `Protein group`), shape = 1) +
      coord_fixed() + ggtitle(ttl, subtitle = pv) +
      theme_bw() + guides(alpha = FALSE, shape = FALSE) +
      scale_color_gradientn(colors = hcl.colors(length(Exp)), breaks = breaks, labels = labels, guide = "legend")
    if (protlist) {
      if (length(g1)) {
        plot <- plot + geom_point(data = scores2, colour = "red", shape = 2, aes(x = PC1, y = PC2, text = `Protein group`))
      }
    }
    ggsave(paste0(wd, "/PCA plots/", ttl, ".jpeg"), plot, dpi = 150)
    ggsave(paste0(wd, "/PCA plots/", ttl, ".pdf"), plot)
    if ("PC3" %in% colnames(scores)) {
      plot_lyPCAProt <- plot_ly(scores, x = ~PC1, y = ~PC2, z = ~PC3, color = ~Direction, text = ~`Protein group`,
                                type = "scatter3d", mode = "markers", showlegend = FALSE, marker = list(size = 1))
      plot_lyPCAProt <- layout(plot_lyPCAProt, title = ttl)
      if (length(g1)) {
        scores3 <- scores2
        scores3$"Protein group" <- gsub(" - | ?, ?", "<br>", scores3$"Protein group")
        plot_lyPCAProt <- add_trace(plot_lyPCAProt, data = scores3, x = ~PC1, y = ~PC2, z = ~PC3,
                                    type = "scatter3d", mode = "markers+text", color = I("red"), marker = list(size = 5),
                                    showlegend = FALSE, textposition = "bottom right")
      }
    } else {
      plot_lyPCAProt <- ggplotly(plot, tooltip = "text")
    }
    saveWidget(partial_bundle(plot_lyPCAProt), paste0(wd, "/PCA plots/", ttl, ".html"))
    system(paste0("open \"", wd, "/PCA plots/", ttl, ".html"))
    plot <- plot + geom_text_repel(data = scores, aes(x = PC1, y = PC2, label = `Protein group`, alpha = Alpha), size = 2.5)
    if (protlist) {
      if (length(g1)) {
        plot <- plot + geom_text_repel(data = scores2, colour = "red", size = 2.5, aes(x = PC1, y = PC2, label = `Protein group`))
      }
    }
    poplot(plot, width = 18)
    ggsave(paste0(wd, "/PCA plots/", ttl, " (labels).jpeg"), plot, dpi = 150)
    ggsave(paste0(wd, "/PCA plots/", ttl, " (labels).pdf"), plot)
  } else { warning("No PCA plots drawn: samples are too similar!") }
}

# Heatmap
## Heatmap with hierarchical clustering:
kol <- grep(topattern(PG.int.col), colnames(PG), value = TRUE)
if (length(kol) > 1) {
  require(gplots)
  require(pdftools)
  setwd(wd); suppressWarnings(dir.create("Heatmaps"))
  temp <- set_colnames(PG[, kol], gsub(topattern(PG.int.col), "", kol))
  av <- apply(temp, 1, function(x) { x <- is.all.good(x); return(sum(x)/length(x)) })
  temp <- sweep(temp, 1, av, "-")
  temp <- as.matrix(temp)
  clust <- hclust(dist(data.matrix(t(temp))))
  nm <- "Heatmap - Samples (PG-level)"
  par(cex.main = 0.8)
  grDevices::windows(width = 18, height = 10)
  heatmap.2(temp, Colv = as.dendrogram(clust), Rowv = FALSE, main = nm,
            xlab = NULL, ylab = NULL, labRow = FALSE, key = TRUE, keysize = 1,
            trace = "none", density.info = c("none"),
            margins = c(12, 9), col = colorRampPalette(colors = c("blue","gray","red")),
            sepcolor = "blue", dendrogram = "column",
            cexRow = 0.5,
            cexCol = 1)
  dev.copy(pdf, paste0("Heatmaps/", nm, ".pdf"))
  dev.off()
  pdf_convert(paste0("Heatmaps/", nm, ".pdf"), format = "jpeg", filenames = paste0("Heatmaps/", nm, ".jpeg"), dpi = 600)
  #system(paste0("open \"", wd, "/Heatmaps/", nm, ".jpeg", "\""))
}
#
## Heatmap with double hierarchical clustering and highlighting proteins of interest
kol <- grep(topattern(PG.int.col), colnames(PG), value = TRUE)
HKlust <- FALSE
if (length(kol) > 1) {
  library(ggplot2)
  library(ggdendro)
  library(reshape2)
  library(gridExtra)
  library(ggpubr)
  setwd(wd); suppressWarnings(dir.create("Heatmaps"))
  temp <- set_rownames(set_colnames(PG[, kol], gsub(topattern(PG.int.col), "", kol)), PG$Label)
  w <- which(apply(temp, 1, function(x) { length(is.all.good(x)) }) == length(kol))
  temp <- temp[w,]
  av <- apply(temp, 1, function(x) { x <- is.all.good(x); return(sum(x)/length(x)) })
  temp <- sweep(temp, 1, av, "-")
  temp <- as.matrix(temp)
  hcluster <- hclust(dist(t(temp)))
  hdendro <- as.dendrogram(hcluster)
  vcluster <- hclust(dist(temp))
  vdendro <- as.dendrogram(vcluster)
  hord <- order.dendrogram(hdendro)
  vord <- order.dendrogram(vdendro)
  hord <- colnames(temp)[hord]
  vord <- rownames(temp)[vord]
  # Create dendrogram plot
  hdendro.plot <- ggdendrogram(data = hdendro) + theme(axis.text.y = element_text(size = 0.1))
  vdendro.plot <- ggdendrogram(data = vdendro, rotate = TRUE, labels = FALSE) +
    theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
          panel.background = element_rect(fill = "transparent", colour = NA), 
          plot.background = element_rect(fill = "transparent", colour = NA))
  # Heatmap
  # Data wrangling
  temp2 <- set_colnames(melt(temp), c("Label", "Sample", "value"))
  temp2$Label <- as.character(temp2$Label)
  temp2$Sample <- as.character(temp2$Sample)
  temp2$"Leading protein IDs" <- PG$"Leading protein IDs"[match(temp2$Label, PG$Label)]
  gr <- grsep(unlist(IDs.list), x = temp2$"Leading protein IDs", mode = "grepl")
  temp2$Colour <- ifelse(gr, "red", "grey")
  temp2$Size <- ifelse(gr, 3, 1)
  # Extract the order of the tips in the dendrograms
  # Order the levels according to their position in the clusters
  temp2$Xmin <- match(temp2$Sample, hord)-1
  temp2$Ymin <- match(temp2$Label, vord)-1
  Xscale <- length(unique(temp2$Sample))
  temp2$Label2 <- temp2$Label
  w <- which(nchar(temp2$Label2) > 25)
  temp2$Label2[w] <- paste0(substr(temp2$Label2[w], 1, 22), "...")
  # Create heatmap plot
  w1 <- which(temp2$Colour == "red")
  w2 <- which((temp2$Xmin == max(temp2$Xmin))&(temp2$Colour == "red"))
  nm <- "Heatmap with hierarchical clustering - Samples (PG-level)"
  heatmap.plot <- ggplot(temp2) + geom_tile(aes(x = Xmin, y = Ymin, fill = value), width = 1, height = 1) +
    geom_tile(data = temp2[w1,], aes(x = Xmin, y = Ymin, fill = value), width = 1, height = 1, colour = "red") +
    geom_text(data = temp2[w2,], aes(y = Ymin+0.75, label = Label2),
              x = length(unique(temp2$Sample))-0.45, colour = "red", hjust = 0, vjust = 0.5, cex = 2) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
          axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
          panel.background = element_rect(fill = "transparent", colour = NA)) +
    scale_fill_gradient2(low = "blue", high = "orange") + xlab("Sample") + ylab("Protein group")
  leg <- get_legend(heatmap.plot)
  heatmap.plot <- heatmap.plot + theme(legend.position = "none")
  # Combine
  htmp <- arrangeGrob(grobs = list(hdendro.plot, leg, heatmap.plot, vdendro.plot),
                      widths = c(1, 2*Xscale-2, 1, 1), heights = c(1, 5),
                      layout_matrix = rbind(c(NA, 1, NA, 2), c(3, 3, 3, 4)),
                      padding = 5)
  setwd(paste0(wd, "/Heatmaps"))
  ggsave(paste0(nm, ".jpeg"), htmp, dpi = 300, width = 10, height = 10, units = "in")
  ggsave(paste0(nm, ".pdf"), htmp)
  system(paste0("open \"", wd, "/Heatmaps/", nm, ".jpeg", "\""))
  msg <- "Look at the hierarchical dendrogram on the right of the heatmap: how many clusters do you want to extract?"
  KutOff <- as.numeric(dlg_input(msg)$res)
  while (is.na(KutOff)) {   KutOff <- as.numeric(dlg_input(paste0("I expect a number, so I'll ask again:\n", msg))$res) }
  clusters <- cutree(vcluster, KutOff)
  PG$"Hierarchical cluster" <- NA
  w <- which(names(clusters) %in% names(clusters))
  PG$"Hierarchical cluster"[w] <- clusters[match(PG$Label[w], names(clusters))]
  temp2$"Hierarchical cluster" <- clusters[match(temp2$Label, names(clusters))]
  temp2$"Hierarchical cluster"[w1] <- "Protein of interest"
  myColors <- setNames(c(colorRampPalette(c("blue", "orange"))(KutOff), "red"),
                       c(1:KutOff, "Protein of interest"))
  colScale <- scale_colour_manual(name = "colour", values = myColors, na.value = "black")
  nm2 <- "Heatmap with hierarchical clustering - Samples (PG-level, coloured clusters)"
  w3 <- which(temp2$Xmin == max(temp2$Xmin))
  heatmap.plot2 <- ggplot(temp2) + geom_tile(aes(x = Xmin, y = Ymin, fill = value), width = 1, height = 1) +
    geom_tile(data = temp2[w1,], aes(x = Xmin, y = Ymin, fill = value), width = 1, height = 1, colour = "red") +
    geom_text(data = temp2[w3,], aes(y = Ymin+0.75, label = Label2, colour = `Hierarchical cluster`),
              x = length(unique(temp2$Sample))-0.45, hjust = 0, vjust = 0.5, cex = 2) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
          axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
          panel.background = element_rect(fill = "transparent", colour = NA)) +
    scale_fill_gradient2(low = "blue", high = "orange") + colScale + xlab("Sample") + ylab("Protein group") +
    guides(colour = guide_legend(title = "Hierarchical cluster"))
  htmp2 <- arrangeGrob(grobs = list(hdendro.plot, leg, heatmap.plot2, vdendro.plot),
                       widths = c(1, 2*Xscale-2, 1, 1), heights = c(1, 5),
                       layout_matrix = rbind(c(NA, 1, NA, 2), c(3, 3, 3, 4)),
                       padding = 5)
  setwd(paste0(wd, "/Heatmaps"))
  ggsave(paste0(nm2, ".jpeg"), htmp2, dpi = 300, width = 10, height = 10, units = "in")
  ggsave(paste0(nm2, ".pdf"), htmp2)
  #system(paste0("open \"", wd, "/Heatmaps/", nm, ".jpeg", "\""))
  system(paste0("open \"", wd, "/Heatmaps/", nm2, ".jpeg", "\""))
  HKlust <- TRUE
  setwd(wd)
}

# Protein expression profiles:
if (length(Exp) > 1) {
  setwd(wd); suppressWarnings(dir.create("Profile plots"))
  g <- grep(topattern(PG.int.col), colnames(PG), value = TRUE)
  temp <- PG[, c("Leading protein IDs", "Protein IDs", "Label", g)]
  colnames(temp) <- gsub(topattern(PG.int.col), "", colnames(temp))
  temp <- set_colnames(melt(temp, id.vars = c("Leading protein IDs", "Protein IDs", "Label")),
                       c("Leading protein IDs", "Protein IDs", "Protein group", "Sample", "Top3"))
  temp <- temp[which(is.all.good(temp$Top3, 2)),]
  m <- median(temp$Top3)
  temp2 <- data.frame("Protein group" = unique(temp$"Protein group"), check.names = FALSE)
  temp2[, c("Leading protein IDs", "Protein IDs")] <- PG[match(temp2$"Protein group", PG$Label),
                                                         c("Leading protein IDs", "Protein IDs")]
  temp2$Top3 <- sapply(temp2$"Protein group", function(x) {
    x <- temp$Top3[which(temp$"Protein group" == x)]
    return(x[length(x)])
  })
  temp2$Sample <- rev(Exp)[1]
  temp2$Alpha <- (temp2$Top3 - m)^2
  if (protlist) {
    gA1 <- grsep(prot.IDs, x = temp$"Protein IDs") # Better to use "Protein IDs", not "Leading ...", for proteins of interest
    gA2 <- grsep(prot.IDs, x = temp2$"Protein IDs")
    if (length(gA1)) {
      gB1 <- grsep(prot.IDs, x = temp$"Protein IDs", inver = TRUE)
      gB2 <- grsep(prot.IDs, x = temp2$"Protein IDs", inver = TRUE)
      tempA <- temp[gA1,]
      temp <- temp[gB1,]
      temp2A <- temp2[gA2,]
      temp2 <- temp2[gB2,]
    }
  }
  ttl <- "Protein group Top3 expression profiles"
  plot <- ggplot(temp, aes(x = Sample, y = Top3, group = `Leading protein IDs`, colour = `Leading protein IDs`, text = `Protein group`)) +
    geom_line(linetype = "dashed") + geom_point(size = 1) + ggtitle(ttl) +
    theme_bw() +  theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  if (protlist) {
    if (length(gA1)) {
      plot <- plot + geom_line(data = tempA, linetype = "solid", colour = "red") + geom_point(data = tempA, colour = "red")
    }
  }
  setwd(paste0(wd, "/Profile plots/"))
  plotlyProfiles <- ggplotly(plot, tooltip = "text")
  saveWidget(partial_bundle(plotlyProfiles), paste0(ttl, ".html"))
  setwd(wd)
  system(paste0("open \"", wd, "/Profile plots/", ttl, ".html"))
  plot <- plot + geom_text(data = temp2, aes(label = `Protein group`, alpha = Alpha), hjust = 0, cex = 2.5)
  if (protlist) {
    if (length(gA1)) {
      plot <- plot + geom_text(data = temp2A, aes(label = `Protein group`), colour = "red", alpha = 1, hjust = 0, cex = 3)
    }
  }
  poplot(plot, 10, 18)
  ggsave(paste0(wd, "/Profile plots/", ttl, " (labels).jpeg"), plot, dpi = 150)
  ggsave(paste0(wd, "/Profile plots/", ttl, " (labels).pdf"), plot)
}
save.image("Backup.RData")
#load("Backup.RData")

# Export tables:
## Evidences
setwd(wd); suppressWarnings(dir.create("Tables"))
write.csv(file = paste0(wd, "/Tables/evidence.csv"), ev, row.names = FALSE)
## Modified Peptides
write.csv(file = paste0(wd, "/Tables/peptides_full.csv"), pep, row.names = FALSE)
## Full protein groups file
write.csv(file = paste0(wd, "/Tables/proteinGroups_full.csv"), PG, row.names = FALSE)

## Simplified peptides file
d <- pep[, c("Modified sequence", "Sequence", "Proteins" , "id", "Protein group IDs",
             grep("^Intensity\\.", colnames(pep), value = TRUE),
             grep("^Ratio\\.", colnames(pep), value = TRUE),
             "PEP")]
colnames(d) <- gsub("\\.", " ", colnames(d))
write.csv(d, file = paste0(wd, "/Tables/Modified Peptides.csv"), row.names = FALSE)
#
# Also sub-table of peptides matching to "interesting proteins":
if ((exists("prot.names"))&&(length(prot.names) > 0)&&(prot.names != "")) {
  temp <- setNames(strsplit(d$Proteins, ";"), d$`New Peptide ID`)
  temp <- melt(temp)
  w <- which(temp$value %in% db$"Protein ID"[which(db$`Common Name` %in% prot.names)])
  if (length(w) > 0) {
    w <- temp$L1[w]
    d1 <- d[match(w, d$`New Peptide ID`),]
    write.csv(d, file = paste0(wd, "/Tables/Modified Peptides - proteins of interest.csv"), row.names = FALSE)
  } else { warning("Not a single peptide for proteins in the list of interest, no sub-table will be created!") }
}
#
#
#
#
## Multi-tabs protein groups file as xlsx file
# Create formatting styles:
Styles <- list(IDs = createStyle(textDecoration = c("bold", "italic"), numFmt = "TEXT", fontSize = 11, valign = "center"))
# - Peptide and evidence counts and IDs
pepcolours <- colorRampPalette(c("thistle2", "thistle1"))(6)
evcolours <- colorRampPalette(c("slategray2", "slategray1"))(6)
Styles[["Global Pep. IDs"]] <- createStyle(textDecoration = "bold", bgFill = pepcolours[1], fgFill = pepcolours[1],
                                           numFmt = "TEXT", fontSize = 11, valign = "center")
Styles[["Global Pep. counts"]] <- createStyle(textDecoration = "bold", bgFill = pepcolours[2], fgFill = pepcolours[2],
                                              numFmt = "COMMA", fontSize = 11, valign = "center")
Styles[["Pep. IDs"]] <- createStyle(bgFill = pepcolours[3], fgFill = pepcolours[3],
                                    numFmt = "TEXT", fontSize = 11, valign = "center")
Styles[["Pep. counts"]] <- createStyle(bgFill = pepcolours[4], fgFill = pepcolours[4],
                                       numFmt = "COMMA", fontSize = 11, valign = "center")
Styles[["Global Ev. IDs"]] <- createStyle(textDecoration = "bold", bgFill = evcolours[1], fgFill = evcolours[1], numFmt = "TEXT", fontSize = 11, valign = "center")
Styles[["Ev. IDs"]] <- createStyle(bgFill = evcolours[3], fgFill = evcolours[3],
                                   numFmt = "TEXT", fontSize = 11, valign = "center")
Styles[["Ev. counts"]] <- createStyle(bgFill = evcolours[4], fgFill = evcolours[4],
                                      numFmt = "COMMA", fontSize = 11, valign = "center")
Styles[["Biot. Pep. IDs"]] <- createStyle(bgFill = pepcolours[5], fgFill = pepcolours[5],
                                          numFmt = "TEXT", fontSize = 11, valign = "center")
Styles[["Biot. Pep. counts"]] <- createStyle(bgFill = pepcolours[6], fgFill = pepcolours[6],
                                             numFmt = "COMMA", fontSize = 11, valign = "center")
Styles[["Biot. Ev. IDs"]] <- createStyle(bgFill = evcolours[5], fgFill = evcolours[5],
                                         numFmt = "TEXT", fontSize = 11, valign = "center")
Styles[["Biot. Ev. counts"]] <- createStyle(bgFill = evcolours[6], fgFill = evcolours[6],
                                            numFmt = "COMMA", fontSize = 11, valign = "center")
Styles[["Biot. Pep. %"]] <- "Biot. Pep. %"
# - Individual Expr
Styles[["Individual Expr"]] <- "Individual Expr"
# - Summary Expr
Styles[["Summary Expr"]] <- "Summary Expr"
# - Proteome Ruler
if (protrul) { Styles[["Proteome Ruler"]] <- "Proteome Ruler" }
# - Individual Ratios
Styles[["Individual Ratios"]] <- "Individual Ratios"
# - Regulated
Styles[["Regulated"]] <- createStyle(textDecoration = "bold", bgFill = "lightgoldenrod1", fgFill = "lightgoldenrod1",
                                     numFmt = "TEXT", fontSize = 11, valign = "center")
# - Annotations
if (Annotate) {
  annot <- c("InterPro", "Pfam", "PIRSF", "PROSITE")
  annot.col2 <- gsub("_names$", " names", annot.col)
  AnnotTbl <- data.frame(Name = c("GO", "Taxonomy", annot, "EMBL", "Other"),
                         Colour = RColorBrewer::brewer.pal(8, "Blues"), check.names = FALSE)
  AnnotTbl$Columns <- list(c("GO", "GO-ID"), c("Taxonomy", "TaxID"), NA, NA, NA, NA, "EMBL", NA)
  for (i in annot) { AnnotTbl$Columns[match(i, AnnotTbl$Name)] <- list(c(i, paste0(i, " names"))) }
  AnnotTbl$Columns[match("Other", AnnotTbl$Name)] <- list(annot.col2[which(!annot.col2 %in% unlist(AnnotTbl$Columns))])
  for (i in 1:nrow(AnnotTbl)) {
    Styles[[paste0(AnnotTbl$Name[i], " annotations")]] <- createStyle(bgFill = AnnotTbl$Colour[i], fgFill = AnnotTbl$Colour[i],
                                                                      numFmt = "TEXT", fontSize = 11, valign = "center")
  }
}
# - PEP
Styles[["PEP"]] <- createStyle(bgFill = "aliceblue", fgFill = "aliceblue",
                               numFmt = "SCIENTIFIC", fontSize = 11, valign = "center")
Styles[["Header 1"]] <- createStyle(textDecoration = "bold", halign = "left", valign = "bottom", wrapText = TRUE,
                                    border = "TopBottomLeftRight", borderColour = "black", borderStyle = "medium", textRotation = 50,
                                    numFmt = "TEXT", fontSize = 12)
Styles[["Header 2"]] <- createStyle(textDecoration = c("bold", "underline"), halign = "left", valign = "bottom", wrapText = TRUE,
                                    border = "TopBottomLeftRight", borderColour = "black", borderStyle = "medium", textRotation = 50,
                                    numFmt = "TEXT", fontSize = 12)
# - Proteins in list
Styles[["Protein list"]] <- createStyle(textDecoration = "bold", fontColour = "brown")
# - Summary tab
Styles[["Interlines"]] <- createStyle(fontSize = 12, halign = "left", valign = "center", textDecoration = "underline2", wrapText = FALSE)
Styles[["Small table - TEXT"]] <- createStyle(fontSize = 10, halign = "center", valign = "center", wrapText = FALSE, numFmt = "TEXT")
Styles[["Small table - NUMBER"]] <- createStyle(fontSize = 10, halign = "center", valign = "center", wrapText = FALSE, numFmt = "NUMBER")
Styles[["Header 3"]] <- createStyle(fontSize = 11, valign = "center", halign = "center", wrapText = TRUE,
                                    border = "TopBottomLeftRight", borderColour = "black", borderStyle = "thick", numFmt = "TEXT")
# Also colour scales and text decorations for conditional formatting
ColScaleList <- list(`Individual Expr` = c("darkred", "green"),
                     `Summary Expr` = c("darkred", "green"),
                     `Proteome Ruler` = c("darkslateblue", "goldenrod1"),
                     `Individual Ratios` = c("blue", "white", "red"))
DecoList <- list(`Summary Expr` = "bold")
#
# Columns to include, their final form and their order
corecol <- c("id", "Common Names", "Leading protein IDs", "Protein IDs", "Names", "Genes", "Gene names", "Mol. weight [kDa]")
pepevcol1 <- c("Peptides count", "Peptide IDs", "Razor peptide IDs", "Unique peptide IDs")
if (IsBioID) { pepevcol1 <- c(pepevcol1, "Biot. peptide IDs", "Biot. peptides count", "Biot. peptides [%]") }
pepevcol2 <- unique(c(pepevcol1, "Evidence IDs", grep("([Pp]eptide|[Ee]vidence)(s count| IDs)", colnames(PG), value = TRUE)))
pepcol2 <- grep("[Pp]eptide", pepevcol2, value = TRUE)
evcol2 <- grep("[Ee]vidence", pepevcol2, value = TRUE)
pepevcol2 <- c(grep("^Biot\\. ", pepcol2, value = TRUE, invert = TRUE), grep("^Biot\\. ", pepcol2, value = TRUE),
               grep("^Biot\\. ", evcol2, value = TRUE, invert = TRUE), grep("^Biot\\. ", evcol2, value = TRUE))
kol <- c(corecol, pepevcol2)
kol <- c(kol, "Mol. weight [kDa]", "Sequence coverage [%]", "PEP",
         grep(topattern(gsub(" - $", "", PG.int.col), start = FALSE), colnames(PG), value = TRUE)) # (so as to also capture the mean intensity column)
if (protrul) { kol <- c(kol, grep(topattern("log10(est. copies/cell) - "), colnames(PG), value = TRUE)) }
if (MakeRatios) {
  kol <- c(kol, grep(topattern(PG.rat.col), colnames(PG), value = TRUE),
           grep(topattern("Sequence coverage [%] -"), colnames(PG), value = TRUE))
  temp <- paste0(c("Enriched", "Regulated")[RatiosThresh_2sided+1], " - ")
  kol <- c(kol, grep(topattern(temp), colnames(PG), value = TRUE))
}
if (HKlust) { kol <- c(kol, "Hierarchical cluster") }
kol <- unique(kol[which(kol %in% colnames(PG))])
PG2 <- PG[,kol]
if (exists("annot.col")) { PG2[, gsub("_names$", " names", annot.col)] <- PG[,annot.col] }
if ("TaxID" %in% colnames(PG)) { PG2$TaxID <- as.numeric(PG2$TaxID) }
# (Insert empty column for visual sequence coverage)
w <- which(colnames(PG2) == "Sequence coverage [%]")
PG2 <- add_column(PG2, `First accession sequence coverage (peptides)` = "", .after = w-1)
# 
# Which columns are affected by each style
# - IDs
Which_col <- list(IDs = c("Leading protein IDs", "Protein IDs", "Names", "Genes", "Gene names", "Common Names"))
# - Evidences and peptides
Which_col[["Global Pep. IDs"]] <- grep("IDs$", pepevcol1, value = TRUE)
Which_col[["Global Pep. counts"]] <- grep("count$", pepevcol1, value = TRUE)
Which_col[["Pep. IDs"]] <- grep("[P,p]eptide IDs - ", pepevcol2, value = TRUE)
Which_col[["Pep. counts"]] <- grep("[P,p]eptides count - ", pepevcol2, value = TRUE)
Which_col[["Global Ev. IDs"]] <- "Evidence IDs"
Which_col[["Ev. IDs"]] <- grep("[Ee]vidence IDs - ", pepevcol2, value = TRUE)
Which_col[["Ev. counts"]] <- grep("[Ee]vidences count - ", pepevcol2, value = TRUE)
if (IsBioID) {
  Which_col[["Biot. pep. IDs"]] <- grep("^Biot\\.", Which_col[["Pep. IDs"]], value = TRUE)
  Which_col[["Pep. IDs"]] <- grep("^Biot\\.", Which_col[["Pep. IDs"]], value = TRUE, invert = TRUE)
  Which_col[["Biot. pep. counts"]] <- grep("^Biot\\.", Which_col[["Pep. counts"]], value = TRUE)
  Which_col[["Pep. counts"]] <- grep("^Biot\\.", Which_col[["Pep. counts"]], value = TRUE, invert = TRUE)
  Which_col[["Biot. ev. IDs"]] <- grep("^Biot\\.",  Which_col[["Ev. IDs"]], value = TRUE)
  Which_col[["Ev. IDs"]] <- grep("^Biot\\.",  Which_col[["Ev. IDs"]], value = TRUE, invert = TRUE)
  Which_col[["Biot. ev. counts"]] <- grep("^Biot\\.", Which_col[["Ev. counts"]], value = TRUE)
  Which_col[["Ev. counts"]] <- grep("^Biot\\.", Which_col[["Ev. counts"]], value = TRUE, invert = TRUE)
  Which_col[["Biot. pep. %"]] <- "Biot. Pep. %"
}
# - Individual Expr
Which_col[["Individual Expr"]] <- grep(topattern(PG.int.col), colnames(PG2), value = TRUE)
# - Summary Expr
Which_col[["Summary Expr"]] <- paste0("Mean ", gsub(" - $", "", PG.int.col))
# - Proteome Ruler
if (protrul) { Which_col[["Proteome Ruler"]] <- grep(topattern("log10(est. copies/cell) - "), colnames(PG2), value = TRUE) }
if (MakeRatios) {
  # - Individual Ratios
  Which_col[["Individual Ratios"]] <- grep(topattern(PG.rat.col), colnames(PG2), value = TRUE)
  # - Regulated
  Which_col[["Regulated"]] <- grep(paste0(c("Enriched", "Regulated")[RatiosThresh_2sided+1], " - "),
                                   colnames(PG2), value = TRUE)
}
# - Annotations
if (Annotate) { for (i in 1:nrow(AnnotTbl)) { Which_col[[paste0(AnnotTbl$Name[i], " annotations")]] <- AnnotTbl$Columns[[i]] } }
# - PEP
Which_col[["PEP"]] <- "PEP"
# Melt
Which_col <- melt(Which_col); colnames(Which_col) <- c("Col", "Grp")
stopifnot(nrow(Which_col) == length(unique(Which_col$Col)))
#
# Column for ordering data:
sortkol <- c(paste0("Mean ", gsub(" - $", "", PG.int.col)), "Biot. peptides count")[1+IsBioID]
AnalysisParam$"PG - Sorting column" <- sortkol
PG2 <- PG2[order(PG2[[sortkol]], decreasing = TRUE),]
# Filter for cases where number of characters is more than the Excel limit:
## Should only ever be character columns (lists of spectra IDs, sequences...)
ExcelMax <- 32767
w <- which(sapply(colnames(PG2), function(x) { class(PG2[[x]]) }) == "character")
if (length(w)) {
  for (i in w) { #i <- w[1]
    w1 <- which(nchar(PG2[[colnames(PG2)[i]]]) > ExcelMax)
    if (length(w1)) {
      PG2[[colnames(PG2)[i]]][w1] <- paste0(substr(PG2[[colnames(PG2)[i]]][w1], 1, ExcelMax-3), "...")
    }
  }
}
# Create workbook:
wb <- createWorkbook()
sheetnms <- c()
wrcol <- list()
sheetrows <- c()
sheetnm <- "PG - all columns"
sheetnms <- unique(c(sheetnms, sheetnm))
if (!sheetnm %in% names(wb)) { addWorksheet(wb, sheetnm) }
wb_tabs <- match(sheetnm, sheetnms)
wrcol[[sheetnm]] <- colnames(PG2)
writeData(wb, sheetnm, PG2)
sheetrows[[sheetnm]] <- nrow(PG2)
# Procedures common to all protein groups tabs:
for (sheetnm in sheetnms) { #sheetnm <- sheetnms[1]
  CS <- sapply(wrcol[[sheetnm]], function(x) {
    x <- match(x, Which_col$Col)
    if (is.na(x)) { return(createStyle()) } else { return(Styles[[Which_col$Grp[x]]]) }
  })
  for (i in 1:length(CS)) {
    if (class(CS[[i]]) == "Style") { addStyle(wb, sheetnm, CS[[i]], rows = (1:sheetrows[[sheetnm]])+1, cols = i, stack = TRUE) }
    if (class(CS[[i]]) == "character") {
      mtch <- Which_col$Grp[match(names(CS)[i], Which_col$Col)]
      if (!is.na(mtch)) {
        if (mtch %in% names(DecoList)) {
          addStyle(wb, sheetnm, createStyle(textDecoration = DecoList[[mtch]]), (1:sheetrows[[sheetnm]])+1, i, stack = TRUE)    
        }
        if (mtch %in% names(ColScaleList)) {
          conditionalFormatting(wb, sheetnm, rows = (1:sheetrows[[sheetnm]])+1, cols = i, type = "colourScale",
                                rule = NULL, style = ColScaleList[[mtch]], stack = TRUE)
        }
      }
    }
  }
  addStyle(wb, sheetnm, Styles[["Header 1"]], 1, 1:length(wrcol[[sheetnm]]))
  addStyle(wb, sheetnm, Styles[["Header 2"]], 1, which(wrcol[[sheetnm]] %in% corecol))
  setColWidths(wb, sheetnm, 1:length(wrcol[[sheetnm]]), 7) # Basic width
  w <- unique(c(which(wrcol[[sheetnm]] %in% c("PEP", corecol[which(corecol != "id")])),
                grep("Regulated - ", wrcol[[sheetnm]]))) # Specific columns
  if (length(w)) { setColWidths(wb, sheetnm, w, 15) }
  w <- which(wrcol[[sheetnm]] == "Common Names")  # Common name column
  if (length(w)) { setColWidths(wb, sheetnm, w, 50) }
  w <- grep("([Ee]vidence|[Pp]eptide) IDs", pepevcol2) # Peptide and protein IDs column
  if (length(w)) { setColWidths(wb, sheetnm, w, 10) }
  setRowHeights(wb, sheetnm, 1:(sheetrows[[sheetnm]]+1), c(200, rep(15, sheetrows[[sheetnm]])))
  freezePane(wb, sheetnm, firstRow = TRUE)
}
if (protlist) {
  prot.list <- unlist(IDs.list)
  if (length(prot.list)) {
    w <- grsep(prot.list, x = PG2$"Protein IDs")
    if (length(w)) {
      for (sheetnm in grep("^PG - ", sheetnms, value = TRUE)) { #sheetnm <- grep("^PG - ", sheetnms, value = TRUE)[1]
        addStyle(wb, sheetnm, Styles[["Protein list"]], w+1, 1:length(wrcol[[sheetnm]]), gridExpand = TRUE, stack = TRUE)
      }
    }
  }
}
if (HKlust) {
  sheetnm <- sheetnms[1]
  w <- which(wrcol[[sheetnm]] == "Hierarchical cluster")
  conditionalFormatting(wb, sheetnms[1], w, 1:nrow(PG2),
                        style = c("blue", "orange"),
                        type = "colourScale")
}
# Special summary tab
nk <- c(ncol(Exp_summary), ncol(AA_biases))+2
nr <- c(nrow(Exp_summary), nrow(AA_biases))+1
temp <- AA_biases
for (i in grep("^% - |^Ratio$", colnames(temp))) { temp[[i]] <- signif(temp[[i]], 3) }
temp1 <- as.matrix(rbind(colnames(Exp_summary), Exp_summary))
temp2 <- as.matrix(rbind(colnames(temp), temp))
temp1 <- cbind(NA, NA, temp1)
temp2 <- cbind(NA, NA, temp2)
if (nk[1] < nk[2]) { temp1 <- cbind(temp1, matrix("", ncol = nk[2]-nk[1], nrow = nr[1])) }
if (nk[1] > nk[2]) { temp2 <- cbind(temp2, matrix("", ncol = nk[1]-nk[2], nrow = nr[2])) }
mnk <- max(nk)
tags <- c("Summary", "Amino acid compositional biases: observed peptides /vs/ parent proteome)")
temp <- rbind(rep(NA, mnk),
              c(NA, tags[1], rep(NA, mnk-2)),
              temp1,
              rep(NA, mnk),
              rep(NA, mnk),
              c(NA, tags[2], rep(NA, mnk-2)),
              temp2,
              rep(NA, mnk))
colnames(temp) <- NULL
w <- sapply(tags, function(x) { which(temp[,2] == x) })
sheetnm <- "Summary"
sheetnms <- unique(c(sheetnms, sheetnm))
if (!sheetnm %in% names(wb)) { addWorksheet(wb, sheetnm) }
wb_tabs <- match(sheetnm, sheetnms)
wrcol[[sheetnm]] <- NA
writeData(wb, sheetnm, temp, colNames = FALSE, rowNames = FALSE, keepNA = FALSE)
sheetrows[[sheetnm]] <- nrow(temp)
setColWidths(wb, sheetnm, 1, 5)
setColWidths(wb, sheetnm, 2:mnk, 15)
for (i in 1:2) { #i <- 1
  r <- w[i]+1
  k <- nk[i]
  #int <- max(c(1, c(0, w)[i])):w[i]
  addStyle(wb, sheetnm, Styles[["Interlines"]], w[i], 2, stack = TRUE)
  addStyle(wb, sheetnm, Styles[["Header 3"]], r, 3:k, stack = TRUE)
  for (j in 3:k) {
    test <- sum(is.na(suppressWarnings(as.numeric(as.character(temp[(r+1):(r+nr[i]-1), j]))))) == 0
    addStyle(wb, sheetnm, Styles[[paste0("Small table - ", c("TEXT", "NUMBER")[test+1])]], (r+1):(r+nr[i]-1), j, stack = TRUE)
  }
  # Borders:
  addStyle(wb, sheetnm, createStyle(border = "right", borderColour = "black", borderStyle = "thin"), (r+1):(r+nr[i]-1), 2, stack = TRUE)
  addStyle(wb, sheetnm, createStyle(border = "left", borderColour = "black", borderStyle = "thin"), (r+1):(r+nr[i]-1), k+1, stack = TRUE)
  addStyle(wb, sheetnm, createStyle(border = "top", borderColour = "black", borderStyle = "thin"), r+nr[i], 3:nk[i], stack = TRUE)
  setRowHeights(wb, sheetnm, r, 50)
}
saveWorkbook(wb, paste0(wd, "/Tables/proteinGroups.xlsx"), TRUE)
#shell(paste0("\"", wd, "/Tables/proteinGroups.xlsx\""))
#
# Create custom conditional formatting with Niklas' first macro
# In case Excel complained about the workbook needing repairing, save first then run the lines below:
temp1 <- data.frame(Accession = sapply(strsplit(PG2$"Leading protein IDs", ";"), function(x) { x[[1]] }),
                    `Peptide IDs` = PG2$"Peptide IDs", check.names = FALSE)
temp1$Sequence <- db$Sequence[match(temp1$Accession, db$"Protein ID")]
temp2 <- data.frame(ID = pep$id, Sequence = pep$Sequence)
setwd(paste0(wd, "/Tables"))
write.csv(file = "temp1.csv", temp1, row.names = FALSE)
write.csv(file = "temp2.csv", temp2, row.names = FALSE)
path <- as.data.frame(library()$results)
path <- normalizePath(path$LibPath[match("proteoCraft", path$Package)])
wdn <- normalizePath(wd)
for (sheetnm in "PG - all columns") { #sheetnm <- "PG - all columns"
  cmd <- paste0("WScript \"", path, "\\proteoCraft\\data\\sequ_color_only_VBS.vbs\" ", #full path of vbaFromVbs.vbs
                "/fpathXlsm:\"", path, "\\proteoCraft\\data\\sequ_color_only_VBA.xlsm\" ", #full path of SequColorVBA.xlsm
                "/fpathPeptides:\"", wdn, "\\Tables\\temp2.csv\" ", #fullpath of peptides.csv
                "/fpathProt:\"", wdn, "\\Tables\\temp1.csv\" ", #fullpath of proteins.csv
                "/fpathGrps:\"", wdn, "\\Tables\\proteinGroups.xlsx\" ", #fullpath of proteinGroups.xlsx
                "/shGrps:\"", sheetnm, "\" ", #sheet name of proteinsGroups.xlsx
                "/colPepID:\"ID\" ", #col name of peptide ID
                "/colPepSequ:\"Sequence\" ", #col name of peptide sequence
                "/colProtAcc:\"Accession\" ", #col name of Accession No.
                "/colProtPepIDs:\"Peptide IDs\" ", #col name of list of peptide IDs
                "/colProtSequ:\"Sequence\" ", #col name of protein sequence
                "/colGrpsProtIDs:\"Leading protein IDs\" ",
                "/colGrpsSequ:\"First accession sequence coverage (peptides)\"")
  #message(cmd) #(for if you need to inspect the command)
  system(command = cmd, wait = TRUE)
}
unlink(paste0(wd, "/Tables/temp1.csv"))
unlink(paste0(wd, "/Tables/temp2.csv"))
setwd(wd)
#
# Format specific columns as progress bars with Niklas' second macro:
seqcolkol <- paste(grep("^Biot\\. peptides \\[%\\]|^Sequence coverage \\[%\\]", colnames(PG2), value = TRUE), collapse = ";")
for (sheetnm in "PG - all columns") { #sheetnm <- "PG - all columns"
  cmd <- paste0("WScript \"", path, "\\proteoCraft\\data\\Format_as_bar_VBS.vbs\" ",
                "/fpathXlsm:\"", path, "\\proteoCraft\\data\\Format_as_bar_VBA.xlsm\" ", #fullpath of SequColorVBA.xlsm
                "/fpathGrps:\"", wdn, "\\Tables\\proteinGroups.xlsx\" ", #fullpath of proteinGroups.xlsx
                "/shGrps:\"", sheetnm, "\" ", #sheet name of proteinsGroups.xlsx
                "/colFormatBar:\"", seqcolkol, "\"" #col names for bar-format separated with ";"
  )
  system(command = cmd, wait = TRUE)
}
# Open Excel file:
#shell(paste0("\"", wd, "/Tables/proteinGroups.xlsx\""))
#
# Save session info
setwd(wd)
if (!dir.exists("Workflow control")) { dir.create("Workflow control")}
writeLines(capture.output(sessionInfo()), "Workflow control/sessionInfo.txt")
#
# Save copy of this script to local work directory
script.dir <- getActiveDocumentContext()
script.dir <- script.dir$path
file.copy(script.dir, wd, overwrite = TRUE)

save.image("Backup.RData")
#load("Backup.RData")
openwd()
