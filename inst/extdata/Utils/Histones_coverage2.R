wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)
require(openxlsx2)
require(svDialogs)
require(parallel)
require(proteoCraft) # Update me!

# Create parallel processing cluster
N.clust <- detectCores()-1
a <- 1
tst <- try(clusterExport(parClust, "a", envir = environment()), silent = TRUE)
if ("try-error" %in% class(tst)) {
  try(stopCluster(parClust), silent = TRUE)
  parClust <- makeCluster(N.clust, type = "SOCK")
}

#
FP_search <- rstudioapi::selectDirectory("Select a FragPipe search folder", path = wd)
fls <- list.files(FP_search, full.names = TRUE, recursive = FALSE)
FP_WorkflowFl <- grep("/fragpipe\\.workflow$", fls, value = TRUE)
FP_ManifestFl <- grep("/fragpipe-files\\.fp-manifest$", fls, value = TRUE)
if (!length(FP_WorkflowFl)) { stop("Could not locate FragPipe workflow file in folder!") }
if (length(FP_WorkflowFl) > 1) { stop("Too many FragPipe workflow files in folder, only 1 expected!") }
if (!length(FP_ManifestFl)) { stop("Could not locate FragPipe files manifest in folder!") }
if (length(FP_ManifestFl) > 1) { stop("Too many FragPipe files manifests in folder, only 1 expected!") }
FP2MQ <- FP_to_MQ(FP_WorkflowFl, FP_ManifestFl, cl = parClust)
pep <- FP2MQ$Evidence
Exp <- unique(pep$Experiment)
mods <- FP2MQ$PTMs
map <- FP2MQ$FracMap
dbFl <- gsub("database.db-path=", "", grep("^database.db-path=", FP2MQ$WorkFlow, value = TRUE))
dbFl <- gsub("\\\\\\\\", "/", dbFl)
dbFl <- gsub("\\\\", "", dbFl)
file.exists(dbFl)

RenameExp <- c(TRUE, FALSE)[match(dlg_message("Should we rename individual samples?", "yesno")$res,
                                  c("yes", "no"))]
if (RenameExp) {
  smplsMapFl <- paste0(FP_search, "/Samples_map.csv")
  ok <- FALSE
  kount <- 0
  while (!ok) {
    if (file.exists(smplsMapFl)) {
      cat(paste0("Samples map template detected at:\n", smplsMapFl, "\n"))
    } else {
      tmp <- data.frame(Old = Exp, New = "?")
      write.csv(tmp, smplsMapFl, row.names = FALSE)
      cat(paste0("Samples map template created at:\n", smplsMapFl, "\n"))
    }
    cmd <- paste0("open \"", smplsMapFl, "\"")
    system(cmd)
    msg <- c("Make your edits in the table then click ok",
             "The last version did not have the expected columns (\"Old\" and \"New\"),\n
             please make your edits in the table then click ok")[(kount > 0)+1]
    dlg_message("Make your edits in the table then click ok", "ok")
    smplsMap <- read.csv(smplsMapFl)
    ok <- sum(!c("Old", "New") %in% colnames(smplsMap)) == 0
    kount <- kount + 1
  }
  pep$Old_Experiment <- pep$Experiment
  pep$Experiment <- smplsMap$New[match(pep$Old_Experiment, smplsMap$Old)]
  Exp <- unique(pep$Experiment)
}

Remove0Int <- c(TRUE, FALSE)[match(dlg_message("Should we remove peptides with 0 intensity values?", "yesno")$res,
                                   c("yes", "no"))]
if (Remove0Int) { pep <- pep[which(pep$Intensity > 0),] }

prsDBfl <- paste0(wd, "/Parsed DB.csv")
#if (!file.exists(prsDBfl)) {
db <- Format.DB(dbFl, cl = parClust)
write.csv(db, prsDBfl, row.names = FALSE)
#} else {
#  db2 <- read.csv(prsDBfl, check.names = FALSE)
#}
db <- db[grep("^>rev_", db$Header, invert = TRUE),]

Hist <- grep("histone", db$Header, ignore.case = TRUE, value = TRUE)
Hist <- grep("ase(,|\\.| )", Hist, ignore.case = TRUE, value = TRUE, invert = TRUE)
Hist <- grep("-(binding|like)", Hist, ignore.case = TRUE, value = TRUE, invert = TRUE)
Hist <- grep("chaperone", Hist, ignore.case = TRUE, value = TRUE, invert = TRUE)
Hist <- grep("ethyl", Hist, ignore.case = TRUE, value = TRUE, invert = TRUE)
Hist <- grep("factor", Hist, ignore.case = TRUE, value = TRUE, invert = TRUE)
Hist <- grep("non-histone", Hist, ignore.case = TRUE, value = TRUE, invert = TRUE)
Hist <- grep("\\(fragment\\)", Hist, ignore.case = TRUE, value = TRUE, invert = TRUE)
HistIDs <- db$`Protein ID`[match(Hist, db$Header)]

#unique(db$`Common Name`[match(Hist, db$Header)])
histDB <- db[match(HistIDs, db$`Protein ID`),]
writeFasta(histDB, paste0(wd, "/H_sapiens_all_histones.fasta"))
histDB$ID_Nm <- apply(histDB[, c("Protein ID", "Common Name")], 1, paste, collapse = " - ")
histDB$ID_Nm <- gsub("/", "_", histDB$ID_Nm)
clusterExport(parClust, "histDB", envir = environment())

# With FP, we need to guess modification names...
# But here, we have a home-made table available:
modTblFl <- "...My_Dataset/histone_modification_masses.xlsx"
modTbl <- read_xlsx(modTblFl)
colnames(modTbl)[1] <- "Priority"
modTbl <- modTbl[which(!is.na(modTbl$"molecular weight")),]
modTbl$`molecular weight` <- as.numeric(modTbl$`molecular weight`)
modTbl$Priority[which(is.na(modTbl$Priority))] <- ""
modTbl <- modTbl[which(modTbl$Priority != "Not analysed"),]
tmp <- listMelt(strsplit(modTbl$`modified residues`, ", *"), 1:nrow(modTbl))
tmp$Mass_Shift <- modTbl$`molecular weight`[tmp$L1]
tmp$Name <- modTbl$Modification[tmp$L1]
tmp <- aggregate(tmp$Name, list(tmp$value, tmp$Mass_Shift), function(x) { paste(sort(unique(x)), collapse = "/") })
colnames(tmp) <- c("value", "Mass_Shift", "Name")
tmp2 <- listMelt(mods$AA, 1:nrow(mods))
tmp2$Mass_Shift <- mods$`Mass delta`[tmp2$L1]
tmp2$Name <- mods$`Full name`[tmp2$L1]
tmp2$new_Name <- tmp2$Name
w <- which(tmp2$Mass_Shift %in% tmp$Mass_Shift)
tmp2$new_Name[w] <- tmp$Name[match(tmp2$Mass_Shift[w], tmp$Mass_Shift)]
# Now we want to rename our modifications
#w <- which(!1:nrow(mods) %in% tmp2$L1)
w <- which(1:nrow(mods) %in% tmp2$L1)
mods$`new Full name` <- mods$`Full name`
mods$`new Full name`[w] <- tmp2$new_Name[match(w, tmp2$L1)]
g <- grep("\\(", pep$`Modified sequence_verbose`)
if (length(g)) {
  tmp <- data.frame(Seq = pep$`Modified sequence_verbose`[g])
  tmp$Seq <- strsplit(tmp$Seq, "\\(|\\)")
  tmp$L <- sapply(tmp$Seq, length)
  l <- max(tmp$L)
  tmp2 <- as.data.frame(t(apply(tmp[, c("Seq", "L")], 1, function(x) {
    x <- c(x[[1]], rep("", l-x[[2]]))
  })))
  k <- (1:((l-1)/2))*2
  for (i in k) {
    w <- which(tmp2[[i]] == "Carbamidomethyl") # No need to show that one!
    tmp2[w, i] <- ""
    w <- which(tmp2[[i]] %in% mods$`Full name`)
    tmp2[w, i] <- paste0("(", mods$`new Full name`[match(tmp2[w, i], mods$`Full name`)], ")")
  }
  tmp2 <- do.call(paste, c(tmp2, sep = ""))
  pep$`Modified sequence_verbose`[g] <- tmp2
}

#
tst <- gsub("^_[A-Z]+_$|^_[A-Z]*\\(|\\)[A-Z]*_$", "", gsub("\\)[A-Z]+\\(", "___", pep$`Modified sequence_verbose`))
tst <- unlist(strsplit(tst, "___"))
tst <- aggregate(tst, list(tst), length)
nrow(tst) # Number of identified PTMs
length(unique(mods$`new Full name`)) # Number of searched PTMs
stopifnot(sum(!tst$Group.1 %in% mods$`new Full name`) == 0)# Checking that the names match

# Check peptide-to-protein assignment, and isolate histones peptides
### NB: This doesn't take into consideration N-terminal methionines!!!
msg <- "Update FragPipe's original protein-to-peptides assignments?"
Update_Prot_matches <- c(TRUE, FALSE)[match(dlg_message(msg, "yesno")$res, c("yes", "no"))]
if (Update_Prot_matches) {
  uSeq <- unique(pep$Sequence)
  evmatch <- ProtMatch2(uSeq, db, cl = parClust)
  evmatch2 <- ProtMatch2(uSeq, histDB, cl = parClust)
  pep$Proteins_FP <- pep$Proteins
  pep$Proteins <- evmatch$Proteins[match(pep$Sequence, evmatch$Sequence)]
  pep$"Histone(s)" <- evmatch2$Proteins[match(pep$Sequence, evmatch2$Sequence)]
  #View(pep[, c("Proteins_FP", "Proteins", "Histone(s)", "Modified sequence_verbose")])
  #sum(pep$Proteins_FP != pep$Proteins)
  #sum(pep$Proteins_FP == pep$Proteins)
} else {
  tmp <- strsplit(pep$Proteins, ";")
  tmp <- parSapply(parClust, tmp, function(x) { x[which(x %in% histDB$"Protein ID")] })
  pep$"Histone(s)" <- parSapply(parClust, tmp, paste, collapse = ";")
}
tmp <- pep[, c("Histone(s)", "Proteins")]
tmp$"Proteins" <- strsplit(tmp$"Proteins", ";")
tmp$"Histone(s)" <- strsplit(tmp$"Histone(s)", ";")
tst <- apply(tmp[, c("Histone(s)", "Proteins")], 1, function(x) {
  sum(!x[[1]] %in% x[[2]])
})
w <- which(tst > 0)
#View(pep[w, c("Modified sequence_verbose", "Histone(s)", "Proteins")])
#View(pep[grep("phospho", pep$`Modified sequence_verbose`, ignore.case = TRUE), c("Modified sequence_verbose", "Histone(s)", "Proteins")])

tst <- aggregate(pep$"Histone(s)", list(pep$"Histone(s)"), length)
View(tst)
w <- which(sapply(strsplit(pep$`Histone(s)`, ";"), function(x) {
  x <- unlist(x)
  x <- x[which(x != "")]
  (length(x) > 1)||(!is.na(x))
}))
histPep <- pep[w,]
tst <- gsub("^_[A-Z]+_$|^_[A-Z]*\\(|\\)[A-Z]*_$", "", gsub("\\)[A-Z]+\\(", "___", histPep$`Modified sequence_verbose`))
tst <- unlist(strsplit(tst, "___"))
tst <- aggregate(tst, list(tst), length)
nrow(tst) # Number of identified PTMs
stopifnot(sum(!unique(unlist(strsplit(histPep$`Histone(s)`, ";"))) %in% histDB$`Protein ID`) == 0)

histDir <- paste0(FP_search, "/Histones_coverage_maps")
if (!dir.exists(histDir)) { dir.create(histDir) }
Exp <- unique(pep$Experiment)

# Check PTMs
tst <- gsub("_[A-Z]+_|_[A-Z]*\\(|\\)[A-Z]*_", "", gsub("\\)[A-Z]*\\(", "___", histPep$`Modified sequence_verbose`))
tst <- unlist(strsplit(tst, "___"))
tst <- aggregate(tst, list(tst), length)
tst <- tst[order(tst$x, decreasing = TRUE),]
colnames(tst) <- c("PTM", "Nb. of peptides")
View(tst)

tmpPep <- histPep[, c("Experiment", "Histone(s)", "Intensity", "Modified sequence", "Modified sequence_verbose")]
clusterExport(parClust, list("tmpPep", "histDir", "Hist", "Coverage"), envir = environment())
#clusterExport(parClust, "Coverage", envir = environment())
for (e in Exp) { #e <- Exp[1]
  w <- which(pep$Experiment == e)
  #w <- which(tmpPep$Experiment %in% Exp)
  if (length(w)) {
    pp <- tmpPep[w,]
    pp <- aggregate(pp$Intensity, list(pp$`Modified sequence`, pp$`Modified sequence_verbose`), sum, na.rm = TRUE)
    colnames(pp) <- c("Modified sequence", "Modified sequence_verbose", "Intensity")
    pp$"Histone(s)" <- tmpPep$"Histone(s)"[match(pp$`Modified sequence`, tmpPep$`Modified sequence`)]
    tmpPP <- listMelt(strsplit(pp$"Histone(s)", ";"), 1:nrow(pp), c("Histone", "row"))
    tmpPP <- aggregate(tmpPP$row, list(tmpPP$"Histone"), c)
    tmpPP <- setNames(lapply(tmpPP$x, function(x) {
      pp[x, c("Modified sequence", "Modified sequence_verbose", "Intensity")]
    }), tmpPP$Group.1)
    clusterExport(parClust, list("e", "tmpPP"), envir = environment())
    # tst <- sapply(names(tmpPP), function(x) {
    #   g <- grep("nyl", tmpPP[[x]]$`Modified sequence_verbose`, value = TRUE)
    #   if (length(g)) { stop(cat(x, "\n", paste(g, collapse = "\n"), "\n")) }
    # })
    # View(tmpPP[["A0A2R8Y619"]])
    w <- which(histDB$`Protein ID`[match(Hist, histDB$Header)] %in% names(tmpPP))
    a <- parLapply(parClust, Hist[w], function(hist) { #hist <- Hist[w][1]
    #for (hist in Hist[w]) {
      #grep("Q07133", Hist[w])
      #hist <- grep("A0A2R8Y619", Hist[w], value = TRUE)
      m <- match(hist, histDB$Header)
      nm <- histDB$ID_Nm[m]
      sq <- setNames(histDB$Sequence[m], histDB$ID_Nm[m])
      tmp <- tmpPP[[histDB$`Protein ID`[m]]]
      Coverage(sq, tmp$`Modified sequence_verbose`, "Align2", FALSE, FALSE,
               title = paste0(nm, " coverage - ", e), save.path = paste0(histDir, "/", nm, " - ", e), save = "jpeg")
    })
  }
}
#openwd(histDir)
