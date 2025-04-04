wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)
require(rstudioapi)
require(qs)
require(plyr)
require(data.table)
require(openxlsx2)
require(svDialogs)
require(parallel)
require(limma)

updateMe <- c(TRUE, FALSE)[match(dlg_message("Update the proteoCraft package?", "yesno")$res, c("yes", "no"))]
if (updateMe) {
  dir <- "...Projects_Doc_Folder/Mass_Spec/proteoCraft_package"
  fls <- list.files(dir, "proteoCraft_[0-9]+\\.[0-9]+\\.[0-9]+\\.[0-9]+\\.tar\\.gz$", full.names = TRUE)
  vers <- as.data.frame(t(sapply(strsplit(gsub(".*/proteoCraft_|\\.tar\\.gz$", "", fls), "\\."), function(x) { as.integer(unlist(x)) })))
  vers$File <- fls
  o <- order(vers$V1, vers$V2, vers$V3, vers$V4, decreasing = TRUE)
  vers <- vers[o,]
  dflt <- vers$File[1]
  fl <- rstudioapi::selectFile(path = dflt)
  if ((!is.na(fl))&&(file.exists(fl))&&(length(fl) == 1)&&(is.character(fl))&&(grepl("\\.tar\\.gz$", fl))) {
    try({
      detach(package:proteoCraft, unload = TRUE)
      remove.packages("proteoCraft")
    }, silent = TRUE)
    install.packages(fl, repos = NULL, type = "source")
  }
}
require(proteoCraft)

# This file is hard-coded in here: this is the table of mass-shift-to-PTM-name matches.
# You can edit it but PLEASE stick to the current layout/logic or the script will BREAK!
modTblFl <- "...My_Dataset/histone_modification_masses.xlsx"
stopifnot(file.exists(modTblFl))

# Create parallel processing cluster
N.clust <- detectCores()-1
a <- 1
tst <- try(clusterExport(parClust, "a", envir = environment()), silent = TRUE)
if ("try-error" %in% class(tst)) {
  try(stopCluster(parClust), silent = TRUE)
  parClust <- makeCluster(N.clust, type = "SOCK")
}

saveFun2 <- function(x, file) {
  tmp <- paste0("qs::qsavem(", deparse(substitute(x)),
                ", file = '", file, "', nthreads = max(c(parallel::detectCores()-1, 1)))")
  #cat(tmp)
  eval(parse(text = tmp))
}
saveImgFun2 <- function(file) { # More elegant rewriting - I think
  args2 <- list(file = file)
  obj <- objects(.GlobalEnv)
  obj <- setNames(lapply(obj, get, envir = .GlobalEnv), obj)
  args2$x <- obj
  do.call(qs::qsave, args2)
}
loadFun2 <- function(file) {
  qs::qload(file, env = globalenv(), nthreads = max(c(parallel::detectCores()-1, 1)))
}

parDir <- dirname(rstudioapi::getSourceEditorContext()$path)
# We will use the version of the Coverage function which was saved in the same folder
useLocals <- FALSE
if (useLocals) { # In case I want to use a local version of these functions
  srcFuns <- paste0(parDir, "/", c("Coverage",
                                   "DIANN_to_MQ",
                                   ".DIANN_to_MQ_server2",
                                   ".DIANN_to_MQ_ui2", 
                                   "FP_to_MQ"), ".R")
  for (Src in srcFuns) { if (file.exists(Src)) { source(Src) } }
}

if (exists("inDir")) { dflt <- inDir } else { dflt <- parDir }
inDir <- rstudioapi::selectDirectory("Select folder in which FragPipe result folder(s) are located",
                                     path = dflt)
FP_searches <- list.dirs(inDir, recursive = FALSE)
w <- sapply(FP_searches, function(FP_search) { #FP_search <- FP_searches[1]
  #print(FP_search)
  fls <- list.files(FP_search, full.names = TRUE, recursive = FALSE)
  FP_WorkflowFl <- grep("/fragpipe\\.workflow$", fls, value = TRUE)
  FP_ManifestFl <- grep("/fragpipe-files\\.fp-manifest$", fls, value = TRUE)
  OK <- TRUE
  if (!length(FP_WorkflowFl)) {
    #warning(paste0(FP_search, ": could not locate FragPipe workflow file in folder!"))
    OK <- FALSE
  }
  if (length(FP_WorkflowFl) > 1) {
    #warning(paste0(FP_search, ": too many FragPipe workflow files in folder, only 1 expected!"))
    OK <- FALSE
  }
  if (!length(FP_ManifestFl)) {
    #warning(paste0(FP_search, ": could not locate FragPipe files manifest in folder!"))
    OK <- FALSE
  }
  if (length(FP_ManifestFl) > 1) {
    #warning(paste0(FP_search, ": too many FragPipe files manifests in folder, only 1 expected!"))
    OK <- FALSE
  }
  return(OK)
})
FP_searches <- FP_searches[w]
l <- length(FP_searches)
stopifnot(l > 0)  
if (l > 1) {
  tmp <- gsub(".*/", "", FP_searches)
  tmp2 <- svDialogs::dlg_list(tmp, tmp, TRUE, title = "Select FragPipe output folders to analyze:")$res  
  FP_searches <- FP_searches[match(tmp2, tmp)]
}
l <- length(FP_searches)
stopifnot(l > 0)

dflt <- inDir
if (length(FP_searches) == 1) { dflt <- FP_searches }
dstDir <- rstudioapi::selectDirectory("Select results directory", path = dflt)
backupFl <- paste0(inDir, "/Backup.RData")

FP2MQs <- setNames(lapply(FP_searches, function(FP_search) { #FP_search <- FP_searches[1]
  print(FP_search)
  fls <- list.files(FP_search, full.names = TRUE, recursive = FALSE)
  FP_WorkflowFl <- grep("/fragpipe\\.workflow$", fls, value = TRUE)
  FP_ManifestFl <- grep("/fragpipe-files\\.fp-manifest$", fls, value = TRUE)
  FP2MQ <- try(FP_to_MQ(FP_WorkflowFl, FP_ManifestFl, cl = parClust), silent = TRUE)
  if (!"try-error" %in% class(FP2MQ)) {
    FP2MQ$Outcome <- TRUE
    Exp <- unique(FP2MQ$Evidence$Experiment)
    if (!length(Exp)) {
      tmp <- unique(FP2MQ$Evidence$"Raw file")
      if (length(tmp) != 1) { tmp <- "Exp1" }
      FP2MQ$Evidence$Experiment <- Exp <- tmp
    }
    dbFl <- gsub("database.db-path=", "", grep("^database.db-path=", FP2MQ$WorkFlow, value = TRUE))
    dbFl <- gsub("\\\\\\\\", "/", dbFl)
    dbFl <- gsub("\\\\", "", dbFl)
    FP2MQ$Exp <- Exp
    FP2MQ$dbFl <- dbFl
  } else { FP2MQ <- list(Outcome = FALSE) }
  return(FP2MQ)
}), FP_searches)
tst <- sapply(FP_searches, function(FP_search) { FP2MQs[[FP_search]]$Outcome })
w <- which(!tst)
if (length(w)) {
  warning(paste0("The following director", c("y", "ies")[(length(w) > 1)+1], " could not be parsed (ignore this warning if they are not FP output directories):\n",
                 paste0(" - ", FP_searches[w], "\n", collapse = "")), "\n")
  FP_searches <- FP_searches[which(tst)]
  FP2MQs <- FP2MQs[FP_searches]
}
#
# Check for results consistency
PSMs <- setNames(lapply(FP_searches, function(FP_search) { #FP_search <- FP_searches[1]
  x <- FP2MQs[[FP_search]]$Evidence
  if (!"Modified sequence_verbose" %in% colnames(x)) {
    tmp <- aggregate(1:nrow(x), list(x$"Modified sequence"), list)
    g <- grep("\\(", tmp$Group.1)
    if (length(g)) {
      tmp$ModSeq <- strsplit(tmp$Group.1, "\\(|\\)")
      mods <- FP2MQs[[FP_search]]$PTMs
    }
  }
  return(x)
}), FP_searches)
pep <- plyr::rbind.fill(PSMs)
pep$`Raw file name` <- gsub(".*/", "", pep$`Raw file`)
raw.files <- unique(pep$`Raw file name`)
pep$Bruker_runID <- as.integer(gsub(".*_|\\.d$", "", pep$`Raw file name`))
pep$Scan_ID <- do.call(paste, c(pep[, c("Bruker_runID", "MS/MS scan number")], sep = "___"))
pep$tmp <- round(pep$"m/z", 3)
pep$Scan_ID_MZ <- do.call(paste, c(pep[, c("Scan_ID", "tmp")], sep = "___"))
pep$tmp <- NULL

# The part below is commented as it doesn't work!!!
# The idea was to ensure consistency when combining different searches...
# But when analyzing a single search, it was still removing PSMs!
# The reason: one MS2 can of course yield several PSMs, with close or even identical m/z...
# and all can be valid!!!
# One would need to also check which peaks were matched:
# - get all sequences and all MS2 spectra involved
# - calculate new scores
# - take the best match
# This would be very work intensive to write so has been dropped for now.
#
# When combining different searches, we want to make sure the same spectrum is not assigned
# different peptidoforms.
# Note: this can only work in DDA mode!
# tmp <- as.data.table(pep[, c("Scan_ID_MZ", "Modified sequence_verbose")])
# tmp <- tmp[, list(ModSeq = list(`Modified sequence_verbose`)), by = list(Scan_ID_MZ = Scan_ID_MZ)]
# tmp <- as.data.frame(tmp)
# tmp$L <- sapply(tmp$ModSeq, length)
# mx <- max(tmp$L)
# if (mx > 1) {
#   tstAmbig <- aggregate(tmp$L, list(tmp$L), length)
#   colnames(tstAmbig) <- c("Matches reported per m/z", "Scans")
#   View(tstAmbig)
#   Ambig <- tmp[which(tmp$L > 1),]
#   wAmbig <- which(pep$Scan_ID_MZ %in% Ambig$Scan_ID_MZ)
#   wClear <- which(!pep$Scan_ID_MZ %in% Ambig$Scan_ID_MZ)
#   kol <- c("Scan_ID", "Modified sequence_verbose", "Theoretical m/z", "m/z", "Score", "PEP")
#   stopifnot(sum(!kol %in% colnames(pep)) == 0)
#   tmp2 <- as.data.table(pep[wAmbig, kol])
#   tmp2 <- tmp2[, list(Score = list(Score),
#                       PEP = list(PEP),
#                       `m/z` = list(`m/z`),
#                       `theor m/z` = list(`Theoretical m/z`),
#                       ModSeq = list(`Modified sequence_verbose`)), by = list(Scan_ID = Scan_ID)]
#   tmp2 <- as.data.frame(tmp2)
#   tmp2$BestModSeq <- sapply(1:nrow(tmp2), function(x) {
#     Sc <- unlist(tmp2$Score[x])
#     mdSq <- unlist(tmp2$ModSeq[x])
#     PEP <- unlist(tmp2$PEP[x])
#     delta <- unlist(tmp2$`m/z`[x])-unlist(tmp2$`theor m/z`[x])
#     w <- which(Sc == max(Sc))
#     if (length(w) > 1) {
#       w <- w[which(delta[w] == min(delta[w]))]
#       if (length(w) > 1) {
#         w <- w[which(PEP[w] == min(PEP[w]))]
#       }
#     }
#     mdSq[w]
#   })
#   #class(tmp2$BestModSeq)
#   wU <- which(sapply(tmp2$BestModSeq, length) == 1)
#   wNU <- which(sapply(tmp2$BestModSeq, length) > 1)
#   if (length(wU)) {
#     tmpA <- do.call(paste, c(pep[wAmbig, c("Scan_ID", "Modified sequence_verbose")], sep = "___"))
#     tmpB <- do.call(paste, c(tmp2[wU, c("Scan_ID", "BestModSeq")], sep = "___"))
#     w <- which(tmpA %in% tmpB)
#     rmv <- nrow(pep)-length(c(wClear, w))
#     cat(paste0("Removing ", rmv, " PSMs in favour of a higher quality PSM to the same spectrum but from another search.\n"))
#     pep <- pep[c(wClear, w),]
#   }
#   if (length(wNU)) {
#     warning(paste0(length(wNU), " ambiguous spectra could not be resolved and were removed."))
#     pep <- pep[which(!pep$Scan_ID %in% tmp2$Scan_ID[wNU]),]
#   }
# }
Exp <- unique(pep$Experiment)
mods <- lapply(FP_searches, function(FP_search) { FP2MQs[[FP_search]]$PTMs })
mods <- plyr::rbind.fill(mods)
if (!"Delta" %in% colnames(mods)) {
  # Let's allow for my inconsistencies in naming the delta mass columns:
  if ("Mass delta" %in% colnames(mods)) { mods$Delta <- mods$"Mass delta" } else {
    if ("Mass shift" %in% colnames(mods)) { mods$Delta <- mods$"Mass shift" } else {
      if ("Delta mass" %in% colnames(mods)) { mods$Delta <- mods$"Delta mass" }
    }
  }
}
mods <- mods[c(which(is.na(mods$`Old mark`)), which(!is.na(mods$`Old mark`))),]
mods$Delta <- as.numeric(mods$Delta)
mods <- mods[order(mods$Delta, decreasing = FALSE),]
kol <- colnames(mods)
kol <- kol[which(kol != "Old mark")]
modsTst <- do.call(paste, c(mods[, kol], sep = ""))
modsTst <- aggregate(1:nrow(mods), list(modsTst), min)
mods <- mods[modsTst$x,]
mods <- mods[c(which(is.na(mods$`Old mark`)), which(!is.na(mods$`Old mark`))),]
mods <- mods[order(mods$Delta, decreasing = FALSE),]
saveImgFun2(backupFl)
#loadFun2(backupFl)

# We now need to resolve cases where we got the PTM's name wrong
# - not that we made a mistake, but automation can only go so far.
# Here the issue is that FP takes as input just a delta mass, not a PTM,
# and to convert to MQ-like format we got a name from UniMod which may be incorrect.
# We also need to resolve ambiguous cases, where a single "mark" is assigned to 2 or more distinct PTMs.
# Let's do all of this at once.
modTbl <- read_xlsx(modTblFl)
colnames(modTbl)[1] <- "Priority"
modTbl <- modTbl[which(!is.na(modTbl$"molecular weight")),]
modTbl$`molecular weight` <- as.numeric(modTbl$`molecular weight`)
modTbl$Priority[which(is.na(modTbl$Priority))] <- ""
modTbl <- modTbl[which(modTbl$Priority != "Not analysed"),]
tmp <- listMelt(strsplit(modTbl$`modified residues`, ","), 1:nrow(modTbl), c("Pos", "row"))
tmp[, c("MW", "Name")] <- modTbl[tmp$row, c("molecular weight", "Modification")]
tst <- aggregate(tmp$Name, list(tmp$Pos, tmp$MW), function(x) { length(unique(x)) })
stopifnot(max(tst$x) == 1) # We assume that there is only one mod name for each combination of position and mass shift 
mods2 <- listMelt(mods$Site, 1:nrow(mods), c("Pos", "row"))
mods2[, c("MW", "Name", "Mark")] <- mods[mods2$row, c("Mass delta", "Full name", "Mark")]
kol <- c("Pos", "MW")
tmp$aggr <- do.call(paste, c(tmp[, kol], sep = "___"))
mods2$aggr <- do.call(paste, c(mods2[, kol], sep = "___"))
mods2$newName <- ""
mods2$newMark <- mods2$Mark
w <- which(mods2$aggr %in% tmp$aggr)
mods2$newName[w] <- tmp$Name[match(mods2$aggr[w], tmp$aggr)]
g1 <- grep("-", mods2$newName[w])
g2 <- grep("-", mods2$newName[w], invert = TRUE)
if (length(g1)) {
  mods2$newMark[w[g1]] <- sapply(strsplit(tolower(mods2$newName[w[g1]]), "-"), function(x) {
    paste0(substr(x[[1]], 1, 1), substr(x[[2]], 1, 1))
  })
}
if (length(g2)) {
  mods2$newMark[w[g2]] <- substr(tolower(mods2$newName[w[g2]]), 1, 2)
}
tst <- aggregate(mods2$MW, list(mods2$newMark), unique)
wMult <- which(sapply(tst$x, length) > 1)
if (length(wMult)) {
  mods2$newMark_tmp <- mods2$newMark
  for (i in wMult) {
    #i <- wMult[1]
    w <- which(mods2$newMark_tmp == tst$Group.1[i])
    m <- mods2[w,]
    m$Pos[which(sapply(m$Pos, length) == 0)] <- "X"
    r <- 1
    s <- 1:nrow(m); s <- s[which(s != r)]
    tst <- lapply(s, function(x) {
      paste0(tolower(m$Pos[[x]]), substr(m$newMark[[x]], 1, 1))
    })
    tst <- lapply(1:length(tst), function(x) {
      rs <- tst[[x]]
      rs[which(!rs %in% mods2$newMark)]
    })
    l <- length(tst)
    if (l > 1) {
      for (i in 2:l) {
        tst[[i]] <- tst[[i]][which(!tst[[i]] %in% unlist(tst[1:(i-1)]))]
      }
    }
    tst <- sapply(tst, function(x) {
      x <- unlist(x)
      if (length(x)) { x <- x[1] } else { x <- "That didnae work, did it?" }
      return(x)
    })
    tst2 <- ((tst %in% mods2$newMark)|(tst == "That didnae work, did it?"))
    w2 <- which(!tst2)
    m$newMark[s][w2] <- tst[w2]
    w1 <- which(tst2)
    if (length(w1)) {
      # not tested
      s1 <- s[w1]
      rs <- c()
      kount <- 1
      char <- c(0:9, letters)
      taken <- unique(c(mods2$newMark, m$newMark))
      for (j in s1) {
        tst <- paste0(tolower(m$Pos[s1]), char[kount])
        while (((tst) %in% taken)&&(kount < length(char))) {
          kount <- kount+1
          tst <- paste0(tolower(m$Pos[s1]), char[kount])
        }
        if (kount == length(char)) {
          stop("I am really out of options here, never thought this would go this far! Check the code just in case, this should never happen.")
        } else {
          rs <- c(rs, tst)
        }
      }
      m$newMark[[s1]] <- rs
    }
    mods2[w,] <- m
  }
}
tst <- aggregate(mods2$MW, list(mods2$newMark), unique)
stopifnot(is.numeric(tst$x))
mods <- mods2; rm(mods2)
w <- which(mods$newName == "")
mods$newName[w] <- mods$Name[w]
mods$newName <- gsub("\\)", ">", gsub("\\(", "<", mods$newName))
unq <- grep("\\(", unique(pep$`Modified sequence_verbose`), value = TRUE)
mods2 <- mods
mods2$newMark <- paste0("(", mods2$newMark, ")")
mods2$newName <- paste0("(", mods2$newName, ")")
clusterExport(parClust, "mods2", envir = environment())
unq2 <- as.data.frame(t(parSapply(parClust, unq, function(x) {
  #x <- unq[1]
  x <- y <- proteoCraft::annot_to_tabl(x)[[1]]
  w <- which(x$Annotations != "")
  m <- match(x$Annotations[w], mods2$Name)
  x$Annotations[w] <- mods2$newMark[m]
  y$Annotations[w] <- mods2$newName[m]
  x <- do.call(paste, c(x, sep = "", collapse = ""))
  y <- do.call(paste, c(y, sep = "", collapse = ""))
  return(c(x, y))
})))
tst1 <- gsub("\\)[A-Z]*\\(", "___", gsub("_[A-Z]+_|_[A-Z]*\\(|\\)[A-Z]*_", "", unq2[[1]]))
tst2 <- gsub("\\)[A-Z]*\\(", "___", gsub("_[A-Z]+_|_[A-Z]*\\(|\\)[A-Z]*_", "", unq2[[2]]))
tst1 <- unique(unlist(strsplit(tst1, "___")))
tst2 <- unique(unlist(strsplit(tst2, "___")))
stopifnot(sum(!tst1 %in% mods$newMark) == 0, sum(!tst2 %in% mods$newName) == 0)
w <- which(pep$`Modified sequence_verbose` %in% unq)
m <- match(pep$`Modified sequence_verbose`[w], unq)
pep$`Modified sequence`[w] <- unq2[m, 1]
pep$`Modified sequence_verbose`[w] <- unq2[m, 2]
# Done!!!
# Pfeewwwww...
map <- setNames(lapply(FP_searches, function(FP_search) {
  x <- FP2MQs[[FP_search]]$FracMap
  x$Path <- NULL
  x
}), FP_searches)
tst <- sapply(map, function(x) { do.call(paste, c(x, collapse = "")) })
tst <- aggregate(1:length(tst), list(tst), list)
tst$x <- sapply(tst$x, min)
map <- plyr::rbind.fill(map[tst$x])
#
dbFl <- unique(unlist(sapply(FP_searches, function(FP_search) {
  WorkFlow <- FP2MQs[[FP_search]]$WorkFlow
  dbFl <- gsub("database.db-path=", "", grep("^database.db-path=", WorkFlow, value = TRUE))
  dbFl <- gsub("\\\\\\\\", "/", dbFl)
  gsub("\\\\", "", dbFl)
  #file.exists(dbFl)
})))
# Hard fix
dbFl <- "D:/Fasta_databases/Mus_musculus/Mus_musculus_(C57BL-6J)_UP_20230201_Iso_noDupl_cont.fasta"

RenameExp <- c(TRUE, FALSE)[match(dlg_message("Should we rename individual samples?", "yesno")$res,
                                  c("yes", "no"))]
if (RenameExp) {
  smplsMapFl <- paste0(dstDir, "/Samples_map.csv")
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

Remove0Int <- FALSE
#Remove0Int <- c(FALSE, TRUE)[match(dlg_message("Should we keep peptides with 0 intensity values?", "yesno")$res, c("yes", "no"))]
if (Remove0Int) { pep <- pep[which(pep$Intensity > 0),] }

prsDBfl <- paste0(dstDir, "/Parsed DB.csv")
if (!file.exists(prsDBfl)) {
  db <- Format.DB(dbFl, cl = parClust)
  write.csv(db, prsDBfl, row.names = FALSE)
} else {
  db <- data.table::fread(prsDBfl, integer64 = "numeric", check.names = FALSE, data.table = FALSE) 
}
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
# Default organism may not be correct:
warning("This script assumes Mus musculus samples: if this is not the case, edit the file name in the next line!")
writeFasta(histDB, paste0(dstDir, "/M_musculus_all_histones.fasta"))
histDB$ID_Nm <- apply(histDB[, c("Protein ID", "Common Name")], 1, paste, collapse = " - ")
histDB$ID_Nm <- gsub("/", "_", histDB$ID_Nm)
clusterExport(parClust, "histDB", envir = environment())

#
tmp <- data.frame(ModSeq =  pep$"Modified sequence_verbose")
tmp$tst <- gsub("^_[A-Z]+_$|^_[A-Z]*\\(|\\)[A-Z]*_$", "",
                gsub("\\)[A-Z]+\\(", "___", tmp$ModSeq))
tst <- unlist(strsplit(tmp$tst, "___"))
tst <- aggregate(tst, list(tst), length)
nrow(tst) # Number of identified PTMs
length(unique(mods$newName)) # Number of searched PTMs
stopifnot(sum(!tst$Group.1 %in% mods$newName) == 0) # Checking that the names match

# Check peptide-to-protein assignment, and isolate histones peptides
### NB: This doesn't take into consideration N-terminal methionines!!!
msg <- "Update FragPipe's original protein-to-peptides assignments?"
Update_Prot_matches <- c(TRUE, FALSE)[match(dlg_message(msg, "yesno")$res, c("yes", "no"))]
if (Update_Prot_matches) {
  fl <- paste0(dstDir, "/evmatch.RData")
  if (file.exists(fl)) { loadFun2(fl) } else {
    uSeq <- unique(pep$Sequence)
    evmatch <- ProtMatch2(uSeq, db, cl = parClust)
    saveFun2(evmatch, fl)
  }
  pep$Proteins_FP <- pep$Proteins
  pep$Proteins <- evmatch$Proteins[match(pep$Sequence, evmatch$Sequence)]
  #View(pep[, c("Proteins_FP", "Proteins", "Histone(s)", "Modified sequence_verbose")])
  #sum(pep$Proteins_FP != pep$Proteins)
  #sum(pep$Proteins_FP == pep$Proteins)
}

tmp <- strsplit(pep$Proteins, ";")
tmp <- proteoCraft::listMelt(tmp, 1:nrow(pep), c("ID", "row"))
tmp <- tmp[which(tmp$ID %in% histDB$`Protein ID`),]
tmp <- data.table::as.data.table(tmp)
tmp <- tmp[, list(IDs = list(ID)), by = list(row = row)]
tmp <- as.data.frame(tmp)
tmp$IDs <- parSapply(parClust, tmp$IDs, paste, collapse = ";")
pep$"Histone(s)" <- ""
pep$"Histone(s)"[tmp$row] <- tmp$IDs
#View(pep[tmp$row, c("Proteins", "Histone(s)")])
sum(pep$Proteins != pep$`Histone(s)`)
sum(pep$Proteins == pep$`Histone(s)`)

tmp <- pep[, c("Histone(s)", "Proteins")]
tmp$"Proteins" <- strsplit(tmp$"Proteins", ";")
tmp$"Histone(s)" <- strsplit(tmp$"Histone(s)", ";")
tst <- apply(tmp[, c("Histone(s)", "Proteins")], 1, function(x) {
  sum(!x[[1]] %in% x[[2]])
})
w <- which(tst > 0)
#View(pep[w, c("Modified sequence_verbose", "Histone(s)", "Proteins")])
#View(pep[grep("phospho", pep$`Modified sequence_verbose`, ignore.case = TRUE), c("Modified sequence_verbose", "Histone(s)", "Proteins")])

saveImgFun2(backupFl)
#loadFun2(backupFl)

# Summarize over charge
if (!"Raw file path" %in% colnames(pep)) {
  stopifnot("Raw file" %in% colnames(pep))
  pep$"Raw file path" <- pep$`Raw file`
}
kol1 <- c("Modified sequence", "Raw file path", "Experiment")
kol2 <- c("Intensity", "Charge", "m/z", "Score")
kol3 <- c("Modified sequence_verbose", "Modifications", "Proteins", "Histone(s)", "Mass",
          "Missed cleavages", "Potential contaminant", "Reverse")
pepSum <- as.data.table(pep[, c(kol1, kol2)])
pepSum <- pepSum[,
                 list(Intensity = sum(Intensity),
                      Charge = paste(Charge, collapse = ";"),
                      `m/z` = paste(`m/z`, collapse = ";"),
                      Score = max(Score, na.rm = TRUE)),
                 by = list(`Modified sequence` = `Modified sequence`,
                           `Raw file path` = `Raw file path`,
                           Experiment = Experiment)]
pepSum <- as.data.frame(pepSum)
pepSum[, kol3] <- pep[match(pepSum[["Modified sequence"]], pep[["Modified sequence"]]), kol3]
pepSum <- pepSum[, c("Modified sequence", "Modified sequence_verbose", "Modifications", "Intensity",
                     "Raw file path", "Experiment", "Proteins", "Histone(s)", "Mass", "Charge", "m/z", "Score",
                     "Missed cleavages", "Potential contaminant", "Reverse")]

tst <- aggregate(pepSum$"Histone(s)", list(pepSum$"Histone(s)"), length)
colnames(tst) <- c("Histone matches", "Peptides")
tst <- tst[which(tst$`Histone matches` != ""),]
tst <- tst[order(tst$Peptides, decreasing = TRUE),]
View(tst)
w <- which(parSapply(parClust, strsplit(pepSum$`Histone(s)`, ";"), function(x) {
  x <- unlist(x)
  x <- x[which(x != "")]
  (length(x) > 1)||(!is.na(x))
}))
histPep <- pepSum[w,]
tst <- gsub("^_[A-Z]+_$|^_[A-Z]*\\(|\\)[A-Z]*_$", "", gsub("\\)[A-Z]+\\(", "___",
                                                           histPep$`Modified sequence_verbose`))
tst <- unlist(strsplit(tst, "___"))
tst <- aggregate(tst, list(tst), length)
nrow(tst) # Number of identified PTMs
stopifnot(sum(!unique(unlist(strsplit(histPep$`Histone(s)`, ";"))) %in% histDB$`Protein ID`) == 0)
histPep$Sequence <- pep$Sequence[match(histPep$`Modified sequence`, pep$`Modified sequence`)]

histDir <- paste0(dstDir, "/Histones_coverage_maps")
if (!dir.exists(histDir)) { dir.create(histDir) } else {
  ls <- list.files(histDir, full.names = TRUE)
  if (length(ls)) {
    msg <- "The coverage maps folder is not empty, clear it?"
    clearDir <- c(TRUE, FALSE)[match(dlg_message(msg, "yesno")$res, c("yes", "no"))]
    if (clearDir) {
      for (fl in ls) { unlink(fl) }
    }
  }
}
Exp <- unique(pep$Experiment)

# Check PTMs
wMod <- grep("\\(", histPep$`Modified sequence_verbose`)
tmpMds <- gsub("_[A-Z]+_|_[A-Z]*\\(|\\)[A-Z]*_", "",
               gsub("\\)[A-Z]*\\(", "___", histPep$`Modified sequence_verbose`[wMod]))
tst <- unlist(strsplit(tmpMds, "___"))
tst <- aggregate(tst, list(tst), length)
tst <- tst[order(tst$x, decreasing = TRUE),]
colnames(tst) <- c("PTM", "Nb. of peptides")
View(tst)
data.table::fwrite(tst, paste0(histDir, "/PTMs summary.csv"), sep = ",", row.names = FALSE, na = "NA")

# Count PTM sites
kol <- c("Sequence", "Modified sequence_verbose", "Proteins")
tmp <- histPep[wMod, kol]
tst <- do.call(paste, c(tmp, sep = "___"))
tst <- unique(tst)
tmp <- Isapply(strsplit(tst, split = "___"), unlist)
colnames(tmp) <- kol
tmp2 <- unlist(strsplit(tmp$Proteins, ";"))
tmp2 <- aggregate(tmp2, list(tmp2), length)
tmp2 <- tmp2[order(tmp2$x, decreasing = TRUE),]
clusterExport(parClust, "tmp2", envir = environment())
tmp$Proteins_by_pep <- parSapply(parClust, strsplit(tmp$Proteins, ";"), function(x) { #x <- strsplit(tmp$Proteins[5], ";")
  x <- unlist(x)
  if (length(x) > 1) {
    o <- tmp2$x[match(x, tmp2$Group.1)]
    x <- x[order(o, decreasing = TRUE)]
    x <- paste(sort(unlist(x)), collapse = ";")
  }
  return(x)
})
tmp$First_Protein <- gsub(";.*", "", tmp$Proteins_by_pep)
tmpSq <- paste0("_", db$Sequence[match(tmp$First_Protein, db$`Protein ID`)])
tmp$Match <- sapply(1:nrow(tmp), function(x) {
  nchar(unlist(strsplit(tmpSq[x], tmp$Sequence[x]))[1])
})
tmp$Tbl <- lapply(tmp$`Modified sequence_verbose`, proteoCraft::annot_to_tabl)
tmp$Sites <- lapply(1:nrow(tmp), function(x) { #x <- 1
  mtch <- tmp$Match[x]
  pr <- tmp$First_Protein[x]
  tbl <- tmp$Tbl[[x]]
  tbl <- tbl[[1]]
  nr <- nrow(tbl)-2
  tbl$Dist <- c(1, 1:nr, nr)
  tbl$Pos <- tbl$Dist + mtch - 1
  tbl <- tbl[which(tbl$Annotations != ""), , drop = FALSE]
  tbl$Site <- apply(tbl[, c("Sequence", "Pos", "Annotations")], 1, function(y) {
    paste0(pr, " ", y[[1]], y[[2]], "(", y[[3]], ")")
  })
  return(tbl$Site)
})
Sites <- unique(unlist(tmp$Sites))
Sites <- data.frame(Site = Sites)
Sites$Type <- gsub(".*\\(|\\)", "", Sites$Site)
nSites <- aggregate(Sites$Type, list(Sites$Type), length)
nSites <- nSites[order(nSites$x, decreasing = TRUE),]
colnames(nSites) <- c("PTM", "Nb. of sites")
View(nSites)
data.table::fwrite(nSites, paste0(histDir, "/PTM sites summary.csv"), sep = ",", row.names = FALSE, na = "NA")
sum(nSites$`Nb. of sites`[which(!nSites$PTM %in% c("Carbamidomethyl", "Oxidation"))])

saveImgFun2(backupFl)
#loadFun2(backupFl)

kol <- c("Experiment", "Histone(s)", "Intensity", "Modified sequence", "Modified sequence_verbose")
kol %in% colnames(histPep)
tmpPep <- histPep[, kol]
clusterExport(parClust, list("tmpPep", "histDir", "Hist", "Coverage", "histDB"), envir = environment())
#clusterExport(parClust, "Coverage", envir = environment())
for (e in Exp) { #e <- Exp[1]
  w <- which(tmpPep$Experiment == e)
  #w <- which(tmpPep$Experiment %in% Exp)
  if (length(w)) {
    pp <- tmpPep[w,]
    pp <- aggregate(pp$Intensity,
                    list(pp$`Modified sequence`, pp$`Modified sequence_verbose`),
                    sum, na.rm = TRUE)
    colnames(pp) <- c("Modified sequence", "Modified sequence_verbose", "Intensity")
    pp$"Histone(s)" <- tmpPep$"Histone(s)"[match(pp$`Modified sequence`, tmpPep$`Modified sequence`)]
    tmpPP <- listMelt(strsplit(pp$"Histone(s)", ";"), 1:nrow(pp), c("Histone", "row"))
    tmpPP <- aggregate(tmpPP$row, list(tmpPP$"Histone"), c)
    tmpPP <- setNames(lapply(tmpPP$x, function(x) {
      pp[x, c("Modified sequence_verbose", "Intensity")]
    }), tmpPP$Group.1)
    e <- gsub("[\\\\/:\\*\\?\"<>\\|]", "_", e)
    clusterExport(parClust, list("e", "tmpPP"), envir = environment())
    # tst <- sapply(names(tmpPP), function(x) {
    #   g <- grep("nyl", tmpPP[[x]]$`Modified sequence_verbose`, value = TRUE)
    #   if (length(g)) { stop(cat(x, "\n", paste(g, collapse = "\n"), "\n")) }
    # })
    # View(tmpPP[["A0A2R8Y619"]])
    w <- which(histDB$`Protein ID`[match(Hist, histDB$Header)] %in% names(tmpPP))
    if (length(w)) {
      a <- parSapply(parClust, Hist[w], function(hist) {
        #hist <- Hist[w][1]
        #hist <- grep("G3UX40", Hist[w], value = TRUE)
        #hist <- grep("G3UWL7", Hist[w], value = TRUE)
        #hist <- grep("A0A0N4SV66", Hist[w], value = TRUE)
        #hist <- grep("A0A8I4SYN6", Hist[w], value = TRUE)
        #for (hist in Hist[w]) {
        #grep("Q07133", Hist[w])
        m <- match(hist, histDB$Header)
        if (length(m)) {
          nm <- histDB$ID_Nm[m]
          sq <- setNames(histDB$Sequence[m], histDB$ID_Nm[m])
          tmpSq <- tmpPP[[histDB$`Protein ID`[m]]]
          Coverage(sq, tmpSq$`Modified sequence_verbose`, "Heat", FALSE, FALSE,
                   title = paste0(nm, " coverage - ", e),
                   intensities = tmpSq$Intensity,
                   save.path = paste0(histDir, "/", nm, " - ", e),
                   save = c("pdf", "jpeg"))
          fl <- paste0(histDir, "/", nm, " - ", e)
          if (file.exists(fl)) { unlink(fl) }
        } else { stop() }
      }) 
    }
  }
}
#openwd(histDir)
data.table::fwrite(histPep, paste0(histDir, "/Hist_peptidoforms.tsv"),
                   sep = "\t", row.names = FALSE, na = "NA")
