# Source to control how we use the splitByOrg argument and the lower level PG_Assemble source
#
pgSrc <- paste0(libPath, "/extdata/Sources/PG_Assemble.R")
#rstudioapi::documentOpen(pgSrc)
if (!exists("splitByOrg")) { splitByOrg <- FALSE }
if (splitByOrg) {
  uOrgs <- unique(db[[dbOrgKol]])
  uOrgsNorm <- gsub("_+", "_", gsub("[<>,:\"/|?.\\\\]", "_", uOrgs))
  stopifnot(length(uOrgsNorm) == length(uOrgs)) # Unhandled exception
  names(uOrgsNorm) <- uOrgs
  PGs_list <- list()
  myTmpPepFilt <- listMelt(strsplit(pep$Proteins, ";"), ColNames = c("ID", "Row"))
  myTmpPepFilt$Org <- db[match(myTmpPepFilt$ID, db$`Protein ID`), dbOrgKol]
  for (currOrg in uOrgs) {
    DB <- db[db[[dbOrgKol]] == currOrg,]
    Pep <- pep[1L:nrow(pep) %in% myTmpPepFilt$Row[myTmpPepFilt$Org == currOrg],]
    Ev <- ev[ev$`Modified sequence` %in% Pep$`Modified sequence`,]
    source(pgSrc)
    PGs_list[[currOrg]] <- PG_assembly
  }
  # Stitch database back together
  db <- plyr::rbind.fill(lapply(PGs_list, \(x) { x$Database }))
  # Protein groups
  PGs <- lapply(PGs_list, \(x) { x$Protein.groups })
  PG <- plyr::rbind.fill(PGs)
  # Peptides and PSMs
  # NB: For now we are allowing one razor protein group per organism!
  #     Mainly because we are lazy and it doesn't sound too important.
  #     True.
  #     We will probably change this in the future.
  kol <- c("Protein group IDs", "Razor protein group ID",
           "Leading proteins", "Leading razor proteins",
           "Protein names", "Protein names (all)", "Common protein names",
           "Gene names", "Gene names (all)")
  tmpPeps <- lapply(PGs_list, \(x) { x$Peptides[, c("Modified sequence", kol)] })
  tmpPeps <- plyr::rbind.fill(tmpPeps)
  tmpPeps <- as.data.table(tmpPeps)
  for (k in kol) {
    if (is.character(tmpPeps[[k]])) { tmpPeps[[k]] <- strsplit(tmpPeps[[k]], ";") }
  }
  tmpPeps <- tmpPeps[,
                     lapply(.SD, \(x) { paste(unique(unlist(x)), collapse = ";") }),
                     keyby = `Modified sequence`]
  tmpPeps <- as.data.frame(tmpPeps)
  pep[, kol] <- tmpPeps[match(pep$`Modified sequence`, tmpPeps$`Modified sequence`), kol]
  ev[, kol] <- tmpPeps[match(ev$`Modified sequence`, tmpPeps$`Modified sequence`), kol]
} else {
  #rstudioapi::documentOpen(pgSrc)
  source(pgSrc)
  #loadFun("PG_assembly.RDS")
  PG <- PG_assembly$Protein.groups
  pep <- PG_assembly$Peptides
  db <- PG_assembly$Database
  if ("Evidences" %in% names(PG_assembly)) { ev <- PG_assembly$Evidences }
}
DatAnalysisTxt$PG_assemble <- paste0("Observed peptides were grouped and assembled into ",
                                     c("", "species-specific ")[splitByOrg + 1L],
                                     "Protein Groups using a heuristic maximum parsimony algorithm.")
PG %<o% PG
