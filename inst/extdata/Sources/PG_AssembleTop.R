#### Assemble protein groups
#
# Top-level source to control how we use the splitByOrg argument and the lower level PG_Assemble source.
cat("Starting protein groups assembly...\n")
if (!exists("N_Pep")) { # Defined in the parameter sources, but we need a default
  N_Pep %<o% 1L
}
pgSrc <- paste0(libPath, "/extdata/Sources/PG_Assemble.R")
#rstudioapi::documentOpen(pgSrc)
if (!exists("splitByOrg")) { splitByOrg <- FALSE }
PG_assembly %<o% list()
if (splitByOrg) {
  uOrgs <- unique(db[[dbOrgKol]])
  uOrgsNorm <- gsub("_+", "_", gsub("[<>,:\"/|?.\\\\]", "_", uOrgs))
  stopifnot(length(uOrgsNorm) == length(uOrgs)) # Unhandled exception
  names(uOrgsNorm) <- uOrgs
  myTmpPepFilt <- listMelt(strsplit(pep$Proteins, ";"), ColNames = c("ID", "Row"))
  myTmpPepFilt$Org <- db[match(myTmpPepFilt$ID, db$`Protein ID`), dbOrgKol]
  for (currOrg in uOrgs) {
    DB <- db[db[[dbOrgKol]] == currOrg,]
    Pep <- pep[1L:nrow(pep) %in% myTmpPepFilt$Row[myTmpPepFilt$Org == currOrg],]
    Ev <- ev[ev$`Modified sequence` %in% Pep$`Modified sequence`,]
    source(pgSrc)
    PG_assembly[[currOrg]] <- PG_assembl
  }
  # Stitch database back together
  db <- plyr::rbind.fill(lapply(PG_assembly, \(x) { x$Database }))
  # Protein groups
  PGs <- lapply(PG_assembly, \(x) { x$Protein.groups })
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
  tmpPeps <- lapply(PG_assembly, \(x) { x$Peptides[, c("Modified sequence", kol)] })
  tmpPeps <- plyr::rbind.fill(tmpPeps)
  requireNamespace("data.table")
  data.table::setDT(tmpPeps)
  for (k in kol) {
    if (is.character(tmpPeps[[k]])) { tmpPeps[[k]] <- strsplit(tmpPeps[[k]], ";") }
  }
  tmpPeps <- tmpPeps[, lapply(.SD, \(x) { paste(unique(unlist(x)), collapse = ";") }),
                     keyby = `Modified sequence`]
  tmpPeps <- as.data.frame(tmpPeps)
  pep[, kol] <- tmpPeps[match(pep$`Modified sequence`, tmpPeps$`Modified sequence`), kol]
  ev[, kol] <- tmpPeps[match(ev$`Modified sequence`, tmpPeps$`Modified sequence`), kol]
} else {
  #rstudioapi::documentOpen(pgSrc)
  source(pgSrc)
  PG_assembly[[1L]] <- PG_assembl
  PG <- PG_assembl$Protein.groups
  pep <- PG_assembl$Peptides
  db <- PG_assembl$Database
  if ("Evidences" %in% names(PG_assembl)) { ev <- PG_assembl$Evidences }
}
DatAnalysisTxt$PG_assemble <- paste0("Observed peptides were grouped and assembled into ",
                                     c("", "species-specific ")[splitByOrg + 1L],
                                     "Protein Groups using a heuristic maximum parsimony algorithm.")
PG %<o% PG
saveFun(PG_assembly, backupFls["PG_assembly"])
