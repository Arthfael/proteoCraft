#' PG_assemble
#'
#' @description 
#' This former function is now a source. And faster for it!
#' It takes a peptides and a data base data.frames as input and returns a list with:
#'  - a data.frame of parsimoniously-inferred protein groups,
#'  - the peptides file with a few extra columns,
#'  - the database file with a few extra columns; this is also cleaned up to remove any NAs which can make for poor protein names/annotations.
#'  
#'  Note:
#'  Pep expects a wide (peptides), rather than long (PSMs/evidences), table! Ev is only used for mapping PSMs, not to assemble protein groups.
#'  We haven't explored the potential implications for the function's logic of supplying to Pep a long table.
#'  But when this was accidentally tested on a large-ish dataset, a powerful PC quickly froze and did not manage to resolve the task overnight.
#'  Just do not do it!
#' 
#' @param Pep The peptides sequences data frame. Should contain "Modified.sequence", "Sequence" and "PEP" columns, as well as IDs, protein and evidence IDs columns (see next two arguments).
#' @param Peptide.IDs Name of the column with peptide IDs. Default = "id"
#' @param Proteins.col Name of the column with matched protein accessions. Default = "Proteins"
#' @param DB A search database parsed as a table by the Format.DB function. Will also act as a filter: only protein accessions present in its "Protein ID" column will be part of the assembly!!!
#' @param Evidence.IDs Name of the column with evidence IDs. Default = "Evidence IDs"
#' @param Ev The evidences data frame. Used for a few annotations columns. Not required.
#' @param N.clust A limit on the number of vCPUs to use. If left as NULL (default), uses the number of available clusters - 1, to a minimum of 1.
#' @param N.reserved Default = 1. Number of reserved vCPUs the function is not to use. Note that for obvious reasons it will always use at least one.
#' @param Custom_PGs Data frame of protein groups known a priori (default = NA). One column required, called "Leading protein IDs" with semicolon-separated IDs; optional 2nd column: a numeric/integer called "Priority level". Any peptide mapping to these groups will be assigned to them to the exclusion of other protein groups. Use the priority column to resolve ambiguities (1 = highest priority). NB: Use with caution, only for known, priority proteins (e.g. artificial constructs).
#' @param ContCol Name of column in the database specifying whether a protein is a contaminant ("+", otherwise ""). Default = "Potential contaminant"
#' @param Npep Minimum number of razor or unique peptides to tag a protein group as a true discovery (others will still be reported). Default = N_Pep if it exists, otherwise 2.
#' @param cl Already have a cluster handy? Why spend time making a new one, which on top of that may invalidate the old one. Just pass it along!
#' @param dbOrgCol Default = "Organism_Full". Name of the column in DB corresponding to the organism annotation.
#' 
#' @examples
#' temp <- PG_assemble(Pep = pep, DB = db)
#' 
#' @import data.table
#' @export

source(parSrc)
# Below: expected values of old arguments.
# These can be superseded by providing a value before calling the source.
if (!exists("Peptide.IDs")) { Peptide.IDs <- "id" }
if (!exists("Proteins.col")) { Proteins.col <- "Proteins" }
if (!exists("Evidence.IDs")) { Evidence.IDs <- "Evidence IDs" }
if (exists("custPGs")) { Custom_PGs <- custPGs }
if (!exists("Custom_PGs")) { Custom_PGs <- NA }
if (!exists("ContCol")) { ContCol <- "Potential contaminant" }
Npep <- if (exists("N_Pep")) { N_Pep } else { 2L }
if (!exists("Npep")) { Npep <- 2L }
if (!exists("DB")) { DB <- db }
if (!exists("Pep")) { Pep <- pep }
if (!exists("Ev")) { Ev <- ev }
if ((!exists("dbOrgCol")) || (!is.character(dbOrgCol)) || (length(dbOrgCol) != 1L) || (!nchar(dbOrgCol)) || (dbOrgCol %in% colnames(DB))) {
  dbOrgCol <- "Organism_Full"
}
if (!dbOrgCol %in% colnames(DB)) { dbOrgCol <- "Organism" }
if (!dbOrgCol %in% colnames(DB)) { dbOrgCol <- NA }
if ((!exists("splitByOrg")) || (!is.logical(splitByOrg)) || (length(splitByOrg) != 1L) || is.na(splitByOrg)) {
  splitByOrg <- FALSE
}
if (splitByOrg && (!exists("n_of_PGs"))) { n_of_PGs <- 0L }

if ((!splitByOrg) || (currOrg == uOrgs[1L])) {
  cat("Protein groups inference step:\nAssembling peptides into the minimum number of protein groups required to explain the data...\n")
}
if (splitByOrg) {
  cat(currOrg, "...\n")
}

#DefArg(PG_assemble); Npep = N_Pep; N.clust <- 55L; Custom_PGs <- custPGs; Ev <- ev; DB <- db
#Pep <- pep
#Pep <- pep[1L:1000L,]
#
tm1 <- Sys.time()
# Note:
# This is not a perfect alternative to missing but will work in most cases, unless x matches a function imported by a package 
misFun <- \(x) { return(!exists(deparse(substitute(x)))) }
#
requireNamespace("data.table")
invisible(parallel::clusterCall(parClust, \() {
  requireNamespace("data.table")
  return()
}))
#
if (!"MW [kDa]" %in% colnames(DB)) {
  DB$"MW [kDa]" <- round(Peptides::mw(DB$Sequence)/1000, 3L)
}
# We do not want any NAs in our protein and gene names!
tmpDB <- as.matrix(DB)
w <- which(is.na(tmpDB), arr.ind = TRUE)
if (nrow(w)) { DB[w] <- "" }
# Peptides (a new data frame in casTRUE# Peptides (a new data frame in case I screw up, but this is essentially the same as the original one)
seq <- data.frame(Proteins = Pep[[Proteins.col]], "Modified sequence" = Pep$"Modified sequence",
                  Sequence = Pep$Sequence, id = Pep[[Peptide.IDs]], PEP = Pep$PEP, check.names = FALSE)
if (!misFun(Ev)) { seq$"Evidence IDs" <- Pep[[Evidence.IDs]] }
tmpids <- DB$"Protein ID"
seq$.Proteins <- strsplit(seq$Proteins, ";")
# Rarely we can have peptides without protein assignments in the input.
# Obviously those are useless for the purpose of assembling protein groups!
# There may also be protein IDs not from the data base provided, which should be removed.
rws <- 1L:nrow(seq)
tmp1 <- seq$.Proteins
tmp1 <- listMelt(tmp1, asDT = TRUE)
tmp1 <- tmp1[tmp1$value %in% DB$"Protein ID",]
tmp1 <- tmp1[, .(value = list(value)),
             by = .(L1)]
tmp1 <- tmp1[order(tmp1$L1),]
seq$.Proteins2 <- lapply(rws, c)
w <- which(rws %in% tmp1$L1)
seq$.Proteins[w] <- tmp1$value
#
nPrts <- lengths(seq$.Proteins)
wY <- which(nPrts > 0L)
wN <- which(nPrts == 0L)
if (length(wN)) { warning("Not all of the protein IDs in the peptides file are in the data base!") }
seq <- seq[wY,]
seq$Proteins <- vapply(seq$Proteins, paste, "", collapse = ";")
# (Creating temporary IDs is convenient as I don't even need to use "match" for this.)
proffset <- pgoffset <- 0L
CustPG <- !is.na(list(Custom = Custom_PGs)) # Allows testing despite differences in dimensions
if (CustPG) {
  #Custom_PGs <- Custom_PGs[, "Leading protein IDs", drop = FALSE]
  stopifnot("Leading protein IDs" %in% colnames(Custom_PGs))
  if (!"Priority level" %in% colnames(Custom_PGs)) { Custom_PGs$"Priority level" <- 1L }
  Custom_PGs <- Custom_PGs[, c("Leading protein IDs", "Priority level")]
  Custom_PGs$"Priority level" <- as.integer(Custom_PGs$"Priority level")
  stopifnot(sum(is.na(Custom_PGs$"Priority level")) == 0L)
  # Create proteins table for custom PGs
  # Expand custom PGs columns
  tmp <- unique(unlist(strsplit(Custom_PGs$`Leading protein IDs`, ";")))
  g <- grsep(tmp, x = seq$Proteins)
  if (length(g)) {
    seq2 <- seq[g,]
    seq <- seq[!(1L:nrow(seq)) %in% g,]
    Custom_PGs$.pep.ids <- list(NA)
    priorities <- sort(unique(Custom_PGs$"Priority level"), decreasing = TRUE) # Sorted by decreasing value (increasing priority)
    for (i in max(priorities):1L) { #i <- max(priorities)
      wcp <- which(Custom_PGs$"Priority level" == i)
      Custom_PGs$.pep.ids[wcp] <- lapply(strsplit(Custom_PGs$`Leading protein IDs`[wcp], ";"), \(x) {
        seq2$id[grsep(unlist(x), x = seq2$Proteins)]
      })
      wlp <- which(Custom_PGs$"Priority level" > i)
      if (length(wlp)) {
        Custom_PGs$.pep.ids[wlp] <- lapply(Custom_PGs$.pep.ids[wlp], \(x) {
          base::setdiff(x, unlist(Custom_PGs$.pep.ids[wcp]))
        })
      }
    }
    temp <- listMelt(strsplit(seq2$Proteins, ";") , seq2$id, c("Protein", "pep.id"), TRUE)
    w <- which(!as.character(temp$Protein) %in% c("", " ", "NA"))
    if (length(w) != nrow(temp)) {
      warning("Some protein IDs seem dodgy, we will remove them but you should check what happened.")
      temp <- temp[w,]
    }
    prot2 <- temp[, .(pep.ids = paste(pep.id, collapse = ";"),
                      .pep.ids = list(pep.id),
                      `Peptides count` = length(pep.id)),
                  keyby = .(Protein = Protein)]
    prot2 <- as.data.frame(prot2)
    row.names(prot2) <- 1L:nrow(prot2)
    # Sort by N of peptides
    prot2 <- prot2[order(prot2$"Peptides count", decreasing = TRUE),]
    prot2$Prot.id <- as.character(1L:nrow(prot2))
    proffset <- nrow(prot2)
    Custom_PGs$.Leading.protein.IDs <- strsplit(Custom_PGs$`Leading protein IDs`, ";")
    Custom_PGs$"Peptide IDs" <- vapply(Custom_PGs$.pep.ids, paste, "", collapse = ";")  
    Custom_PGs$.lead.protein.ids <- lapply(Custom_PGs$.Leading.protein.IDs, \(x) {
      prot2$Prot.id[match(x, prot2$Protein)]
    })
    Custom_PGs$lead.protein.ids <- vapply(Custom_PGs$.lead.protein.ids, paste, "", collapse = ";")
    Custom_PGs$"Peptides count" <- lengths(Custom_PGs$.pep.ids)
    Custom_PGs$temp.pg.id <- as.character(1L:nrow(Custom_PGs))
    pgoffset <- nrow(Custom_PGs)
    seq2$.Proteins <- strsplit(seq2$Proteins, ";")
    seq2$.Prot.ids <- lapply(seq2$.Proteins, \(x) { prot2$Prot.id[match(unlist(x), prot2$Protein)] })
    temp <- listMelt(Custom_PGs$.pep.ids, Custom_PGs$temp.pg.id, c("pep.id", "temp.pg.id"), TRUE)
    temp <- temp[, .(temp.pg.ids = paste(temp.pg.id, collapse = ";")),
                 by = .(pep.id = pep.id)]
    seq2$temp.pg.ids <- temp$temp.pg.ids[match(seq2$id, temp$pep.id)]
    seq2$.temp.pg.ids <- strsplit(seq2$temp.pg.ids, ";")
    w1 <- which(Custom_PGs$`Peptides count` > 0L)
    w0 <- which(Custom_PGs$`Peptides count` == 0L)
    l <- length(w0)
    if (l) {
      warning(paste0("The following custom protein group ", c("was", "were")[(l>1L)+1L], " not found:",
                     paste0("\n - ", Custom_PGs$`Leading protein IDs`[w0]), "\n"))
    }
    Custom_PGs <- Custom_PGs[w1,]
  } else {
    warning("No peptide matches to provided custom protein group(s) found, argument \"Custom_PGs\" will be ignored.")
    CustPG <- FALSE
  }
}
# Proteins to peptides table:
# Unique proteins...
temp <- listMelt(strsplit(seq$Proteins, ";"), seq$id, c("Protein", "pep.id"), TRUE)
w <- which(!as.character(temp$Protein) %in% c("", " ", "NA"))
if (length(w) != nrow(temp)) {
  warning("Some protein IDs seem dodgy, we will remove them but you should check what happened.")
  temp <- temp[w,]
}
prot <- temp[, .(pep.ids = paste(pep.id, collapse = ";"),
                 .pep.ids = list(pep.id),
                 `Peptides count` = length(pep.id)),
             keyby = .(Protein = Protein)]
prot <- as.data.frame(prot)
row.names(prot) <- 1L:nrow(prot)
# Sort by N of peptides
prot <- prot[order(prot$"Peptides count", decreasing = TRUE),]
# Temporary protein IDs
prot$Prot.id <- as.character(1L:nrow(prot) + proffset)
if (CustPG) {
  prot <- rbind(prot2, prot)
  stopifnot(length(unique(prot$Prot.id)) == nrow(prot))
}
rownames(prot) <- stringi::stri_join("Pr_", prot$Prot.id)
# assign temp prot IDs to peptides
tmp1 <- seq$.Proteins
tmp1 <- listMelt(tmp1, ColNames = c("Protein", "row"), asDT = TRUE)
tmp1$ProtID <- prot$Prot.id[match(tmp1$Protein, prot$Protein)]
tmp1 <- tmp1[, .(.Prot.ids = list(ProtID)),
             by = .(row = row)]
tmp1 <- as.data.frame(tmp1)
tmp1 <- tmp1[order(tmp1$row),]
rws <- 1L:nrow(seq)
seq$.Prot.ids <- lapply(rws, c) 
w <- which(rws %in% tmp1$row)
seq$.Prot.ids <- tmp1$.Prot.ids
#
# For each protein P, find all other proteins Q it shares peptides with...
#a1 <- Sys.time()
tmp1 <- listMelt(prot$.pep.ids)
tmp1$prIDs <- seq$.Prot.ids[match(tmp1$value, seq$id)]
tmp1 <- listMelt(tmp1$prIDs, tmp1$L1, asDT = TRUE)
tmp1 <- tmp1[, .(prIDs = list(sort(unique(value)))),
             by = (row = L1)]
prot$pep_to_P <- tmp1$prIDs[match(1L:nrow(prot), tmp1$row)]
#
# Identify minimal set of proteins required to explain the identified peptides:
cat(" - Identifying minimal set of protein groups required to explain observed peptides.\n")
## Checking which proteins are contained within another?
## A protein is Leading (at this stage) if it cannot be subsumed into another.
cat("   - Identifying, for each protein ID, which other protein IDs are subsumable within it.\n")
# We actually start in reverse, by checking for each protein which other proteins it can explain
tmp1 <- prot[, c(".pep.ids", "pep_to_P", "Prot.id")]
tmpFl <- tempfile(fileext = ".RDS")
parallel::clusterExport(parClust, "tmpFl", envir = environment())
readr::write_rds(tmp1, tmpFl)
invisible(parallel::clusterCall(parClust, \() {
  assign("tmp1", readr::read_rds(tmpFl), .GlobalEnv)
  return()
}))
unlink(tmpFl)
f0 <- \(x) {
  #x <- tmp1[1L,]
  pep.ids <- unlist(x[[1L]]) #This protein's peptides
  pep_to_P <- unlist(x[[2L]])
  pep_to_P_to_pep <- tmp1$.pep.ids[match(pep_to_P, tmp1$Prot.id)] # list, all peptides from other proteins this protein shares some peptides with
  #Prot.id <- Prot.id[[1L]] #(for testing)
  x3l <- lengths(pep_to_P_to_pep)
  test <- vapply(pep_to_P_to_pep, \(y) { sum(y %in% pep.ids) }, 1L) == x3l # Can the 2nd protein be subsumed in the 1st?
  res <- pep_to_P[test]
  res <- base::setdiff(res, unlist(x[[3L]]))
  names(res) <- NULL
  return(list(res))
}
prot$Contains <- parallel::parApply(parClust, tmp1, 1L, f0)
## Contained proteins and their temp peptide ids:
contained <- data.frame(Contained = sort(as.integer(unique(unlist(prot$Contains)))))
contained$.pep.ids <- prot$.pep.ids[match(contained$Contained, prot$Prot.id)]
## Assign a temporary Leading status...
### A protein may be contained by (an)other protein(s) with exactly the same peptides.
### All non-contained proteins get Leading status.
cat("   - Identifying temporary Leading protein IDs - i.e. those which are not subsumable into a single other protein ID.\n")
prot$Is.Leading <- !prot$Prot.id %in% contained$Contained
## Identify containers of each contained protein:
cat("   - Identifying Leading protein ID containing each subsumable protein ID.\n")
w <- which(lengths(prot$Contains) > 0L)
a <- setNames(prot$Contains[w], prot$Prot.id[w])
a <- a[lengths(a) > 0L]
if (length(a)) {
  temp <- listMelt(a, asDT = TRUE)
  temp <- temp[, .(byWhom = list(L1)),
               by = .(Contained = value)]
  temp <- as.data.frame(temp)
  contained$.byWhom <- temp$byWhom[match(contained$Contained, temp$Contained)]
  ## Container peptide IDs:
  tmp3 <- listMelt(contained$.byWhom)
  tmp3$IDs <- prot$.pep.ids[match(tmp3$value, prot$Prot.id)]
  tmp4 <- listMelt(tmp3$IDs, ColNames = c("value", "row3"), asDT = TRUE)
  tmp4[, c("byWhom", "row")] <- tmp3[tmp4$row3, c("value", "L1")]
  tmp4 <- tmp4[, .(value = list(value)),
               by = .(row = row,
                      byWhom = byWhom)]
  tmp4 <- tmp4[, .(value = list(value)),
               by = .(row = row)]
  contained$Container.pep.ids <- tmp4$value[match(1L:nrow(contained), tmp4$row)]
  ## Now these containers, are they "larger" (N of peptides)?
  l <- lengths(contained$.pep.ids)
  contained$Container.is.larger <- lapply(1L:nrow(contained), \(i) {
    lengths(contained$Container.pep.ids[[i]]) > l[[i]]
  })
  #
  tst <- listMelt(contained$Container.is.larger, asDT = TRUE)
  tst <- tst[, .(sum = sum(value)),
             by = .(row = L1)]
  contained$test <- tst$sum[match(1L:nrow(contained), tst$row)]
  ## Contained proteins get (temporary) Leading status here if they are equivalent to all of their containers,
  ## i.e. if they have no container larger than them:
  contained$Is.Leading <- contained$test == 0L
  prot$Is.Leading[match(contained$Contained[contained$Is.Leading], prot$Prot.id)] <- TRUE
}
## Create protein groups:
cat("   - Creating temporary protein groups.\n")
pg <- prot
if (CustPG) { pg <- pg[!pg$Prot.id %in% prot2$Prot.id,] }
data.table::setDT(pg)
pg <- pg[pg$Is.Leading, c("Protein", "Prot.id", "pep.ids")]
pg <- pg[, .(`Leading protein IDs` = paste(Protein, collapse = ";"),
             lead.protein.ids = paste(Prot.id, collapse = ";")),
         keyby = .(`Peptide IDs` = pep.ids)]
pg <- as.data.frame(pg)
pg$temp.pg.id <- as.character(1L:nrow(pg) + pgoffset)
pg$.pep.ids <- lapply(strsplit(pg$"Peptide IDs", ";"), as.integer)
pg$"Peptides count" <- lengths(pg$.pep.ids)
pg$.Leading.protein.IDs <- strsplit(pg$"Leading protein IDs", ";")
pg$.lead.protein.ids <- strsplit(pg$lead.protein.ids, ";")
#temp <- prot[prot$Is.Leading,]
tmp <- listMelt(pg$.lead.protein.ids, 1L:nrow(pg), ColNames = c("protID", "rowPG"))
tmp$Contains <- prot$Contains[match(tmp$protID, prot$Prot.id)]
tmp2 <- listMelt(tmp$Contains, 1L:nrow(tmp), ColNames = c("Contained", "row1"), asDT = TRUE)
tmp2[, c("rowPG", "protID")] <- tmp[tmp2$row1, c("rowPG", "protID")]
tmp2 <- tmp2[, .(Contained = list(base::setdiff(Contained, protID))),
             by = .(rowPG = rowPG)]
pg$"Also contains" <- tmp2$Contained[match(1L:nrow(pg), tmp2$rowPG)]
#
tmp1 <- pg[, c(".lead.protein.ids", "Also contains")]
tmp2 <- prot[, c("Protein", "Prot.id")]
tmpFls <- vapply(1L:2L, \(x) { tempfile(fileext = ".RDS") }, "")
tmp1 <- rbind(listMelt(pg$".lead.protein.ids", ColNames = c("protID", "row"), asDT = TRUE),
              listMelt(pg$"Also contains", ColNames = c("protID", "row"), asDT = TRUE))
tmp1$Protein <- prot$Protein[match(tmp1$protID, tmp2$Prot.id)]
tmp1 <- tmp1[, .(.Protein.IDs = list(Protein),
                 `Protein IDs` = paste(Protein, collapse = ";")),
             by = .(row = row)]
pg[, c(".Protein.IDs", "Protein IDs")] <- tmp1[match(1L:nrow(pg), tmp1$row), c(".Protein.IDs", "Protein IDs")]
if (CustPG) {
  Custom_PGs$"Also contains" <- lapply(Custom_PGs$.lead.protein.ids, \(x) {
    res <- unlist(prot$Contains[match(unlist(x), prot$Prot.id)])
    res <- base::setdiff(res, x)
    return(res)
  })
  Custom_PGs$.Protein.IDs <- apply(Custom_PGs[,c(".lead.protein.ids", "Also contains")], 1L, \(x) {
    list(prot$Protein[match(unlist(x), prot$Prot.id)])
  })
  Custom_PGs$"Protein IDs" <- vapply(Custom_PGs$.Protein.IDs, \(x) { paste(unlist(x), collapse = ";") }, "")
}
# Assign temporary PG IDs to seq
temp <- listMelt(pg$.pep.ids, pg$temp.pg.id, c("pep.id", "temp.pg.id"), asDT = TRUE)
temp <- temp[, .(temp.pg.ids = paste(temp.pg.id, collapse = ";")),
             by = .(pep.id = pep.id)] 
temp <- as.data.frame(temp)
seq$temp.pg.ids <- temp$temp.pg.ids[match(seq$id, temp$pep.id)]
seq$.temp.pg.ids <- strsplit(seq$temp.pg.ids, ";")
cat(paste0(" ---> ", nrow(pg), " protein groups are required at this stage.\n"))
## Test whether there are further protein groups we can remove?
## In the previous steps, we identified proteins which are not contained within another single protein.
## These then became the founding stones of protein groups.
## However, we should also remove PGs which can be explained by a combination of other PGs.
cat("   - Identifying protein groups which can be subsumed into a combination of other protein groups.\n")
tmp1 <- listMelt(pg$.pep.ids, ColNames = c("pepID", "row"), asDT = TRUE)
# Are all peptides for that PG contained in the others?
# A simpler way to answer this is to ask": are all peptides in the PG assigned to 2+ PGs?
tmp1$wPGs <- seq$.temp.pg.ids[match(tmp1$pepID, seq$id)]
tmp1$tst <- lengths(tmp1$wPGs) > 1L
tmp1 <- tmp1[, .(tst = as.logical(min(tst))),
             by = .(row = row)]
pg$Removable <- tmp1$tst[match(1L:nrow(pg), tmp1$row)]
#
# Calculate individual and global PEPs - useful for later
cat(" - Calculating individual and global protein group PEPs (Posterior Error Probabilities).\n")
temp <- listMelt(strsplit(seq$Proteins, ";"), seq$PEP, asDT = TRUE)
temp <- temp[, .(x = prod(L1)),
             by = .(Group.1 = value)] 
temp <- as.data.frame(temp)
prot$PEP <- NA_real_
w <- which(prot$Protein %in% temp$Group.1)
prot$PEP[w] <- temp$x[match(prot$Protein[w], temp$Group.1)]
if (CustPG) {
  temp <- listMelt(strsplit(seq2$Proteins, ";"), seq2$PEP, asDT = TRUE)
  temp <- temp[, .(x = prod(L1)),
               by = .(Group.1 = value)] 
  temp <- as.data.frame(temp)
  w <- which(prot$Protein %in% temp$Group.1)
  prot$PEP[w] <- temp$x[match(prot$Protein[w], temp$Group.1)]
}
#
tmp1 <- pg$.Protein.IDs
tmp2 <- prot[, c("PEP", "Protein")]
f0 <- \(x) { tmp2$PEP[match(unlist(x), tmp2$Protein)] }
tmp1 <- listMelt(pg$.Protein.IDs, asDT = TRUE)
tmp1$PEP <- prot$PEP[match(tmp1$value, prot$Protein)]
tmp1 <- tmp1[, .(.PEPs = list(PEP),
                 PEPs = paste(PEP, collapse = ";"),
                 PEP = min(PEP)),
             by = .(row = L1)]
pg[, c(".PEPs", "PEPs", "PEP")] <- tmp1[match(1L:nrow(pg), tmp1$row), c(".PEPs", "PEPs", "PEP")]
if (CustPG) {
  Custom_PGs$.Protein.IDs
  Custom_PGs$.PEPs <- lapply(Custom_PGs$.Protein.IDs, \(x) { prot$PEP[match(x, prot$Protein)] })
  Custom_PGs$PEPs <- vapply(Custom_PGs$.PEPs, paste, "", collapse = ";")
  Custom_PGs$PEP <- vapply(Custom_PGs$.PEPs, min, 1)
}
##
## So there are a few PGs which we can subsume into others. Oh joy. Er. Ok. How do we deal with this?
## The first step is to identify a closed group of overlapping PGs and their peptides.
## We then look at these collections of overlapping PGs, and determine at least in approximation the minimum set of protein groups
## required to explain their shared peptides.
## Now doing this the strictly correct way would go like this: for a group of P protein groups,
## "removability" is defined as a property such that a removable group may be explained by a combination of other protein groups.
## - identify the N non-removable PGs in the group: this is the "core"; there are M removable PGs,
## - from m = 1 to m = M-1, pick all possible combinations of m removables: does adding these to the core suffice to explain the peptides?
## Stop as soon as the answer is yes, and keep any of the m removables which figure in a winning combination.
##
## Now this is all fine, but can be terribly ineffective as soon as there are more than a few removables. In fact, there can be several tens of them!
## I thus have to use a simpler, potentially imperfect approximation:
## - after identifying the core, sort pgs by PEP (I am now calculating PEP before, not after) and number of peptides
## - remove all protein groups which do not add to the count
## - also remove one by one (starting from the largest) any group which, if individually removed, does not affect the final peptide count
##
W <- which(pg$Removable)
pg$Remove <- FALSE
if (length(W)) {
  # Get peptides of removable protein groups and sort them by current number of protein groups
  pt <- seq[seq$id %in% unlist(pg$.pep.ids[W]),]
  pt <- pt[order(lengths(pt$.temp.pg.ids), decreasing = TRUE),]
  # For now none of these peptides have been dealt with:
  pt$Done <- FALSE
  while (sum(!pt$Done)) {
    i <- which(!pt$Done)[1L]
    cat(paste0(" ---> Still ", length(unique(unlist(pt$.temp.pg.ids[!pt$Done]))), " unresolved protein groups...\n"))
    #cat(paste0("   i = ", i, "\n"))
    # Get all protein groups for the current peptide:
    ptpg <- pg[pg$temp.pg.id %in% unlist(pt$.temp.pg.ids[i]),]
    # Peptides from those PGs:
    p <- unlist(ptpg$.pep.ids)
    # overlappers = PGs containing at least one "p" peptide:
    tst0 <- vapply(pg$.pep.ids, \(x) { sum(unlist(x) %in% p) }, 1L) > 0L
    overlappers <- which(tst0) # overlappers and all derived variables (overlappers2, goners, keepers, removables) are row indices, not protein group IDs
    # All of the peptides of overlapping PGs:
    P <- unique(unlist(pg$.pep.ids[overlappers])); L <- length(P)
    # Extend overlappers:
    tst2 <- vapply(pg$.pep.ids, \(x) { sum(unlist(x) %in% P) }, 1L) > 0L # Proteins with at least one "P" peptide
    overlappers2 <- which(tst2)
    # Get extended overlappers2 peptides:
    P2 <- unique(unlist(pg$.pep.ids[overlappers2])); L2 <- length(P2)
    # Are they the same?
    #kount <- 1
    #cat(paste0("   Nb. of affected peptides = ", L, "\n"))
    #cat(paste0("   Nb. of affected protein groups = ", length(overlappers), "\n"))
    # Run loop until we have a closed ensemble of peptides and protein groups
    while (L2 > L) {
      #kount <- kount + 1
      P <- P2; L <- L2
      tst2 <- vapply(pg$.pep.ids, \(x) { sum(unlist(x) %in% P) }, 1L) > 0L
      overlappers2 <- which(tst2)
      P2 <- unique(unlist(pg$.pep.ids[overlappers2])); L2 <- length(P2)
      #cat(paste0("   ", kount, "\n"))
      #cat(paste0("   Nb. of affected peptides = ", L2, "\n"))
      #cat(paste0("   Nb. of affected protein groups = ", length(overlappers2), "\n"))
    }
    overlappers <- overlappers2 # There should be at least 3 of these
    P <- P2; L <- length(P)
    n <- length(overlappers)
    # NB: presumably n = 2 is impossible, if we did our job correctly defining Leading proteins above.
    if (n < 2L) { warning("   Check code: I expected at least 3, not 2 or less, overlapping groups here.") } # Always good to check that code is not buggy.
    # Table of peptides and PGs
    peptopg <- as.data.frame(sapply(overlappers, \(x) {
      P %in% unlist(pg$.pep.ids[x])
    }))
    rownames(peptopg) <- P
    colnames(peptopg) <- overlappers
    # Heuristic: we will sort protein groups by number of peptides then PEP and remove those which do not increase coverage
    peep <- setNames(pg$PEP[overlappers], overlappers) # PEP is used to tie-break protein groups of equal ranking
    peep[is.na(peep)] <- 1
    peptopg <- peptopg[, order(peep, decreasing = FALSE)]
    smcl <- setNames(colSums(peptopg), colnames(peptopg)) # Then the main ranking is applied: number of peptides
    peptopg <- peptopg[, order(smcl, decreasing = TRUE)]
    #
    aim <- nrow(peptopg)
    incl <- setNames(rep(FALSE, ncol(peptopg)), colnames(peptopg))
    incl[1L] <- TRUE
    covpep <- peptopg[, 1L]
    curr <- sum(covpep)
    # First filter
    # Heuristic? we take the straighter path to complete peptides coverage
    # There may be cases where the final number of protein groups is not the smallest? Not sure.
    while (curr < aim) {
      wnc <- which(!covpep)
      wni <- names(incl)[!incl]
      tst <- apply(peptopg[wnc, wni, drop = FALSE], 2L, sum)
      bst <- wni[tst == max(tst)]
      if (length(bst) > 1L) { bst <- bst[smcl[bst] == max(smcl[bst])][1L] }
      covpep[wnc] <- peptopg[wnc, bst]
      incl[bst] <- TRUE
      curr <- sum(covpep)
    }
    goners <- as.integer(names(incl)[!incl])
    keepers <- as.integer(names(incl)[incl])
    removables <- keepers[pg$Removable[keepers]] # self explanatory
    P2 <- unique(unlist(pg$.pep.ids[keepers])); L2 <- length(P2)
    stopifnot(L == L2)
    # Second filter: trying to remove some more, one by one
    tst1a <- lengths(pg$.pep.ids[removables]) # Their N of peptides
    tst1b <- pg$PEP[removables] # Their PEPs
    tst1b[is.na(tst1b)] <- 1
    go <- length(removables) > 0L
    while (go) {
      tst2 <- vapply(removables, \(x) {
        length(unique(unlist(pg$.pep.ids[base::setdiff(keepers, x)])))
      }, 1L) == L # Can any of these be dropped without affecting the number of peptides?
      witch <- which(tst2)
      if (length(witch)) {
        witch <- witch[tst1a[witch] == max(tst1a[witch])]
        if (length(witch) > 1L) { witch <- witch[tst1b[witch] == min(tst1b[witch])][1L] }
        wutch <- which(!c(1L:length(removables)) %in% witch)
        goners <- c(goners, removables[witch])
        keepers <- base::setdiff(keepers, goners)
        tst1a <- tst1a[wutch]
        tst1b <- tst1b[wutch]
        removables <- removables[wutch]
        if (!length(removables)) { go <- FALSE }
      } else { go <- FALSE }
    }
    P2 <- unique(unlist(pg$.pep.ids[keepers])); L2 <- length(P2)
    stopifnot(L == L2)
    cat(paste0(" ---> Removing redundant protein group", c("", "s")[(length(goners)>1L)+1L],
               " #", paste(sort(goners), collapse = "-"), "\n"))
    #
    # Kept for reference: alternative, logically correct (notwithstanding unknown potential bugs),
    # but takes centuries to complete for some protein groups.
    # Now how many do we want to remove?
    # We know how many we want to keep at least: all those not flagged pg$Removable == TRUE
    #core <- overlappers[!pg$Removable[overlappers]]; lc <- length(core)
    #removables <- overlappers[pg$Removable[overlappers]]; lr <- length(removables)
    # We want to try all ways to choose between 0 and lr removables, going up, such that we get all peptides.
    #cp <- unique(unlist(pg$.pep.ids[core]))
    #tst <- which(!P %in% cp)
    #if (length(tst)) {
    #  done <- FALSE
    #  r <- 0L
    #  while (!done) {
    #    r <- r+1L
    #    Com <- gtools::combinations(n = lr, r = r, v = removables)
    #    tst <- apply(Com, 1L, \(x) { sum(!P %in% unique(c(cp, unlist(pg$.pep.ids[unlist(x)])))) }) == 0L
    #    if (sum(tst)) { done <- TRUE }
    #  }
    #  if (r == lr) { stop("Something went wrong!") }
    #  keepers <- sort(c(core, unique(as.integer(Com[tst,]))))
    #} else { keepers <- core }
    #goners <- base::setdiff(overlappers, keepers)
    pg$Remove[goners] <- TRUE
    pt$Done[pt$id %in% unlist(pg$.pep.ids[overlappers])] <- TRUE
  }
}
cat("   All protein groups successfully resolved!\n")
## Remove those protein groups we don't need anymore:
pg1 <- pg[!pg$Remove,]
p1 <- unique(unlist(pg1$.pep.ids))
p <- unique(unlist(pg$.pep.ids))
if (sum(!p %in% p1)) {
  # (Sanity check)
  stop("There may be a bug, check!")
}
pg <- pg1; rm(pg1)
pg$Removable <- NULL
pg$Remove <- NULL
pg$Origin <- "Inferred"
if (CustPG) {
  Custom_PGs$Origin <- "Custom"
  #Custom_PGs$"Priority level" <- NULL
  custpgpep <- unique(unlist(Custom_PGs$.pep.ids))
  pg$.pep.ids <- lapply(pg$.pep.ids, \(x) { base::setdiff(x, custpgpep) })
  pg$"Priority level" <- Inf
  pg <- rbind(Custom_PGs, pg)
  seq <- rbind(seq2, seq)
  seq <- seq[match(Pep$"Modified sequence", seq$"Modified sequence"),]
}
doCont <- ((!is.null(ContCol))&&(ContCol %in% colnames(DB)))
if (doCont) {
  tmp <- listMelt(strsplit(pg$`Protein IDs`, ";"), pg$temp.pg.id)
  tmp2 <- DB$`Protein ID`[DB[[ContCol]] == "+"]
  tmp <- tmp[sub("^CON_+", "", tmp$value) %in% sub("^CON_+", "", tmp2),]
  pg$"Potential contaminant" <- c("", "+")[(pg$temp.pg.id %in% tmp$L1)+1L]
}
cat(paste0("   Final number of protein groups: ", nrow(pg), "\n"))
# Protein group IDs
pg <- pg[order(pg$"Peptides count", decreasing = TRUE),]
if (doCont) {
  pg <- pg[order(pg$"Potential contaminant", decreasing = FALSE),]
}
#
pg$id <- 1L:nrow(pg)
if (splitByOrg) {
  pg$id <- pg$id + n_of_PGs
}
# Peptide IDs
pepcolnm <- if (grepl("peptide", Peptide.IDs, ignore.case = TRUE)) {
  sub("ss$", "s", paste0(Peptide.IDs, "s"))
} else { "Peptide IDs" }
pg[[pepcolnm]] <- vapply(pg$.pep.ids, paste, "", collapse = ";")
# Assign final PG IDs to seq
tmp1 <- listMelt(seq$.temp.pg.ids, asDT = TRUE)
tmp1$id <- pg$id[match(tmp1$value, pg$temp.pg.id)]
tmp1 <- tmp1[!is.na(tmp1$id),]
tmp1 <- tmp1[order(tmp1$id),]
tmp1 <- tmp1[, .(.IDs = list(id),
                 IDs = paste(id, collapse = ";")),
             by = .(row = L1)]
seq[, c(".PG.ids", "Protein group IDs")] <- tmp1[match(1L:nrow(seq), tmp1$row), c(".IDs", "IDs")]
pg$temp.pg.id <- NULL
# Unique peptides: peptides unique to a protein group
cat(" - Identifying unique peptides.\n")
l1 <- !grepl(";", seq$"Protein group IDs")
temp <- data.table::data.table(id = seq$id[l1],
                               pgID = seq$"Protein group IDs"[l1])
temp <- temp[, .(pep.ids = list(id)),
             keyby = .(pgID = pgID)]
pg$.unique.pep.ids <- list(NULL)
w <- which(pg$id %in% temp$pgID)
pg$.unique.pep.ids[w] <- temp$pep.ids[match(pg$id[w], temp$pgID)]
pg$"Unique peptide IDs" <- vapply(pg$.unique.pep.ids, paste, "", collapse = ";")
pg$"Unique peptides" <- lengths(pg$.unique.pep.ids)
# Apply Occam's razor:
# If a peptide is shared, it should go to the protein group with the most identified peptides in total.
# (Here, "peptide is razor" means unique OR razor.)
# In cases of ties, we will go for lowest PEP.
# Note that these ties would mean that some few protein groups which do get reported may not have any razor peptides.
cat(" - Applying Occam's razor,\n     i.e. assigning \"razor\" status to shared peptides for the parent group with the highest overall peptides count,\n     with protein priority (if available) and PEP as tie-breakers (in that order).\n")
seq$.Protein.group.IDs <- strsplit(seq$"Protein group IDs", ";")
## Note here: the code below allows for multiple razor PGs.
# This can happen if all have the same priority, number of peptides and PEPs (can happen with NA PEP values).
tmp <- data.frame(PG_IDs = unique(seq$"Protein group IDs"))
tmp$.PG_IDs <- strsplit(tmp$PG_IDs, ";")
kol <- c("id", "Peptides count", "PEP", "Priority level")
kol <- base::intersect(kol, colnames(pg))
tmp1 <- tmp$.PG_IDs
tmp2 <- pg[, kol]
tmpFls <- vapply(1L:3L, \(x) { tempfile(fileext = ".RDS") }, "")
parallel::clusterExport(parClust, "tmpFls", envir = environment())
readr::write_rds(tmp1, tmpFls[1L])
readr::write_rds(tmp2, tmpFls[2L])
readr::write_rds(CustPG, tmpFls[3L])
invisible(parallel::clusterCall(parClust, \() {
  assign("tmp1", readr::read_rds(tmpFls[1L]), .GlobalEnv)
  assign("tmp2", readr::read_rds(tmpFls[2L]), .GlobalEnv)
  assign("CustPG", readr::read_rds(tmpFls[3L]), .GlobalEnv)
  return()
}))
unlink(tmpFls)
f0 <- \(x) {
  #sapply(tmp1, \(x) {
  #x <- unlist(tmp1[1L])
  x <- unlist(x)
  w1 <- match(x, tmp2$id) # ~ x
  if (CustPG) {
    p <- tmp2$`Priority level`[w1] # ~ x
    mp <- min(p) # ~ x
    w2 <- which(tmp2$`Priority level`[w1] == mp) # ~ x[w2]  
  } else { w2 <- 1L:length(w1) }
  k <- tmp2$"Peptides count"[w1[w2]]
  Mk <- max(k) # ~ x[w2]
  w3 <- which(k == Mk) # ~ x[w2][w3]
  if (length(w3) > 1) {
    z <- tmp2$PEP[w1[w2]] # ~ x[w2]
    z[is.na(z)] <- 1
    z[-w3] <- 1
    mP <- min(z, na.rm = TRUE)
    w3 <- which((k == Mk)&(z == mP)) # ~ x[w2][w3]
  }
  return(x[w2][w3])
}
#environment(f0) <- .GlobalEnv
tmp$.Razor_PG_IDs <- parallel::parSapply(parClust, tmp1, f0)
tmp$Razor_PG_IDs <- vapply(tmp$.Razor_PG_IDs, paste, "", collapse = ";")
seq[, c(".Razor.protein.group.ID", "Razor protein group ID")] <- tmp[match(seq$"Protein group IDs", tmp$PG_IDs),
                                                                     c(".Razor_PG_IDs", "Razor_PG_IDs")]
#
c1 <- c("Leading proteins", "Leading razor proteins")
c2 <- c(".Protein.group.IDs", ".Razor.protein.group.ID")
for (i in 1L:length(c1)) { #i <- 1L
  tmp1 <- listMelt(seq[[c2[i]]], ColNames = c("pgID", "row"))
  tmp1 <- listMelt(pg$.Leading.protein.IDs[match(tmp1$pgID, pg$id)], tmp1$row, c("leadPrID", "row"), TRUE)
  tmp1$PEP <- prot$PEP[match(tmp1$leadPrID, prot$Protein)]
  tmp1 <- tmp1[order(tmp1$PEP, decreasing = FALSE), c("leadPrID", "row")]
  tmp1 <- tmp1[, .(leadPrIDs = paste(leadPrID, collapse = ";")),
               by = .(row = row)]
  seq[[c1[i]]] <- tmp1$leadPrIDs[match(1L:nrow(seq), tmp1$row)]
}
# Sort by PEP
if (CustPG) { # Priority-based razor peptides decision for custom protein groups
  for (i in priorities) { #i <- 1
    wcp <- which(pg$"Priority level" == i)
    p <- unique(unlist(pg$.pep.ids[wcp]))
    whp <- which(pg$"Priority level" < i)
    if (length(whp)) { p <- base::setdiff(!p, unique(unlist(pg$.pep.ids[whp]))) } # Remove peptides matching protein groups of higher priority
    m <- match(p, seq$id)
    seq$.Razor.protein.group.ID[m] <- lapply(seq$.Protein.group.IDs[m], \(x) { base::intersect(x %in% pg$id[wcp]) })
  }
}
# Update columns
seq$"Razor protein group ID" <- vapply(seq$.Razor.protein.group.ID, paste, "", collapse = ";")
seq$"Protein group IDs" <- vapply(seq$.Protein.group.IDs, paste, "", collapse = ";")
#
seq$"Leading razor proteins" <- vapply(strsplit(seq$"Leading razor proteins", ";"), \(x) { unlist(x)[1L] }, "")
if (length(pg$.lead.protein.ids) != length(unique(pg$.lead.protein.ids))) {
  stop("A \"leading\" protein cannot be in several protein groups!")
}
# Assign razor peptides to protein groups
temp <- listMelt(as.list(seq$.Razor.protein.group.ID), seq$id, asDT = TRUE)
temp <- temp[, .(x = paste(sort(as.integer(unique(L1))), collapse = ";")), 
             keyby = .(Group.1 = value)]
pg$"Razor peptide IDs" <- temp$x[match(pg$id, as.integer(temp$Group.1))]
pg$.razor.pep.ids <- lapply(strsplit(pg$`Razor peptide IDs`, ";"), as.integer)
pg$"Razor + unique peptides" <- lengths(pg$.razor.pep.ids)
pg$.Peptide.is.razor <- apply(pg[, c(".pep.ids", ".razor.pep.ids")], 1L, \(x) { #x <- pg[1, c(".pep.ids", ".razor.pep.ids")]
  c("False", "True")[(x[[1L]] %in% x[[2L]])+1L]
})
# 
pg$"Peptide is razor" <- vapply(pg$.Peptide.is.razor, paste, "", collapse = ";")
# Rewrite below using data.table was done by Claude: mine was much more cumbersome.
c1  <- paste0("Peptide counts (", c("all", "unique", "razor+unique"), ")")
c2  <- c(".pep.ids", ".unique.pep.ids", ".razor.pep.ids")
cnt <- c("n_all", "n_unique", "n_razor")
ids <- listMelt(pg$.Protein.IDs, ColNames = c("Protein", "id"), asDT = TRUE)
ids[, ord := seq_len(.N),
    by = id]
keep <- !duplicated(prot$Protein)
prot_pep <- data.table::data.table(Protein = rep.int(prot$Protein[keep], lengths(prot$.pep.ids[keep])),
                                   pep = unlist(prot$.pep.ids[keep], use.names = FALSE))
pairs <- prot_pep[ids, on = "Protein", allow.cartesian = TRUE]
for (k in 1L:3L) {
  hit <- listMelt(pg[[c2[k]]], ColNames = c("pep", "id"), asDT = TRUE)[!is.na(pep)]
  pairs[, (cnt[k]) := 0L]
  pairs[hit, on = .(id, pep), (cnt[k]) := 1L]
}
res <- pairs[, lapply(.SD, sum),
             by = .(id, ord),
             .SDcols = cnt][order(id, ord),
                            lapply(.SD, paste, collapse = ";"),
                            by = id,
                            .SDcols = cnt]
pg[, c1] <- ""
pg[res$id, c1] <- as.data.frame(res[, ..cnt])
#
test1 <- as.integer(unlist(strsplit(pg$"Peptide counts (all)", ";")))
test2 <- as.integer(unlist(strsplit(pg$"Peptide counts (unique)", ";")))
test3 <- as.integer(unlist(strsplit(pg$"Peptide counts (razor+unique)", ";")))
temp <- listMelt(strsplit(pg$`Protein IDs`, ";"), pg$id, c("Protein ID", "Protein Group ID"))
temp$"All peptides" <- test1
temp$"Unique peptides" <- test2
temp$"Razor+unique peptides" <- test3
w <- which(test2 > test1)
if (length(w)) {
  stop("No way a protein group can have more unique peptides than total peptides! Check your code!")
  #View(temp[w,])
  #m <- match(unique(temp$`Protein Group ID`[w]), pg$id)
  #View(pg[m, c(c1, c2)])
}
w <- which(test2 > test3)
if (length(w)) {
  stop("No way a protein group can have more unique peptides than razor peptides! Check your code!")
  #View(temp[w,])
  #m <- match(unique(temp$`Protein Group ID`[w]), pg$id)
  #View(pg[m, c(c1, c2)])
}
w <- which(test3 > test1)
if (length(w)) {
  stop("No way a protein group can have more razor peptides than total peptides! Check your code!")
  #View(temp[w,])
  #m <- match(unique(temp$`Protein Group ID`[w]), pg$id)
  #View(pg[m, c(c1, c2)])
}
# Mark protein groups which the user should want to keep if they want to be stringent
# (at least Npep unique + razor peptides)
if (Npep > 1L) {
  pg[[paste0("Quality filter: min ", Npep, " razor&unique pep.")]] <- c("", "Keep")[
    vapply(strsplit(pg$`Peptide counts (razor+unique)`, ";"), \(x) {
      as.integer(x[[1L]]) >= Npep
    }, TRUE) + 1L]
}
#
cat(" - Getting sequence coverage.\n")
c1 <- c("Sequence coverage [%]", "Unique + razor sequence coverage [%]", "Unique sequence coverage [%]")
c2 <- c(".Leading.protein.IDs", ".pep.ids", ".razor.pep.ids", ".unique.pep.ids")
tmp1 <- DB[, c("Protein ID", "Sequence")]
tmp1$`Protein ID` <- sub("^CON_+", "", tmp1$`Protein ID`)
tmp2 <- pg[, c2]
tmp3 <- seq[, c("id", "Sequence")]
tmpFls <- vapply(1L:4L, \(x) { tempfile(fileext = ".RDS") }, "")
parallel::clusterExport(parClust, "tmpFls", envir = environment())
readr::write_rds(tmp1, tmpFls[1L])
readr::write_rds(tmp2, tmpFls[2L])
readr::write_rds(tmp3, tmpFls[3L])
readr::write_rds(c2, tmpFls[4L])
invisible(parallel::clusterCall(parClust, \() {
  assign("tmp1", readr::read_rds(tmpFls[1L]), .GlobalEnv)
  assign("tmp2", readr::read_rds(tmpFls[2L]), .GlobalEnv)
  assign("tmp3", readr::read_rds(tmpFls[3L]), .GlobalEnv)
  assign("c2", readr::read_rds(tmpFls[4L]), .GlobalEnv)
  return()
}))
unlink(tmpFls)
f0 <- \(x) {
  #x <- tmp2[1, c2]
  s <- sub("^CON_+", "", unlist(x[[1L]])[1L])
  w <- match(s, tmp1$"Protein ID")
  ptids <- unlist(x[[2L]])
  prids <- unlist(x[[3L]])
  puids <- unlist(x[[4L]])
  pt <- tmp3$Sequence[match(ptids, tmp3$id)]
  pr <- tmp3$Sequence[match(prids, tmp3$id)]
  pu <- tmp3$Sequence[match(puids, tmp3$id)]
  print(pt)
  print(pr)
  print(pu)
  res <- rep(NA, 3L)
  if (length(w)) {
    s <- setNames(tmp1$Sequence[w], s)
    stopifnot(length(pt) > 0L)
    res[1L] <- round(100*Coverage(s, pt), 1L)
    if (length(pr)) { res[2L] <- round(100*Coverage(s, pr), 1L) } 
    if (length(pu)) { res[3L] <- round(100*Coverage(s, pu), 1L) } 
  } else { warning(paste0("Protein accession ", s, " was not found in the database!")) }
  #}
  return(res)
}
source(parSrc)
pg[, c1] <- t(parallel::parApply(parClust, tmp2[, c2], 1L, f0))
# Note:
#   The "Gene names" and "Protein names" columns will not match the MQ ones...
#   but after checking that seems to be because the MQ list is missing some!
cat(" - Getting Gene and Protein names + annotations...\n")
ca <- c("Common Names", "Common Name")
if (sum(ca %in% colnames(DB))) {
  ca <- base::intersect(ca, colnames(DB))
  test <- vapply(ca, \(x) { length(is.na(DB[[x]])) }, 1L) > 0L
  ca <- ca[test][1L]
}
cb <- c("Genes", "Gene", "Gene IDs", "Gene ID")
if (sum(cb %in% colnames(DB))) {
  cb <- base::intersect(cb, colnames(DB))
  test <- vapply(cb, \(x) { length(is.na(DB[[x]])) }, 1L) > 0L
  cb <- cb[test][1L]
}
c1 <- c("Number of proteins", "Fasta headers")
c2 <- "temp"
if (length(ca) == 1L) {
  c1 <- c(c1, "Protein names")
  c2 <- c(c2, ca)
}
if (length(cb) == 1L) {
  c1 <- c(c1, "Gene names")
  c2 <- c(c2, cb)
}
kol <- c("Protein ID", "Header", ca, cb)
kol <- base::intersect(kol, colnames(DB))
temp <- DB[, kol]
temp$temp <- sub("^>", "", temp$Header)
temp$"Protein ID" <- sub("^CON_+", "", temp$"Protein ID")
temp2b <- sub("^CON_+", "", gsub(";CON_+", ";", pg$"Protein IDs"))
tmp2 <- strsplit(temp2b, ";")
tmp1 <- temp[, c("Protein ID", c2)]
tmpFls <- vapply(1L:6L, \(x) { tempfile(fileext = ".RDS") }, "")
parallel::clusterExport(parClust, "tmpFls", envir = environment())
readr::write_rds(tmp1, tmpFls[1L])
readr::write_rds(tmp2, tmpFls[2L])
readr::write_rds(c1, tmpFls[3L])
readr::write_rds(c2, tmpFls[4L])
readr::write_rds(ca, tmpFls[5L])
readr::write_rds(cb, tmpFls[6L])
invisible(parallel::clusterCall(parClust, \() {
  assign("tmp1", readr::read_rds(tmpFls[1L]), .GlobalEnv)
  assign("tmp2", readr::read_rds(tmpFls[2L]), .GlobalEnv)
  assign("c1", readr::read_rds(tmpFls[3L]), .GlobalEnv)
  assign("c2", readr::read_rds(tmpFls[4L]), .GlobalEnv)
  assign("ca", readr::read_rds(tmpFls[5L]), .GlobalEnv)
  assign("cb", readr::read_rds(tmpFls[6L]), .GlobalEnv)
  return()
}))
unlink(tmpFls)
f0 <- \(x) { #x <- tmp2[1L]
  x <- unlist(x)
  w <- match(x, tmp1$"Protein ID")
  res <- tmp1[w, c2]
  if (length(w) > 1L) {
    res2 <- paste(res[[c2[1L]]], collapse = ";")
    if (length(ca) == 1L) { res2 <- c(res2, paste(res[[ca]], collapse = ";")) }
    if (length(cb) == 1L) {
      y <- data.table::data.table(a = res[[cb]])
      y <- y[, .(x = length(a)),
             by = .(Group.1 = a)]
      y <- as.data.frame(y)
      y <- y[order(y$x, decreasing = TRUE),]
      res2 <- c(res2, paste(y$Group.1, collapse = ";"))
    }
  } else { res2 <- unlist(res) }
  if (length(res2) != (length(c1)-1L)) { stop(paste(x, collapse = ";")) }
  return(setNames(c(length(x), res2), c1))
}
#environment(f0) <- .GlobalEnv
tst2 <- parallel::parLapply(parClust, tmp2, f0)
tst2 <- as.data.frame(do.call(rbind, tst2))
pg[, c1] <- tst2
#vapply(1L:4L, \(x) { sum(tst1[[x]] != tst2[[x]]) }, 1L)
pg$"Number of proteins" <- as.integer(pg$"Number of proteins")
DB$temp <- NULL
DB$"Sequence length" <- nchar(DB$Sequence)
pg[, c("Mol. weight [kDa]", "Sequence length")] <- DB[match(sub(";.*", "", pg$"Leading protein IDs"), DB$"Protein ID"),
                                                      c("MW [kDa]", "Sequence length")]
tmp1 <- listMelt(strsplit(pg$"Leading protein IDs", ";"), asDT = TRUE)
tmp1$L <- as.character(DB$`Sequence length`[match(tmp1$value, DB$`Protein ID`)])
tmp1 <- tmp1[, .(lengths = paste(L, collapse = ";")),
             by = (row = L1)]
pg$"Sequence lengths2" <- tmp1$lengths[match(1L:nrow(pg), tmp1$row)]
tmp1 <- listMelt(pg$.pep.ids, ColNames = c("pepID", "row"), asDT = TRUE)
tmp2 <- nchar(gsub(paste(c(AA, "_"), collapse = "|"), "", seq$"Modified sequence"))
tmp1$tmp2 <- tmp2[match(tmp1$pepID, seq$id)]
f0 <- \(x) {
  c("", "+")[min(x > 0L)+1L]
}
tmp1 <- tmp1[, .(tst = f0(tmp2)),
             by = .(row = row)]
temp <- tmp1$tst[match(1L:nrow(pg), tmp1$row)]
u <- unique(temp)
if (identical(u, "+")) {
  warning("Do we really expect all of our peptides to be modified?")
}
pg$"Only identified by site" <- temp
if (!misFun(Ev)) {
  seq$.Evidence.IDs <- strsplit(seq$"Evidence IDs", ";")
  tmp2 <- tmp1 <- listMelt(pg$.pep.ids, ColNames = c("pepID", "row"))
  tmp2$evIDs <- seq$.Evidence.IDs[match(tmp2$pepID, seq$id)]
  tmp2 <- listMelt(tmp2$evIDs, ColNames = c("evID", "tmpRow"), asDT = TRUE)
  tmp2[, c("pepID", "row")] <- tmp1[tmp2$tmpRow, c("pepID", "row")]
  tmp2 <- tmp2[, .(evIDs = paste(sort(as.integer(unique(evID))), collapse = ";")),
               by = .(row = row)]
  pg$"Evidence IDs" <- tmp2$evIDs[match(1L:nrow(pg), tmp2$row)]
  Ev$"Protein group IDs" <- seq$"Protein group IDs"[match(Ev$"Modified sequence", seq$"Modified sequence")]
  if ("MS/MS count" %in% colnames(Ev)) {
    temp <- data.table::setDT(Ev[, c("MS/MS count", "Modified sequence")])
    temp <- temp[, .(n_MS2s = sum(`MS/MS count`)),
                 keyby = .(modSeq = `Modified sequence`)]
    seq$"MS/MS count" <- temp$n_MS2s[match(seq$"Modified sequence", temp$modSeq)]
    temp <- listMelt(seq$.Protein.group.IDs, seq$"MS/MS count", c("pgIDs", "n_MS2s"), TRUE)
    temp$n_MS2s <- as.integer(temp$n_MS2s)
    temp <- temp[, .(n_MS2s = sum(n_MS2s)),
                 keyby = .(pgIDs = pgIDs)]
    pg$"MS/MS count" <- temp$n_MS2s[match(pg$id, temp$pgIDs)]
  }
  # Note:
  # I don't understand how MaxQuant determines identification type.
  # There are protein groups with plenty of direct evidences with their own MS/MS, PEP etc...
  # but the protein group identification type is "By matching"?!?!
}
pg <- pg[, which(!grepl("^\\.", colnames(pg)))]
pg$"Also contains" <- NULL
pg$lead.protein.ids <- NULL
pg$unique.pep.ids <- NULL
pg$"Unique peptide IDs" <- ""
temp <- data.table::setDT(seq[!grepl(";", seq$"Protein group IDs"),
                              c("id", "Protein group IDs")])
temp <- temp[, .(ids = paste(sort(id), collapse = ";")),
             keyby = .(pgIDs = `Protein group IDs`)]
w <- which(pg$id %in% temp$pgIDs)
pg$"Unique peptide IDs"[w] <- temp$ids[match(pg$id[w], temp$pgIDs)]
# Some more annotations
if ("Common Name" %in% colnames(DB)) {
  tmp1 <- listMelt(strsplit(pg$"Leading protein IDs", ";"), ColNames = c("leadPrIDs", "row"), asDT = TRUE)
  tmp1$name <- DB$"Common Name"[match(tmp1$leadPrIDs, db$"Protein ID")]
  tmp1 <- tmp1[(!is.na(tmp1$name)) & (tmp1$name != ""),]
  tmp1 <- tmp1[, .(name = paste(name, collapse = ";")),
               by = .(row = row)]
  pg$"Common Names" <- tmp1$name[match(1L:nrow(pg), tmp1$row)]
  w <- which(is.na(pg$"Common Names"))
  pg$"Common Names"[w] <- ""
  #
  w <- which(pg$"Common Names" == "")
  if (length(w)) { if ("Protein names" %in% colnames(pg)) {
    pg$"Common Names"[w] <- pg$"Protein names"[w]
  } }
  w <- which(pg$"Common Names" == "")
  if ((length(w))&&("Names" %in% colnames(pg))) { pg$"Common Names"[w] <- pg$Names[w] }
  pg$"Common Name (short)" <- vapply(strsplit(pg$"Common Names", ";"), \(x) {
    x <- unlist(x)[1L]
    if (!length(x)) { x <- "" }
    return(x)
  }, "")
  id.col <- base::intersect(c("Names", "Protein IDs", "Gene names"), colnames(pg))
  w <- which(pg$"Common Name (short)" == "")
  if ((length(id.col))&&(length(w))) {
    pg$"Common Name (short)"[w] <- apply(pg[w, id.col], 1L, \(x) {
      ww <- which(!x %in% c("", " ", "NA", NA))
      x <- if (length(ww)) { x[ww[1L]] } else { "" }
      return(x)
    })
  }
}
# Taxonomy
if (!is.na(dbOrgCol)) {
  tmp <- listMelt(strsplit(pg$`Protein IDs`, ";"), ColNames = c("ID", "Row"), asDT = TRUE)
  tmp$Org <- DB[match(tmp$ID, DB$`Protein ID`), dbOrgCol]
  tmp <- tmp[, .(Org = paste(sort(unique(Org)), collapse = ";")),
             by = .(Row = Row)]
  tmp <- as.data.frame(tmp)
  pgOrgCol %<o% c("Organism", "Organism(s)")[(sum(grepl(";", tmp$Org)) > 0L)+1L]
  pg[[pgOrgCol]] <- tmp$Org[match(1L:nrow(pg), tmp$Row)]
}
# PG label column
pg$temp <- gsub(";.+", ";...", pg$"Leading protein IDs")
tmp <- pg[, c("temp", "Common Name (short)")]
tmp$Label <- apply(tmp, 1L, \(x) {
  x <- { if ((is.na(x[[2L]]))||(x[[2L]] == "NA")) { x[[1L]] } else { paste0(x[[1L]], " - ", x[[2L]]) } }
  return(x)
})
pg$Label <- tmp$Label
pg$temp <- NULL
if (!CustPG) {
  pg$"Priority level" <- NULL
  pg$Origin <- NULL
}
seq <- seq[, which(!grepl("^\\.|^temp", colnames(seq)))]
colnames(seq)[match("id", colnames(seq))] <- Peptide.IDs
colnames(seq)[match("Evidence IDs", colnames(seq))] <- Evidence.IDs
colnames(seq)[match("Proteins", colnames(seq))] <- Proteins.col
Pep[, colnames(seq)] <- seq[match(Pep$"Modified sequence", seq$"Modified sequence"),]
#
PG_assembl <- list(Protein.groups = pg,
                   Peptides = Pep,
                   Database = DB)
if (!misFun(Ev)) { PG_assembl$Evidences <- Ev }
invisible(parallel::clusterCall(parClust, \(x) {
  try(rm(tmp1, tmp2, tmp3, tmp4, CustPG, c1, c2, ca, cb), silent = TRUE)
  return()
}))
tm2 <- Sys.time()
try({
  cat(paste0(nrow(pg), " protein groups assembled in ", sub("^Time difference of ", "", capture.output(tm2-tm1)), "\n"))
}, silent = TRUE)
#
if (splitByOrg) {
  n_of_PGs <- n_of_PGs + nrow(pg)
}
cat("\nDone!\n\n")
