#' PG_assemble
#'
#' @description 
#' This former function is now a source. And faster.
#' It takes a peptides and a data base data frames and returns a list with:
#'  - a data.frame of parsimoniously-inferred protein groups,
#'  - the peptides file with a few extra columns,
#'  - the database file with a few extra columns.
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
#' @param Npep Minimum number of razor or unique peptides to tag a protein group as a true discovery (others will still be reported). Default = 2
#' @param cl Already have a cluster handy? Why spend time making a new one, which on top of that may invalidate the old one. Just pass it along!
#' 
#' @examples
#' temp <- PG_assemble(Pep = pep, DB = db)
#' 
#' @import data.table
#' @export

cat("Protein groups inference step:\nAssembling peptides into the minimum number of protein groups required to explain the data...\n")
Peptide.IDs <- "id"
Proteins.col <- "Proteins"
Evidence.IDs <- "Evidence IDs"
N.reserved <- 1
Custom_PGs <- NA
ContCol <- "Potential contaminant"
Npep <- 2
cl <- parClust
Pep <- pep
Ev <- ev
DB <- db
N.clust <- 55
Custom_PGs <- custPGs
Npep = NPep
TESTING <- TRUE

# An immense thank you goes to Steve Weston from
# https://stackoverflow.com/questions/19467133/performance-of-clusterapply-deteriorates-when-called-inside-a-function?rq=4
# for this thread explaining how to use parallel within a function!
#
TESTING <- FALSE
#proteoCraft::DefArg(proteoCraft::PG_assemble)
#Pep <- pep ;Ev <- ev ;DB <- db ;N.clust <- 55; Custom_PGs <- custPGs; Npep = NPep; TESTING <- TRUE
#Pep <- pep[1:1000,] ;Ev <- ev ;DB <- db ;N.clust <- 55; Custom_PGs <- custPGs; Npep = NPep; TESTING <- TRUE
#
if (TESTING) {
  tm1 <<- Sys.time()
  # Note:
  # This is not a perfect alternative to missing but will work in most cases, unless x matches a function imported by a package 
  misFun <- function(x) { return(!exists(deparse(substitute(x)))) }
} else { misFun <- missing }
#
dc <- parallel::detectCores()
if (misFun(N.reserved)) { N.reserved <- 1 }
if (misFun(N.clust)) {
  N.clust <- max(c(dc-N.reserved, 1))
} else {
  if (N.clust > max(c(dc-N.reserved, 1))) {
    warning("More cores specified than allowed, I will ignore the specified number! You should always leave at least one free for other processes, see the \"N.reserved\" argument.")
    N.clust <- max(c(dc-N.reserved, 1))
  }
}
# Create cluster (some steps are slow otherwise)
cleanUp <- FALSE
if (misFun(cl)) {
  dc <- parallel::detectCores()
  if (misFun(N.reserved)) { N.reserved <- 1 }
  if (misFun(N.clust)) {
    N.clust <- max(c(dc-N.reserved, 1))
  } else {
    if (N.clust > max(c(dc-N.reserved, 1))) {
      warning("More cores specified than allowed, I will ignore the specified number! You should always leave at least one free for other processes, see the \"N.reserved\" argument.")
      N.clust <- max(c(dc-N.reserved, 1))
    }
  }
  cl <- parallel::makeCluster(N.clust, type = "SOCK")
  cleanUp <- TRUE
}
#
invisible(parallel::clusterCall(cl, function() {
  library(data.table)
  return()
}))
#
if (!"MW [kDa]" %in% colnames(DB)) {
  DB$"MW [kDa]" <- round(suppressWarnings(parallel::parSapply(cl, DB$Sequence, Peptides::mw))/1000, 3)
}
# Peptides (a new data frame in case I screw up, but this is essentially the same as the original one)
seq <- data.frame(Proteins = Pep[[Proteins.col]], "Modified sequence" = Pep$"Modified sequence",
                  Sequence = Pep$Sequence, id = Pep[[Peptide.IDs]], PEP = Pep$PEP, check.names = FALSE)
if (!misFun(Ev)) { seq$"Evidence IDs" <- Pep[[Evidence.IDs]] }
tmpids <- DB$"Protein ID"
seq$.Proteins <- strsplit(seq$Proteins, ";")
# Rarely we can have peptides without protein assignments in the input.
# Obviously those are useless for the purpose of assembling protein groups!
# There may also be protein IDs not from the data base provided, which should be removed.
rws <- 1:nrow(seq)
tmp1 <- seq$.Proteins
tmp1 <- proteoCraft::listMelt(tmp1)
tmp1 <- tmp1[which(tmp1$value %in% DB$"Protein ID"),]
tmp1 <- data.table::as.data.table(tmp1)
tmp1 <- tmp1[, list(value = list(value)), by = list(L1)]
tmp1 <- tmp1[order(tmp1$L1),]
seq$.Proteins2 <- lapply(rws, c)
w <- which(rws %in% tmp1$L1)
seq$.Proteins[w] <- tmp1$value
#
wY <- which(vapply(seq$.Proteins, length, 1) > 0)
wN <- which(vapply(seq$.Proteins, length, 1) == 0)
if (length(wN) > 0) { warning("Not all of the protein IDs in the peptides file are in the data base!") }
seq <- seq[wY,]
seq$Proteins <- vapply(seq$Proteins, paste, "", collapse = ";")
# (Creating temporary IDs is convenient as I don't even need to use "match" for this.)
proffset <- pgoffset <- 0
CustPG <- !is.na(list(Custom = Custom_PGs)) # Allows testing despite differences in dimensions
if (CustPG) {
  #Custom_PGs <- Custom_PGs[, "Leading protein IDs", drop = FALSE]
  stopifnot("Leading protein IDs" %in% colnames(Custom_PGs))
  if (!"Priority level" %in% colnames(Custom_PGs)) { Custom_PGs$"Priority level" <- 1 }
  Custom_PGs <- Custom_PGs[,c("Leading protein IDs", "Priority level")]
  Custom_PGs$"Priority level" <- as.integer(Custom_PGs$"Priority level")
  stopifnot(sum(is.na(Custom_PGs$"Priority level")) == 0)
  # Create proteins table for custom PGs
  # Expand custom PGs columns
  tmp <- unique(unlist(strsplit(Custom_PGs$`Leading protein IDs`, ";")))
  g <- proteoCraft::grsep(tmp, x = seq$Proteins)
  if (length(g)) {
    seq2 <- seq[g,]
    seq <- seq[which(! (1:nrow(seq)) %in% g),]
    Custom_PGs$.pep.ids <- list(NA)
    priorities <- sort(unique(Custom_PGs$"Priority level"), decreasing = TRUE) # Sorted by decreasing value (increasing priority)
    for (i in max(priorities):1) { #i <- max(priorities)
      wcp <- which(Custom_PGs$"Priority level" == i)
      Custom_PGs$.pep.ids[wcp] <- lapply(strsplit(Custom_PGs$`Leading protein IDs`[wcp], ";"), function(x) {
        seq2$id[proteoCraft::grsep(unlist(x), x = seq2$Proteins)]
      })
      wlp <- which(Custom_PGs$"Priority level" > i)
      if (length(wlp)) {
        Custom_PGs$.pep.ids[wlp] <- lapply(Custom_PGs$.pep.ids[wlp], function(x) {
          x[which(!x %in% unlist(Custom_PGs$.pep.ids[wcp]))]
        })
      }
    }
    temp <- proteoCraft::listMelt(strsplit(seq2$Proteins, ";") , seq2$id, c("Protein", "pep.id"))
    w <- which(!as.character(temp$Protein) %in% c("", " ", "NA"))
    if (length(w) != nrow(temp)) {
      warning("Some protein IDs seem dodgy, we will remove them but you should check what happened.")
      temp <- temp[w,]
    }
    prot2 <- magrittr::set_colnames(aggregate(temp$pep.id, list(temp$Protein),
                                              paste,
                                              collapse = ";"),
                                    c("Protein", "pep.ids"))
    prot2$.pep.ids <- lapply(strsplit(prot2$pep.ids, ";"), as.integer)
    # Sort by N of peptides
    prot2$"Peptides count" <- vapply(prot2$.pep.ids, length, 1)
    prot2 <- prot2[order(prot2$"Peptides count", decreasing = TRUE),]
    prot2$Prot.id <- as.character(c(1:nrow(prot2)))
    proffset <- nrow(prot2)
    Custom_PGs$.Leading.protein.IDs <- strsplit(Custom_PGs$`Leading protein IDs`, ";")
    Custom_PGs$"Peptide IDs" <- vapply(Custom_PGs$.pep.ids, paste, "", collapse = ";")  
    Custom_PGs$.lead.protein.ids <- lapply(Custom_PGs$.Leading.protein.IDs, function(x) {
      prot2$Prot.id[match(x, prot2$Protein)]
    })
    Custom_PGs$lead.protein.ids <- vapply(Custom_PGs$.lead.protein.ids, paste, "", collapse = ";")
    Custom_PGs$"Peptides count" <- vapply(Custom_PGs$.pep.ids, length, 1)
    Custom_PGs$temp.pg.id <- as.character(1:nrow(Custom_PGs))
    pgoffset <- nrow(Custom_PGs)
    seq2$.Proteins <- strsplit(seq2$Proteins, ";")
    seq2$.Prot.ids <- lapply(seq2$.Proteins, function(x) { prot2$Prot.id[match(unlist(x), prot2$Protein)] })
    temp <- proteoCraft::listMelt(Custom_PGs$.pep.ids, Custom_PGs$temp.pg.id, c("pep.id", "temp.pg.id"))
    temp <- magrittr::set_colnames(aggregate(temp$temp.pg.id,
                                             list(temp$pep.id),
                                             paste,
                                             collapse = ";"),
                                   c("pep.id", "temp.pg.ids"))
    seq2$temp.pg.ids <- temp$temp.pg.ids[match(seq2$id, temp$pep.id)]
    seq2$.temp.pg.ids <- strsplit(seq2$temp.pg.ids, ";")
    w1 <- which(Custom_PGs$`Peptides count` > 0)
    w0 <- which(Custom_PGs$`Peptides count` == 0)
    l <- length(w0)
    if (l) {
      warning(paste0("The following custom protein group ", c("was", "were")[(l>1)+1], " not found:",
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
temp <- proteoCraft::listMelt(strsplit(seq$Proteins, ";"), seq$id, c("Protein", "pep.id"))
w <- which(!as.character(temp$Protein) %in% c("", " ", "NA"))
if (length(w) != nrow(temp)) {
  warning("Some protein IDs seem dodgy, we will remove them but you should check what happened.")
  temp <- temp[w,]
}
prot <- magrittr::set_colnames(aggregate(temp$pep.id,
                                         list(temp$Protein),
                                         paste,
                                         collapse = ";"),
                               c("Protein", "pep.ids"))
# data.table rewrite of above is actually slightly slower
# tmp1 <- data.table::data.table(pep.id = temp$pep.id,
#                                Protein = temp$Protein)
# tmp1 <- tmp1[, list(pep.ids = paste(pep.id, collapse = ";")), by = list(Protein)]
# prot <- as.data.frame(tmp1)
# ... done! With their temporary peptide IDs.
prot$.pep.ids <- lapply(strsplit(prot$pep.ids, ";"), as.integer)
# Sort by N of peptides
prot$"Peptides count" <- vapply(prot$.pep.ids, length, 1)
prot <- prot[order(prot$"Peptides count", decreasing = TRUE),]
# Temporary protein IDs
prot$Prot.id <- as.character(c(1:nrow(prot)) + proffset)
if (CustPG) {
  prot <- rbind(prot2, prot)
  stopifnot(length(unique(prot$Prot.id)) == nrow(prot))
}
rownames(prot) <- stringi::stri_join("Pr_", prot$Prot.id)
# assign temp prot IDs to peptides
tmp1 <- seq$.Proteins
tmp1 <- proteoCraft::listMelt(tmp1, ColNames = c("Protein", "row"))
tmp1$Prot.id <- prot$Prot.id[match(tmp1$Protein, prot$Protein)]
tmp1 <- data.table::data.table(ProtID = tmp1$Prot.id, row = tmp1$row)
tmp1 <- tmp1[, list(.Prot.ids = list(ProtID)), by = list(row = row)]
tmp1 <- as.data.frame(tmp1)
tmp1 <- tmp1[order(tmp1$row),]
rws <- 1:nrow(seq)
seq$.Prot.ids <- lapply(rws, c) 
w <- which(rws %in% tmp1$row)
seq$.Prot.ids <- tmp1$.Prot.ids
# For each protein P, find all other proteins Q it shares peptides with...
#a1 <- Sys.time()
tmp1 <- prot$.pep.ids
tmp2 <- seq[, c("id", ".Prot.ids")]
saveRDS(tmp1, paste0(wd, "/tmp1.RDS"))
saveRDS(tmp2, paste0(wd, "/tmp2.RDS"))
parallel::clusterExport(cl, "wd", envir = environment())
invisible(parallel::clusterCall(cl, function() {
  tmp1 <<- readRDS(paste0(wd, "/tmp1.RDS"))
  tmp2 <<- readRDS(paste0(wd, "/tmp2.RDS"))
  return()
}))
unlink(paste0(wd, "/tmp1.RDS"))
unlink(paste0(wd, "/tmp2.RDS"))
f0 <- function(x) { sort(unique(unlist(tmp2$.Prot.ids[match(unlist(x), tmp2$id)]))) }
#environment(f0) <- .GlobalEnv
prot$pep_to_P <- parallel::parSapply(cl, tmp1, f0)
# Alternative slower rewrite
#a2 <- Sys.time()
# b1 <- Sys.time()
# tmp1 <- prot$.pep.ids
# tmp1 <- proteoCraft::listMelt(tmp1, ColNames = c("pep.id", "row"))
# tmp1$.Prot.ids <- seq$.Prot.ids[match(tmp1$pep.id, seq$id)]
# tmp1 <- proteoCraft::listMelt(tmp1$.Prot.ids, tmp1$row, c("Prot.id", "row"))
# tmp1 <- data.table::data.table(Prot.id = tmp1$Prot.id, row = tmp1$row)
# tmp1 <- tmp1[, list(.Prot.ids = sort(unique((Prot.id)))), by = list(row = row)]
# tmp1 <- as.data.frame(tmp1)
# tmp1 <- tmp1[order(tmp1$row),]
# rws <- 1:nrow(prot)
# prot$pep_to_P <- "" 
# w <- which(rws %in% tmp1$row)
# prot$pep_to_P[w] <- tmp1$.Prot.ids
# b2 <- Sys.time()
# 
#... then get those proteins' peptides
# Skipping, this takes too long
# tmp1 <- prot[, c("pep_to_P", ".pep.ids", "Prot.id")]
# parallel::clusterExport(cl, "tmp1", envir = environment())
# f0 <- function(x) { lapply(unlist(x), function(y) { tmp1$.pep.ids[match(y, tmp1$Prot.id)] }) }
# #environment(f0) <- .GlobalEnv
# prot$pep_to_P_to_pep <- parallel::parSapply(cl, tmp1$pep_to_P, f0)
# Identify minimal set of proteins required to explain the identified peptides:
cat(" - Identifying minimal set of proteins required to explain observed peptides.\n")
## Checking which proteins are contained within another?
## A protein is Leading (at this stage) if it cannot be subsumed into another:
cat(" - Identifying, for each protein ID, whether its peptides are contained by a single other protein ID.\n")
tmp1 <- prot[, c(".pep.ids", "pep_to_P", "Prot.id")]
saveRDS(tmp1, paste0(wd, "/tmp1.RDS"))
invisible(parallel::clusterCall(cl, function() {
  tmp1 <<- readRDS(paste0(wd, "/tmp1.RDS"))
  return()
}))
unlink(paste0(wd, "/tmp1.RDS"))
f0 <- function(x) {
  #x <- tmp1[1,]
  pep.ids <- unlist(x[[1]]) #This protein's peptides
  pep_to_P <- unlist(x[[2]])
  pep_to_P_to_pep <- tmp1$.pep.ids[match(pep_to_P, tmp1$Prot.id)] # list, all peptides from other proteins this protein shares some peptides with
  #Prot.id <- Prot.id[[1]] #(for testing)
  x3l <- vapply(pep_to_P_to_pep, length, 1)
  test <- vapply(pep_to_P_to_pep, function(y) { sum(y %in% pep.ids) }, 1) == x3l # Can the 2nd protein be subsumed in the 1st?
  res <- pep_to_P[which(test)]
  return(res)
}
#environment(f0) <- .GlobalEnv
protTemp <- parallel::parApply(cl, tmp1[, c(".pep.ids", "pep_to_P")], 1, f0)
#
## For each protein ID, which proteins does it contain?
cat(" - Identifying, for each protein ID, which other IDs are subsumable within it.\n")
tmp1 <- prot[, "Prot.id", drop = FALSE]
tmp1$temp <- protTemp
saveRDS(tmp1, paste0(wd, "/tmp1.RDS"))
invisible(parallel::clusterCall(cl, function() {
  tmp1 <<- readRDS(paste0(wd, "/tmp1.RDS"))
  return()
}))
unlink(paste0(wd, "/tmp1.RDS"))
f0 <- function(x) {
  x2 <- unlist(x[[2]])
  return(x2[which(x2 != x[[1]])])
}
#environment(f0) <- .GlobalEnv
prot$Contains <- parallel::parApply(cl, tmp1[, c("Prot.id", "temp")], 1, f0)
## Contained proteins and their temp peptide ids:
contained <- data.frame(Contained = sort(as.integer(unique(unlist(prot$Contains)))))
contained$.pep.ids <- prot$.pep.ids[match(contained$Contained, prot$Prot.id)]
## Assign a temporary Leading status...
### A protein may be contained by (an)other protein(s) with exactly the same peptides.
### All non-contained proteins get Leading status.
cat(" - Identifying temporary Leading protein IDs - i.e. those which are not subsumable into a single other protein ID.\n")
prot$Is.Leading <- !prot$Prot.id %in% contained$Contained
## Identify containers of each contained protein:
cat(" - Identifying Leading protein ID containing each subsumable protein ID.\n")
w <- which(vapply(prot$Contains, length, 1) > 0)
a <- setNames(prot$Contains[w], prot$Prot.id[w])
temp <- proteoCraft::listMelt(a)
temp <- data.table::data.table(value = temp$value, L1 = temp$L1)
temp <- temp[, list(byWhom = list(L1)), by = list(Contained = value)]
temp <- as.data.frame(temp)
contained$.byWhom <- temp$byWhom[match(contained$Contained, temp$Contained)]
## Container peptide IDs:
tmp3 <- contained$.byWhom
tmp4 <- prot[, c(".pep.ids", "Prot.id")]
saveRDS(tmp3, paste0(wd, "/tmp3.RDS"))
saveRDS(tmp4, paste0(wd, "/tmp4.RDS"))
invisible(parallel::clusterCall(cl, function() {
  tmp3 <<- readRDS(paste0(wd, "/tmp3.RDS"))
  tmp4 <<- readRDS(paste0(wd, "/tmp4.RDS"))
  return()
}))
unlink(paste0(wd, "/tmp3.RDS"))
unlink(paste0(wd, "/tmp4.RDS"))
f0 <- function(x) { tmp4$.pep.ids[match(unlist(x), tmp4$Prot.id)] }
#environment(f0) <- .GlobalEnv
contained$Container.pep.ids <- parallel::parLapply(cl, tmp3, f0)
## Now these containers, are they "larger" (N of peptides)?
tmp1 <- contained[, c(".pep.ids", "Container.pep.ids")]
saveRDS(tmp1, paste0(wd, "/tmp1.RDS"))
invisible(parallel::clusterCall(cl, function() {
  tmp1 <<- readRDS(paste0(wd, "/tmp1.RDS"))
  return()
}))
unlink(paste0(wd, "/tmp1.RDS"))
f0 <- function(x) {
  x1 <- unlist(x[[1]])
  x1l <- length(x1)
  x2 <- x[[2]]
  x2l <- vapply(x2, length, 1)
  res <- x2l > x1l
  return(res)
}
#environment(f0) <- .GlobalEnv
contained$Container.is.larger <- parallel::parApply(cl, tmp1, 1, f0)
#
tmp1 <- contained$Container.is.larger
saveRDS(tmp1, paste0(wd, "/tmp1.RDS"))
invisible(parallel::clusterCall(cl, function() {
  tmp1 <<- readRDS(paste0(wd, "/tmp1.RDS"))
  return()
}))
unlink(paste0(wd, "/tmp1.RDS"))
f0 <- function(x) { sum(x) }
#environment(f0) <- .GlobalEnv
contained$test <- parallel::parSapply(cl, tmp1, f0)
## Contained proteins get (temporary) Leading status here if they are equivalent to all of their containers,
## i.e. if they have no container larger than them:
contained$Is.Leading <- contained$test == 0
prot$Is.Leading[match(contained$Contained[which(contained$Is.Leading)], prot$Prot.id)] <- TRUE
## Create protein groups:
cat(" - Creating temporary protein groups.\n")
pg <- prot
if (CustPG) { pg <- pg[which(!pg$Prot.id %in% prot2$Prot.id),] }
pg <- pg[which(pg$Is.Leading), c("Protein", "Prot.id", "pep.ids")]
pg <- magrittr::set_colnames(aggregate(pg[,c("Protein", "Prot.id")],
                                       list(pg$pep.ids),
                                       paste,
                                       collapse = ";"),
                             c("Peptide IDs", "Leading protein IDs", "lead.protein.ids"))
pg$temp.pg.id <- as.character(1:nrow(pg) + pgoffset)
pg$.pep.ids <- lapply(strsplit(pg$"Peptide IDs", ";"), as.integer)
pg$"Peptides count" <- vapply(pg$.pep.ids, length, 1)
pg$.Leading.protein.IDs <- strsplit(pg$"Leading protein IDs", ";")
pg$.lead.protein.ids <- strsplit(pg$lead.protein.ids, ";")
#temp <- prot[which(prot$Is.Leading),]
tmp1 <- pg$.lead.protein.ids
tmp2 <- prot[, c("Contains", "Prot.id")]
saveRDS(tmp1, paste0(wd, "/tmp1.RDS"))
saveRDS(tmp2, paste0(wd, "/tmp2.RDS"))
invisible(parallel::clusterCall(cl, function() {
  tmp1 <<- readRDS(paste0(wd, "/tmp1.RDS"))
  tmp2 <<- readRDS(paste0(wd, "/tmp2.RDS"))
  return()
}))
unlink(paste0(wd, "/tmp1.RDS"))
unlink(paste0(wd, "/tmp2.RDS"))
f0 <- function(x) {
  x <- unlist(x)
  res <- unlist(tmp2$Contains[match(x, tmp2$Prot.id)])
  res <- res[which(!res %in% x)]
  return(res)
}
#environment(f0) <- .GlobalEnv
pg$"Also contains" <- parallel::parSapply(cl, tmp1, f0)
#
tmp1 <- pg[,c(".lead.protein.ids", "Also contains")]
tmp2 <- prot[, c("Protein", "Prot.id")]
saveRDS(tmp1, paste0(wd, "/tmp1.RDS"))
saveRDS(tmp2, paste0(wd, "/tmp2.RDS"))
invisible(parallel::clusterCall(cl, function() {
  tmp1 <<- readRDS(paste0(wd, "/tmp1.RDS"))
  tmp2 <<- readRDS(paste0(wd, "/tmp2.RDS"))
  return()
}))
unlink(paste0(wd, "/tmp1.RDS"))
unlink(paste0(wd, "/tmp2.RDS"))
f0 <- function(x) { list(tmp2$Protein[match(unlist(x), tmp2$Prot.id)]) }
#environment(f0) <- .GlobalEnv
pg$.Protein.IDs <- parallel::parApply(cl, tmp1, 1, f0)
pg$"Protein IDs" <- vapply(pg$.Protein.IDs, function(x) { paste(unlist(x), collapse = ";") }, "")
if (CustPG) {
  Custom_PGs$"Also contains" <- lapply(Custom_PGs$.lead.protein.ids, function(x) {
    res <- unlist(prot$Contains[match(unlist(x), prot$Prot.id)])
    res <- res[which(!res %in% x)]
    return(res)
  })
  Custom_PGs$.Protein.IDs <- apply(Custom_PGs[,c(".lead.protein.ids", "Also contains")], 1, function(x) {
    list(prot$Protein[match(unlist(x), prot$Prot.id)])
  })
  Custom_PGs$"Protein IDs" <- vapply(Custom_PGs$.Protein.IDs, function(x) { paste(unlist(x), collapse = ";") }, "")
}
# Assign temporary PG IDs to seq
temp <- proteoCraft::listMelt(pg$.pep.ids, pg$temp.pg.id, c("pep.id", "temp.pg.id"))
temp <- data.table::data.table(pep.id = temp$pep.id, temp.pg.id = temp$temp.pg.id)
temp <- temp[, list(temp.pg.ids = paste(temp.pg.id, collapse = ";")), by = list(pep.id = pep.id)] 
temp <- as.data.frame(temp)
seq$temp.pg.ids <- temp$temp.pg.ids[match(seq$id, temp$pep.id)]
seq$.temp.pg.ids <- strsplit(seq$temp.pg.ids, ";")
cat(paste0(" ---> ", nrow(pg), " protein groups are required at this stage.\n"))
## Test whether there are further protein groups we can remove?
## In the previous steps, we identified proteins which are not contained within another single protein.
## These then became the founding stones of protein groups.
## However, we should also remove PGs which can be explained by a combination of other PGs.
cat(" - Identifying protein groups which can be subsumed into a combination of other protein groups.\n")
tmp1 <- 1:nrow(pg)
tmp2 <- pg$.pep.ids
saveRDS(tmp1, paste0(wd, "/tmp1.RDS"))
saveRDS(tmp2, paste0(wd, "/tmp2.RDS"))
invisible(parallel::clusterCall(cl, function() {
  tmp1 <<- readRDS(paste0(wd, "/tmp1.RDS"))
  tmp2 <<- readRDS(paste0(wd, "/tmp2.RDS"))
  return()
}))
unlink(paste0(wd, "/tmp1.RDS"))
unlink(paste0(wd, "/tmp2.RDS"))
f0 <- function(x) {
  x <- unlist(x)
  a1 <- tmp1[which(tmp1 != x)]
  x1 <- unlist(tmp2[x]) # pep ids for that PG
  y <- unique(unlist(tmp2[a1])) # pep ids for all other PGs
  sum(x1 %in% y) == length(x1) # Are all peptides for that PG contained in the others?
}
#environment(f0) <- .GlobalEnv
pg$Removable <- parallel::parSapply(cl, tmp1, f0)
# Calculate individual and global PEPs - useful for later
cat(" - Calculating individual and global protein group PEPs (Posterior Error Probabilities).\n")
temp <- proteoCraft::listMelt(strsplit(seq$Proteins, ";"), seq$PEP)
temp <- data.table::data.table(value = temp$value, L1 = temp$L1)
temp <- temp[, list(x = prod(L1)), by = list(Group.1 = value)] 
temp <- as.data.frame(temp)
prot$PEP <- NA
w <- which(prot$Protein %in% temp$Group.1)
prot$PEP[w] <- temp$x[match(prot$Protein[w], temp$Group.1)]
if (CustPG) {
  temp <- proteoCraft::listMelt(strsplit(seq2$Proteins, ";"), seq2$PEP)
  temp <- data.table::data.table(value = temp$value, L1 = temp$L1)
  temp <- temp[, list(x = prod(L1)), by = list(Group.1 = value)] 
  temp <- as.data.frame(temp)
  w <- which(prot$Protein %in% temp$Group.1)
  prot$PEP[w] <- temp$x[match(prot$Protein[w], temp$Group.1)]
}
#
tmp1 <- pg$.Protein.IDs
tmp2 <- prot[, c("PEP", "Protein")]
saveRDS(tmp1, paste0(wd, "/tmp1.RDS"))
saveRDS(tmp2, paste0(wd, "/tmp2.RDS"))
invisible(parallel::clusterCall(cl, function() {
  tmp1 <<- readRDS(paste0(wd, "/tmp1.RDS"))
  tmp2 <<- readRDS(paste0(wd, "/tmp2.RDS"))
  return()
}))
unlink(paste0(wd, "/tmp1.RDS"))
unlink(paste0(wd, "/tmp2.RDS"))
f0 <- function(x) { tmp2$PEP[match(unlist(x), tmp2$Protein)] }
#environment(f0) <- .GlobalEnv
pg$.PEPs <- parallel::parSapply(cl, tmp1, f0)
#
tmp1 <- pg$.PEPs
saveRDS(tmp1, paste0(wd, "/tmp1.RDS"))
invisible(parallel::clusterCall(cl, function() {
  tmp1 <<- readRDS(paste0(wd, "/tmp1.RDS"))
  return()
}))
unlink(paste0(wd, "/tmp1.RDS"))
f0 <- function(x) { paste(x, collapse = ";") }
#environment(f0) <- .GlobalEnv
pg$PEPs <- parallel::parSapply(cl, tmp1, f0)
f0 <- function(x) { min(x) }
#environment(f0) <- .GlobalEnv
pg$PEP <- parallel::parSapply(cl, tmp1, f0)
if (CustPG) {
  Custom_PGs$.Protein.IDs
  Custom_PGs$.PEPs <- lapply(Custom_PGs$.Protein.IDs, function(x) { prot$PEP[match(x, prot$Protein)] })
  Custom_PGs$PEPs <- vapply(Custom_PGs$.PEPs, paste, "", collapse = ";")
  Custom_PGs$PEP <- vapply(Custom_PGs$.PEPs, min, 1)
}
#Removable <- aggregate(pg$Removable, list(pg$Removable), length)
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
  pt <- seq[which(seq$id %in% unlist(pg$.pep.ids[W])),]
  pt <- pt[order(vapply(pt$.temp.pg.ids, length, 1), decreasing = TRUE),]
  # For now none of these peptides have been dealt with:
  pt$Done <- FALSE
  while (sum(!pt$Done)) {
    i <- which(!pt$Done)[1]
    cat(paste0(" ---> Still ", length(unique(unlist(pt$.temp.pg.ids[which(!pt$Done)]))), " unresolved protein groups...\n"))
    #cat(paste0("   i = ", i, "\n"))
    # Get all protein groups for the current peptide:
    ptpg <- pg[which(pg$temp.pg.id %in% unlist(pt$.temp.pg.ids[i])),]
    # Peptides from those PGs:
    p <- unlist(ptpg$.pep.ids)
    # overlappers = PGs containing at least one "p" peptide:
    tst0 <- vapply(pg$.pep.ids, function(x) { sum(unlist(x) %in% p) }, 1) > 0
    overlappers <- which(tst0) # overlappers and all derived variables (overlappers2, goners, keepers, removables) are row indices, not protein group IDs
    # All of the peptides of overlapping PGs:
    P <- unique(unlist(pg$.pep.ids[overlappers])); L <- length(P)
    # Extend overlappers:
    tst2 <- vapply(pg$.pep.ids, function(x) { sum(unlist(x) %in% P) }, 1) > 0 # Proteins with at least one "P" peptide
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
      tst2 <- vapply(pg$.pep.ids, function(x) { sum(unlist(x) %in% P) }, 1) > 0
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
    if (n < 2) { warning("   Check code: I expected at least 3, not 2 or less, overlapping groups here.") } # Always good to check that code is not buggy.
    # Table of peptides and PGs
    peptopg <- set_rownames(magrittr::set_colnames(as.data.frame(sapply(overlappers, function(x) {
      P %in% unlist(pg$.pep.ids[x])
    })), overlappers), P)
    # Heuristic: we will sort protein groups by number of peptides then PEP and remove those which do not increase coverage
    peep <- setNames(pg$PEP[overlappers], overlappers) # PEP is used to tie-break protein groups of equal ranking
    peep[which(is.na(peep))] <- 1
    peptopg <- peptopg[, order(peep, decreasing = FALSE)]
    smcl <- setNames(colSums(peptopg), colnames(peptopg)) # Then the main ranking is applied: number of peptides
    peptopg <- peptopg[, order(smcl, decreasing = TRUE)]
    #
    aim <- nrow(peptopg)
    incl <- setNames(rep(FALSE, ncol(peptopg)), colnames(peptopg))
    incl[1] <- TRUE
    covpep <- peptopg[,1]
    curr <- sum(covpep)
    # First filter
    # Heuristic? we take the straighter path to complete peptides coverage
    # There may be cases where the final number of protein groups is not the smallest? Not sure.
    while (curr < aim) {
      wnc <- which(!covpep)
      wni <- names(incl)[which(!incl)]
      tst <- apply(peptopg[wnc, wni, drop = FALSE], 2, sum)
      bst <- wni[which(tst == max(tst))]
      if (length(bst) > 1) { bst <- bst[which(smcl[bst] == max(smcl[bst]))[1]] }
      covpep[wnc] <- peptopg[wnc, bst]
      incl[bst] <- TRUE
      curr <- sum(covpep)
    }
    goners <- as.integer(names(incl)[which(!incl)])
    keepers <- as.integer(names(incl)[which(incl)])
    removables <- keepers[which(pg$Removable[keepers])] # self explanatory
    P2 <- unique(unlist(pg$.pep.ids[keepers])); L2 <- length(P2)
    stopifnot(L == L2)
    # Second filter: trying to remove some more, one by one
    tst1a <- vapply(pg$.pep.ids[removables], length, 1) # Their N of peptides
    tst1b <- pg$PEP[removables] # Their PEPs
    tst1b[which(is.na(tst1b))] <- 1
    go <- length(removables) > 0
    while (go) {
      tst2 <- vapply(removables, function(x) {
        length(unique(unlist(pg$.pep.ids[keepers[which(keepers != x)]])))
      }, 1) == L # Can any of these be dropped without affecting the number of peptides?
      witch <- which(tst2)
      if (length(witch)) {
        witch <- witch[which(tst1a[witch] == max(tst1a[witch]))]
        if (length(witch) > 1) { witch <- witch[which(tst1b[witch] == min(tst1b[witch]))][1] }
        wutch <- which(!c(1:length(removables)) %in% witch)
        goners <- c(goners, removables[witch])
        keepers <- keepers[which(!keepers %in% goners)]
        tst1a <- tst1a[wutch]
        tst1b <- tst1b[wutch]
        removables <- removables[wutch]
        if (!length(removables)) { go <- FALSE }
      } else { go <- FALSE }
    }
    P2 <- unique(unlist(pg$.pep.ids[keepers])); L2 <- length(P2)
    stopifnot(L == L2)
    cat(paste0(" ---> Removing redundant protein group", c("", "s")[(length(goners)>1)+1],
               " #", paste(sort(goners), collapse = "-"), "\n"))
    #
    # Kept for reference: alternative, logically correct (notwithstanding unknown potential bugs),
    # but takes centuries to complete for some protein groups.
    # Now how many do we want to remove?
    # We know how many we want to keep at least: all those not flagged pg$Removable == TRUE
    #core <- overlappers[which(!pg$Removable[overlappers])]; lc <- length(core)
    #removables <- overlappers[which(pg$Removable[overlappers])]; lr <- length(removables)
    # We want to try all ways to choose between 0 and lr removables, going up, such that we get all peptides.
    #cp <- unique(unlist(pg$.pep.ids[core]))
    #tst <- which(!P %in% cp)
    #if (length(tst) > 0) {
    #  done <- FALSE
    #  r <- 0
    #  while (!done) {
    #    r <- r+1
    #    Com <- gtools::combinations(n = lr, r = r, v = removables)
    #    tst <- apply(Com, 1, function(x) { sum(!P %in% unique(c(cp, unlist(pg$.pep.ids[unlist(x)])))) }) == 0
    #    if (length(which(tst)) > 0) { done <- TRUE }
    #  }
    #  if (r == lr) { stop("Something went wrong!") }
    #  keepers <- sort(c(core, unique(as.integer(Com[which(tst),]))))
    #} else { keepers <- core }
    #goners <- overlappers[which(!overlappers %in% keepers)]
    pg$Remove[goners] <- TRUE
    pt$Done[which(pt$id %in% unlist(pg$.pep.ids[overlappers]))] <- TRUE
  }
}
cat("   All protein groups successfully resolved!\n")
## Remove those protein groups we don't need anymore:
pg1 <- pg[which(!pg$Remove),]
p1 <- unique(unlist(pg1$.pep.ids))
p <- unique(unlist(pg$.pep.ids))
if (sum(!p %in% p1)) {
  # (Sanity check)
  if (cleanUp) { stopCluster(cl) }
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
  pg$.pep.ids <- lapply(pg$.pep.ids, function(x) { x[which(!x %in% custpgpep)] })
  pg$"Priority level" <- Inf
  pg <- rbind(Custom_PGs, pg)
  seq <- rbind(seq2, seq)
  seq <- seq[match(Pep$"Modified sequence", seq$"Modified sequence"),]
}
if ((!is.null(ContCol))&&(ContCol %in% colnames(DB))) {
  tmp <- proteoCraft::listMelt(strsplit(pg$`Protein IDs`, ";"), pg$temp.pg.id)
  tmp2 <- DB$`Protein ID`[which(DB[[ContCol]] == "+")]
  tmp <- tmp[which(gsub("^CON__", "", tmp$value) %in% gsub("^CON__", "", tmp2)),]
  pg$"Potential contaminant" <- c("", "+")[(pg$temp.pg.id %in% tmp$L1)+1]
}
cat(paste0("   Final number of protein groups: ", nrow(pg), "\n"))
# Protein group IDs
pg <- pg[order(pg$"Peptides count", decreasing = TRUE),]
if ((!is.null(ContCol))&&(ContCol %in% colnames(DB))) {
  pg <- pg[order(pg$"Potential contaminant", decreasing = FALSE),]
}
pg$id <- 1:nrow(pg)
# Peptide IDs
if (grepl("peptide", Peptide.IDs, ignore.case = TRUE)) { pepcolnm <- gsub("ss$", "s", paste0(Peptide.IDs, "s")) } else {
  pepcolnm <- "Peptide IDs"
}
pg[[pepcolnm]] <- vapply(pg$.pep.ids, paste, "", collapse = ";")
# Assign final PG IDs to seq
tmp1 <- seq$.temp.pg.ids
tmp2 <- pg[, c("id", "temp.pg.id")]
saveRDS(tmp1, paste0(wd, "/tmp1.RDS"))
saveRDS(tmp2, paste0(wd, "/tmp2.RDS"))
invisible(parallel::clusterCall(cl, function() {
  tmp1 <<- readRDS(paste0(wd, "/tmp1.RDS"))
  tmp2 <<- readRDS(paste0(wd, "/tmp2.RDS"))
  return()
}))
unlink(paste0(wd, "/tmp1.RDS"))
unlink(paste0(wd, "/tmp2.RDS"))
f0 <- function(x) { tmp2$id[which(tmp2$temp.pg.id %in% unlist(x))] }
#environment(f0) <- .GlobalEnv
tmp1 <- seq$.PG.ids <- parallel::parSapply(cl, tmp1, f0)
saveRDS(tmp1, paste0(wd, "/tmp1.RDS"))
invisible(parallel::clusterCall(cl, function() {
  tmp1 <<- readRDS(paste0(wd, "/tmp1.RDS"))
  return()
}))
unlink(paste0(wd, "/tmp1.RDS"))
f0 <- function(x) { paste(x, collapse = ";") }
#environment(f0) <- .GlobalEnv
seq$"Protein group IDs" <- parallel::parSapply(cl, tmp1, f0)
pg$temp.pg.id <- NULL
# Unique peptides: peptides unique to a protein group
cat(" - Identifying unique peptides.\n")
temp <- seq[which(!grepl(";", seq$"Protein group IDs")), c("id", "Protein group IDs")]
temp <- magrittr::set_colnames(aggregate(temp$id, list(temp$"Protein group IDs"), function(x) {
  x
}), c("id", "pep.ids"))
pg$.unique.pep.ids <- list(NULL)
w <- which(pg$id %in% temp$id)
pg$.unique.pep.ids[w] <- temp$pep.ids[match(pg$id[w], temp$id)]
pg$"Unique peptide IDs" <- vapply(pg$.unique.pep.ids, paste, "", collapse = ";")
pg$"Unique peptides" <- vapply(pg$.unique.pep.ids, length, 1)
# Apply Occam's razor:
# If a peptide is shared, it should go to the protein group with the most identified peptides in total.
# (Here, "peptide is razor" means unique OR razor.)
# In cases of ties, we will go for lowest PEP.
# Note that these ties would mean that some few protein groups which do get reported may not have any razor peptides.
cat(" - Applying Occam's razor,\n     i.e. assigning \"razor\" status to shared peptides for the parent group with the highest overall peptides count,\n     with protein priority (if available) and PEP as tie-breakers (in that order).")
seq$.Protein.group.IDs <- strsplit(seq$"Protein group IDs", ";")
## Note here: the code below allows for multiple razor PGs.
# This can happen if all have the same priority, number of peptides and PEPs (can happen with NA PEP values).
tmp <- data.frame(PG_IDs = unique(seq$"Protein group IDs"))
tmp$.PG_IDs <- strsplit(tmp$PG_IDs, ";")
kol <- c("id", "Peptides count", "PEP", "Priority level")
kol <- kol[which(kol %in% colnames(pg))]
tmp1 <- tmp$.PG_IDs
tmp2 <- pg[, kol]
saveRDS(tmp1, paste0(wd, "/tmp1.RDS"))
saveRDS(tmp2, paste0(wd, "/tmp2.RDS"))
saveRDS(CustPG, paste0(wd, "/CustPG.RDS"))
invisible(parallel::clusterCall(cl, function() {
  tmp1 <<- readRDS(paste0(wd, "/tmp1.RDS"))
  tmp2 <<- readRDS(paste0(wd, "/tmp2.RDS"))
  CustPG <<- readRDS(paste0(wd, "/CustPG.RDS"))
  return()
}))
unlink(paste0(wd, "/tmp1.RDS"))
unlink(paste0(wd, "/tmp2.RDS"))
unlink(paste0(wd, "/CustPG.RDS"))
f0 <- function(x) {
  #sapply(tmp1, function(x) {
  #x <- unlist(tmp1[1])
  x <- unlist(x)
  w1 <- match(x, tmp2$id) # ~ x
  if (CustPG) {
    p <- tmp2$`Priority level`[w1] # ~ x
    mp <- min(p) # ~ x
    w2 <- which(tmp2$`Priority level`[w1] == mp) # ~ x[w2]  
  } else { w2 <- 1:length(w1) }
  k <- tmp2$"Peptides count"[w1[w2]]
  Mk <- max(k) # ~ x[w2]
  w3 <- which(k == Mk) # ~ x[w2][w3]
  if (length(w3) > 1) {
    z <- tmp2$PEP[w1[w2]] # ~ x[w2]
    z[which(is.na(z))] <- 1
    z[-w3] <- 1
    mP <- min(z, na.rm = TRUE)
    w3 <- which((k == Mk)&(z == mP)) # ~ x[w2][w3]
  }
  return(x[w2][w3])
}
#environment(f0) <- .GlobalEnv
tmp$.Razor_PG_IDs <- parallel::parSapply(cl, tmp1, f0)
tmp$Razor_PG_IDs <- vapply(tmp$.Razor_PG_IDs, paste, "", collapse = ";")
seq[, c(".Razor.protein.group.ID", "Razor protein group ID")] <- tmp[match(seq$"Protein group IDs", tmp$PG_IDs), c(".Razor_PG_IDs", "Razor_PG_IDs")]
c1 <- c("Leading proteins", "Leading razor proteins")
c2 <- c(".Protein.group.IDs", ".Razor.protein.group.ID")
tmp1 <- pg[, c("id", ".Leading.protein.IDs")]
tmp2 <- prot[, c("PEP", "Protein")]
tmp3 <- seq[, c2]
saveRDS(tmp1, paste0(wd, "/tmp1.RDS"))
saveRDS(tmp2, paste0(wd, "/tmp2.RDS"))
saveRDS(tmp3, paste0(wd, "/tmp3.RDS"))
saveRDS(c2, paste0(wd, "/c2.RDS"))
invisible(parallel::clusterCall(cl, function() {
  tmp1 <<- readRDS(paste0(wd, "/tmp1.RDS"))
  tmp2 <<- readRDS(paste0(wd, "/tmp2.RDS"))
  tmp3 <<- readRDS(paste0(wd, "/tmp3.RDS"))
  c2 <<- readRDS(paste0(wd, "/c2.RDS"))
  return()
}))
unlink(paste0(wd, "/tmp1.RDS"))
unlink(paste0(wd, "/tmp2.RDS"))
unlink(paste0(wd, "/tmp3.RDS"))
unlink(paste0(wd, "/c2.RDS"))
f0 <- function(x) {
  x <- unlist(tmp1$.Leading.protein.IDs[match(unlist(x), tmp1$id)])
  y <- tmp2$PEP[match(x, tmp2$Protein)]
  x <- x[order(y, decreasing = FALSE)]
  return(paste(x, collapse = ";"))
}
#environment(f0) <- .GlobalEnv
for (i in 1:length(c1)) {
  seq[[c1[i]]] <- parallel::parSapply(cl, tmp3[[c2[i]]], f0)
}
# Sort by PEP
if (CustPG) { # Priority-based razor peptides decision for custom protein groups
  for (i in priorities) { #i <- 1
    wcp <- which(pg$"Priority level" == i)
    p <- unique(unlist(pg$.pep.ids[wcp]))
    whp <- which(pg$"Priority level" < i)
    if (length(whp)) { p <- p[which(!p %in% unique(unlist(pg$.pep.ids[whp])))] } # Remove peptides matching protein groups of higher priority
    m <- match(p, seq$id)
    seq$.Razor.protein.group.ID[m] <- lapply(seq$.Protein.group.IDs[m], function(x) { x[which(x %in% pg$id[wcp])] })
  }
}
# Update columns
seq$"Razor protein group ID" <- vapply(seq$.Razor.protein.group.ID, paste, "", collapse = ";")
seq$"Protein group ID" <- vapply(seq$.Protein.group.IDs, paste, "", collapse = ";")
#
seq$"Leading razor proteins" <- vapply(strsplit(seq$"Leading razor proteins", ";"), function(x) { unlist(x)[1] }, "")
if (length(pg$.lead.protein.ids) != length(unique(pg$.lead.protein.ids))) {
  stop("A \"leading\" protein cannot be in several protein groups!")
}
# Assign razor peptides to protein groups
temp <- proteoCraft::listMelt(as.list(seq$.Razor.protein.group.ID), seq$id)
temp <- aggregate(temp$L1, list(temp$value), function(x) { paste(sort(as.integer(unique(x))), collapse = ";") })
pg$"Razor peptide IDs" <- temp$x[match(pg$id, as.integer(temp$Group.1))]
pg$.razor.pep.ids <- lapply(strsplit(pg$`Razor peptide IDs`, ";"), as.integer)
pg$"Razor + unique peptides" <- vapply(pg$.razor.pep.ids, length, 1)
pg$.Peptide.is.razor <- apply(pg[, c(".pep.ids", ".razor.pep.ids")], 1, function(x) { #x <- pg[1, c(".pep.ids", ".razor.pep.ids")]
  c("False", "True")[(x[[1]] %in% x[[2]])+1]
})
# 
pg$"Peptide is razor" <- vapply(pg$.Peptide.is.razor, paste, "", collapse = ";")
c1 <- paste0("Peptide counts (", c("all", "unique", "razor+unique"), ")")
c2 <- c(".pep.ids", ".unique.pep.ids", ".razor.pep.ids", ".Protein.IDs")
tmp1 <- pg[, c2]
tmp2 <- prot[, c(".pep.ids", "Protein")]
saveRDS(tmp1, paste0(wd, "/tmp1.RDS"))
saveRDS(tmp2, paste0(wd, "/tmp2.RDS"))
saveRDS(c2, paste0(wd, "/c2.RDS"))
invisible(parallel::clusterCall(cl, function() {
  tmp1 <<- readRDS(paste0(wd, "/tmp1.RDS"))
  tmp2 <<- readRDS(paste0(wd, "/tmp2.RDS"))
  c2 <<- readRDS(paste0(wd, "/c2.RDS"))
  return()
}))
unlink(paste0(wd, "/tmp1.RDS"))
unlink(paste0(wd, "/tmp2.RDS"))
unlink(paste0(wd, "/c2.RDS"))
f0 <- function(x) {
  x1 <- unlist(x[[1]]); x2 <- unlist(x[[2]]); x3 <- unlist(x[[3]]); x4 <- unlist(x[[4]])
  x4 <- tmp2$.pep.ids[match(x4, tmp2$Protein)]
  x1 <- paste(vapply(x4, function(y) { sum(y %in% x1) }, 1), collapse = ";")
  x2 <- paste(vapply(x4, function(y) { sum(y %in% x2) }, 1), collapse = ";")
  x3 <- paste(vapply(x4, function(y) { sum(y %in% x3) }, 1), collapse = ";")
  return(c(x1, x2, x3))
}
#environment(f0) <- .GlobalEnv
pg[, c1] <- as.data.frame(t(parallel::parApply(cl, tmp1[, c2], 1, f0)))
test1 <- as.integer(unlist(strsplit(pg$"Peptide counts (all)", ";")))
test2 <- as.integer(unlist(strsplit(pg$"Peptide counts (unique)", ";")))
test3 <- as.integer(unlist(strsplit(pg$"Peptide counts (razor+unique)", ";")))
temp <- proteoCraft::listMelt(strsplit(pg$`Protein IDs`, ";"), pg$id, c("Protein ID", "Protein Group ID"))
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
if (Npep > 1) {
  pg[[paste0("Quality filter: min ", Npep, " razor&unique pep.")]] <- c("", "Keep")[
    vapply(strsplit(pg$`Peptide counts (razor+unique)`, ";"), function(x) {
      as.integer(x[[1]]) >= Npep
    }, TRUE) + 1]
}
#
cat(" - Getting sequence coverage.\n")
c1 <- c("Sequence coverage [%]", "Unique + razor sequence coverage [%]", "Unique sequence coverage [%]")
c2 <- c(".Leading.protein.IDs", ".pep.ids", ".razor.pep.ids", ".unique.pep.ids")
tmp1 <- DB[, c("Protein ID", "Sequence")]
tmp1$`Protein ID` <- gsub("^CON__", "", tmp1$`Protein ID`)
tmp2 <- pg[, c2]
tmp3 <- seq[, c("id", "Sequence")]
saveRDS(tmp1, paste0(wd, "/tmp1.RDS"))
saveRDS(tmp2, paste0(wd, "/tmp2.RDS"))
saveRDS(tmp3, paste0(wd, "/tmp3.RDS"))
saveRDS(c2, paste0(wd, "/c2.RDS"))
invisible(parallel::clusterCall(cl, function() {
  tmp1 <<- readRDS(paste0(wd, "/tmp1.RDS"))
  tmp2 <<- readRDS(paste0(wd, "/tmp2.RDS"))
  tmp3 <<- readRDS(paste0(wd, "/tmp3.RDS"))
  c2 <<- readRDS(paste0(wd, "/c2.RDS"))
  return()
}))
unlink(paste0(wd, "/tmp1.RDS"))
unlink(paste0(wd, "/tmp2.RDS"))
unlink(paste0(wd, "/tmp3.RDS"))
unlink(paste0(wd, "/c2.RDS"))
f0 <- function(x) {
  #x <- tmp2[1, c2]
  s <- gsub("^CON__", "", unlist(x[[1]])[1])
  w <- match(s, tmp1$"Protein ID")
  ptids <- unlist(x[[2]])
  prids <- unlist(x[[3]])
  puids <- unlist(x[[4]])
  pt <- tmp3$Sequence[match(ptids, tmp3$id)]
  pr <- tmp3$Sequence[match(prids, tmp3$id)]
  pu <- tmp3$Sequence[match(puids, tmp3$id)]
  print(pt)
  print(pr)
  print(pu)
  res <- rep(NA, 3)
  if (length(w)) {
    s <- setNames(tmp1$Sequence[w], s)
    stopifnot(length(pt) > 0)
    res[1] <- round(100*proteoCraft::Coverage(proteins = s, peptides = pt), 1)
    if (length(pr)) { res[2] <- round(100*proteoCraft::Coverage(proteins = s, peptides = pr), 1) } 
    if (length(pu)) { res[3] <- round(100*proteoCraft::Coverage(proteins = s, peptides = pu), 1) } 
  } else { warning(paste0("Protein accession ", s, " was not found in the database!")) }
  #}
  return(res)
}
#environment(f0) <- .GlobalEnv
#pg[, c1] <- as.data.frame(t(parallel::parApply(cl, pg[, c2], 1, function(x) {
pg[, c1] <- t(parallel::parApply(cl, tmp2[, c2], 1, f0))
# Note:
#   The "Gene names" and "Protein names" columns will not match the MQ ones...
#   but after checking that seems to be because the MQ list is missing some!
cat(" - Getting Gene and Protein names, miscellaneous annotations...\n")
ca <- c("Common Names", "Common Name")
if (sum(ca %in% colnames(DB))) {
  ca <- ca[which(ca %in% colnames(DB))]
  test <- vapply(ca, function(x) { length(is.na(DB[[x]])) }, 1) > 0 
  ca <- ca[which(test)[1]]
}
cb <- c("Genes", "Gene", "Gene IDs", "Gene ID")
if (sum(cb %in% colnames(DB))) {
  cb <- cb[which(cb %in% colnames(DB))]
  test <- vapply(cb, function(x) { length(is.na(DB[[x]])) }, 1) > 0
  cb <- cb[which(test)[1]]
}
c1 <- c("Number of proteins", "Fasta headers")
c2 <- "temp"
if (length(ca) == 1) {
  c1 <- c(c1, "Protein names")
  c2 <- c(c2, ca)
}
if (length(cb) == 1) {
  c1 <- c(c1, "Gene names")
  c2 <- c(c2, cb)
}
kol <- c("Protein ID", "Header", ca, cb)
kol <- kol[which(kol %in% colnames(DB))]
temp <- DB[, kol]
temp$temp <- gsub("^>", "", temp$Header)
temp$"Protein ID" <- gsub("^CON__", "", temp$"Protein ID")
temp2b <- gsub("^CON__", "", gsub(";CON__", ";", pg$"Protein IDs"))
tmp2 <- strsplit(temp2b, ";")
tmp1 <- temp[, c("Protein ID", c2)]
saveRDS(tmp1, paste0(wd, "/tmp1.RDS"))
saveRDS(tmp2, paste0(wd, "/tmp2.RDS"))
saveRDS(c1, paste0(wd, "/c1.RDS"))
saveRDS(c2, paste0(wd, "/c2.RDS"))
saveRDS(ca, paste0(wd, "/ca.RDS"))
saveRDS(cb, paste0(wd, "/cb.RDS"))
invisible(parallel::clusterCall(cl, function() {
  tmp1 <<- readRDS(paste0(wd, "/tmp1.RDS"))
  tmp2 <<- readRDS(paste0(wd, "/tmp2.RDS"))
  c1 <<- readRDS(paste0(wd, "/c1.RDS"))
  c2 <<- readRDS(paste0(wd, "/c2.RDS"))
  ca <<- readRDS(paste0(wd, "/ca.RDS"))
  cb <<- readRDS(paste0(wd, "/cb.RDS"))
  return()
}))
unlink(paste0(wd, "/tmp1.RDS"))
unlink(paste0(wd, "/tmp2.RDS"))
unlink(paste0(wd, "/c1.RDS"))
unlink(paste0(wd, "/c2.RDS"))
unlink(paste0(wd, "/ca.RDS"))
unlink(paste0(wd, "/cb.RDS"))
f0 <- function(x) { #x <- tmp2[1]
  x <- unlist(x)
  w <- match(x, tmp1$"Protein ID")
  res <- tmp1[w, c2]
  if (length(w) > 1) {
    res2 <- paste(res[[c2[1]]], collapse = ";")
    if (length(ca) == 1) { res2 <- c(res2, paste(res[[ca]], collapse = ";")) }
    if (length(cb) == 1) {
      #y <- aggregate(res[[cb]], list(res[[cb]]), length)
      y <- data.table(a = res[[cb]])
      y <- y[, list(x = length(a)), by = list(Group.1 = a)]
      y <- as.data.frame(y)
      y <- y[order(y$x, decreasing = TRUE),]
      res2 <- c(res2, paste(y$Group.1, collapse = ";"))
    }
  } else { res2 <- unlist(res) }
  if (length(res2) != (length(c1)-1)) { stop(paste(x, collapse = ";")) }
  return(setNames(c(length(x), res2), c1))
}
#environment(f0) <- .GlobalEnv
tst2 <- parallel::parLapply(cl, tmp2, f0)
tst2 <- as.data.frame(do.call(rbind, tst2))
pg[, c1] <- tst2
#vapply(1:4, function(x) { sum(tst1[[x]] != tst2[[x]]) }, 1)
pg$"Number of proteins" <- as.integer(pg$"Number of proteins")
DB$temp <- NULL
DB$"Sequence length" <- nchar(DB$Sequence)
tmp1 <- DB[, c("Protein ID", "MW [kDa]", "Sequence length")]
tmp2 <- vapply(strsplit(pg$"Leading protein IDs", ";"), function(x) { x[1] }, "")
saveRDS(tmp1, paste0(wd, "/tmp1.RDS"))
saveRDS(tmp2, paste0(wd, "/tmp2.RDS"))
invisible(parallel::clusterCall(cl, function() {
  tmp1 <<- readRDS(paste0(wd, "/tmp1.RDS"))
  tmp2 <<- readRDS(paste0(wd, "/tmp2.RDS"))
  return()
}))
unlink(paste0(wd, "/tmp1.RDS"))
unlink(paste0(wd, "/tmp2.RDS"))
f0 <- function(x) { unlist(tmp1[match(x, tmp1$"Protein ID"), c("MW [kDa]", "Sequence length")]) }
#environment(f0) <- .GlobalEnv
pg[, c("Mol. weight [kDa]", "Sequence length")] <- as.data.frame(t(parallel::parSapply(cl, tmp2, f0)))
f0 <- function(x) { paste(tmp1$"Sequence length"[match(x, tmp1$"Protein ID")], collapse = ";") } # (I checked, the order is fine) 
#environment(f0) <- .GlobalEnv
pg$"Sequence lengths" <- parallel::parSapply(cl, tmp2, f0)
tmp1 <- pg$.pep.ids
tmp2 <- seq$id
tmp3 <- nchar(gsub(paste(c(AA, "_"), collapse = "|"), "", seq$"Modified sequence"))
saveRDS(tmp1, paste0(wd, "/tmp1.RDS"))
saveRDS(tmp2, paste0(wd, "/tmp2.RDS"))
saveRDS(tmp3, paste0(wd, "/tmp3.RDS"))
invisible(parallel::clusterCall(cl, function() {
  tmp1 <<- readRDS(paste0(wd, "/tmp1.RDS"))
  tmp2 <<- readRDS(paste0(wd, "/tmp2.RDS"))
  tmp3 <<- readRDS(paste0(wd, "/tmp3.RDS"))
  return()
}))
unlink(paste0(wd, "/tmp1.RDS"))
unlink(paste0(wd, "/tmp2.RDS"))
unlink(paste0(wd, "/tmp3.RDS"))
f0 <- function(x) { #x <- tmp1[1]
  x <- unlist(x)
  x <- tmp3[match(x, tmp2)]
  x <- c("", "+")[min(x > 0)+1]
  return(x)
}
#environment(f0) <- .GlobalEnv
temp <- parallel::parSapply(cl, tmp1, f0)
u <- unique(temp)
if (length(u) > 1) {
  pg$"Only identified by site" <- temp
} else {
  if (u == "+") { warning("Do we really expect all of our peptides to be modified?") }
  pg$"Only identified by site" <- temp
}
if (!misFun(Ev)) {
  seq$.Evidence.IDs <- strsplit(seq$"Evidence IDs", ";")
  tmp1 <- pg$.pep.ids
  tmp2 <- seq[, c("id", ".Evidence.IDs")]
  saveRDS(tmp1, paste0(wd, "/tmp1.RDS"))
  saveRDS(tmp2, paste0(wd, "/tmp2.RDS"))
  invisible(parallel::clusterCall(cl, function() {
    tmp1 <<- readRDS(paste0(wd, "/tmp1.RDS"))
    tmp2 <<- readRDS(paste0(wd, "/tmp2.RDS"))
    return()
  }))
  unlink(paste0(wd, "/tmp1.RDS"))
  unlink(paste0(wd, "/tmp2.RDS"))
  f0 <- function(x) {
    x <- unlist(x)
    paste(sort(as.integer(unique(unlist(tmp2$.Evidence.IDs[match(x, tmp2$id)])))),
          collapse = ";")
  }
  #environment(f0) <- .GlobalEnv
  pg$"Evidence IDs" <- parallel::parSapply(cl, tmp1, f0)
  if ("MS/MS count" %in% colnames(Ev)) {
    temp <- aggregate(Ev$"MS/MS count", list(Ev$"Modified sequence"), sum)
    seq$"MS/MS count" <- temp$x[match(seq$"Modified sequence", temp$Group.1)]
    temp <- proteoCraft::listMelt(seq$.Protein.group.IDs, seq$"MS/MS count", c("Protein group ID", "MS/MS count"))
    temp <- magrittr::set_colnames(aggregate(as.integer(temp$"MS/MS count"),
                                             list(temp$"Protein group ID"),
                                             sum),
                                   c("Protein group ID", "MS/MS count"))
    pg$"MS/MS count" <- temp$"MS/MS count"[match(pg$id, temp$"Protein group ID")]
    Ev$"Protein group IDs" <- seq$"Protein group IDs"[match(Ev$"Modified sequence", seq$"Modified sequence")]
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
temp <- seq[which(!grepl(";", seq$"Protein group IDs")), c("id", "Protein group IDs")]
temp <- aggregate(temp$id, list(temp$"Protein group IDs"), function(x) {
  paste(sort(x), collapse = ";")
})
w <- which(pg$id %in% temp$Group.1)
pg$"Unique peptide IDs"[w] <- temp$x[match(pg$id[w], temp$Group.1)]
# Some more annotations
if ("Common Name" %in% colnames(DB)) {
  tmp1 <- DB[, c("Protein ID", "Common Name")]
  tmp2 <- strsplit(pg$"Leading protein IDs", ";")
  saveRDS(tmp1, paste0(wd, "/tmp1.RDS"))
  saveRDS(tmp2, paste0(wd, "/tmp2.RDS"))
  invisible(parallel::clusterCall(cl, function() {
    tmp1 <<- readRDS(paste0(wd, "/tmp1.RDS"))
    tmp2 <<- readRDS(paste0(wd, "/tmp2.RDS"))
    return()
  }))
  unlink(paste0(wd, "/tmp1.RDS"))
  unlink(paste0(wd, "/tmp2.RDS"))      
  f0 <- function(x) {
    x <- tmp1$"Common Name"[match(unlist(x), tmp1$"Protein ID")]
    x <- x[which((!is.na(x))&(x != ""))]
    paste(x, collapse = ";")
  }
  #environment(f0) <- .GlobalEnv
  pg$"Common Names" <- parallel::parSapply(cl, tmp2, f0)
  w <- which(pg$"Common Names" == "")
  if (length(w)) { if ("Protein names" %in% colnames(pg)) { pg$"Common Names"[w] <- pg$"Protein names"[w] } }
  w <- which(pg$"Common Names" == "")
  if (length(w)) { if ("Names" %in% colnames(pg)) { pg$"Common Names"[w] <- pg$Names[w] } }
  pg$"Common Name (short)" <- vapply(strsplit(pg$"Common Names", ";"), function(x) {
    x <- unlist(x)[1]
    if (!length(x)) { x <- "" }
    return(x)
  }, "")
  id.col <- c("Names", "Protein IDs", "Gene names")
  id.col <- id.col[which(id.col %in% colnames(pg))]
  w <- which(pg$"Common Name (short)" == "")
  if ((length(id.col))&&(length(w) > 0)) {
    pg$"Common Name (short)"[w] <- apply(pg[w, id.col], 1, function(x) {
      ww <- which(!x %in% c("", " ", "NA", NA))
      if (length(ww) > 0) { x <- x[ww[1]] } else { x <- "" }
      return(x)
    })
  }
}
# PG label column
pg$temp <- gsub(";.+", ";...", pg$"Leading protein IDs")
pg$Label <- apply(pg[, c("temp", "Common Name (short)")], 1, function(x) {
  if (is.na(x[[2]])) { x <- x[[1]] } else { x <- paste0(x[[1]], " - ", x[[2]]) }
  return(x)
})
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
PG_assembly <- list(Protein.groups = pg, Peptides = Pep, Database = DB)
if (!misFun(Ev)) { PG_assembly$Evidences <- Ev }
invisible(parallel::clusterCall(cl, function(x) {
  try(rm(tmp1, tmp2, tmp3, tmp4, CustPG, c1, c2, ca, cb), silent = TRUE)
  return()
}))
if (cleanUp) { parallel::stopCluster(cl) }
if (TESTING) {
  tm2 <<- Sys.time()
  print(tm2-tm1)
}
saveFun(PG_assembly, file = "PG_assembly.RData")
