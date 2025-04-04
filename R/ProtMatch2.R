#' ProtMatch2
#'
#' @description
#' 
#' Once again, the function has been radically re-written, though without renaming it this time!
#' The indexing from the previous version is gone, replaced by better parallelisation and data.tables use, yielding a small improvement.
#' A small, rare bug is also corrected in the process.
#'
#' @param Seq Character vector of peptide sequences.
#' @param DB The formated protein database, with proteins IDs and sequences.
#' @param Cut For digestion. The character vector of amino acid(s) after which to cut. For each, provide one letter followed or preceded by underscore indicating whether the cut is C- or N-terminal, respectively (usually C-terminal). Default (for Trypsin) = c("K_", "R_")
#' @param N.clust A limit on the number of vCPUs to use. If left as NULL (default), uses the number of available clusters - 1, to a minimum of 1. Only used if missed > 0
#' @param N.reserved Default = 1. Number of reserved vCPUs the function is not to use. Note that for obvious reasons it will always use at least one. Only used if missed > 0
#' @param cl Already have a cluster handy? Why spend time making a new one, which on top of that may invalidate the old one. Just pass it along!
#'
#' @import data.table
#' @export

ProtMatch2 <- function(Seq,
                       DB,
                       Cut = c("K_", "R_"),
                       N.clust = get("N.clust"),
                       N.reserved = 1,
                       cl) {
  TESTING <- FALSE
  #proteoCraft::DefArg(proteoCraft::ProtMatch2)
  #Seq = unique(ev$Sequence);DB = db
  #Seq = unique(pep$Sequence);DB = db
  #w <- 1:nrow(ev); Seq = unique(ev$Sequence[w]); DB = db
  #TESTING <- TRUE
  if (TESTING) {
    # Note:
    # This is not a perfect alternative to missing but will work in most cases, unless x matches a function imported by a package 
    misFun <- function(x) { return(!exists(deparse(substitute(x)))) }
    if (!require(tictoc)) { install.packages("tictoc") }
  } else { misFun <- missing }
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
  fF <- function(x) { paste(sort(x), collapse = ";") } # Final protein aggregate function
  environment(fF) <- .GlobalEnv
  # Pre-process observed sequences
  Seq2 <- data.frame(Sequence = unique(Seq)) # Unique peptide sequences
  Seq2$DigSeq <- Seq2$Sequence
  for (C in Cut) {
    Seq2$DigSeq <- gsub(gsub("_", "", C), C, Seq2$DigSeq) # Insert cut sites before conversion of Ile to Leu!
  }
  Seq2$DigSeq <- gsub("^_|_$", "", Seq2$DigSeq)
  Seq2$DigSeq <- gsub("I", "L", Seq2$DigSeq) # Convert Ile to Leu
  Seq2$myDigest <- strsplit(Seq2$DigSeq, "_") # Digest
  Seq2$L <- sapply(Seq2$myDigest, length)
  Seq2$Parent_ID <- paste0("lP", 1:nrow(Seq2)) # Here lP stands for "long peptide" (may have missed cleavages)
  #
  # Digest database and match to predicted peptides:
  Nms <- DB$`Protein ID`
  dbSeq <- DB$Sequence
  for (C in Cut) {
    repl <- gsub("_", "", C)
    dbSeq <- stringi::stri_replace_all_regex(dbSeq, repl, C) # Insert cut sites before conversion of Ile to Leu!
  }
  dbSeq <- stringi::stri_replace_all_regex(dbSeq, "^_|_$", "")
  dbSeq <- stringi::stri_replace_all_regex(dbSeq, "I", "L") # Convert Ile to Leu
  dbSeq <- setNames(dbSeq, Nms)
  Dig <- setNames(strsplit(dbSeq, "_"), Nms) # Digest
  g <- grep("^M", DB$Sequence)
  frstPep <- setNames(sapply(Dig[g], function(x) { x[[1]] }), Nms[g])
  frstPep <- gsub("^M", "", frstPep)
  parallel::clusterExport(cl, c("Dig", "frstPep", "Nms", "g"), envir = environment())
  f0 <- function(x) {
    seq <- Dig[[Nms[x]]]
    pos <- 1:length(seq)
    if (x %in% g) { # Deal with N-term Meth
      m <- match(x, g)
      seq <- c(seq, frstPep[m])
      pos <- c(pos, 1)
    }
    data.frame(Pos = pos,
               Seq = seq,
               Prot = Nms[x])
  }
  environment(f0) <- .GlobalEnv
  Dig <- setNames(parallel::parLapply(cl, 1:length(Dig), f0), Nms)
  # Table mapping database proteins to their fragments
  Dig2 <- plyr::rbind.fill(Dig)
  #
  # Filter proteins using observed fragments
  allProt <- unique(Dig2$Prot[which(Dig2$Seq %in% unlist(Seq2$myDigest))])
  Dig2 <- Dig2[which(Dig2$Prot %in% allProt),]
  Dig <- Dig[which(names(Dig) %in% allProt)]
  #
  # Match fragments to proteins which could have generated them
  Frag2Prot <- data.table::as.data.table(Dig2[, c("Seq", "Prot")])
  Frag2Prot <- as.data.frame(Frag2Prot[, list(Prot = list(unique(Prot))), by = list(Seq = Seq)]) 
  #
  Res <- list()
  #
  w1 <- which(Seq2$L == 1)
  wM <- which(Seq2$L > 1)
  stopifnot(sum(!Seq %in% Seq2$Sequence) == 0)
  stopifnot(sum(!unique(c(w1, wM) %in% 1:nrow(Seq2))) == 0)
  if (length(w1)) {
    m1 <- match(Seq2$DigSeq[w1], Frag2Prot$Seq) # Ok to use match here because unicity
    Res$Full <- data.frame(Proteins = parallel::parSapply(cl, Frag2Prot$Prot[m1], fF),
                           Sequence = Seq2$Sequence[w1])
    #View(Res[[1]])
  }
  # if (TESTING) {
  #   kount <<- 0
  #     wM1 <- which(Seq %in% Seq2$Sequence[wM])
  #     tstFun <- function() {
  #     kount <<- kount + 1
  #     kol <- c("Parent", "Sequence")
  #     kol <- kol[which(kol %in% colnames(Par2Prot))][1]
  #     a <- length(which(!Seq[wM1] %in% Par2Prot[[kol]])) == 0
  #     cat("Step ", kount, "\n  -> ", a, "\n\n")
  #   }
  # }
  if (length(wM)) {
    Seq2flt <- Seq2[wM,]
    #
    # Table of matches between observed peptides and fragments
    Frags <- proteoCraft::listMelt(Seq2flt$myDigest, Seq2flt$Parent_ID, c("Fragment", "Parent_ID"))
    tmp <- c(0, cumsum(Seq2flt$L[1:(nrow(Seq2flt)-1)]))
    tmp <- tmp[match(Frags$Parent_ID, Seq2flt$Parent_ID)]
    Frags$Frag_ID <- 1
    Frags$Frag_ID <- cumsum(Frags$Frag_ID)-tmp
    Frags <- Frags[which(Frags$Fragment %in% Dig2$Seq),] # Filter
    Frags$fullFrag_ID <- do.call(paste, c(Frags[, c("Parent_ID", "Fragment", "Frag_ID")], sep = "_"))
    # Reciprocally filter Dig and Dig2 for observed fragments
    Dig2 <- Dig2[which(Dig2$Seq %in% Frags$Fragment),]
    Dig <- Dig[which(names(Dig) %in% Dig2$Prot)]
    # Table of matches between database proteins and fragments
    Frags$allProts <- Frag2Prot$Prot[match(Frags$Fragment, Frag2Prot$Seq)]
    # Check which proteins match the requirements in terms of matching at least once all fragments
    Par2Prot <- proteoCraft::listMelt(Frags$allProts, Frags$fullFrag_ID, c("Prot", "fullFrag_ID"))
    Par2Prot$Parent_ID <- Frags$Parent_ID[match(Par2Prot$fullFrag_ID, Frags$fullFrag_ID)]
    #which(Par2Prot$Parent_ID != gsub("_.*", "", Par2Prot$fullFrag_ID))
    Par2Prot <- data.table::as.data.table(Par2Prot)
    Par2Prot <- Par2Prot[, list(nFrags = length(fullFrag_ID)), by = list(Parent_ID = Parent_ID, Prot = Prot)]
    Par2Prot[, c("Parent", "Parent_L")] <- Seq2flt[match(Par2Prot$Parent_ID, Seq2flt$Parent_ID), c("Sequence", "L")]
    if (TESTING) { tstFun() }
    #
    w <- which(Par2Prot$nFrags >= Par2Prot$Parent_L)
    Par2Prot <- Par2Prot[w,]
    #if (TESTING) { tstFun() }
    Par2Prot <- as.data.frame(Par2Prot)
    #if (TESTING) { tstFun() }
    allProt <- unique(Par2Prot$Prot)
    # Filter proteins using fragments
    Dig2 <- Dig2[which(Dig2$Prot %in% allProt),]
    Dig <- Dig[which(names(Dig) %in% allProt)]
    # We will test all possible Fragment/Protein combinations in Par2Prot
    tmpFrags <- setNames(Seq2flt$myDigest, Seq2flt$Parent_ID)
    parallel::clusterExport(cl, list("tmpFrags", "Dig"), envir = environment())
    fm <- function(x) {
      #x <- Par2Prot[1, c("Parent_ID", "Prot")]
      fr1 <- tmpFrags[[x[[1]]]]
      fr2 <- Dig[[x[[2]]]]
      l <- length(fr1)
      tst <- unlist(lapply(1:l, function(x) { fr2$Pos[which(fr2$Seq == fr1[x])]-x+1 }))
      rs <- FALSE
      if (length(tst)) { rs <- max(aggregate(tst, list(tst), length)$x) == l }
      rs
    }
    environment(fm) <- .GlobalEnv
    Par2Prot$Test <- parallel::parApply(cl, Par2Prot[, c("Parent_ID", "Prot")], 1, fm)
    if (TESTING) { tstFun() }
    #View(Par2Prot[which(Par2Prot$Parent_ID %in% Seq2flt$Parent_ID[match(Aha, Seq2flt$Sequence)]),])
    #
    # (These are fully digested parent peptides)
    Par2Prot <- Par2Prot[which(Par2Prot$Test),]
    #if (TESTING) { tstFun() }
    Par2Prot <- data.table::data.table(Par2Prot[, c("Parent", "Prot")])
    #if (TESTING) { tstFun() }
    Par2Prot <- Par2Prot[, list(Proteins = fF(Prot)), by = list(Sequence = Parent)]
    #if (TESTING) { tstFun() }
    Par2Prot <- as.data.frame(Par2Prot)
    #if (TESTING) { tstFun() }
    #Par2Prot <- aggregate(Par2Prot$Prot, list(Par2Prot$Parent), fF)
    #colnames(Par2Prot) <- c("Sequence", "Proteins")
    Res$Misses <- Par2Prot
  }
  Res <- plyr::rbind.fill(Res)
  Res <- Res[, c("Sequence", "Proteins")]
  #stopifnot(sum(!Seq %in% Res$Sequence) == 0)
  w <- which(!Seq %in% Res$Sequence)
  if (length(w)) {
    tmp <- data.frame(Sequence = Seq[w], Proteins = "")
    Res <- rbind(Res, tmp)
  }
  row.names(Res) <- 1:nrow(Res)
  if (cleanUp) { parallel::stopCluster(cl) }
  return(Res)
}
