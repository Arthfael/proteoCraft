#' ProtMatch
#'
#' @description
#' A source to match observed peptide sequences to a protein digest:
#'   - Some search software do not report all proteins from the database which match a specific sequence.
#'   - Others, such as MaxQuant, do report weird matches for a minority of PSMs, or miss a few perfectly valid matches.
#'  Whilst this may be correct from their point of view, this function allows re-checking matches for consistency between software.
#'  This will return a data frame of matches between observed peptides and protein IDs from a formatted protein sequences database.
#' 
#' 
#' 17/07/2025
#' ----------
#' Continuing issues with parallelisation mean I have decided to switch this to a source, as with many such functions.
#' 
#' 
#' 10/06/2025
#' ----------
#' Aaand... once again, this function has been radically re-written. Again, I did not rename it, because as far as I can tell the outcome is the same!
#' 
#' The way the previous version used data.tables was not that great (creating un-necessarily large tables), and the function remained too slow.
#' I have dropped those parts and instead have focused on filtering stuff as early as possible.
#' I have also improved the parallelisation: now instead of using parallel::clusterExport() for large objects, the functions writes (then deletes) temporary .RDS files.
#' If at some point writing those should fail, we could try the following solutions:
#'  - using tmpfile() for those instead of a fixed-name file,
#'  - or allow reverting to clusterExport() if .RDS serialization doesn't work.
#' 
#' These changes have yielded an about 20x gain in speed, turning the function, then its script heir, from a cumbersome monster to a decently fast one.
#'
#'
#' @param Seq Character vector of peptide sequences.
#' @param DB The formated protein database, with proteins IDs and sequences.
#' @param Cut For digestion. The character vector of amino acid(s) after which to cut. For each, provide one letter followed or preceded by underscore indicating whether the cut is C- or N-terminal, respectively (usually C-terminal). Default (for Trypsin) = c("K_", "R_")
#' @param N.clust A limit on the number of vCPUs to use. If left as NULL (default), uses the number of available clusters - 1, to a minimum of 1. Only used if missed > 0
#' @param N.reserved Default = 1. Number of reserved vCPUs the function is not to use. Note that for obvious reasons it will always use at least one. Only used if missed > 0
#' @param cl Already have a cluster handy? Why spend time making a new one, which on top of that may invalidate the old one. Just pass it along!
#' @param I_eq_L Should we consider I and L identical? Currently, by default, TRUE for both DIA and DDA: see https://github.com/vdemichev/DiaNN/discussions/1631
#'
#' @import data.table
#' @export

if (!exists("I_eq_L")) { I_eq_L <- TRUE } # Sensible default
I_eq_L %<o% I_eq_L
if (!exists("Cut")) { Cut <- c("K_", "R_") } # Provide default here,
# but defining Cut before running this allows to override it.
# Should be in the future pulled from the search parameters.
misFun <- function(x) { return(!exists(deparse(substitute(x)))) }
#
RPath <- as.data.frame(library()$results)
RPath <- normalizePath(RPath$LibPath[match("proteoCraft", RPath$Package)], winslash = "/")
libPath <- paste0(RPath, "/proteoCraft")
parSrc <- paste0(libPath, "/extdata/R scripts/Sources/make_check_Cluster.R")
#rstudioapi::documentOpen(parSrc)
source(parSrc, local = FALSE)
#
# Define final protein aggregate function
fAggr0 <- function(x) { paste(sort(unique(x)), collapse = ";") }
#
# Pre-process observed sequences
Seq2 <- data.frame(Sequence = unique(Seq)) # Unique observed peptide sequences
Seq2$DigSeq <- Seq2$Sequence
for (C in Cut) {
  repl <- gsub("_", "", C)
  Seq2$DigSeq <- stringi::stri_replace_all_regex(Seq2$DigSeq, repl, C,
                                                 opts_regex = stringi::stri_opts_regex(dotall = TRUE)) # -> insert cut sites before conversion of Ile to Leu!
}
Seq2$DigSeq <- stringi::stri_replace_all_regex(Seq2$DigSeq, "^_|_$", "",
                                               opts_regex = stringi::stri_opts_regex(dotall = TRUE))
if (I_eq_L) { # Convert Ile to Leu:
  # Should be after we introduce the cuts: those could be I/L sensitive
  Seq2$DigSeq <- gsub("I", "L", Seq2$DigSeq)
}
Seq2$myDigest <- strsplit(Seq2$DigSeq, "_") # -> digest
Seq2$L <- vapply(Seq2$myDigest, length, 1)
Seq2$Parent_ID <- paste0("lP", 1:nrow(Seq2)) # -> here "lP" stands for "long peptide" (may have missed cleavages)
#
# Fully digest database and match to predicted peptides:
Nms <- DB$`Protein ID`
DBSeq <- DB$Sequence # -> all protein sequences in the database
for (C in Cut) {
  repl <- gsub("_", "", C)
  DBSeq <- stringi::stri_replace_all_regex(DBSeq, repl, C) # -> insert cut sites before conversion of Ile to Leu!
}
DBSeq <- stringi::stri_replace_all_regex(DBSeq, "^_|_$", "",
                                         opts_regex = stringi::stri_opts_regex(dotall = TRUE))
if (I_eq_L) { # Convert Ile to Leu:
  # Should be after we introduce the cuts: those could be I/L sensitive
  DBSeq <- stringi::stri_replace_all_regex(DBSeq, "I", "L",
                                           opts_regex = stringi::stri_opts_regex(dotall = TRUE)) 
}
DBSeq <- setNames(DBSeq, Nms)
Dig <- setNames(strsplit(DBSeq, "_"), Nms) # -> our full digest
#
# Allow for N-term methionine loss
g <- grep("^M", DB$Sequence)
frstPep <- setNames(vapply(Dig[g], function(x) { x[[1]] }, ""), Nms[g])
frstPepNoMeth <- gsub("^M", "", frstPep)
# - Cluster export small objects
myWD <- getwd()
if ((exists("wd"))&&(dir.exists(wd))) { myWD <- wd }
parallel::clusterExport(parClust, c("Nms", "g", "myWD", "fAggr0"), envir = environment())
# - Use serialization to export efficiently large objects
readr::write_rds(Dig, paste0(myWD, "/tmpDig.RDS"))
readr::write_rds(frstPepNoMeth, paste0(myWD, "/1stPepNoMeth.RDS"))
invisible(parallel::clusterCall(parClust, function(x) {
  Dig <<- readr::read_rds(paste0(myWD, "/tmpDig.RDS")) # So it stays in cluster for next call!
  frstPepNoMeth <<- readr::read_rds(paste0(myWD, "/1stPepNoMeth.RDS")) # Same as above
  return()
}))
#f0 <- function(x) { exists("Dig") & exists("frstPepNoMeth") }
#tst <- parallel::clusterCall(parClust, f0)
f0 <- function(x) {
  seq <- Dig[[Nms[x]]]
  pos <- 1:length(seq)
  if (x %in% g) { # Deal with N-term Meth
    m <- match(x, g)
    seq <- c(seq, frstPepNoMeth[m])
    pos <- c(pos, 1)
  }
  #rm(Dig)
  return(data.frame(Pos = pos,
                    Seq = seq,
                    Prot = Nms[x]))
}
Dig <- setNames(parallel::parLapply(parClust, 1:length(Dig), f0), Nms)
# - Cleanup
unlink(paste0(myWD, "/tmpDig.RDS"))
unlink(paste0(myWD, "/1stPepNoMeth.RDS"))
#
# Table mapping database proteins to their fragments
Dig2 <- data.table::rbindlist(Dig)
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
wM <- which(!(1:nrow(Seq2) %in% w1))
stopifnot(sum(!Seq %in% Seq2$Sequence) == 0)
stopifnot(sum(!unique(c(w1, wM) %in% 1:nrow(Seq2))) == 0)
if (length(w1)) { # Full-digest peptides: easy peasy...
  m1 <- match(Seq2$DigSeq[w1], Frag2Prot$Seq) # Ok to use match here because unicity
  Res$Full <- data.frame(Sequence = Seq2$Sequence[w1],
                         Proteins = parallel::parSapply(parClust, Frag2Prot$Prot[m1], fAggr0))
  #View(Res$Full)
}
if (length(wM)) { # Peptides with missed cleavages: the real fun begins...
  #
  Seq2flt <- Seq2[wM,]
  #
  f0Mtch <- function(j, nFrag = 2, pepRange = myRng) { #j <- 1 #j <- j+1 #j <- j-1
    #
    res <- c()
    stopifnot(nFrag >= 1, nFrag == as.integer(nFrag))
    #
    fragRg <- 1:nFrag
    Seq2flt_i <- readr::read_rds(paste0(myWD, "/tmpA.RDS"))
    Frag2Prot_i <- readr::read_rds(paste0(myWD, "/tmpB.RDS"))
    Dig2_i <- readr::read_rds(paste0(myWD, "/tmpC.RDS"))
    a1 <- (c(0, pepRange)+1)[j]
    a2 <- pepRange[j]
    if (a2 >= a1) {
      rg_ij <- a1:a2
      Seq2flt_ij <- Seq2flt_i[rg_ij,]
      n_ij <- nrow(Seq2flt_ij) # = how many ij-peptides there are (size of current range)
      Seq2flt_ij_Dig <- Seq2flt_ij$myDigest # = full ij-digest list
      allFr_ij <- unique(unlist(Seq2flt_ij_Dig)) # = all ij-fragments
      w_ij <- which(Frag2Prot_i$Seq %in% allFr_ij) # = ij-filter for our frag-2-prot table
      if (length(w_ij)) {
        Frag2Prot_ij <- Frag2Prot_i[w_ij,] # = corresponding extract of the frag-2-prot table
        Frag2Prot_ij <- proteoCraft::listMelt(Frag2Prot_ij$Prot, Frag2Prot_ij$Seq, c("Protein", "Fragment"))
        wSinA <- which(Dig2_i$Seq %in% allFr_ij)
        if (length(wSinA)) {
          Dig2_ij <- Dig2_i[wSinA,] # = all frag-2-prot with potential positions
          #Dig2_j <- Dig2[which(Dig2$Seq %in% Frag2Prot_j$Seq),] # equivalent to above
          Dig2_ij$ID <- do.call(paste, c(Dig2_ij[, c("Prot", "Seq", "Pos")], sep = "_")) # All protein/peptide/pos fragments which are covered
          # If our missed-cleavage peptide exists, then its full digest sub-peptides are in there!
          fr_ij <- do.call(rbind, Seq2flt_ij_Dig)
          fr_ij_L <- do.call(rbind, lapply(Seq2flt_ij_Dig, nchar))
          # Identify the longest peptide fragment:
          # (We use the longest, most-specific peptide, so as to get a shorter list of candidate parent proteins!!!)
          wLongest <- apply(fr_ij_L, 1, function(x) {
            which(x == max(x))[1]
          })
          Longest <- vapply(1:n_ij, function(x) { fr_ij[x, wLongest[x]]}, "")
          # Other fragments
          wOthers <- do.call(rbind, lapply(wLongest, function(x) { fragRg[which(!fragRg %in% x)] }))
          Others <- do.call(rbind, lapply(1:n_ij, function(x) { fr_ij[x, wOthers[x,]] }))
          # Offset to apply to the position of those vs the longest peptide's own position, when searchging for matches
          posOffsets <- sweep(wOthers, 1, wLongest, "-")
          # All candidate proteins for each longest peptide
          candidates <- lapply(Longest, function(x) { which(Dig2_ij$Seq == x) })
          wLnon0 <- which(vapply(candidates, length, 1) > 0)
          if (length(wLnon0)) {
            candidates <- proteoCraft::listMelt(candidates[wLnon0], wLnon0, c("Row", "missedPeptide"))
            Dig2_ij <- as.data.frame(Dig2_ij) # Important for next step!!! At this stage, it's a data.table and does not respond well to the next call, somehow.
            candidates[, colnames(Dig2_ij)] <- Dig2_ij[candidates$Row, colnames(Dig2_ij)]
            candidates$Others <- lapply(candidates$missedPeptide, function(x) { Others[x,] })
            candidates$posOffsets <- lapply(candidates$missedPeptide, function(x) { posOffsets[x,] })
            candidates$Others_expect <- apply(candidates[, c("Prot", "Others", "Pos", "posOffsets")], 1, function(x) {
              list(paste(x[[1]], x[[2]], as.character(x[[3]]+x[[4]]), sep = "_"))
            }, simplify = FALSE)
            Others_expect <- as.data.frame(t(as.data.frame(candidates$Others_expect)))
            Others_notFound <- do.call(cbind, lapply(1:(i-1), function(x) {
              !Others_expect[[x]] %in% Dig2_ij$ID
            }))
            wAllFound <- which(rowSums(Others_notFound) == 0)
            if (length(wAllFound)) {
              candidates <- candidates[wAllFound,]
              candidates$missedPeptide <- Seq2flt_ij$Sequence[candidates$missedPeptide]
              res <- aggregate(candidates$Prot, list(candidates$missedPeptide), fAggr0)
              colnames(res) <- c("Sequence", "Proteins")
            }
          }
        }
      }
    }
    return(res)
  }
  #
  # For each number of missed cleavages (missed = i-1)
  m <- max(Seq2flt$L)
  for (i in 2:m) { #i <- 2 #i <- 3
    Seq2flt_i <- Seq2flt[which(Seq2flt$L == i),] # = peptides with i-1 missed cleavages
    n_i <- nrow(Seq2flt_i) # = how many we have
    cat(paste0(" -> Processing ", n_i, " peptides with ", i-1, " missed cleavage", c("", "s")[((i-1) > 1)+1], "...\n"))
    allFr_i <- unique(unlist(Seq2flt_i$myDigest)) # = all their fragments
    w_i <- which(Frag2Prot$Seq %in% allFr_i) # = i-filter for our frag-2-prot table
    if (length(w_i)) { # This check is necessary because of the rare - but not impossible - case where we find no peptide with N misses but some with N+1!
      Frag2Prot_i <- Frag2Prot[w_i,]
      Dig2_i <- Dig2[which(Dig2$Seq %in% allFr_i),] # = all i-frag-2-prot with potential positions
      myRng <- as.integer(round(as.numeric(1:N.clust)*n_i/N.clust)) # = ranges assigning i-peptides to cluster cores: we want to distribute our fragments over the cluster for efficient parallelisation
      myRng <- unique(myRng)
      myRng <- myRng[which(myRng > 0)]
      # Export temporary objects which will be read
      readr::write_rds(Seq2flt_i, paste0(myWD, "/tmpA.RDS"))
      readr::write_rds(Frag2Prot_i, paste0(myWD, "/tmpB.RDS"))
      readr::write_rds(Dig2_i, paste0(myWD, "/tmpC.RDS"))
      parallel::clusterExport(parClust, list("i", "myWD"), envir = environment())
      tmpResI <- parallel::parLapply(parClust, 1:length(myRng), f0Mtch, nFrag = i, pepRange = myRng)
      tmpResI <- do.call(rbind, tmpResI)
      Res[[paste0(i, "-missed cleavages")]] <- tmpResI
      unlink(paste0(myWD, "/tmpA.RDS"))
      unlink(paste0(myWD, "/tmpB.RDS"))
      unlink(paste0(myWD, "/tmpC.RDS"))
    }
  }
}
Res <- as.data.frame(data.table::rbindlist(Res))[, c("Sequence", "Proteins")]
#stopifnot(sum(!Seq %in% Res$Sequence) == 0)
w <- which(!Seq %in% Res$Sequence)
if (length(w)) {
  tmp <- data.frame(Sequence = Seq[w], Proteins = "")
  Res <- rbind(Res, tmp)
}
row.names(Res) <- 1:nrow(Res)
#
invisible(parallel::clusterCall(parClust, function(x) {
  try(rm(Dig, frstPepNoMeth, Seq2flt_i, Frag2Prot_i, Dig2_i), silent = TRUE)
  return()
}))
#
evmatch <- Res
