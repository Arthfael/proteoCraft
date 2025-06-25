#' Digest
#'
#' @description
#' A function to digest unmodified protein Sequence into peptides.
#' 
#' Will by default allow for loss of N-terminal methionine.
#' Default parameters digest like trypsin.
#' Proline after the R or K is assumed to block some but not all cleavages, hence both missed and cleaved peptides are included.
#' 
#' Uses parallel processing by default for now.
#' Still needs some speed boosts. The whole thing should be parallelized in one big go, as opposed to by bits.
#' The regex approach currently used may not be ideal. It provided a small speed bump when introduced to replace list-based approaches, but may more less safe.
#' 
#' Used by ProtMatch, but not by the much faster/better ProtMatch2.
#' 
#' @param Seq The character vector of protein Sequence(s) to digest.
#' @param Cut The character vector of amino acid(s) after which to cut. For each, provide one letter followed or preceded by underscore indicating whether the cut is C- or N-terminal, respectively. Default (for Trypsin) = c("K_", "R_")
#' @param strict.avoid The pattern(s) that will strictly block cleavage when overlapping with cleavage sites. Write cleavage site as underscore. Default = ""
#' @param loose.avoid The pattern(s) that will partially block cleavage when overlapping with cleavage sites. Write cleavage site as underscore. Default =  c("K_P", "R_P")
#' @param missed How many missed cleavages should be considered? Default = 2.
#' @param min The minimum length of peptides to report. Default = 7
#' @param max The maximum length of peptides to report. Set to FALSE (default) to allow any peptide length above min.
#' @param characters.test Should the compatibility of Sequences with the vector of valid amino acids be tested (True by default)?
#' @param AA Character vector of which amino Acids letters are considered valid? Default is the 20 standard ones plus S for selenocysteine: c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","U","V","W","Y")
#' @param RemoveNtermMet The default "loose" allows for both retention and loss of each N-terminal methionine. Set to "strict" to force loss. Set to "predict" to apply the rules as summarized in Wingfield PT. N-Terminal Methionine Processing. Curr Protoc Protein Sci. 2017 Apr 3;88:6.14.1-6.14.3. doi: 10.1002/cpps.29. PMID: 28369664; PMCID: PMC5663234. 
#' @param N.clust A limit on the number of vCPUs to use. If left as NULL (default), uses the number of available clusters - 1, to a minimum of 1. Only used if missed > 0
#' @param N.reserved Default = 1. Number of reserved vCPUs the function is not to use. Note that for obvious reasons it will always use at least one. Only used if missed > 0
#' @param collapse For compatibility with older usage: option to specify a character (e.g. ";") to collapse the peptides list for each protein into a single character vector.
#' @param ChnkSz Deprecated
#' @param cl Already have a cluster handy? Why spend time making a new one, which on top of that may invalidate the old one. Just pass it along!
#' 
#' @examples
#' Digest(Seq = "MAAAAAAAAAKCCCCCCCRD")
#' 
#' @import data.table
#' @export

Digest <- function(Seq,
                   Cut = c("K_", "R_"),
                   strict.avoid = "",
                   loose.avoid = c("K_P", "R_P"),
                   missed = 2,
                   min = 7,
                   max = FALSE,
                   characters.test = TRUE,
                   AA = proteoCraft::AA,
                   RemoveNtermMet = "loose",
                   N.clust,
                   N.reserved = 1,
                   collapse,
                   ChnkSz = Inf,
                   cl) {
  parThresh <- 500
  TESTING <- FALSE
  #proteoCraft::DefArg(proteoCraft::Digest)
  
  #DB <- db; Seq = setNames(DB$Sequence, DB$`Protein ID`)
  #characters.test = TRUE; N.clust = 5
  #TESTING <- TRUE
  #Cut = c("Y_", "W_", "F_")
  #Cut = c("E_")
  if (TESTING) {
    # Note:
    # This is not a perfect alternative to missing but will work in most cases, unless x matches a function imported by a package 
    misFun <- function(x) { return(!exists(deparse(substitute(x)))) }
  } else { misFun <- missing }
  Cut <- toupper(Cut)
  stopifnot(class(missed) %in% c("numeric", "integer", "integer64") & missed == round(missed),
            missed >= 0,
            min(nchar(Cut)) == 2,
            max(nchar(Cut)) == 2)
  usePar <- FALSE
  if ((length(Seq) > parThresh)&&((missed > 0)||(length(loose.avoid)))) {
    usePar <- require(parallel)
    if (usePar) {
      #
      # Create cluster
      tstCl <- misFun(cl)
      if (!misFun(cl)) {
        tstCl <- suppressWarnings(try({
          a <- 1
          clusterExport(cl, "a", envir = environment())
        }, silent = TRUE))
        tstCl <- !"try-error" %in% class(tstCl)
      }
      if ((misFun(cl))||(!tstCl)) {
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
      }
      N.clust <- length(cl)
      #
    }
  }
  RemoveNtermMet <- tolower(RemoveNtermMet)
  strict.avoid <- toupper(strict.avoid)
  strict.avoid <- strict.avoid[which(strict.avoid != "")]
  loose.avoid <- toupper(loose.avoid)
  loose.avoid <- loose.avoid[which(loose.avoid != "")]
  l <- which(strict.avoid %in% loose.avoid)
  if (length(l)) {
    warning("There is overlap between which patterns to avoid strictly vs loosely, we will be conservative and set those to loose!")
    strict.avoid <- strict.avoid[which(!strict.avoid %in% loose.avoid)]
  }
  SEQ <- toupper(gsub("\\*.*", "", Seq)) # Remove stop codons (*)
  SEQ <- SEQ[which(nchar(SEQ) > 0)]
  # The little alchemy below is to deal with cases where the sequences include duplicates of the same accession
  lSEQ <- length(SEQ)
  if (is.null(names(SEQ))) { names(SEQ) <- paste0("Protein ", 1:lSEQ) }
  SEQ <- data.table::data.table(SEQ = SEQ, Name = names(SEQ))
  SEQ <- SEQ[, list(x = list(SEQ)), by = list(Group.1 = Name)]
  SEQ <- as.data.frame(SEQ)
  stopifnot(min(vapply(SEQ$x, length, 1)) == 1)
  stopifnot(max(vapply(SEQ$x, length, 1)) == 1) # This indicates that the same accession has been provided with different sequences,
  # which is not acceptable!
  # Create non-redundant protein names
  SEQ <- setNames(unlist(SEQ$x), SEQ$Group.1)
  if (usePar) {
    ChnkSz <- lSEQ
    ChnkSz <- ceiling(ChnkSz/N.clust)
    Chnks <- data.frame(Start = 1, End = lSEQ)
    n <- ceiling(lSEQ/ChnkSz)
    Chnks <- round((1:(n-1))*lSEQ/n)
    Chnks <- c(Chnks, lSEQ)
    Chnks <- data.frame(Start = c(1, Chnks[1:(n-1)]+1),
                        End = Chnks)
    stopifnot(Chnks$Start[1] == 1,
              rev(Chnks$End)[1] == lSEQ,
              sum(vapply(1:(nrow(Chnks)-1), function(x) {
                Chnks$Start[x+1]-Chnks$End[x]-1
              }, 1)) == 0)
    Chnks <- apply(Chnks, 1, function(x) { SEQ[x[[1]]:x[[2]]] })
  }
  # Optional characters test
  if (characters.test) {
    if (usePar) {
      test <- unlist(parallel::parLapply(cl, Chnks, function(x) {
        nchar(gsub(paste(proteoCraft::AA, collapse = "|"), "", x))
      }))
    } else {
      test <- nchar(gsub(paste(AA, collapse = "|"), "", SEQ))
    }
    if (max(test)) {
      if (lSEQ > 1) {
        w <- which(test > 0)
        l <- length(w)
        if (l > 1) {
          if (l > 100) { w <- paste0(paste(w[1:100], collapse = ", "), "...") } else {
            w <- paste0(paste(w[1:(l-1)], collapse = ", "), " and ", w[l])
          }
        }
        msg <- paste0("There are illegal characters in input Sequence", c("", "s")[(l > 1)+1], " number ", w, "!")
      } else { msg <- "There are illegal characters in the input Sequence!" }
      warning(msg)
    }
  }
  cllpsTst <- FALSE
  cllps <- ""
  if ((!misFun(collapse))&&(length(collapse) == 1)) {
    cllpsTst <- TRUE
    cllps <- as.character(collapse)
  }
  F0 <- function(Sq) { #Sq <- SEQ #Sq <- Chnks[[1]]
    # N-terminal Methionine:
    if (RemoveNtermMet == "strict") { Sq <- gsub("^M", "", Sq) }
    if (RemoveNtermMet == "loose") { Sq <- gsub("^M", "@M", Sq) }
    # For loose removal, the above method:
    # - ensures that we know which M are protein N-terminal and that M is still available normally for potential cleavage sites
    # - and avoids duplicating the list - we can deal with these peptides later (before applying size filters!!!)
    # It's a bit more complicated for the prediction method:
    # - "@" will still denote N-terminal Methionine, but only for ambiguous amino acids (T and V)
    if (RemoveNtermMet == "predict") {
      # Identify cases where Proline blocks cleavage
      gPr <- grep("^M[GASCPTV]P", Sq)
      # Identify N-terminal Methionines which should be systematically removed
      gSmll <- grep("^M[GASCP]", Sq)
      gSmll <- gSmll[which(!gSmll %in% gPr)] # Apply proline block
      # Identify ambiguous cases
      gAmb <- grep("^M[VT]", Sq)
      gAmb <- gAmb[which(!gAmb %in% gPr)] # Apply proline block
      # Apply
      Sq[gSmll] <- gsub("^M", "", Sq[gSmll])
      Sq[gAmb] <- paste0("@", Sq[gAmb])
    }
    if (RemoveNtermMet %in% c("loose", "predict")) { AT <- grepl("^@M", Sq) }
    # Introduce breaks in the sequence
    if (TESTING) { print("Identifying cleavage sites...") }
    for (C in Cut) { Sq <- gsub("_$", "", gsub(gsub("_", "", C), C, Sq)) }
    # Deal with patterns strictly blocking digest 
    if (length(strict.avoid)) { #strict.avoid <- loose.avoid # (for TESTING)
      for (sa in strict.avoid) { Sq <- gsub(sa, gsub("_", "", sa), Sq) }
    }
    # It is too inefficient to deal with combinatorics of loose avoidance sites at this stage
    # I will first assume that they are strictly blocking, and allow the maximum of missed cleavage sites,
    # then I will select those which have at least one and apply combinatorics on them.
    if (length(loose.avoid)) {
      for (la in loose.avoid) { Sq <- gsub(la, gsub("_", "", la), Sq) }
    }
    #
    # Break Sequence into peptides
    if (TESTING) { print(paste0("Creating fully cleaved peptides")) }
    Sq <- strsplit(Sq, "_")
    # Apply missed cleavages
    Rs <- Sq
    if (missed) {
      L <- vapply(Sq, length, 1)
      # Seeing the size of the input, a good old "for" loop here is actually probably a good idea
      # so as not to overly tax memory...
      for (ms in 1:missed) { #ms <- 1
        if (TESTING) { print(paste0("Creating peptides with ", ms, "-missed cleavages...")) }
        w <- which(L >= ms+1)
        if (length(w)) {
          temp <- setNames(lapply(w, function(x) { #x <- 1
            vapply(1:(L[x]-ms), function(p) { paste(Sq[[x]][p:(p+ms)], collapse = "") }, "")
          }), names(Sq)[w])
          Rs[w] <- mapply(c, Rs[w], temp, SIMPLIFY = FALSE)
        }
      }
    }
    # Now deal with loose.avoid
    if (length(loose.avoid)) {
      if (TESTING) { print(paste0("Creating effect of patterns partially but not fully blocking digest...")) }
      temp <- proteoCraft::listMelt(Rs)
      gr <- grep(paste(gsub("_", "", loose.avoid), collapse = "|"), temp$value)
      temp <- temp[gr,]
      # Restore cuts
      for (C in Cut) { temp$value <- gsub("_$", "", gsub(gsub("_", "", C), C, temp$value)) }
      temp$value <- strsplit(temp$value, "_")
      L <- vapply(temp$value, length, 1)
      if (length(L)) {
        temp$value <- vapply(1:length(L), function(x) { #x <- 1
          paste(unlist(sapply(1:L[x], function(p1) {
            vapply(p1:L[x], function(p2) {
              paste(temp$value[[x]][p1:p2], collapse = "")
            }, "")
          })), collapse = ";")
          # Collapsing seems actually more efficient here - the vectors end up shorter
        }, "")
        temp <- aggregate(temp$value, list(temp$L1), paste, collapse = ";") # Not accelerated by data.table!
        temp <- setNames(strsplit(temp$x, ";"), temp$Group.1)
        w <- which(names(Rs) %in% names(temp))
        Rs[w] <- mapply(c, Rs[w], temp, SIMPLIFY = FALSE)
      }
    }
    if (RemoveNtermMet %in% c("loose", "predict")) {
      Rs[AT] <- lapply(Rs[AT], function(x) { unique(c(gsub("^@M", "M", x),
                                                      gsub("^@M", "", x)))
      })
    }
    Rs <- lapply(Rs, function(x) { unique(x[which(nchar(x) >= min)]) }) # Filter for min and remove rare duplicates
    if (max) { Rs <- lapply(Rs, function(x) { x[which(nchar(x) <= max)] }) } # Optional: filter for max
    #tst <- unique(gsub("[A-Z]", "", unlist(Rs)))
    if (cllpsTst) { Rs <- vapply(Rs, paste, "", collapse = cllps) }
    return(Rs)
  }
  if (usePar) {
    environment(F0) <- .GlobalEnv
    parallel::clusterExport(cl, list("TESTING", "misFun", "RemoveNtermMet", "Cut", "strict.avoid",
                                     "loose.avoid", "missed", "min", "max", "cllpsTst", "cllps"),
                            envir = environment())
    RES <- parallel::parLapply(cl, Chnks, F0)
    #
    parallel::stopCluster(cl)
    RES <- do.call(c, RES)
  } else {
    RES <- F0(SEQ)
  }
  if (!TESTING) { return(RES) }
}
