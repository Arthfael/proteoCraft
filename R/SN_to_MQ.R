#' SN_to_MQ
#'
#' @description 
#' Converts a Spectronaut PSMs table to a MaxQuant evidence.txt-like table.
#' As with the other functions of this family, the idea is not necessarily to get perfect conversion, but close enough that the data can be fed into this package's analysis scripts.
#' 
#' @returns
#' A list of 2 items:\cr
#' - Evidence: a data.frame similar to the evidence.txt table MaxQuant creates\cr
#' - PTMs: a modifications table\cr
#' 
#' @param SN_fls Path to (a) Spectronaut PSMs table file (.tsv, .parquet, or both).
#' @param Fixed.mods A vector of those Fixed modifications which will not be reported in the results (reported by diaNN, but MQ ignores them). Default = "Carbamidomethyl"
#' @param digPattern Default = "KR" for trypsin. Only used for the "Missed cleavages" column. Vector of amino acids after which the enzyme cleaves. 
#' @param N.clust A limit on the number of vCPUs to use. If left as NULL (default), uses the number of available clusters - 1, to a minimum of 1.
#' @param N.reserved Default = 1. Number of reserved vCPUs the function is not to use. Note that for obvious reasons it will always use at least one.
#' @param cl Already have a cluster handy? Why spend time making a new one, which on top of that may invalidate the old one. Just pass it along!
#' 
#' @examples
#' library(MSstatsPTM)
#' SN_fls <- "C:/Users/MYUSER/desktop/testSN.tsv"
#' write.table(MSstatsPTM::spectronaut_input, SN_fls, sep = "\t", row.names = FALSE)
#' temp <- SN_to_MQ(SN_fls)
#' ev <- temp$Evidences
#' Modifs <- temp$PTMs
#' 
#' @import data.table
#' @export

SN_to_MQ <- function(SN_fls,
                     Fixed.mods = c("Carbamidomethyl"),
                     digPattern = "KR",
                     N.clust,
                     N.reserved = 1L,
                     cl) {
  TESTING <- FALSE
  #DefArg(SN_to_MQ);TESTING <- TRUE;cl <- parClust
  misFun <- if (TESTING) {
    \(x) {
      return(!exists(deparse(substitute(x))))
    }
  } else { missing }
  #
  # Create cluster
  stopCl <- FALSE
  if (misFun(cl) || is.null(cl) || (!inherits(cl, "cluster"))) {
    dc <- parallel::detectCores()
    if (misFun(N.reserved)) { N.reserved <- 1L }
    nMax <- max(c(dc - N.reserved, 1L))
    if (misFun(N.clust)) { N.clust <- nMax } else {
      if (N.clust > nMax) {
        warning("More cores specified than allowed, I will ignore the specified number! You should always leave at least one free for other processes, see the \"N.reserved\" argument.")
        N.clust <- nMax
      }
    }
    cl <- parallel::makeCluster(N.clust, type = "SOCK")
    stopCl <- TRUE
  }
  N.clust <- length(cl)
  #
  utils::data(modifications, package = "PTMods", envir = environment())
  UniMod <- modifications
  rm(modifications)
  # Remove:
  # - substitutions
  gS <- grepl("^[A-Z][a-z]{2}->[A-Z][a-z]{2} substitution$", UniMod$Description)
  UniMod_S <- UniMod[gS,]
  UniMod <- UniMod[!gS,]
  # - deltas
  gD <- grepl("^Delta:", UniMod$Name)
  UniMod_D <- UniMod[gD,]
  UniMod <- UniMod[!gD,]
  #
  # Read PSMs (report file)
  wTsv <- grep("\\.tsv$", SN_fls)
  wPrqt <- grep("\\.parquet$", SN_fls)
  lTsv <- length(wTsv)
  lPrqt <- length(wPrqt)
  stopifnot(lTsv <= 1L,
            lPrqt <= 1L,
            sum(lTsv, lPrqt) > 0L)
  if (lTsv) {
    SN <- data.table::fread(SN_fls[wTsv], integer64 = "numeric", check.names = FALSE, data.table = FALSE)
  }
  if (lPrqt) {
    SN_prqt <- arrow::read_parquet(SN_fls[wPrqt])
    if (lTsv) {
      # In case we have both, we merge them, with priority to the parquet file (newer format)...
      # ... but first we check our assumptions:
      stopifnot(sum(SN$R.FileName != SN_prqt$R.FileName) == 0L,
                sum(SN$EG.PrecursorId != SN_prqt$EG.PrecursorId) == 0L)
      # If ok, proceed:
      k2 <- colnames(SN_prqt)
      SN[, k2] <- SN[, k2]
    } else {
      SN <- SN_prqt
    }
  }
  #
  # We only want precursors, not fragments!
  SN <- SN[SN$F.FrgLossType == "noloss",]
  #
  colnames(SN) <- gsub("[^A-Za-z0-9]", ".", colnames(SN)) # To be more resilient, in case the tables have been saved/processed in different manners
  #
  uMdSeq <- unique(SN$EG.ModifiedSequence)
  m1 <- match(SN$EG.ModifiedSequence, uMdSeq)
  uSeq <- gsub("^_|_$|\\[[^]]+\\]", "", uMdSeq)
  SN$Sequence <- uSeq[m1]
  #AA <- proteoCraft::AA # for testing
  test <- gsub(paste(AA, collapse = "|"), "", uSeq)
  test <- test[match(SN$Sequence, uSeq)]
  w <- which(test == "")
  lX <- nrow(SN) - length(w)
  if (lX) {
    u <- unique(unlist(strsplit(SN$Sequence, "")))
    u <- paste(u[!u %in% AA], collapse = "-")
    warning(paste0("Removing ", lX, " peptide", c("", "s")[(lX > 
                                                              1L) + 1L], " with unknown amino acids ", u, "!"))
    SN <- SN[w,]
  }
  #
  # Process PTMs
  tmp <- unique(sub("\\].*", "", unlist(strsplit(sub("^_[^[]+", "", grep("\\[", uMdSeq, value = TRUE)), "\\["))))
  tmp <- base::setdiff(tmp, "")
  allPTMs <- data.frame("Full mark" = tmp,
                     "Full name" = sub(" *\\(.*", "", tmp),
                     check.names = FALSE)
  #
  allPTMs$UniMod_ID <- UniMod$UnimodId[match(allPTMs$"Full name", UniMod$Description)]
  w <- which(is.na(allPTMs$UniMod_ID))
  allPTMs$UniMod_ID[w] <- UniMod$UnimodId[match(allPTMs$"Full name"[w], UniMod$Name)]
  #
  allPTMs$Site_long <- lapply(1L:nrow(allPTMs), \(x) { NULL })
  g <- grep("\\(", tmp)
  if (length(g)) {
    allPTMs$Site_long[g] <- gsub(".*\\(|\\)$", "", tmp[g])
    w <- g[which(!vapply(allPTMs$Site_long[g], \(x) { gsub("[ -]", "", tolower(x)) %in% c("nterm", "cterm", "proteinnterm", "proteincterm") }, TRUE))]
    if (length(w)) {
      allPTMs$Site_long[w] <- lapply(allPTMs$Site_long[w], strsplit, split = "")
    }
  }
  #
  w <- which((lengths(allPTMs$Site_long) == 0L) & (!is.na(allPTMs$UniMod_ID)))
  if (length(w)) {
    allPTMs$Site_long[w] <- lapply(allPTMs$UniMod_ID[w], \(x) { UniMod$Site[which(UniMod$UnimodId == x)] })
  }
  # Mass shift
  w <- which(!is.na(allPTMs$UniMod_ID))
  allPTMs$UniMod[w] <- paste0("UniMod:", allPTMs$UniMod_ID[w])
  allPTMs$"Mass shift" <- NA
  allPTMs$"Mass shift"[w] <- UniMod$MonoMass[match(allPTMs$UniMod_ID[w], UniMod$UnimodId)]
  w <- which(is.na(allPTMs$"Mass shift"))
  SN$"m/z if unmodified" <- parallel::parSapply(cl, SN$Sequence, Peptides::mw, monoisotopic = TRUE)
  if (length(w) || (!"FG.PrecMz" %in% colnames(SN))) {
    tmp <- strsplit(gsub("\\][^[]*\\[", "]><[", gsub("^_[^[]*\\[?|\\]?[^]]*_$", "", SN$EG.ModifiedSequence)), "\\]><\\[")
    uPTM <- unique(unlist(tmp))
    ptmMatr <- listMelt(tmp, ColNames = c("PTM", "w"), asDT = TRUE)
    ptmMatr$count <- 1L
    ptmMatr <- data.table::dcast(ptmMatr, w~PTM, sum)
    wPTMs <- ptmMatr$w
    ptmMatr <- as.data.frame(ptmMatr)
    ptmMatr <- as.matrix(ptmMatr[match(1L:nrow(SN), ptmMatr$w), allPTMs$`Full mark`])
    wNA <- which(is.na(ptmMatr), arr.ind = TRUE)
    ptmMatr[wNA] <- 0L
  }
  if (length(w)) {
    if (sum(!c("FG.PrecMz", "F.FrgLossType") %in% colnames(SN))) {
      stop("We could not map all PTMs to UniMod. We will need to map them to empirical delta masses, but for this we need columns \"FG.PrecMz\" and \"F.FrgLossType\" in the Spectronaut report!")
    }
    SN$PTM_effect <- (SN$FG.PrecMz-1.007276466879)*SN$FG.Charge - SN$"m/z if unmodified"
    fit <- MASS::rlm(x = ptmMatr[wPTMs,],
                     y = SN$PTM_effect[wPTMs],
                     weights = SN$"m/z if unmodified"[wPTMs],
                     wt.method = "inv.var",
                     psi = MASS::psi.huber,
                     maxit = 100L)
    coefFit <- coef(fit)
    allPTMs$"Mass shift - empirical" <- coefFit[allPTMs$`Full mark`]
    allPTMs$"Mass shift"[w] <- allPTMs$"Mass shift - empirical"[w]
  }
  EV <- data.frame(`Raw file` = SN$R.FileName, 
                   Proteins = SN$PG.ProteinGroups,
                   `Full Protein ID` = SN$PG.ProteinNames,
                   `Gene names` = SN$PG.Genes, 
                   Sequence = SN$Sequence,
                   `Mod. seq. (Spectronaut format)` = SN$EG.ModifiedSequence, 
                   Length = nchar(SN$Sequence),
                   Charge = SN$FG.Charge,
                   Intensity = SN$EG.TotalQuantity, # Left intentionally vague in case it goes further (was SN$EG.TotalQuantity..Settings. in my example )
                   Reverse = SN$EG.IsDecoy,
                   Modifications = "Unmodified",
                   check.names = FALSE)
  EV$"m/z" <- if ("FG.PrecMz" %in% colnames(SN)) {
    SN$FG.PrecMz
  } else {
    SN$"m/z if unmodified" + rowSums(sweep(ptmMatr, 2L, allPTMs$`Mass shift`, "*"))/SN$FG.Charge
  }
  EV$Mass <- (EV$"m/z"-1.007276466879)*EV$Charge
  EV$PEP <- if ("EG.PEP" %in% colnames(SN)) { SN$EG.PEP } else {
    if  ("EG.Qvalue" %in% colnames(SN)) { SN$EG.Qvalue } else {
      0 # We need some value because some of the downstream processing may expect a PEP; better an overly optimistic value than none!
    }
  }
  #
  allPTMs$Type <- c("Variable", "Fixed")[(allPTMs$"Full name" %in% Fixed.mods) + 1L]
  allPTMs$Position <- "Anywhere"
  allPTMs$Position[grep("C-term", allPTMs$"Full name", ignore.case = TRUE)] <- "Any C-term"
  allPTMs$Position[grep("N-term", allPTMs$"Full name", ignore.case = TRUE)] <- "Any N-term"
  allPTMs$Position[grep("protein C-term", allPTMs$"Full name", ignore.case = TRUE)] <- "Protein C-term"
  allPTMs$Position[grep("protein N-term", allPTMs$"Full name", ignore.case = TRUE)] <- "Protein N-term"
  #
  # PTM marks
  allPTMs$Mark <- tolower(substr(allPTMs$"Full name", 1L, 2L))
  ## Sometimes, some marks are duplicates, e.g. if you searched for "Acetyl (Protein N-term)" and "Acetyl (K)" together!
  ## We want to fix this so that each modification has a unique mark:
  tstMark <- aggregate(allPTMs$Mark, list(allPTMs$Mark), length)
  W <- which(tstMark$x > 1L)
  #allPTMs$Mark <- allPTMs$"Old mark"
  if (length(W)) {
    allPTMs$"Old mark" <- allPTMs$Mark
    for (i in W) {
      #i <- W[1L]
      w <- which(allPTMs$Mark == tstMark$Group.1[i])
      m <- allPTMs[w,]
      m$AA[lengths(m$AA) == 0L] <- "X"
      r <- if ("Acetyl" %in% m$"Full name") { which(m$"Full name" == "Acetyl") } else { 1L }
      s <- c(1L:nrow(m)); s <- s[s != r]
      test <- apply(m[s, c("AA", "Mark")], 1L, \(x) { paste0(tolower(x[[1L]]), substr(x[[2L]], 1L, 1L))[1L] })
      w0 <- which(!test %in% allPTMs$Mark)
      if (length(w0)) {
        m$Mark[s][w0] <- test[w0[1L]]
      }
      w1 <- which(test %in% allPTMs$Mark)
      if (length(w1)) {
        # not tested
        s <- s[w1]
        test <- sapply(s, list)
        kount <- 1L
        char <- c(as.character(0L:9L), letters)
        taken <- unique(c(allPTMs$Mark, m$Mark))
        for (j in s) {
          tst <- paste0(tolower(m$AA[s]), char[kount])
          while ((tst %in% taken) && (kount < length(char))) {
            kount <- kount+1L
            tst <- paste0(tolower(m$AA[s]), char[kount])
          }
          if (kount == length(char)) {
            stop("I am really out of options here, never thought this would go this far! Check the code just in case, this should never happen.")
          } else {
            test[[j]] <- tst
          }
        }
        m$Mark[[s]] <- unlist(test)
      }
      allPTMs[w,] <- m
    }
  }
  #
  wMod <- grep("\\[", SN$EG.ModifiedSequence)
  uMdSq <- unique(SN$EG.ModifiedSequence[wMod])
  tmpLst <- listMelt(strsplit(uMdSq, "\\[|\\]"), ColNames = c("modSeq", "wh"))
  tmpLst$Verbose <- tmpLst$modSeq
  wV <- which(tmpLst$modSeq %in% allPTMs$`Full mark`[allPTMs$Type == "Variable"])
  wF <- which(tmpLst$modSeq %in% allPTMs$`Full mark`[allPTMs$Type == "Fixed"])
  if (length(wV)) {
    mV <- match(tmpLst$modSeq[wV], allPTMs$`Full mark`)
    tmpLst$modSeq[wV] <- paste0("(", allPTMs$Mark[mV], ")")
    tmpLst$Verbose[wV] <- paste0("(", allPTMs$`Full name`[mV], ")")
  }
  if (length(wF)) {
    mF <- match(tmpLst$modSeq[wF], allPTMs$`Full mark`)
    tmpLst$modSeq[wF] <- ""
    tmpLst$Verbose[wF] <- paste0("(", allPTMs$`Full name`[mF], ")")
  }
  # Important: create data.tables individually, or use data.table::copy()
  tmp1 <- data.table::setDT(tmpLst[, c("modSeq", "wh")])
  tmp2 <- data.table::setDT(tmpLst[, c("Verbose", "wh")])
  tmp3 <- tmpLst
  tmp3$Mods <- gsub("^\\(|\\)$", "", tmp3$Verbose)
  tmp3 <- tmp3[tmp3$Mods %in% allPTMs$"Full name", c("Mods", "wh")]
  tmp3$Mods <- allPTMs$Mark[match(tmp3$Mods, allPTMs$"Full name")]
  data.table::setDT(tmp3)
  tmp1 <- tmp1[, .(modSeq = paste0(modSeq, collapse = "")),
               by = .(wh = wh)]
  tmp2 <- tmp2[, .(Verbose = paste0(Verbose, collapse = "")),
               by = .(wh = wh)]
  tmp3 <- tmp3[, .(n = .N),
               by = .(wh = wh, Mods = Mods)]
  tmp3$n <- as.character(tmp3$n)
  tmp3$modTag <- tmp3$Mods
  w2p <- which(tmp3$n != "1")
  tmp3$modTag[w2p] <- do.call(paste, c(tmp3[w2p, c("n", "Mods")], sep = " "))
  tmp3 <- tmp3[order(tmp3$Mods),]
  tmp3 <- tmp3[, .(Mods = paste(modTag, collapse = ",")), by = .(wh = wh),]
  tmp3 <- as.data.frame(tmp3)
  tmp3 <- tmp3[order(tmp3$wh),]
  EV$`Modified sequence` <- EV$`Modified sequence_verbose` <- EV$`Mod. seq. (Spectronaut format)`
  tmp1$Verbose <- tmp2$Verbose[match(tmp1$wh, tmp2$wh)]
  tmp1$Mods <- tmp3$Mods[match(tmp1$wh, tmp3$wh)]
  tmp1$Orig <- uMdSq[tmp1$wh]
  uM <- match(SN$EG.ModifiedSequence[wMod], tmp1$Orig)
  EV[wMod, c("Modified sequence", "Modified sequence_verbose", "Modifications")] <- tmp1[uM, c("modSeq", "Verbose", "Mods")]
  #
  # Retention time and IM:
  # For now I do not have a good example to map those, but will eventually want to...
  #
  digPat <- paste(c("[", digPattern, "]"), collapse = "")
  tmp <- gsub(paste0(digPat, "$"), "", EV$Sequence)
  EV$"Missed cleavages" <- nchar(tmp) - nchar(gsub(digPat, "", tmp))
  EV$"Potential contaminant" <- ""
  EV$Search_ID <- basename(SN_fls[c(wPrqt, wTsv)[1L]])
  #
  EV <- cbind(data.frame(id = 1L:nrow(EV)), EV)
  Res <- list(Evidence = EV,
              PTMs = allPTMs)
  if (stopCl) { parallel::stopCluster(cl) }
  return(Res)
}
