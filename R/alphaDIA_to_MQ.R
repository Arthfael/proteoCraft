#' alphaDIA_to_MQ
#'
#' @description 
#' Converts an alphaDIA precursors .tsv output table to a MaxQuant evidence.txt-like table.
#' As with the other functions of this type, the idea is not to get perfect conversion, but close enough that the data can be fed into this package's analysis scripts.
#' 
#' This will evolve as the Mann lab update https://alphadia.readthedocs.io/en/latest/methods/output-format.html
#' which as of 2025/07/16 is incomplete.
#' 
#' The output is a list of 2 items:
#' - Evidence: a data.frame similar to the evidence.txt table MaxQuant creates
#' - PTMs: a modifications table
#' 
#' @param alphaDIA_fl External alphaDIA precursors.tsv file.
#' @param Fixed.mods A vector of those Fixed modifications which will not be reported in the results (reported by alphaDIA, but MQ ignores them). Default = "Carbamidomethyl"
#' @param N.clust A limit on the number of vCPUs to use. If left as NULL (default), uses the number of available clusters - 1, to a minimum of 1.
#' @param N.reserved Default = 1. Number of reserved vCPUs the function is not to use. Note that for obvious reasons it will always use at least one.
#' @param cl Already have a cluster handy? Why spend time making a new one, which on top of that may invalidate the old one. Just pass it along!
#' @param shinyOpt One of "popup" (default), "pane", "dialog", or "browser".
#' 
#' @examples
#' alphaDIA_fl <- rstudioapi::selectFile()
#' alphaDIA_fl <- Param$PSMs
#' temp <- alphaDIA_to_MQ(alphaDIA_fl)
#' ev <- temp$Evidences
#' Modifs <- temp$PTMs
#' 
#' @import data.table
#' @export

alphaDIA_to_MQ <- function(alphaDIA_fl,
                           Fixed.mods = c("Carbamidomethyl"),
                           N.clust,
                           N.reserved = 1,
                           cl,
                           digPattern = "KR") {
  TESTING <- FALSE
  #proteoCraft::DefArg(proteoCraft::alphaDIA_to_MQ);TESTING <- TRUE; cl <- parClust
  #alphaDIA_fl <- rstudioapi::selectFile(path = paste0(getwd(), "/*.tsv"))
  #
  # NB:
  # This function deliberately uses a different object name for the PTMs table to other related functions,
  # to simplify environment-related debugging.
  # Do not change that!
  # Here, use Mods, elsewhere Modifs (or anything else).
  if (TESTING) {
    # Note:
    # This is not a perfect alternative to missing but will work in most cases, unless x matches a function imported by a package 
    misFun <- function(x) { return(!exists(deparse(substitute(x)))) }
  } else { misFun <- missing }
  #
  # Create cluster
  tstCl <- stopCl <- misFun(cl)
  if (!misFun(cl)) {
    tstCl <- suppressWarnings(try({
      a <- 1
      parallel::clusterExport(cl, "a", envir = environment())
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
  UniMod <- unimod::modifications
  # Remove:
  # - substitutions
  UniMod <- UniMod[grep("^[A-Z][a-z]{2}->[A-Z][a-z]{2} substitution$", UniMod$Description, invert = TRUE),]
  # - deltas
  g <- grepl("^Delta:", UniMod$Name)
  UniMod_D <- UniMod[which(g),]
  UniMod <- UniMod[which(!g),]
  #
  # Read PSMs ("precursors" file)
  xt <- gsub(".*\\.", "", alphaDIA_fl)
  if (xt == "tsv") {
    alphaDIA <- data.table::fread(alphaDIA_fl, integer64 = "numeric", check.names = FALSE, data.table = FALSE)
  }
  if (xt == "parquet") {
    alphaDIA <- arrow::read_parquet(alphaDIA_fl)
  }
  # Filter out peptides with ambiguous/non-classical amino acids
  uSeq <- unique(alphaDIA$sequence)
  test <- gsub(paste(proteoCraft::AA, collapse = "|"), "", uSeq)
  test <- test[match(alphaDIA$sequence, uSeq)]
  w <- which(test == "")
  lX <- nrow(alphaDIA)-length(w) 
  if (lX) {
    u <- unique(unlist(strsplit(alphaDIA$sequence, "")))
    u <- paste(u[which(!u %in% AA)], collapse = "-")
    warning(paste0("Removing ", lX, " peptide", c("", "s")[(lX > 1)+1], " with unknown amino acids ", u, "!"))
    alphaDIA <- alphaDIA[w,]
  }
  #
  stopifnot("run" %in% colnames(alphaDIA))
  log_Fl <- paste0(dirname(alphaDIA_fl), "/log.txt")
  if (file.exists(log_Fl)) { 
    log <- readLines(log_Fl)
    gP <- grep("^[0-9]+:[0-9]+:[0-9]+(\\.[0-9]+)? +INFO: +├──", log) # All parameter groups
    gA <- grep("^[0-9]+:[0-9]+:[0-9]+(\\.[0-9]+)? +INFO: +├──raw_paths:", log)
    gB <- vapply(gA, function(x) { min(gP[which(gP > x)]) }, 1)
    g0 <- unlist(lapply(1:length(gA), function(x) { (gA[x]+1):(gB[x]-1) }))
    rwFlPths <- gsub("\\\\", "/", unique(gsub(" +\\[.*", "", gsub("^[0-9]+:[0-9]+:[0-9]+(\\.[0-9]+)? +INFO: │   \033\\[[^ ]+ +", "", log[g0]))))
    rwFls <- gsub(".*/|\\.[^\\.]+$", "", rwFlPths)
  } else {
    msg <- "   log.txt file missing"
    rwFls <- unique(alphaDIA$run)
    rwFlPths <- paste0(dirname(alphaDIA_fl), "/", rwFls)
    if (sum(!file.exists(rwFlPths))) {
      rwFlPths <- paste0(gsub("/[^/]+$", "", dirname(alphaDIA_fl)), "/", rwFls)
    }
    if (sum(!file.exists(rwFlPths))) {
      msg <- paste0(", cannot get full paths of input raw files!")
      rwFlPths <- rwFls
    } else {
      msg <- paste0("... but we were to detect them automatically ;)")
    }
    warning(msg)
  }
  #
  # Quick check: it looks like "proba" is a PEP-like value:
  # High log(score) means high -log10(PEP)
  # Thus, for now we will treat it as such... until a guide is released or we are proven wrong.
  # require(ggplot2)
  # plot <- ggplot(alphaDIA) + geom_point(aes(x = log(score), y = -log10(proba))) + theme_bw() +
  #   facet_wrap(~run)
  # proteoCraft::poplot(plot)
  #
  # Create MQ-like table
  # Note: it seems that at least in some cases columns "pg", "pg_master" and "genes" are identical,
  # but that might be parameters- or fasta header-dependent.
  #w <- which(alphaDIA$pg != alphaDIA$pg_master)
  #length(w)
  #w <- which(alphaDIA$pg != alphaDIA$genes)
  #length(w)
  kol <- c("run", "pg", "sequence", "charge", "proba")
  #print(kol[which(!kol %in% colnames(alphaDIA))])
  stopifnot(sum(!kol %in% colnames(alphaDIA)) == 0)
  EV <- data.frame("Raw file" = alphaDIA$run,
                   "Raw file path" = rwFlPths[match(alphaDIA$run, rwFls)],
                   "Proteins" = alphaDIA$pg, # Good idea to check how reliable this is using ProtMatch
                   "Sequence" = alphaDIA$sequence,
                   "Length" = nchar(alphaDIA$sequence),
                   "Charge" = alphaDIA$charge,
                   "PEP" = alphaDIA$proba, 
                   check.names = FALSE)
  #
  # Process modifications
  # Do not detect them from the log.txt only, because they can also come with the library!
  # Better to just get what is in the file.
  mods <- unique(unlist(strsplit(alphaDIA$mods, ";")))
  Mods <- as.data.frame(t(sapply(strsplit(mods, "@"), unlist)))
  colnames(Mods) <- c("Full name", "Position")
  Mods$"Old mark" <- mods
  Mods$AA <- Mods$Site <- Mods$Position
  w <- which(Mods$Position %in% paste0("Any", c("", " ", "_", "-"), "N-term"))
  if (length(w)) {
    Mods$Site[w] <- "n"
    Mods$AA[w] <- "_"
    Mods$Position[w] <- "Any N-term"
  }
  w <- which(Mods$Position %in% paste0("Any", c("", " ", "_", "-"), "C-term"))
  if (length(w)) {
    Mods$Site[w] <- "c"
    Mods$AA[w] <- "_"
    Mods$Position[w] <- "Any C-term"
  }
  w <- which(Mods$Position %in% paste0("Protein", c("", " ", "_", "-"), "N-term"))
  if (length(w)) {
    Mods$Site[w] <- "[^"
    Mods$AA[w] <- "_"
    Mods$Position[w] <- "Protein N-term"
  }
  w <- which(Mods$Position %in% paste0("Protein", c("", " ", "_", "-"), "C-term"))
  if (length(w)) {
    Mods$Site[w] <- "^]"
    Mods$AA[w] <- "_"
    Mods$Position[w] <- "Protein C-term"
  }
  allPTMs <- aggregate(Mods[, c("Position", "AA", "Site", "Old mark")], list(Mods$"Full name"), list)
  colnames(allPTMs)[1] <- "Full name"
  allPTMs$Type <- c("Variable", "Fixed")[(allPTMs$`Full name` %in% Fixed.mods)+1]
  allPTMs$UniMod <- paste0("UniMod:", UniMod$UnimodId[match(allPTMs$`Full name`, UniMod$Name)])
  allPTMs$Mark <- tolower(substr(allPTMs$"Full name", 1, 2))
  # Sometimes, some marks are duplicates, e.g. if you searched for "Acetyl (Protein N-term)" and "Acetyl (K)" together!
  # We want to fix this so that each modification has a unique mark:
  tstMark <- aggregate(allPTMs$Mark, list(allPTMs$Mark), length)
  W <- which(tstMark$x > 1)
  #allPTMs$Mark <- allPTMs$"tmp mark"
  if (length(W)) {
    allPTMs$"tmp mark" <- allPTMs$Mark
    for (i in W) {
      #i <- W[1]
      w <- which(allPTMs$Mark == tstMark$Group.1[i])
      m <- allPTMs[w,]
      m$AA[which(sapply(m$AA, length) == 0)] <- "X"
      if ("Acetyl" %in% m$"Full name") { r <- which(m$"Full name" == "Acetyl") } else { r <- 1 }
      s <- c(1:nrow(m)); s <- s[which(s != r)]
      test <- apply(m[s, c("AA", "Mark")], 1, function(x) { paste0(tolower(x[[1]]), substr(x[[2]], 1, 1))[1] })
      w0 <- which(!test %in% allPTMs$Mark)
      if (length(w0)) {
        m$Mark[s][w0] <- test[w0[1]]
      }
      w1 <- which(test %in% allPTMs$Mark)
      if (length(w1)) {
        # not tested
        s <- s[w1]
        test <- sapply(s, list)
        kount <- 1
        char <- c(0:9, letters)
        taken <- unique(c(allPTMs$Mark, m$Mark))
        for (j in s) {
          tst <- paste0(tolower(m$AA[s]), char[kount])
          while (((tst) %in% taken)&&(kount < length(char))) {
            kount <- kount+1
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
  allPTMs$`tmp mark` <- NULL
  #
  allPTMs$"Mass shift" <- UniMod$MonoMass[match(gsub("^UniMod:", "", allPTMs$UniMod), UniMod$UnimodId)]
  if ("Mass delta" %in% colnames(allPTMs)) {
    tst <- allPTMs$"Mass shift" - allPTMs$"Mass delta"
    stopifnot(max(tst, na.rm = TRUE) < 0.1)
  }
  # Modified sequence
  tmpMds <- proteoCraft::listMelt(setNames(allPTMs$`Old mark`, allPTMs$Mark))
  tmpMds$Name <- Mods$`Full name`[match(tmpMds$value, Mods$`Old mark`)]
  parallel::clusterExport(cl, "tmpMds", envir = environment())
  wMod <- which(alphaDIA$mods != "")
  tmp1 <- alphaDIA[wMod, c("sequence", "mods", "mod_sites")]
  tmp1$sequence <- strsplit(tmp1$sequence, "")
  tmp1$mods <- strsplit(tmp1$mods, ";")
  tmp1$mod_sites <- strsplit(tmp1$mod_sites, ";")
  tmp1$mod_sites <- lapply(tmp1$mod_sites, as.integer)
  saveRDS(tmp1, "tmp1.RDS")
  invisible(parallel::clusterCall(cl, function(x) {
    tmp1 <<- readRDS("tmp1.RDS")
  }))
  unlink("tmp1.RDS")
  f0 <- function(x) {
    #x <- tmp1[1,]
    sq1 <- sq2 <- unlist(x[[1]])
    mds <- unlist(x[[2]])
    pos <- unlist(x[[3]])
    m <- match(mds, tmpMds$value)
    mds1 <- tmpMds$L1[m]
    mds2 <- tmpMds$Name[m]
    sq1[pos] <- paste0(sq1[pos], "(", mds1, ")")
    sq2[pos] <- paste0(sq2[pos], "(", mds2, ")")
    sq1 <- paste(sq1, collapse = "")
    sq2 <- paste(sq2, collapse = "")
    return(c(sq1, sq2))
  }
  environment(f0) <- .GlobalEnv
  EV$"Modified sequence_verbose" <- EV$"Modified sequence" <- EV$Sequence
  tmp2 <- as.data.frame(t(parallel::parApply(cl, tmp1, 1, f0)))
  #View(cbind(tmp1, tmp2))
  EV$"Modified sequence"[wMod] <- tmp2[, 1]
  EV$"Modified sequence_verbose"[wMod] <- tmp2[, 2]
  EV$"Modified sequence" <- paste0("_", EV$"Modified sequence", "_")
  EV$"Modified sequence_verbose" <- paste0("_", EV$"Modified sequence_verbose", "_")
  #
  EV$Modifications <- alphaDIA$mods
  EV$Modification_sites <- alphaDIA$mod_sites
  #
  # M/Z
  #cat("   Calculating theoretical m/z values...\n")
  EV$"m/z" <- alphaDIA$mz_calibrated
  EV$"Uncalibrated m/z" <- alphaDIA$mz_observed
  EV$"m/z (library)" <- alphaDIA$mz_library
  parallel::clusterExport(cl, "allPTMs", envir = environment())
  tmp0 <- strsplit(unique(EV$Modifications[wMod]), ";")
  mdShft0 <- vapply(tmp0, function(x) {
    sum(allPTMs$`Mass shift`[match(Mods$`Full name`[match(x, Mods$`Old mark`)], allPTMs$`Full name`)]) # No na.rm = TRUE: if we fail to calculate, I want to see it
  }, 1)
  mdShft <- mdShft0[match(EV$Modifications[wMod], tmp0)]
  tmp1 <- do.call(paste, c(EV[, c("Sequence", "Charge")], sep = "___"))
  tmp2 <- unique(tmp1)
  tmp3 <- EV[match(tmp2, tmp1), c("Sequence", "Charge")]                              
  f0 <- function(x) { Peptides::mz(x[[1]], as.integer(x[[2]], cysteins = 0)) }
  environment(f0) <- .GlobalEnv
  tmp3 <- parallel::parApply(cl, tmp3, 1, f0) # Peptides::mz doesn't like vectorization clearly
  EV$"Theoretical m/z" <- tmp3[match(tmp1, tmp2)] 
  EV$"Theoretical m/z"[wMod] <- EV$"Theoretical m/z"[wMod] + mdShft
  #
  EV$Mass <- (EV$"m/z"-1.007276466879)*EV$Charge
  #
  # Retention time (in minutes, as in MaxQuant)
  EV$"Retention time" <- as.numeric(alphaDIA$rt_observed)/60
  EV$"Calibrated retention time" <- as.numeric(alphaDIA$rt_calibrated)/60
  EV$"Retention length" <- as.numeric(alphaDIA$base_width_rt)/60
  EV$"Retention time (library)" <- as.numeric(alphaDIA$rt_library)/60
  #
  # Ion mobility
  EV$IM <- as.numeric(alphaDIA$mobility_observed)
  EV$"Calibrated IM" <- as.numeric(alphaDIA$mobility_calibrated)
  EV$"IM width" <- as.numeric(alphaDIA$base_width_mobility)
  EV$"IM (library)" <- as.numeric(alphaDIA$mobility_library)
  #
  # Decoys
  if ("decoy" %in% colnames(alphaDIA)) {
    EV$Reverse <- c("", "+")[alphaDIA$decoy+1]
  } else {
    EV$Reverse <- ""
  }
  #
  # Score
  EV$Score <- as.numeric(alphaDIA$score)
  #
  # Q-values
  EV$"Q-value" <- as.numeric(alphaDIA$qval)
  EV$"PG q-value" <- as.numeric(alphaDIA$pg_qval)
  #
  # Quantitative values
  EV$"Intensity" <- as.numeric(alphaDIA$intensity)
  #
  EV$Type <- "LIB-DIA" # This may need to change if important.
  #
  digPat <- paste(c("[", digPattern, "]"), collapse = "")
  tmp <- gsub(paste0(digPat, "$"), "", EV$Sequence)
  EV$"Missed cleavages" <- nchar(tmp) - nchar(gsub(digPat, "", tmp))
  EV$"Potential contaminant" <- ""
  #
  if ("Search_ID" %in% colnames(alphaDIA)) { # Keeping this for now eventhough we do not have such a script for alphaDIA yet 
    EV$Search_ID <- alphaDIA$PSMs_file # In cases where the input is a hybrid report created by the merging script in .../Utils
  } else {
    EV$Search_ID <- sort(alphaDIA_fl, decreasing = TRUE)[1] # Provide default: there can be one or 2 files, the .tsv or .parquet, we want just one (ideally the .tsv)
  }
  #
  EV$id <- NULL
  EV <- cbind(data.frame(id = 1:nrow(EV)),
              EV)
  Res <- list(Evidence = EV,
              PTMs = allPTMs)
  #
  invisible(parallel::clusterCall(cl, function(x) {
    try(rm(tmp1), silent = TRUE)
    return()
  }))
  if (stopCl) { parallel::stopCluster(cl) }
  return(Res)
}
