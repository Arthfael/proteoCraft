#' DIANN_to_MQ
#'
#' @description 
#' Converts a diaNN PSMs table to a MaxQuant evidence.txt-like table.
#' As with the other functions of this type, the idea is not to get perfect conversion, but close enough that the data can be fed into this package's analysis scripts.
#' 
#' The output is a list of 3 items:
#' - Evidence: a data.frame similar to the evidence.txt table MaxQuant creates
#' - PTMs: a modifications table
#' - QuantUMS: was QuantUMS quantification on?
#' 
#' Surprisingly, the diaNN PSMs table does not include MZ values! Is there a way to do this by adding a command line argument?
#' For now they are either detected from the library file (same folder) or computed.
#' Cysteine alkylation is assumed to be Carbamidomethyl
#' 
#' @param DIANN_fl External diaNN PSMs table file (.tsv, .parquet format, or both).
#' @param Fixed.mods A vector of those Fixed modifications which will not be reported in the results (reported by diaNN, but MQ ignores them). Default = "Carbamidomethyl"
#' @param PG.Genes If TRUE, will include per Protein Group/Gene quantitative and Q-values columns. Default = FALSE, because we will mostly be using this function in contexts where we want to perform protein groups inference re-calculate quantitative values ourselves.
#' @param Detect.Lib diaNN does not output M/Z values, but those are in the library generated at the end of some runs. If this is set to TRUE, it will try to detect the library and if available get M/Z values from there. Otherwise they will be computed.
#' @param Min.Delta.Score From which minimum delta score should we report a PSM? Default = -Inf (we assume that those few peptides with low delta score are because randomly the decoy database will include real-like peptides)
#' @param N.clust A limit on the number of vCPUs to use. If left as NULL (default), uses the number of available clusters - 1, to a minimum of 1.
#' @param N.reserved Default = 1. Number of reserved vCPUs the function is not to use. Note that for obvious reasons it will always use at least one.
#' @param cl Already have a cluster handy? Why spend time making a new one, which on top of that may invalidate the old one. Just pass it along!
#' @param shinyOpt One of "popup" (default), "pane", "dialog", or "browser".
#' 
#' @examples
#' DIANN_fl <- rstudioapi::selectFile()
#' DIANN_fl <- Param$PSMs
#' temp <- DIANN_to_MQ(DIANN_fl)
#' ev <- temp$Evidences
#' Modifs <- temp$PTMs
#' 
#' @import data.table
#' @export

DIANN_to_MQ <- function(DIANN_fl,
                        Fixed.mods = c("Carbamidomethyl"),
                        PG.Genes = FALSE,
                        Detect.Lib = TRUE,
                        Min.Delta.Score = -Inf,
                        N.clust,
                        N.reserved = 1,
                        cl,
                        shinyOpt = "popup") {
  TESTING <- FALSE
  #proteoCraft::DefArg(proteoCraft::DIANN_to_MQ);TESTING <- TRUE
  #DIANN_fl <- PSMsFls; cl <- parClust
  #DIANN_fl <- rstudioapi::selectFile(path = paste0(getwd(), "/*.tsv"))
  #DIANN_fl <- diannRep
  shinyOpt <- tolower(shinyOpt)
  if (!shinyOpt %in% c("popup", "dialog", "pane", "browser")) { shinyOpt <- "popup" }
  #
  # NB:
  # This function deliberately uses a different object name for the PTMs table to other related functions,
  # to simplify environment-related debugging.
  # Do not change that!
  # Here, use allPTMs, elsewhere Modifs (or anything else).
  if (TESTING) {
    tm1 <<- Sys.time()
    # Note:
    # This is not a perfect alternative to missing but will work in most cases, unless x matches a function imported by a package 
    misFun <- function(x) { return(!exists(deparse(substitute(x)))) }
  } else { misFun <- missing }
  cleanUp <- FALSE
  if (misFun(cl)) {
    dc <- parallel::detectCores()
    if (misFun(N.reserved)) { N.reserved <- 1 }
    if (misFun(N.clust)) {
      N.clust <- max(c(dc-N.reserved, 1))
    } else {
      if (N.clust > max(c(dc-N.reserved, 1))) {
        warning("   More cores specified than allowed, I will ignore the specified number! You should always leave at least one free for other processes, see the \"N.reserved\" argument.")
        N.clust <- max(c(dc-N.reserved, 1))
      }
    }
    cl <- parallel::makeCluster(N.clust, type = "SOCK")
    cleanUp <- TRUE
  }
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
  # Read PSMs
  wTsv <- grep("\\.tsv$", DIANN_fl)
  wPrqt <- grep("\\.parquet$", DIANN_fl)
  lTsv <- length(wTsv)
  lPrqt <- length(wPrqt)
  stopifnot(lTsv <= 1,
            lPrqt <= 1,
            sum(lTsv, lPrqt) > 0)
  if (lTsv) {
    DIANN <- data.table::fread(DIANN_fl[wTsv], integer64 = "numeric", check.names = FALSE, data.table = FALSE)
  }
  if (lPrqt) {
    DIANN_prqt <- arrow::read_parquet(DIANN_fl[wPrqt])
    if (lTsv) {
      k2 <- colnames(DIANN_prqt)
      DIANN[, k2] <- DIANN_prqt[, k2]
      # The reason for overwriting the columns from the tsv is:
      # - The order is the same
      # - The data, if a column is present in both forms of the report, is either identical or, where there is a discrepancy, I have 1 more digit with the .parquet file
    } else {
      DIANN <- DIANN_prqt
    }
  }
  # Filter peptides with ambiguous/non-classical amino acids 
  uSeq <- unique(DIANN$Stripped.Sequence)
  test <- gsub(paste(proteoCraft::AA, collapse = "|"), "", uSeq)
  test <- test[match(DIANN$Stripped.Sequence, uSeq)]
  w <- which(test == "")
  lX <- nrow(DIANN)-length(w) 
  if (lX) {
    u <- unique(unlist(strsplit(DIANN$Stripped.Sequence, "")))
    u <- paste(u[which(!u %in% AA)], collapse = "-")
    warning(paste0("Removing ", lX, " peptide", c("", "s")[(lX > 1)+1], " with unknown amino acids ", u, "!"))
    DIANN <- DIANN[w,]
  }
  #
  log_Fl <- unique(gsub("\\.((tsv)|(parquet))$", ".log.txt", DIANN_fl))
  if (file.exists(log_Fl)) { 
    log <- readLines(log_Fl)
    diannCall <- grep(" --[a-zA-Z]", log, value = TRUE)[1]
    diannCall <- unlist(strsplit(diannCall, " +--"))
    diannCall <- diannCall[2:length(diannCall)]
  }
  #
  libTst <- try({
    libKol <- c()
    if (Detect.Lib) {
      if (file.exists(log_Fl)) {
        DIANNLib_fl <- gsub(".+\\] Saving the library to |\\.speclib$", "", grep("\\] Saving the library to ", log, value = TRUE))
        x <- gsub("^lib ", "", grep("^lib ", diannCall, value = TRUE))
        x <- x[which(nchar(x) > 0)]
        if (length(x)) { diannCall <- c(diannCall, x) }
        DIANNLib_fl <- gsub("\\\\", "/", DIANNLib_fl)
        DIANNLib_fl <- unique(DIANNLib_fl[which(DIANNLib_fl != "")])
        DIANNLib_fl <- unique(c(DIANNLib_fl, gsub("\\.skyline$", "", DIANNLib_fl)))
        if (length(DIANNLib_fl)) {
          if (sum(dirname(DIANNLib_fl) != dirname(DIANN_fl))) {
            DIANNLib_fl <- c(DIANNLib_fl, paste0(dirname(DIANN_fl), "/", basename(DIANNLib_fl)))
          }
          DIANNLib_fl <- as.character(sapply(DIANNLib_fl, function(x) { paste0(x, c(".tsv", ".parquet", "")) }))
          DIANNLib_fl <- DIANNLib_fl[which(file.exists(DIANNLib_fl))]
          if (length(DIANNLib_fl)) {
            DIANNLib_fl <- DIANNLib_fl[1]
            if (grepl("\\parquet$", DIANNLib_fl)) {
              tst <- try({ DIANNLib <- arrow::read_parquet(DIANNLib_fl) }, silent = TRUE)
              if ("try-error" %in% class(tst)) { # Should never happen...
                tst <- try({ DIANNLib <- data.table::fread(DIANNLib_fl, integer64 = "numeric", check.names = FALSE, data.table = FALSE) }, silent = TRUE)
              }
            } else {
              tst <- try({ DIANNLib <- data.table::fread(DIANNLib_fl, integer64 = "numeric", check.names = FALSE, data.table = FALSE) }, silent = TRUE)
              if ("try-error" %in% class(tst)) {
                tst <- try({ DIANNLib <- arrow::read_parquet(DIANNLib_fl) }, silent = TRUE)
              }
            }
            if ("transition_group_id" %in% colnames(DIANNLib)) {
              DIANNLib$transition_group_id <- gsub("\\[", "(", gsub("\\]", ")", DIANNLib$transition_group_id))
            } else {
              DIANNLib$ModifiedPeptideSequence <- gsub("\\[", "(", gsub("\\]", ")", DIANNLib$ModifiedPeptideSequence))
              DIANNLib$transition_group_id <- do.call(paste, c(DIANNLib[, c("ModifiedPeptideSequence", "PrecursorCharge")], sep = "+"))
            }
            #a <- unique(unlist(strsplit(gsub("^[A-Z]*\\(|\\)[A-Z]*$|^[A-Z]+$", "", gsub("\\)[A-Z]*\\(", "_", DIANNLib$ModifiedPeptideSequence)), "_")))
            tmp <- do.call(paste, c(DIANN[, c("Modified.Sequence", "Precursor.Charge")], sep = "+"))
            m <- match(tmp, DIANNLib$transition_group_id)
            DIANN$"m/z" <- DIANNLib$PrecursorMz[m]
            # Also library index + rt
            DIANN$"Library index" <- DIANNLib$transition_group_id[m]
            DIANN$"Library rt" <- NA
            if ("RetentionTime" %in% colnames(DIANNLib)) { DIANN$"Library rt" <- DIANNLib$RetentionTime[m] } else {
              if ("AverageExperimentalRetentionTime" %in% colnames(DIANNLib)) { DIANN$"Library rt" <- DIANNLib$AverageExperimentalRetentionTime[m] } else {
                if ("NormalizedRetentionTime" %in% colnames(DIANNLib)) { DIANN$"Library rt" <- DIANNLib$NormalizedRetentionTime[m] }
              }
            }
            libKol <- c("m/z", "Library index", "Library rt")
          } else { cat("   (NB: no diaNN library tsv found in folder)\n") }
        } else { cat("   (NB: no diaNN library tsv found in folder)\n") }
      } else { cat("   (NB: no diaNN log file found in folder)\n") }
    }
  }, silent = TRUE)
  #
  # Create MQ-like file
  if (!"File.Name" %in% colnames(DIANN)) {
    stopifnot("Run.Index" %in% colnames(DIANN),
              file.exists(log_Fl))
    tst <- aggregate(DIANN$Run, list(DIANN$Run.Index), unique)
    fls <- gsub("(\\\\)+", "/", gsub("^f ", "", grep("^f ", diannCall, value = TRUE)))
    flNms <- gsub(".*/|\\.[^\\.]+$", "", fls)
    tst$x2 <- flNms[tst$Group.1+1]
    stopifnot(sum(tst$x != tst$x2) == 0)
    tst$File <- fls[tst$Group.1+1]
    DIANN$File.Name <- tst$File[match(DIANN$Run.Index, tst$Group.1)]
  }
  kol <- c("Run", "File.Name", "Protein.Ids", "Protein.Names", "Proteotypic", "Genes",
           "Stripped.Sequence", "Modified.Sequence", "Precursor.Charge", "PEP")
  #print(kol[which(!kol %in% colnames(DIANN))])
  stopifnot(sum(!kol %in% colnames(DIANN)) == 0)
  EV <- data.frame("Raw file" = DIANN$Run,
                   "Raw file path" = DIANN$File.Name,
                   "Proteins" = DIANN$Protein.Ids, # Good idea to check how reliable this is using ProtMatch
                   "Full Protein ID" = DIANN$Protein.Names, # Same as above
                   "Proteotypic" = as.logical(DIANN$Proteotypic),
                   "Gene names" = DIANN$Genes,
                   "Sequence" = DIANN$Stripped.Sequence,
                   "Mod. seq. (DiaNN format)" = DIANN$Modified.Sequence, # Useful for matching to other DiaNN outputs!
                   "Modified sequence" = paste0("_", DIANN$Modified.Sequence, "_"), # Will need some reprocessing further down
                   "Length" = nchar(DIANN$Stripped.Sequence),
                   "Charge" = DIANN$Precursor.Charge,
                   "PEP" = DIANN$PEP,
                   check.names = FALSE)
  #
  if (!"try-error" %in% class(libTst)) {
    if (length(libKol)) { EV[, libKol] <- DIANN[, libKol] }
  }
  if ("First.Protein.Description" %in% colnames(DIANN)) {
    EV$"First protein name" <- DIANN$First.Protein.Description
  } else {
    cat("   Column \"First.Protein.Description\" not found, no column of human-readable names generated.\n")
  }
  # Process modifications
  wMod <- grep("\\(", EV$`Modified sequence`)
  mWMod <- match(unique(EV$`Modified sequence`[wMod]), EV$`Modified sequence`)
  tmp <- data.frame(Seq = EV$Sequence[mWMod],
                    ModSeq = EV$`Modified sequence`[mWMod],
                    Z = EV$Charge[mWMod])
  f0 <- function(x) { Peptides::mz(x[[1]], as.integer(x[[2]]), cysteins = 0) }
  environment(f0) <- .GlobalEnv
  tmp$Theoretical <- parallel::parApply(cl, tmp[, c("Seq", "Z")], 1, f0)
  tmp$Comp <- gsub("\\)[A-Z]*\\(", "___",
                   gsub("^_[A-Z]*\\(|\\)[A-Z]*_$|^_[A-Z]+_$", "",
                        gsub("\\[", "(",
                             gsub("\\]", ")", tmp$ModSeq)))) 
  tmp$Comp <- strsplit(tmp$Comp, "___")
  comps <- unique(unlist(tmp$Comp))
  gC <- sum(grepl("^UniMod:", comps))
  l <- length(comps)
  modsType <- "No idea"
  if (gC == l) { modsType <- "UniMod" }
  if (modsType == "No idea") {
    if (shinyOpt == "popup") {
      runApp1 <- "print(shiny::shinyApp(proteoCraft::DIANN_to_MQ_ui1(tstMap), proteoCraft::DIANN_to_MQ_server1, options = list(height = screenRes$height, width = screenRes$width)))"
      runApp2 <- "print(shiny::shinyApp(proteoCraft::DIANN_to_MQ_ui2, proteoCraft::DIANN_to_MQ_server2, options = list(height = screenRes$height, width = screenRes$width)))"
    }
    if (shinyOpt %in% c("dialog", "pane")) {
      if (shinyOpt == "pane") {
        myViewer <- shiny::paneViewer("Viewer", minHeight = "maximize")
      }
      if (shinyOpt == "dialog") {
        myViewer <- shiny::dialogViewer("Viewer", width = screenRes$width, height = screenRes$height)
      }
      runApp1 <- "shiny::runGadget(shiny::shinyApp(proteoCraft::DIANN_to_MQ_ui1(tstMap), proteoCraft::DIANN_to_MQ_server1, options = list(height = screenRes$height, width = screenRes$width)), viewer = myViewer, stopOnCancel = FALSE)"
      runApp2 <- "shiny::runGadget(shiny::shinyApp(proteoCraft::DIANN_to_MQ_ui2, proteoCraft::DIANN_to_MQ_server2, options = list(height = screenRes$height, width = screenRes$width)), viewer = myViewer, stopOnCancel = FALSE)"
    }
    if (shinyOpt == "browser") {
      runApp1 <- "print(shiny::shinyApp(proteoCraft::DIANN_to_MQ_ui1(tstMap), proteoCraft::DIANN_to_MQ_server1, options = list(height = \"100%\", width = \"100%\", launch.browser = TRUE)))"
      runApp2 <- "print(shiny::shinyApp(proteoCraft::DIANN_to_MQ_ui2, proteoCraft::DIANN_to_MQ_server2, options = list(height = \"100%\", width = \"100%\", launch.browser = TRUE)))"
    }
    if (!"m/z" %in% colnames(EV)) {
      msg <- "This DIANN run was searched with a library with non standard PTMs.\nCurrently to make sense of this we need PSM m/z, but these are not written by DiaNN, and they could not be obtained from the library since it could not be found..."
      warning(msg)
      modsType <- "StrangerThings"
      tstMap <- FALSE
    } else {
      tmp$MZ <- EV$"m/z"[mWMod]
      tmp$DM <- (tmp$MZ-tmp$Theoretical)*tmp$Z
      screenRes <- rpanel::rp.screenresolution()
      mrks <- paste0("- ", sort(comps))
      # Commented because of issues with shinyApps
      #tstMap <- TRUE
      #eval(parse(text = runApp1))
      #tstMap <- c(FALSE, TRUE)[match(tstMap, c("no", "yes"))]
      msg <- paste0("Non-classical PTMs marks detected (e.g. \"", sort(comps)[1], "\")...")
      opt <- c("Keep them as they are                                                                                                                 ",
               "Manually map each to a more convenient mark                                                                                           ")
      tstMap <- c(FALSE, TRUE)[match(dlg_list(opt, opt[1], title = msg)$res, opt)]
    }
    if (tstMap) {
      # Here we will want to make sense of those mods
      # First, get the precise mass shift for each, so we can narrow down options
      f0 <- function(x) {
        #x <- tmp$ModSeq[1]
        x <- proteoCraft::annot_to_tabl(x)[[1]]
        w <- which(x$Annotations != "")
        x$Sequence[2] <- paste0("_", x$Sequence[2])
        do.call(paste, c(x[w, , drop = FALSE], sep = "mod_"))
      }
      environment(f0) <- .GlobalEnv
      tmp$Mods <- parallel::parLapply(cl, tmp$ModSeq, f0)
      Mods <- sort(unique(unlist(tmp$Mods)))
      parallel::clusterExport(cl, "Mods", envir = environment())
      f0 <- function(x) {
        x <- aggregate(x, list(x), length)
        res <- setNames(rep(0, length(Mods)), Mods)
        w <- which(Mods %in% x$Group.1)
        m <- match(Mods[w], x$Group.1)
        res[w] <- x$x[m]
        return(res)
      }
      environment(f0) <- .GlobalEnv
      d.matr <- t(parallel::parSapply(cl, tmp$Mods, f0))
      fit <- lm.fit(d.matr, tmp$DM)
      Mods <- data.frame(Name = Mods,
                         DM = fit$coefficients)
      Mods$LibDM <- suppressWarnings(as.numeric(gsub(".*mod_", "", Mods$Name)))
      Mods$AA <- gsub("mod_.*", "", Mods$Name)
      Mods$Options <- apply(Mods[, c("DM", "AA")], 1, function(x) {
        dm <- as.numeric(x[[1]])
        if (substr(x[[2]], 1, 1) == "_") {
          w <- which((abs(UniMod$MonoMass - dm) <= 0.02)&(UniMod$Site %in% c(substr(x[[2]], 2, 2), "N-term")))
        } else {
          w <- which((abs(UniMod$MonoMass - dm) <= 0.02)&(UniMod$Site == x[[2]]))
        }
        if (length(w)) {
          Ddm <- dm - UniMod$MonoMass[w]
          w <- w[order(Ddm, decreasing = FALSE)]
        }
        return(w)
      })
      Mods$tag <- do.call(paste, c(data.frame(AA = gsub("^_", "", Mods$AA),
                                              Name = paste0("+", as.character(round(Mods$DM), 4)),
                                              NTerm = c("", "(N-term)")[grepl("^_", Mods$AA)+1]),
                                   sep = " "))
      opt <- sort(unique(unlist(Mods$Options)))
      UniMod2 <- UniMod[opt,]
      UniMod2$Option <- opt
      Mods <- Mods[order(Mods$DM, Mods$tag),]
      commonPTMs <- c("Acetyl:N-term",
                      "Acetyl:K",
                      "Carbamidomethyl:C",
                      "Deamidated:N",
                      "Deamidated:Q",
                      "Gln->pyro-Glu:Q",
                      "Oxidation:M",
                      "Phospho:S",
                      "Phospho:T",
                      "Phospho:Y",
                      "Methyl:K",
                      "Dimethyl:K",
                      "Trimethyl:K",
                      "Methyl:R",
                      "Dimethyl:R",
                      "Trimethyl:R",
                      "GG:K")
      Mods$Options2 <- lapply(1:nrow(Mods), function(i) {
        m <- match(Mods$Options[[i]], UniMod2$Option)
        Ddm <- Mods$DM[i] - UniMod$MonoMass[m]
        opt <- UniMod2$Id[m]
        return(opt[order(Ddm, decreasing = FALSE)])
      })
      Mods$Match <- sapply(1:nrow(Mods), function(i) {
        opt <- Mods$Options2[[i]]
        if (length(opt) >= 1) {
          dflt <- opt[1]
          m <- match(commonPTMs, opt)
          m <- m[which(!is.na(m))]
          if (length(m)) { dflt <- opt[m[1]] }
        } else { dflt <- Mods$tag[i] }
        return(dflt)
      })
      homePath <- paste0(normalizePath(Sys.getenv("HOME"), winslash = "/"), "/R/proteoCraft")
      fl <- paste0(homePath, "/tmpMods.tmp")
      file.create(fl)
      con <- file(fl, "r+")
      serialize(Mods, con)
      suppressWarnings(suppressMessages(try(close(con), silent = TRUE))) # Only necessary in function mode
      eval(parse(text = runApp2))
      #print(shiny::shinyApp(DIANN_to_MQ_ui2, DIANN_to_MQ_server2, options = list(height = screenRes$height, width = screenRes$width)))
      con <- file(fl, "r")
      Mods <- unserialize(con)
      suppressWarnings(suppressMessages(try(close(con), silent = TRUE))) # Only necessary in function mode
      unlink(fl)
      #
      allPTMs <- aggregate(Mods$DM, list(Mods$Match), function(x) {
        rg <- c(min(x), max(x))
        stopifnot(rg[2]-rg[1] < 0.010)
        return(mean(x))
      })
      colnames(allPTMs) <- c("Full name", "Mass delta")
      allPTMs$Position <- "Anywhere"
      allPTMs$Position[grep("C-term", allPTMs$"Full name", ignore.case = TRUE)] <- "Any C-term"
      allPTMs$Position[grep("N-term", allPTMs$"Full name", ignore.case = TRUE)] <- "Any N-term"
      allPTMs$Position[grep("protein C-term", allPTMs$"Full name", ignore.case = TRUE)] <- "Protein C-term"
      allPTMs$Position[grep("protein N-term", allPTMs$"Full name", ignore.case = TRUE)] <- "Protein N-term"
    }
  }
  if (modsType == "UniMod") {
    allPTMs <- data.frame("UniMod" = comps,
                          check.names = FALSE)
    allPTMs$"Full name" <- UniMod$Name[match(gsub("UniMod:", "", allPTMs$UniMod),
                                             UniMod$UnimodId)]
  }
  if (modsType == "StrangerThings") {
    allPTMs <- data.frame("Full name" = comps,
                          check.names = FALSE)
    allPTMs$UniMod <- NA
  }
  allPTMs$Type <- "Variable"
  allPTMs$Type[which(allPTMs$"Full name" %in% Fixed.mods)] <- "Fixed"
  allPTMs$Mark <- tolower(substr(allPTMs$"Full name", 1, 2))  
  ## Identify affected Amino Acids
  if (modsType == "UniMod") {
    allPTMs$AA <- list(NA)
    uMdSq <- unique(EV$`Modified sequence`)
    for (i in 1:nrow(allPTMs)) { #i <- 1
      uMdSq <- unique(EV$`Modified sequence`)
      temp <- gsub(proteoCraft::topattern(paste0("(", allPTMs$UniMod[i], ")"), start = FALSE),
                   ">>>.___", uMdSq)
      temp2 <- grep(">>>$", unlist(strsplit(temp, "\\.___")), value = TRUE)
      nc <- nchar(temp2)
      allPTMs$AA[[i]] <- unique(substr(temp2, nc-3, nc-3))
    }
  }
  if (modsType == "StrangerThings") {
    allPTMs$AA <- list(NA)
    for (i in 1:nrow(allPTMs)) { #i <- 1
      temp <- gsub(proteoCraft::topattern(paste0("(", allPTMs$`Full name`[i], ")"), start = FALSE),
                   ">>>.___", EV$`Modified sequence`)
      temp2 <- grep(">>>$", unlist(strsplit(temp, "\\.___")), value = TRUE)
      nc <- nchar(temp2)
      allPTMs$AA[[i]] <- unique(substr(temp2, nc-3, nc-3))
    }
  }
  if (modsType == "No idea") {
    tmp <- aggregate(gsub("^_", "", Mods$AA), list(Mods$Match), function(x) { sort(unique(x)) })
    allPTMs$Site <- allPTMs$AA <- tmp$x[match(allPTMs$`Full name`, tmp$Group.1)]
    w <- which(allPTMs$Position == "Any N-term")
    if (length(w)) { allPTMs$Site[w] <- sapply(allPTMs$Site[w], function(x) { paste0("n", x) }) }
    w <- which(allPTMs$Position == "Any C-term")
    if (length(w)) { allPTMs$Site[w] <- sapply(allPTMs$Site[w], function(x) { paste0("c", x) }) }
    w <- which(allPTMs$Position == "Protein N-term")
    if (length(w)) { allPTMs$Site[w] <- "[^" }
    w <- which(allPTMs$Position == "Protein C-term")
    if (length(w)) { allPTMs$Site[w] <- "^]" }
    allPTMs$Site_long <- sapply(strsplit(gsub("\\^", "", sapply(allPTMs$Site, paste, collapse = "")), ""), function(x) {
      x[which(x %in% c("[", "n"))] <- "N-term"
      x[which(x %in% c("]", "c"))] <- "C-term"
      return(x)
    })
  }
  ## Sometimes, some marks are duplicates, e.g. if you searched for "Acetyl (Protein N-term)" and "Acetyl (K)" together!
  ## We want to fix this so that each modification has a unique mark:
  tstMark <- aggregate(allPTMs$Mark, list(allPTMs$Mark), length)
  W <- which(tstMark$x > 1)
  #allPTMs$Mark <- allPTMs$"Old mark"
  if (length(W)) {
    allPTMs$"Old mark" <- allPTMs$Mark
    for (i in W) {
      #i <- W[1]
      w <- which(allPTMs$Mark == tstMark$Group.1[i])
      m <- allPTMs[w,]
      m$AA[which(sapply(m$AA, length) == 0)] <- "X"
      if ("Acetyl" %in% m$"Full name") { r <- which(m$"Full name" == "Acetyl") } else { r <- 1 }
      s <- c(1:nrow(m)); s <- s[which(s != r)]
      test <- apply(m[s, c("AA", "Mark")], 1, function(x) { paste0(tolower(x[[1]]), substr(x[[2]], 1, 1)) })
      w1 <- which(!test %in% allPTMs$Mark)
      m$Mark[s][w1] <- test
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
  if (modsType == "No idea") {
    tmp <- aggregate(gsub(".*mod_", "", Mods$Name), list(Mods$Match), function(x) { sort(unique(x)) })
    allPTMs$"Old mark" <- tmp$x[match(allPTMs$`Full name`, tmp$Group.1)] # Could be a list!!!
  }
  if (!"UniMod" %in% colnames(allPTMs)) {
    allPTMs$UniMod <- UniMod$UnimodId[match(allPTMs$`Full name`, UniMod$Name)] # Is not meant to always work
    # Sometimes there will be discrepancies, where UniMod has a compatible isobaric PTM with a different name!
  }
  if (!"Mass shift" %in% colnames(allPTMs)) {
    allPTMs$"Mass shift" <- UniMod$MonoMass[match(gsub("^UniMod:", "", allPTMs$UniMod), UniMod$UnimodId)]
    if ("Mass delta" %in% colnames(allPTMs)) {
      tst <- allPTMs$"Mass shift" - allPTMs$"Mass delta"
      stopifnot(max(tst, na.rm = TRUE) < 0.1)
    }
  }
  #
  if (modsType == "UniMod") {
    tmpMds <- proteoCraft::listMelt(allPTMs$AA, 1:nrow(allPTMs), c("AA", "Row"))
    kol <- c("UniMod", "Full name", "Mark")
    #kol %in% colnames(allPTMs)
    tmpMds[, c("UniMod", "Name", "Mark")] <- allPTMs[tmpMds$Row, kol]
    w <- grep("^UniMod:", tmpMds$UniMod, invert = TRUE)
    if (length(w)) {
      tmpMds$UniMod[w] <- paste0("UniMod:", tmpMds$UniMod[w])
    }
  }
  if (modsType == "StrangerThings") {
    tmpMds <- proteoCraft::listMelt(allPTMs$AA, 1:nrow(allPTMs), c("AA", "Row"))
    kol <- c("Full name", "Mark")
    #kol %in% colnames(allPTMs)
    tmpMds[, c("Name", "Mark")] <- allPTMs[tmpMds$Row, kol]
  }
  if (modsType == "No idea") {
    tmpMds <- proteoCraft::listMelt(allPTMs$`Old mark`, allPTMs$`Full name`, c("Old", "Name"))
    tmpMds$Mark <- allPTMs$Mark[match(tmpMds$Name, allPTMs$`Full name`)]
  }
  parallel::clusterExport(cl, list("tmpMds", "modsType"), envir = environment())
  uMdSq <- unique(EV$`Modified sequence`[wMod])
  temp1 <- strsplit(uMdSq, "\\(|\\)")
  f0 <- function(x) {
    #x <- temp1[[1]]
    x1 <- x2 <- unlist(x)
    l <- length(x1)
    rg <- c(1:((l-1)/2))*2
    nc <- nchar(x1[rg-1])
    aa <- substr(x1[rg-1], nc, nc)
    if (nc[1] == 1) { aa[1] <- substr(x1[3], 1, 1) }
    if (modsType == "UniMod") {
      m <- match(x1[rg], tmpMds$UniMod)
    }
    if (modsType == "StrangerThings") {
      m <- match(x1[rg], tmpMds$Name)
    }
    if (modsType == "No idea") {
      m <- match(x1[rg], tmpMds$Old)
    }
    x1[rg] <- paste0("(", tmpMds$Mark[m], ")")
    x2[rg] <- paste0("(", tmpMds$Name[m], ")")
    x1 <- paste(x1, collapse = "")
    x2 <- paste(x2, collapse = "")
    return(c(x1, x2))
  }
  environment(f0) <- .GlobalEnv
  temp2 <- as.data.frame(t(parSapply(cl, temp1, f0)))
  temp2 <- temp2[match(EV$`Modified sequence`[wMod], uMdSq),]
  EV$`Modified sequence_verbose` <- paste0("_", EV$Sequence, "_")
  EV$`Modified sequence`[wMod] <- temp2[, 1]
  EV$`Modified sequence_verbose`[wMod] <- temp2[, 2]
  #
  EV$Modifications <- "Unmodified"
  uMdSq <- unique(EV$`Modified sequence`[wMod])
  tempMod1 <- gsub("[_|\\)][A-Z]*[_|\\(]", ",", uMdSq)
  tempMod2 <- strsplit(tempMod1, ",")
  parallel::clusterExport(cl, "Fixed.mods", envir = environment())
  f0 <- function(x) { x[which(!x %in% c("", Fixed.mods))] }
  environment(f0) <- .GlobalEnv
  tempMod2 <- parallel::parLapply(cl, tempMod2, f0)
  if (modsType == "No idea") {
    f0 <- function(x) {
      tmpMds$Name[match(x, tmpMds$Mark)]
    } 
    environment(f0) <- .GlobalEnv
    tempMod2 <- parallel::parLapply(cl, tempMod2, f0)
  }
  f0 <- function(x) {
    if (!length(x)) { x <- "Unmodified" } else {
      x <- aggregate(x, list(x), length)
      x <- x[order(x$Group.1, decreasing = FALSE), 2:1]
      x$x <- gsub("1 ", "", paste0(as.character(x$x), " "))
      x <- paste0(apply(x, 1, paste, collapse = ""), collapse = ",") 
    }
    return(x)
  }
  environment(f0) <- .GlobalEnv
  tempMod3 <- parSapply(cl, tempMod2, f0)
  EV$Modifications[wMod] <- tempMod3[match(EV$`Modified sequence`[wMod], uMdSq)]
  # M/Z
  cat("   Calculating theoretical m/z values...\n")
  parallel::clusterExport(cl, "allPTMs", envir = environment())
  f0 <- function(x) {
    x <- x[which(x != "")]
    x <- x[which(x != "UniMod:4")] # Peptides::mz can deal with alkylation
    return(sum(allPTMs$`Mass shift`[match(x, paste0("UniMod:", allPTMs$UniMod))]))
  }
  environment(f0) <- .GlobalEnv
  tempMod4 <- parallel::parSapply(cl, tempMod2, f0)
  tmp1 <- do.call(paste, c(EV[, c("Sequence", "Charge")], sep = "___"))
  tmp2 <- unique(tmp1)
  tmp3 <- EV[match(tmp2, tmp1), c("Sequence", "Charge")]                              
  f0 <- function(x) { Peptides::mz(x[[1]], as.integer(x[[2]])) }
  environment(f0) <- .GlobalEnv
  tmp3 <- parallel::parApply(cl, tmp3, 1, f0) # Peptides::mz doesn't like vectorization clearly
  EV$"Theoretical m/z" <- tmp3[match(tmp1, tmp2)] 
  EV$"Theoretical m/z"[wMod] <- EV$"Theoretical m/z"[wMod] + tempMod4[match(EV$`Modified sequence`[wMod], uMdSq)]/EV$Charge[wMod]
  if (!"m/z" %in% colnames(EV)) { EV$"m/z" <- EV$"Theoretical m/z" }
  w <- which(is.na(EV$"m/z"))
  if (length(w)) { EV$"m/z"[w] <- EV$"Theoretical m/z"[w] }
  #
  EV$Mass <- (EV$"m/z"-1.007276466879)*EV$Charge
  # Retention time (in minutes, as in MaxQuant)
  EV$"Retention time" <- as.numeric(DIANN$RT)
  EV$"Retention time (predicted)" <- as.numeric(DIANN$Predicted.RT)
  EV$iRT <- as.numeric(DIANN$iRT)
  EV$"iRT (predicted)" <- as.numeric(DIANN$Predicted.iRT)
  EV$"Retention time (start)" <- as.numeric(DIANN$RT.Start)
  EV$"Retention time (end)" <- as.numeric(DIANN$RT.Stop)
  EV$"Retention length" <- EV$"Retention time (end)" - EV$"Retention time (start)"
  # 
  if ("Decoy" %in% colnames(DIANN)) {
    EV$Reverse <- c("", "+")[as.logical(DIANN$Decoy)+1]
  } else {
    EV$Reverse <- ""
  }
  if ("Score" %in% colnames(DIANN)) {
    EV$Score <- as.numeric(DIANN$CScore)
  }
  if ("Decoy.CScore" %in% colnames(DIANN)) {
    EV$"Delta score" <- as.numeric(DIANN$CScore) - as.numeric(DIANN$Decoy.CScore)
  }
  #
  if (PG.Genes) {
    EV$"Leading proteins" <- DIANN$Protein.Group # (Will be removed later anyway in the classic workflow, since we perform protein inference again)
    kol <- c("PG.Quantity", "PG.Normalised", "PG.MaxLFQ", "Genes.Quantity", "Genes.Normalised", "Genes.MaxLFQ", "Genes.MaxLFQ.Unique")
    kol <- kol[which(kol %in% colnames(DIANN))]
    EV[, gsub(" Normalised$", " Normalised quantity", gsub("\\.", " ", kol))] <- DIANN[, kol]
    kol <- c("Protein.Q.Value", "PG.Q.Value", "Global.PG.Q.Value", "GG.Q.Value")
    EV[, gsub(" Q Value$", " q-value", gsub("\\.", " ", kol))] <- DIANN[, kol]
  }
  kol <- c("Precursor.Quantity", "Precursor.Normalised", # Should be QuantUMS values, see https://github.com/vdemichev/DiaNN/discussions/764
           # As I understand:
           #  - Precursor.Quantity (and .Normalised?) uses MS2, and also MS1 if the precursor is found in all samples
           #  - Quantity.Quality is a new metric used to determine how good the quantification is
           "Quantity.Quality", "Empirical.Quality", "Normalisation.Noise", "Ms1.Profile.Corr", "FWHM", # Quality metrics
           "Precursor.Translated", "Translated.Quality", "Translated.Q.Value", "Channel.Q.Value", # If I understand correctly these ones are for SILAC-DIA... see https://github.com/vdemichev/DiaNN/discussions/542
           "Ms1.Normalised", "Ms1.Apex.Area", "Ms1.Apex.Mz.Delta", "Normalisation.Factor", "Ms1.Total.Signal.Before", "Ms1.Total.Signal.After",
           "Q.Value", "Global.Q.Value", "Lib.Q.Value",
           # PTM-related values:
           "Peptidoform.Q.Value", "Global.Peptidoform.Q.Value", "Lib.Peptidoform.Q.Value", "PTM.Site.Confidence", "Site.Occupancy.Probabilities",
           "Protein.Sites", "Lib.PTM.Site.Confidence"
  )
  kol <- kol[which(kol %in% colnames(DIANN))]
  EV[, gsub("Ms1", "MS1", gsub("\\.", " ", kol))] <- DIANN[, kol]
  #
  # Quantitative values
  EV$"MS1 Area" <- as.numeric(DIANN$Ms1.Area)
  QuantUMS <- ("Precursor.Quantity" %in% colnames(DIANN))
  if (QuantUMS) {
    cat("   QuantUMS output columns detected, using those as Intensity by default. For your convenience the Ms1.Area column is still present in the output.\n")
    # If the QuantUMS output is present, we use it by default: it should be higher quality than "Ms1.Area" (MS1-based quant only).
    # We do not use the normalized version by default, in case normalization would actually be bad for this particular experiment
    # ("normalization should never be on by default!")
    EV$Precursor.Quantity <- as.numeric(DIANN$Precursor.Quantity)
    EV$Intensity <- EV$Precursor.Quantity
  } else {
    EV$Intensity <- EV$"MS1 Area"
  }
  IntColz <- c("Fragment.Quant.Corrected", "Fragment.Quant.Raw")
  IntColz <- IntColz[which(IntColz %in% colnames(DIANN))]
  if (length(IntColz)) {
    EV$"MS2 intensities" <- DIANN[[IntColz[1]]]
  }
  # Fragment annotations are not provided by default,
  # but each PSM for a given peptidoform*Charge combination always has the same number of values.
  # Assuming the order is then also always the same (which appears to be the case), we have all we need for MS2-based quant.
  #
  # Careful: this does not work where 2 different DiaNN files are combined, because the number of fragments will be different in different runs (best on which are found).
  #
  # Below commented code if you want to check
  #tst <- aggregate(EV$`MS2 intensities`, list(EV$`Modified sequence`, EV$Charge), function(x) {
  #  nchar(x) - nchar(gsub(";", "", x))
  #})
  #tst2 <- sapply(tst$x, function(x) { length(unique(x)) })
  #stopifnot(max(tst2) == 1)
  if ("MS2.Scan" %in% colnames(DIANN)) {
    EV$"MS/MS scan number" <- as.integer(DIANN$MS2.Scan)
  }
  # Ion mobility
  EV$IM <- as.numeric(DIANN$IM)
  EV$"IM (predicted)" <- as.numeric(DIANN$Predicted.IM)
  EV$iIM <- as.numeric(DIANN$iIM)
  EV$"iIM (predicted)" <- as.numeric(DIANN$Predicted.iIM)
  #
  EV$Type <- "LIB-DIA" # This may need to change if important; usually diaNN recommends library-free + MBR, not library-based searches;
  # Not sure if MaxQuant assigns a type for each different way to run DIA
  tmp <- gsub("[K,R]$", "", EV$Sequence)
  EV$"Missed cleavages" <- nchar(tmp) - nchar(gsub("[K,R]", "", tmp))
  EV$"Potential contaminant" <- ""
  if ("Search_ID" %in% colnames(DIANN)) {
    EV$Search_ID <- DIANN$PSMs_file # In cases where the input is a hybrid report created by the merging script in .../Utils
  } else {
    EV$Search_ID <- sort(DIANN_fl, decreasing = TRUE)[1] # Provide default: there can be one or 2 files, the .tsv or .parquet, we want just one (ideally the .tsv)
  }
  # Optional: filter by delta score
  if (("Delta score" %in% colnames(EV))&&(!is.na(Min.Delta.Score))) {
    EV <- EV[which(EV$"Delta score" > Min.Delta.Score),]
  }
  #
  kol <- c("Precursor.Id", "Precursor.Lib.Index", "Decoy")
  #
  EV$id <- NULL
  EV <- cbind(data.frame(id = 1:nrow(EV)),
                         EV)
  Res <- list(Evidence = EV,
              PTMs = allPTMs,
              QuantUMS = QuantUMS)
  if (cleanUp) { parallel::stopCluster(cl) }
  return(Res)
}
