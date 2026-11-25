#' DIANN_to_MQ
#'
#' @description 
#' Converts a diaNN PSMs table to a MaxQuant evidence.txt-like table.
#' As with the other functions of this family, the idea is not necessarily to get perfect conversion, but close enough that the data can be fed into this package's analysis scripts.
#' 
#' @returns
#' A list of 3 items:\cr
#' - Evidence: a data.frame similar to the evidence.txt table MaxQuant creates\cr
#' - PTMs: a modifications table\cr
#' - QuantUMS: was QuantUMS quantification on?
#' 
#' Surprisingly, the diaNN PSMs table does not include MZ values! Is there a way to do this by adding a command line argument?
#' For now they are either detected from the library file (same folder) or computed.
#' Cysteine alkylation is assumed to be Carbamidomethyl
#' 
#' @param DIANN_fl Path to a diaNN PSMs table file (.tsv, .parquet, or both).
#' @param Fixed.mods A vector of those Fixed modifications which will not be reported in the results (reported by diaNN, but MQ ignores them). Default = "Carbamidomethyl"
#' @param PG.Genes If TRUE, will include per Protein Group/Gene quantitative and Q-values columns. Default = FALSE, because we will mostly be using this function in contexts where we want to perform protein groups inference re-calculate quantitative values ourselves.
#' @param Detect.Lib diaNN does not output M/Z values, but those are in the library generated at the end of some runs. If this is set to TRUE, it will try to detect the library and if available get M/Z values from there. Otherwise they will be computed.
#' @param Min.Delta.Score From which minimum delta score should we report a PSM? Default = -Inf (we assume that those few peptides with low delta score are because randomly the decoy database will include real-like peptides)
#' @param N.clust A limit on the number of vCPUs to use. If left as NULL (default), uses the number of available clusters - 1, to a minimum of 1.
#' @param N.reserved Default = 1. Number of reserved vCPUs the function is not to use. Note that for obvious reasons it will always use at least one.
#' @param cl Already have a cluster handy? Why spend time making a new one, which on top of that may invalidate the old one. Just pass it along!
#' @param shinyOpt One of "popup" (default), "pane", "dialog", or "browser".
#' @param digPattern Default = "KR" for trypsin. Only used for the "Missed cleavages" column. Vector of amino acids after which the enzyme cleaves. 
#' @param log_Fl If missing, the function will try to detect the log file based on the name of the report - but this may fail!
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
                        N.reserved = 1L,
                        cl,
                        shinyOpt = "popup",
                        digPattern = "KR",
                        log_Fl) {
  TESTING <- FALSE
  #DefArg(DIANN_to_MQ);TESTING <- TRUE; cl <- parClust
  #DIANN_fl <- PSMsFls
  #DIANN_fl <- rstudioapi::selectFile(path = paste0(getwd(), "/*.tsv"))
  #DIANN_fl <- diannRep
  #DIANN_fl <- psmFls_i; log_Fl = diaNN_logFl_i
  shinyOpt <- tolower(shinyOpt)
  if (!shinyOpt %in% c("popup", "dialog", "pane", "browser")) { shinyOpt <- "popup" }
  #
  # NB:
  # This function deliberately uses a different object name for the PTMs table to other related functions,
  # to simplify environment-related debugging.
  # Do not change that!
  # Here, use allPTMs, elsewhere Modifs (or anything else).
  misFun <- if (TESTING) {
    # Note:
    # This is not a perfect alternative to missing but will work in most cases, unless x matches a function imported by a package 
    \(x) { return(!exists(deparse(substitute(x)))) }
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
  wTsv <- grep("\\.tsv$", DIANN_fl)
  wPrqt <- grep("\\.parquet$", DIANN_fl)
  lTsv <- length(wTsv)
  lPrqt <- length(wPrqt)
  stopifnot(lTsv <= 1L,
            lPrqt <= 1L,
            sum(lTsv, lPrqt) > 0L)
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
  # Filter out peptides with ambiguous/non-classical amino acids
  uSeq <- unique(DIANN$Stripped.Sequence)
  test <- gsub(paste(AA, collapse = "|"), "", uSeq)
  test <- test[match(DIANN$Stripped.Sequence, uSeq)]
  w <- which(test == "")
  lX <- nrow(DIANN)-length(w) 
  if (lX) {
    u <- unique(unlist(strsplit(DIANN$Stripped.Sequence, "")))
    u <- paste(u[!u %in% AA], collapse = "-")
    warning(paste0("Removing ", lX, " peptide", c("", "s")[(lX > 1L)+1L], " with unknown amino acids ", u, "!"))
    DIANN <- DIANN[w,]
  }
  #
  if (misFun(log_Fl)) {
    log_Fl <- unique(sub("\\.(tsv|parquet)$", ".log.txt", DIANN_fl))
  }
  if (file.exists(log_Fl)) { 
    log <- readr::read_lines(log_Fl)
    diannCall <- grep(" --[a-zA-Z]", log, value = TRUE)[1L]
    diannCall <- unlist(strsplit(diannCall, " +--"))
    diannCall <- diannCall[2L:length(diannCall)]
  }
  #
  libTst <- try({
    libKol <- c()
    if (Detect.Lib) {
      if (file.exists(log_Fl)) {
        DIANNLib_fls <- gsub(".+\\] Saving the library to |\\.speclib$", "", grep("\\] Saving the library to ", log, value = TRUE))
        DIANNLib_fls <- c(grep("\\.predicted$", DIANNLib_fls, value = TRUE, invert = TRUE),
                          grep("\\.predicted$", DIANNLib_fls, value = TRUE))
        x <- sub("^lib ", "", grep("^lib ", diannCall, value = TRUE))
        x <- x[nchar(x) > 0L]
        if (length(x)) { diannCall <- c(diannCall, x) }
        DIANNLib_fls <- gsub("\\\\", "/", DIANNLib_fls)
        DIANNLib_fls <- base::setdiff(DIANNLib_fls, "")
        DIANNLib_fls <- base::union(sub("\\.skyline$", "", DIANNLib_fls), DIANNLib_fls)
        if (length(DIANNLib_fls)) {
          if (sum(dirname(DIANNLib_fls) != dirname(DIANN_fl))) {
            DIANNLib_fls <- c(DIANNLib_fls, paste0(dirname(DIANN_fl), "/", basename(DIANNLib_fls)))
          }
          DIANNLib_fls <- DIANNLib_fls[file.exists(DIANNLib_fls)]
          nLibs <- length(DIANNLib_fls)
          if (nLibs) {
            DIANNLibs <- lapply(DIANNLib_fls, \(diannLib_fl) { #diannLib_fl <- DIANNLib_fls[1L]
              if (grepl("\\.parquet$", diannLib_fl)) {
                diannLib <- arrow::read_parquet(diannLib_fl)
              } else {
                diannLib <- data.table::fread(diannLib_fl, integer64 = "numeric", check.names = FALSE, data.table = FALSE)
              }
              return(diannLib)
            })
            DIANNLibs <- plyr::rbind.fill(DIANNLibs)
            if (identical(c("Modified.Sequence", "ModifiedPeptideSequence") %in% colnames(DIANNLib),
                          c("FALSE, TRUE"))) {
              colnames(DIANNLib)[match("ModifiedPeptideSequence", colnames(DIANNLib))] <- "Modified.Sequence"
            }
            if (identical(c("Precursor.Charge", "PrecursorCharge") %in% colnames(DIANNLib),
                          c("FALSE, TRUE"))) {
              colnames(DIANNLib)[match("PrecursorCharge", colnames(DIANNLib))] <- "Precursor.Charge"
            }
            if (identical(c("Precursor.Mz", "PrecursorMz") %in% colnames(DIANNLib),
                          c("FALSE, TRUE"))) {
              colnames(DIANNLib)[match("PrecursorMz", colnames(DIANNLib))] <- "Precursor.Mz"
            }
            if ("transition_group_id" %in% colnames(DIANNLib)) {
              DIANNLib$transition_group_id <- gsub("\\[", "(", gsub("\\]", ")", DIANNLib$transition_group_id))
            } else {
              DIANNLib$Modified.Sequence <- gsub("\\[", "(", gsub("\\]", ")", DIANNLib$Modified.Sequence))
              DIANNLib$transition_group_id <- do.call(paste, c(DIANNLib[, c("Modified.Sequence", "Precursor.Charge")], sep = "+"))
            }
            #a <- unique(unlist(strsplit(gsub("^[A-Z]*\\(|\\)[A-Z]*$|^[A-Z]+$", "", gsub("\\)[A-Z]*\\(", "_", DIANNLib$Modified.Sequence)), "_")))
            tmp <- do.call(paste, c(DIANN[, c("Modified.Sequence", "Precursor.Charge")], sep = "+"))
            m <- match(tmp, DIANNLib$transition_group_id)
            DIANN$"m/z" <- DIANNLib$Precursor.Mz[m]
            # Also library index + rt
            DIANN$"Library index" <- DIANNLib$transition_group_id[m]
            DIANN$"Library rt" <- if ("RetentionTime" %in% colnames(DIANNLib)) { DIANNLib$RetentionTime[m] } else {
              if ("AverageExperimentalRetentionTime" %in% colnames(DIANNLib)) { DIANNLib$AverageExperimentalRetentionTime[m] } else {
                if ("NormalizedRetentionTime" %in% colnames(DIANNLib)) { DIANNLib$NormalizedRetentionTime[m] } else { NA }
              }
            }
            libKol <- c("m/z", "Library index", "Library rt")
          } else { cat("   (NB: no diaNN library tsv found in folder)\n") }
        } else { cat("   (NB: no diaNN library tsv found in folder)\n") }
      } else { cat("   (NB: no diaNN log file found in folder)\n") }
    }
  }, silent = TRUE)
  #
  # Create MQ-like table
  if (!"File.Name" %in% colnames(DIANN)) {
    stopifnot("Run.Index" %in% colnames(DIANN),
              file.exists(log_Fl))
    tst <- aggregate(DIANN$Run, list(DIANN$Run.Index), unique)
    fls <- gsub("(\\\\)+", "/", sub("^f ", "", grep("^f ", diannCall, value = TRUE)))
    flNms <- gsub(".*/|\\.[^.]+$", "", fls)
    tst$x2 <- flNms[tst$Group.1 + 1L]
    stopifnot(sum(tst$x != tst$x2) == 0L)
    tst$File <- fls[tst$Group.1 + 1L]
    DIANN$File.Name <- tst$File[match(DIANN$Run.Index, tst$Group.1)]
  }
  kol <- c("Run", "File.Name", "Protein.Ids", "Protein.Names", "Proteotypic", "Genes",
           "Stripped.Sequence", "Modified.Sequence", "Precursor.Charge", "PEP")
  #print(base::setdiff(kol, colnames(DIANN)))
  stopifnot(sum(!kol %in% colnames(DIANN)) == 0L)
  EV <- data.frame("Raw file" = DIANN$Run,
                   "Raw file path" = DIANN$File.Name,
                   "Proteins" = DIANN$Protein.Ids, # Good idea to check how reliable this is using ProtMatch
                   "Full Protein ID" = DIANN$Protein.Names, # Same as above
                   "Proteotypic" = as.logical(DIANN$Proteotypic),
                   "Gene names" = DIANN$Genes,
                   "Sequence" = DIANN$Stripped.Sequence,
                   "Mod. seq. (DiaNN format)" = DIANN$Modified.Sequence, # Useful for matching to other DiaNN outputs!
                   "Modified sequence" = paste0("_", DIANN$Modified.Sequence, "_"), # Will need some reprocessing further down
                   "Modifications" = "Unmodified",
                   "Length" = nchar(DIANN$Stripped.Sequence),
                   "Charge" = DIANN$Precursor.Charge,
                   "PEP" = DIANN$PEP,
                   check.names = FALSE)
  #
  if (!inherits(libTst, "try-error")) {
    if (length(libKol)) { EV[, libKol] <- DIANN[, libKol] }
  }
  if ("First.Protein.Description" %in% colnames(DIANN)) {
    EV$"First protein name" <- DIANN$First.Protein.Description
  } #else {
  #cat("   (Column \"First.Protein.Description\" not found, no column of human-readable names generated...)\n")
  #}
  # Process modifications
  wMod <- grep("\\(", EV$`Modified sequence`)
  wNotMod <- base::setdiff(1L:nrow(EV), wMod)
  mWMod <- match(unique(EV$`Modified sequence`[wMod]), EV$`Modified sequence`)
  tmp <- data.frame(Seq = EV$Sequence[mWMod],
                    ModSeq = EV$`Modified sequence`[mWMod],
                    Z = EV$Charge[mWMod])
  f0 <- \(x) { Peptides::mz(x[[1L]], as.integer(x[[2L]]), cysteins = 0) }
  environment(f0) <- .GlobalEnv
  tmp$Theoretical <- parallel::parApply(cl, tmp[, c("Seq", "Z")], 1L, f0)
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
    screenRes <- rpanel::rp.screenresolution()
    if (shinyOpt == "popup") {
      run_App1 <- "print(shiny::shinyApp(.DIANN_to_MQ_ui1(tstMap), .DIANN_to_MQ_server1, options = list(height = screenRes$height, width = screenRes$width)))"
      run_App2 <- "print(shiny::shinyApp(.DIANN_to_MQ_ui2, .DIANN_to_MQ_server2, options = list(height = screenRes$height, width = screenRes$width)))"
    }
    if (shinyOpt %in% c("dialog", "pane")) {
      if (shinyOpt == "pane") {
        myViewer <- shiny::paneViewer("Viewer", minHeight = "maximize")
      }
      if (shinyOpt == "dialog") {
        myViewer <- shiny::dialogViewer("Viewer", width = screenRes$width, height = screenRes$height)
      }
      run_App1 <- "shiny::runGadget(shiny::shinyApp(.DIANN_to_MQ_ui1(tstMap), .DIANN_to_MQ_server1, options = list(height = screenRes$height, width = screenRes$width)), viewer = myViewer, stopOnCancel = FALSE)"
      run_App2 <- "shiny::runGadget(shiny::shinyApp(.DIANN_to_MQ_ui2, .DIANN_to_MQ_server2, options = list(height = screenRes$height, width = screenRes$width)), viewer = myViewer, stopOnCancel = FALSE)"
    }
    if (shinyOpt == "browser") {
      run_App1 <- "print(shiny::shinyApp(.DIANN_to_MQ_ui1(tstMap), .DIANN_to_MQ_server1, options = list(height = \"100%\", width = \"100%\", launch.browser = TRUE)))"
      run_App2 <- "print(shiny::shinyApp(.DIANN_to_MQ_ui2, .DIANN_to_MQ_server2, options = list(height = \"100%\", width = \"100%\", launch.browser = TRUE)))"
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
      #g <- shiny:::.globals
      #g$appState <- NULL
      #eval(parse(text = run_App1))
      #shinyCleanup()
      #tstMap <- c(FALSE, TRUE)[match(tstMap, c("no", "yes"))]
      msg <- paste0("Non-classical PTMs marks detected (e.g. \"", sort(comps)[1L], "\")...")
      opt <- c("Keep them as they are                                                                                                                 ",
               "Manually map each to a more convenient mark                                                                                           ")
      tstMap <- c(FALSE, TRUE)[match(svDialogs::dlg_list(opt, opt[1L], title = msg)$res, opt)]
    }
    if (tstMap) {
      # Here we will want to make sense of those mods
      # First, get the precise mass shift for each, so we can narrow down options
      f0 <- \(x) {
        #x <- tmp$ModSeq[1L]
        x <- annot_to_tabl(x)[[1L]]
        w <- which(x$Annotations != "")
        x$Sequence[2L] <- paste0("_", x$Sequence[2L])
        do.call(paste, c(x[w, , drop = FALSE], sep = "mod_"))
      }
      environment(f0) <- .GlobalEnv
      tmp$Mods <- parallel::parLapply(cl, tmp$ModSeq, f0)
      Mods <- sort(unique(unlist(tmp$Mods)))
      parallel::clusterExport(cl, "Mods", envir = environment())
      f0 <- \(x) {
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
      precVal <- 0.02
      Mods$Options <- apply(Mods[, c("DM", "AA")], 1L, \(x) {
        dm <- as.numeric(x[[1L]])
        w <- if (substr(x[[2L]], 1L, 1L) == "_") {
          which((abs(UniMod$MonoMass - dm) <= precVal) & (UniMod$Site %in% c(substr(x[[2L]], 2L, 2L), "N-term")))
        } else {
          which((abs(UniMod$MonoMass - dm) <= precVal) & (UniMod$Site == x[[2L]]))
        }
        if (length(w)) {
          Ddm <- dm - UniMod$MonoMass[w]
          w <- w[order(Ddm, decreasing = FALSE)]
        }
        return(w)
      })
      Mods$tag <- do.call(paste, c(data.frame(AA = sub("^_", "", Mods$AA),
                                              Name = paste0("+", as.character(round(Mods$DM), 4L)),
                                              NTerm = c("", "(N-term)")[grepl("^_", Mods$AA)+1L]),
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
      Mods$Options2 <- lapply(1L:nrow(Mods), \(i) {
        m <- match(Mods$Options[[i]], UniMod2$Option)
        Ddm <- Mods$DM[i] - UniMod$MonoMass[m]
        opt <- UniMod2$Id[m]
        return(opt[order(Ddm, decreasing = FALSE)])
      })
      Mods$Match <- sapply(1L:nrow(Mods), \(i) {
        opt <- Mods$Options2[[i]]
        if (length(opt) >= 1L) {
          dflt <- opt[1L]
          m <- match(commonPTMs, opt)
          m <- m[!is.na(m)]
          if (length(m)) { dflt <- opt[m[1L]] }
        } else { dflt <- Mods$tag[i] }
        return(dflt)
      })
      homePath <- paste0(normalizePath(Sys.getenv("HOME"), winslash = "/"), "/R/proteoCraft")
      fl <- paste0(homePath, "/tmpMods.tmp")
      file.create(fl)
      con <- file(fl, "r+")
      serialize(Mods, con)
      suppressWarnings(suppressMessages(try(base::close(con), silent = TRUE))) # Only necessary in function mode
      #run_App2 <- "print(shiny::shinyApp(.DIANN_to_MQ_ui2, .DIANN_to_MQ_server2, options = list(height = screenRes$height, width = screenRes$width)))"
      g <- shiny:::.globals
      g$appState <- NULL
      eval(parse(text = run_App2))
      shinyCleanup()
      con <- file(fl, "r")
      Mods <- unserialize(con)
      suppressWarnings(suppressMessages(try(base::close(con), silent = TRUE))) # Only necessary in function mode
      unlink(fl)
      #
      allPTMs <- aggregate(Mods$DM, list(Mods$Match), \(x) {
        rg <- c(min(x), max(x))
        stopifnot(rg[2L]-rg[1L] < 0.010)
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
  allPTMs$Type <- c("Variable", "Fixed")[(allPTMs$"Full name" %in% Fixed.mods) + 1L]
  allPTMs$Mark <- tolower(substr(allPTMs$"Full name", 1L, 2L))
  ## Identify affected Amino Acids
  uMdSq <- unique(EV$`Modified sequence`[wMod])
  if (modsType == "UniMod") {
    allPTMs$AA <- list(NA)
    for (i in 1L:nrow(allPTMs)) { #i <- 1L
      temp <- gsub(topattern(paste0("(", allPTMs$UniMod[i], ")"), start = FALSE),
                   ">>>.___", uMdSq)
      temp2 <- grep(">>>$", unlist(strsplit(temp, "\\.___")), value = TRUE)
      nc <- nchar(temp2)
      allPTMs$AA[[i]] <- unique(substr(temp2, nc-3L, nc-3L))
    }
  }
  if (modsType == "StrangerThings") {
    allPTMs$AA <- list(NA)
    for (i in 1L:nrow(allPTMs)) { #i <- 1L
      temp <- gsub(topattern(paste0("(", allPTMs$`Full name`[i], ")"), start = FALSE),
                   ">>>.___", EV$`Modified sequence`)
      temp2 <- grep(">>>$", unlist(strsplit(temp, "\\.___")), value = TRUE)
      nc <- nchar(temp2)
      allPTMs$AA[[i]] <- unique(substr(temp2, nc-3L, nc-3L))
    }
  }
  if (modsType == "No idea") {
    tmp <- aggregate(sub("^_", "", Mods$AA), list(Mods$Match), \(x) { sort(unique(x)) })
    allPTMs$Site <- allPTMs$AA <- tmp$x[match(allPTMs$`Full name`, tmp$Group.1)]
    w <- which(allPTMs$Position == "Any N-term")
    if (length(w)) { allPTMs$Site[w] <- lapply(allPTMs$Site[w], \(x) { paste0("n", x) }) }
    w <- which(allPTMs$Position == "Any C-term")
    if (length(w)) { allPTMs$Site[w] <- lapply(allPTMs$Site[w], \(x) { paste0("c", x) }) }
    w <- which(allPTMs$Position == "Protein N-term")
    if (length(w)) { allPTMs$Site[w] <- "[^" }
    w <- which(allPTMs$Position == "Protein C-term")
    if (length(w)) { allPTMs$Site[w] <- "^]" }
    allPTMs$Site_long <- lapply(strsplit(gsub("\\^", "", lapply(allPTMs$Site, paste, collapse = "")), ""), \(x) {
      x[x %in% c("[", "n")] <- "N-term"
      x[x %in% c("]", "c")] <- "C-term"
      return(x)
    })
  }
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
  if (modsType == "No idea") {
    tmp <- aggregate(gsub(".*mod_", "", Mods$Name), list(Mods$Match), \(x) { sort(unique(x)) })
    allPTMs$"Old mark" <- tmp$x[match(allPTMs$`Full name`, tmp$Group.1)] # Could be a list!!!
  }
  if (!"UniMod" %in% colnames(allPTMs)) {
    allPTMs$UniMod <- UniMod$UnimodId[match(allPTMs$`Full name`, UniMod$Name)] # Is not meant to always work
    # Sometimes there will be discrepancies, where UniMod has a compatible isobaric PTM with a different name!
  }
  if (!"Mass shift" %in% colnames(allPTMs)) {
    allPTMs$"Mass shift" <- UniMod$MonoMass[match(sub("^UniMod:", "", allPTMs$UniMod), UniMod$UnimodId)]
    if ("Mass delta" %in% colnames(allPTMs)) {
      tst <- allPTMs$"Mass shift" - allPTMs$"Mass delta"
      stopifnot(max(tst, na.rm = TRUE) < 0.1)
    }
  }
  #
  tmpLst <- listMelt(strsplit(uMdSq, "\\(|\\)"), ColNames = c("modSeq", "wh"))
  tmpLst$Verbose <- tmpLst$modSeq
  if (modsType == "UniMod") {
    mtchKol <- "UniMod"
  }
  if (modsType == "StrangerThings") {
    mtchKol <- "Name"
  }
  if (modsType == "No idea") {
    mtchKol <- "Old"
  }
  wV <- which(tmpLst$modSeq %in% allPTMs[allPTMs$Type == "Variable", mtchKol])
  wF <- which(tmpLst$modSeq %in% allPTMs[allPTMs$Type == "Fixed", mtchKol])
  if (length(wV)) {
    mV <- match(tmpLst$modSeq[wV], allPTMs[[mtchKol]])
    tmpLst$modSeq[wV] <- paste0("(", allPTMs$Mark[mV], ")")
    tmpLst$Verbose[wV] <- paste0("(", allPTMs$"Full name"[mV], ")")
  }
  if (length(wF)) {
    mF <- match(tmpLst$modSeq[wF], allPTMs[[mtchKol]])
    tmpLst$modSeq[wF] <- ""
    tmpLst$Verbose[wF] <- paste0("(", allPTMs$"Full name"[mF], ")")
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
  EV$`Modified sequence_verbose` <- EV$"Modified sequence"
  tmp1$Verbose <- tmp2$Verbose[match(tmp1$wh, tmp2$wh)]
  tmp1$Mods <- tmp3$Mods[match(tmp1$wh, tmp3$wh)]
  tmp1$Orig <- uMdSq[tmp1$wh]
  uM <- match(EV$"Modified sequence"[wMod], tmp1$Orig)
  stopifnot(sum(is.na(uM)) == 0L)
  EV[wMod, c("Modified sequence", "Modified sequence_verbose", "Modifications")] <- tmp1[uM, c("modSeq", "Verbose", "Mods")]
  #
  # Mass, theoretical mass and M/Z
  #cat("   Calculating theoretical m/z values...\n")
  EV$"Theoretical mass" <- NA_real_
  tmp1 <- tmpLst[!tmpLst$Verbose %in% c("", "_"), c("Verbose", "wh")]
  tmp1$Verbose <- gsub("^[(_]|[_)]$", "", tmp1$Verbose) # Verbose because it includes fixed mods!
  tmp1 <- tmp1[tmp1$Verbose != "", c("Verbose", "wh")]
  wPTM <- which(tmp1$Verbose %in% allPTMs$"Full name")
  wSeq <- base::setdiff(1L:nrow(tmp1), wPTM)
  tmp2 <- data.table::data.table(tmp1[wSeq,]) # We need to re-aggregate those fragments because the mass(PEPTIDER) = mass(PEP) + mass(TIDER)... - mass(H20) -> don't forget the water!!!
  tmp1 <- tmp1[wPTM,]
  tmp2 <- as.data.frame(tmp2[, .(Verbose = paste0(Verbose, collapse = "")), by = .(wh = wh)])
  tmp1 <- rbind(tmp1, tmp2)
  tmp1$Mass <- 0
  # update
  wPTM <- which(tmp1$Verbose %in% allPTMs$"Full name")
  wSeq <- base::setdiff(1L:nrow(tmp1), wPTM)
  #
  tmp1$Mass[wPTM] <- allPTMs$`Mass shift`[match(tmp1$Verbose[wPTM], allPTMs$"Full name")]
  uS <- base::union(tmp1$Verbose[wSeq], EV$Sequence[wNotMod])
  mss <- parallel::parSapply(cl, uS, Peptides::mw, monoisotopic = TRUE)
  tmp1$Mass[wSeq] <- mss[match(tmp1$Verbose[wSeq], uS)]
  EV$"Theoretical mass"[wNotMod] <- mss[match(EV$Sequence[wNotMod], uS)]
  data.table::setDT(tmp1)
  tmp1 <- tmp1[, .(Mass = sum(Mass)),
               by = .(wh = wh)]
  tmp1$Orig <- sub("^_", "", sub("_$", "", uMdSq[tmp1$wh]))
  uM <- match(EV$"Mod. seq. (DiaNN format)"[wMod], tmp1$Orig)
  stopifnot(sum(is.na(uM)) == 0L)
  EV$"Theoretical mass"[wMod] <- tmp1$Mass[uM]
  EV$"Theoretical m/z" <- EV$"Theoretical mass"/EV$Charge + 1.007276466879
  if ("m/z" %in% colnames(EV)) {
    w <- which(is.na(EV$"m/z"))
    if (length(w)) { EV$"m/z"[w] <- EV$"Theoretical m/z"[w] }
    EV$Mass <- (EV$"m/z"-1.007276466879)*EV$Charge
  } else {
    colnames(EV)[match("Theoretical m/z", colnames(EV))] <- "m/z"
    colnames(EV)[match("Theoretical mass", colnames(EV))] <- "Mass"
  }
  #View(EV[, c("Sequence", "Modified sequence", "Modified sequence_verbose", "Modifications")])
  #View(EV[wMod, c("Sequence", "Modified sequence", "Charge", "Mass", "m/z")])
  #
  # Retention time (in minutes, as in MaxQuant)
  EV$"Retention time" <- as.numeric(DIANN$RT)
  EV$"Retention time (predicted)" <- as.numeric(DIANN$Predicted.RT)
  EV$iRT <- as.numeric(DIANN$iRT)
  EV$"iRT (predicted)" <- as.numeric(DIANN$Predicted.iRT)
  EV$"Retention time (start)" <- as.numeric(DIANN$RT.Start)
  EV$"Retention time (end)" <- as.numeric(DIANN$RT.Stop)
  EV$"Retention length" <- EV$"Retention time (end)" - EV$"Retention time (start)"
  #
  EV$Reverse <- if ("Decoy" %in% colnames(DIANN)) {
    c("", "+")[as.logical(DIANN$Decoy)+1L]
  } else { "" }
  if ("Score" %in% colnames(DIANN)) {
    EV$Score <- as.numeric(DIANN$CScore)
  }
  if ("Decoy.CScore" %in% colnames(DIANN)) {
    EV$"Delta score" <- as.numeric(DIANN$CScore) - as.numeric(DIANN$Decoy.CScore)
  }
  #
  if (PG.Genes) {
    EV$"Leading proteins" <- DIANN$Protein.Group # (Will be removed later anyway in the classic workflow, since we perform protein inference again)
    kol <- base::intersect(c("PG.Quantity", "PG.Normalised", "PG.MaxLFQ", "Genes.Quantity",
                             "Genes.Normalised", "Genes.MaxLFQ", "Genes.MaxLFQ.Unique"),
                           colnames(DIANN))
    EV[, sub(" Normalised$", " Normalised quantity", gsub("\\.", " ", kol))] <- DIANN[, kol]
    kol <- c("Protein.Q.Value", "PG.Q.Value", "Global.PG.Q.Value", "GG.Q.Value")
    EV[, sub(" Q Value$", " q-value", gsub("\\.", " ", kol))] <- DIANN[, kol]
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
  kol <- base::intersect(kol, colnames(DIANN))
  EV[, gsub("Ms1", "MS1", gsub("\\.", " ", kol))] <- DIANN[, kol]
  #
  # Quantitative values
  EV$"MS1 Area" <- as.numeric(DIANN$Ms1.Area)
  hasQuantUMS <- ("Precursor.Quantity" %in% colnames(DIANN))
  useQuantUMS <- hasQuantUMS
  # if (hasQuantUMS && (misFun(useQuantUMS) || (!is.logical(useQuantUMS)) || (length(useQuantUMS) != 1L) || is.na(useQuantUMS))) {
  #   useQuantUMS <- c(TRUE, FALSE)[match(svDialogs::dlg_message("QuantUMS output columns detected, using those instead of Ms1.Area for Intensity?",
  #                                                              "yesno")$res, c("yes", "no"))]
  # } else { useQuantUMS <- FALSE }
  if (useQuantUMS) {
    #cat("   QuantUMS output columns detected, using those as Intensity instead of Ms1.Area.\n   (For your convenience the output Evidence (PSMs) table also contains a separate Ms1.Area column.)\n")
    # If the QuantUMS output is present, we use it by default: it should be higher quality than "Ms1.Area" (MS1-based quant only).
    # We do not use the normalized version by default, in case normalization would actually be unsuitable for this particular experiment
    # ("normalization should never be on by default!")
    EV$Precursor.Quantity <- as.numeric(DIANN$Precursor.Quantity)
    EV$Intensity <- EV$Precursor.Quantity
  } else {
    EV$Intensity <- EV$"MS1 Area"
  }
  IntColz <- base::intersect(c("Fragment.Quant.Corrected", "Fragment.Quant.Raw"),
                             colnames(DIANN))
  if (length(IntColz)) {
    EV$"MS2 intensities" <- DIANN[[IntColz[1L]]]
  }
  # Fragment annotations are not provided by default,
  # but each PSM for a given peptidoform*Charge combination always has the same number of values.
  # Assuming the order is then also always the same (which appears to be the case), we have all we need for MS2-based quant.
  #
  # Careful: this does not work where 2 different DiaNN files are combined, because the number of fragments will be different in different runs (best on which are found).
  #
  # Below commented code if you want to check
  #tst <- aggregate(EV$`MS2 intensities`, list(EV$`Modified sequence`, EV$Charge), \(x) {
  #  nchar(x) - nchar(gsub(";", "", x))
  #})
  #tst2 <- sapply(tst$x, \(x) { length(unique(x)) })
  #stopifnot(max(tst2) == 1L)
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
  #
  digPat <- paste(c("[", digPattern, "]"), collapse = "")
  tmp <- gsub(paste0(digPat, "$"), "", EV$Sequence)
  EV$"Missed cleavages" <- nchar(tmp) - nchar(gsub(digPat, "", tmp))
  EV$"Potential contaminant" <- ""
  EV$Search_ID <- if ("Search_ID" %in% colnames(DIANN)) {
    DIANN$PSMs_file # In cases where the input is a hybrid report created by the merging script in .../Utils
  } else {
    sort(DIANN_fl, decreasing = TRUE)[1L] # Provide default: there can be one or 2 files, the .tsv or .parquet, we want just one (ideally the .tsv)
  }
  # Optional: filter by delta score
  if (("Delta score" %in% colnames(EV)) && (!is.na(Min.Delta.Score))) {
    EV <- EV[EV$"Delta score" > Min.Delta.Score,]
  }
  #
  EV$id <- NULL
  EV <- cbind(data.frame(id = 1L:nrow(EV)),
              EV)
  Res <- list(Evidence = EV,
              PTMs = allPTMs,
              QuantUMS = useQuantUMS)
  #
  if (stopCl) { parallel::stopCluster(cl) }
  return(Res)
}
