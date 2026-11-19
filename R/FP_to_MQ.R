#' FP_to_MQ
#'
#' @description 
#' Reads FragPipe's workflow and manifest files.\cr
#' From them, detects and parses psm.tsv tables (one per Experiment) to a single MaxQuant evidence.txt-like table.\cr
#' As with the other functions of this type, the idea is not to get perfect conversion, but close enough that the data can be fed into this package's analysis scripts.
#' 
#' @param FP_Workflow Fragpipe workflow. Parsed to identify FragPipe's output directory, search parameters, incl. modifications.
#' @param FP_Manifest FragPipe's samples manifest, required for mapping raw files to samples ("Experiments" in the MaxQuant sense).
#' @param Min.Delta.Score From which minimum delta score should we report a PSM? Default = -Inf (we assume that those few peptides with low delta score are because randomly the decoy database will include real-like peptides)
#' @param UniProtIDs If TRUE (default), long protein IDs are re-processed to only get the accession.
#' @param FailIfNoQuant If TRUE, the function will fail with an error if it does not detect quantitative values (which, in FragPipe, are optional). It is useful in most case (otherwise a workflow may fail further down the road), but I am setting this to FALSE by default for now to avoid stupid fails for ident-only workflows. 
#' @param N.clust A limit on the number of vCPUs to use. If left as NULL (default), uses the number of available clusters - 1, to a minimum of 1.
#' @param N.reserved Default = 1. Number of reserved vCPUs the function is not to use. Note that for obvious reasons it will always use at least one.
#' @param cl Already have a cluster handy? Why spend time making a new one, which on top of that may invalidate the old one. Just pass it along!
#' @param OpenSearch Is it an open search? If TRUE, will also ignore Open Search-type columns even if they are present. If TRUE, will try to detect and process them - but this is probably still a bit buggy.
#' 
#' @returns
#' A list of (at least) 4 items:\cr
#' - "Evidence": a PSMs table similar (nor identical!) to the evidence table MaxQuant creates\cr
#' - "PTMs": a modifications table\cr
#' - "FracMap": a map of the different MS raw files and their relationship to Experiments\cr
#' - "WorkFlow": the FragPipe workflow
#' 
#' @examples
#' FP.Wrkflw <- paste0(wd, "/fragpipe.workflow")
#' FP.Mnfst <- paste0(wd, "/fragpipe-files.fp-manifest")
#' temp <- FP_to_MQ(FP.Wrkflw, FP.Mnfst)
#' ev <- temp$Evidences
#' Modifs <- temp$PTMs
#' 
#' @import data.table
#' @export

FP_to_MQ <- function(FP_Workflow,
                     FP_Manifest,
                     Min.Delta.Score = -Inf,
                     UniProtIDs = TRUE,
                     FailIfNoQuant = FALSE,
                     N.clust,
                     N.reserved = 1L,
                     cl,
                     OpenSearch = FALSE) {
  TESTING <- FALSE
  #DefArg(FP_to_MQ); cl = parClust;TESTING <- TRUE
  #OpenSearch = TRUE
  #FP_Workflow <- fpWorkflowFl_i; FP_Manifest <- fpManifestFl_i
  #FP_Workflow <- FP_WorkflowFl; FP_Manifest <- FP_ManifestFl
  #FP_Workflow <- wrkflw; FP_Manifest <- mnfst
  #FP_Workflow <- paste0(wd, "/fragpipe.workflow"); FP_Manifest <- paste0(wd, "/fragpipe-files.fp-manifest")
  #
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
    cat("     Making fresh cluster...\n")
    cl <- parallel::makeCluster(N.clust, type = "SOCK")
    stopCl <- TRUE
  }
  N.clust <- length(cl)
  #
  # Load FragPipe parameter file
  if ((!TESTING) && (missing("FP_Workflow") || (!file.exists(FP_Workflow)))) {
    stop("Invalid parameter \"FP_Workflow\": file not found!")
  }
  FP_Wrkflw <- readr::read_lines(FP_Workflow)
  pat <- topattern("diann.run-dia-nn=")
  isActuallyDIANN <- as.logical(toupper(sub(pat, "", grep(pat, FP_Wrkflw, value = TRUE))))
  #
  FP_Dir <- gsub("\\\\", "",
                 gsub("\\\\\\\\", "/",
                      sub("^workdir=", "",
                          grep("^workdir=", FP_Wrkflw, value = TRUE))))
  FixPaths <- FALSE
  if (!dir.exists(FP_Dir)) {
    FixPaths <- TRUE
    warning(paste0("FragPipe's output directory was not found: \"", FP_Dir, "\", assuming the parent folder of the FragPipe Workflow file instead!"))
    FP_Dir <- dirname(FP_Workflow)
  }
  if (isActuallyDIANN) {
    tst1 <- dir.exists(paste0(FP_Dir, "/diann-output"))
    diannRep <- paste0(FP_Dir, "/diann-output/report.tsv")
    tst2 <- file.exists(diannRep)
    if (!tst1) {
      warning("DiaNN quantitation was run by FragPipe but we cannot locate the DiaNN output folder! Ignoring...")
      isActuallyDIANN <- FALSE
    } else {
      if (!tst2) {
        warning("DiaNN quantitation was run by FragPipe but we cannot locate the DiaNN report! Ignoring...")
        isActuallyDIANN <- FALSE
      } else {
        warning("DiaNN quantitation was run, we will deliver DiaNN PSMs by default but also process and output FragPipe PSMs and PTMs as the \"FP_Evidence\" and \"FP_PTMs\" items in the function's output list.")
        #DIANN <- DIANN_to_MQ(diannRep, cl = cl)
      }
    }
  }
  pat <- topattern("tmtintegrator.run-tmtintegrator=")
  isTMT <- as.logical(sub(pat, "", grep(pat, FP_Wrkflw, value = TRUE)))
  if (isTMT) {
    pat <- topattern("tmtintegrator.channel_num=TMT-")
    TMTplex <- as.integer(sub(pat, "", grep(pat, FP_Wrkflw, value = TRUE)))
    TMTtblFl <- paste0(FP_Dir, "/experiment_annotation.tsv")
    if (!file.exists(TMTtblFl)) {
      stop(paste0("This is a TMT dataset but file \"", TMTtblFl, "\" could not be found in the results folder!"))
    }
    TMTtbl <- read.delim(TMTtblFl, check.names = FALSE)
    TMTtbl$plex[is.na(TMTtbl$plex)] <- TMTplex
    Channels <- c(126L, unlist(lapply(127L:134L, \(x) { paste0(as.character(x), c("N", "C")) })), "135N")
    TMTtbl$channel_code <- match(TMTtbl$channel, Channels) - 1L
    tmtKols <- as.character(TMTtbl$sample)
  }
  #
  # Parse FragPipe's samples manifest
  if ((!TESTING) && (missing("FP_Manifest") || (!file.exists(FP_Manifest)))) {
    stop("Invalid parameter \"FP_Manifest\": file not found!")
  }
  FP_Mnfst <- suppressWarnings(read.delim(FP_Manifest, header = FALSE))
  colnames(FP_Mnfst) <- c("Path", "Experiment", "Replicate", "Data type")
  FP_Mnfst$Path <- gsub("\\\\", "/", FP_Mnfst$Path)
  if (FixPaths) { FP_Mnfst$Path <- gsub(".+/", paste0(FP_Dir, "/"), FP_Mnfst$Path) }
  # The bit below is bad and causes inconsistent behaviour! Commented but kept for the wall of shame!
  #if (sum(file.exists(FP_Mnfst$Path)) < nrow(FP_Mnfst)) {
  #  warning("Raw files could not be located, keeping only file names without full path!")
  #  FP_Mnfst$Path <- gsub(".+/", "", FP_Mnfst$Path)
  #}
  FP_Mnfst$"File name" <- gsub(".+/|\\.(raw|mz(X?ML|BIN)|mgf|d)$", "", FP_Mnfst$Path)
  FP_Mnfst$Experiment[is.na(FP_Mnfst$Experiment)] <- ""
  Exp <- Samples <- FP_Mnfst$Experiment
  Rep <- FP_Mnfst$Replicate
  if (length(Rep[!is.na(Rep)])) {
    Samples <- data.frame(Exp = Exp, Rep = Rep)
    # Rep seems to be duplicated somehow, but only sometimes... Why oh why!!!
    Samples1 <- apply(Samples, 1L, \(x) {
      x <- c(x[[1L]], x[[2L]])
      x <- x[!is.na(x)]
      return(paste(x, collapse = "_"))
    })
    Samples2 <- apply(Samples, 1L, \(x) {
      x <- c(x[[1L]], x[[2L]], x[[2L]])
      x <- x[!is.na(x)]
      return(paste(x, collapse = "_"))
    })
    tst1 <- sum(dir.exists(paste0(FP_Dir, "/", Samples1))) == length(Samples1)
    tst2 <- sum(dir.exists(paste0(FP_Dir, "/", Samples2))) == length(Samples2)
    stopifnot(tst1 + tst2 == 1L)
    Samples <- get(paste0("Samples", 1L:2L)[which(c(tst1, tst2))])
  }
  FP_Mnfst$Samples <- Samples
  #
  FP_PSMs <- unique(gsub("//", "/", paste0(FP_Dir, "/", Samples, "/psm.tsv")))
  #FP_Fas <- paste0(FP_Dir, "/", Samples, "/protein.fas")
  tst <- file.exists(FP_PSMs)
  w <- which(!tst)
  l <- length(w)
  FP_Mnfst2 <- FP_Mnfst
  if (l) {
    #SamplesNo <- unique(SamplesNo[w]); l <- length(SamplesNo)
    #if (l > 1) { SamplesNo <- paste0(paste(SamplesNo[1L:(l-1L)], collapse = ", "), " and ", SamplesNo[l]) }
    #plur <- c("", "s")[(l > 1)+1]
    #warning(paste0("The \"psm.tsv\" file", plur, " for sample", plur, " ", SamplesNo,
    #               " could not be found, skipping missing sample", plur, "!"))
    FP_PSMs <- FP_PSMs[tst]
    #FP_Fas <- FP_Fas[tst]
    Samples <- Samples[tst]
    FP_Mnfst2 <- FP_Mnfst2[tst,]
    Exp <- Exp[Exp %in% FP_Mnfst2$Experiment]
    Rep <- Rep[Rep %in% FP_Mnfst2$Replicate]
  }
  isWellBehaved <- TRUE
  if (!length(FP_PSMs)) {
    stopifnot(length(Samples) == 0L,
              length(Exp) == 0L,
              dim(FP_Mnfst2)[1L] == 0L)
    FP_PSMs <- paste0(FP_Dir, "/psm.tsv")
    stopifnot(file.exists(FP_PSMs))
    isWellBehaved <- FALSE
  }
  #tst <- file.exists(FP_Fas)
  #w <- which(!tst); l <- length(w)
  #if (l) {
  #  SamplesNo <- unique(SamplesNo[w])
  #  l <- length(SamplesNo)
  #  SamplesNo <- .commAnd(SamplesNo)
  #  plur <- c("", "s")[(l > 1L)+1L]
  #  warning(paste0("The \"protein.fas\" file", plur, " for sample", plur, " ", SamplesNo,
  #                 " could not be found, skipping missing sample", plur, "!"))
  #  FP_PSMs <- FP_PSMs[tst]
  #  FP_Fas <- FP_Fas[tst]
  #  Samples <- Samples[tst]
  #  FP_Mnfst <- FP_Mnfst[tst,]
  #  Exp <- Exp[Exp %in% FP_Mnfst$Experiment]
  #  Rep <- Rep[Rep %in% FP_Mnfst$Replicate]
  #}
  # Update values
  Exp <- unique(Exp)
  Samples <- unique(Samples)
  Rep <- unique(Rep)
  FP_PSMs <- gsub("//", "/", paste0(FP_Dir, "/", Samples, "/psm.tsv"))
  #
  # Process PTMs
  # - Load and post-process UniMod object
  utils::data(modifications, package = "PTMods", envir = environment())
  UniMod <- modifications
  rm(modifications)
  # As of early 2023 this package did not have the whole of Unimod,
  # in particular TMT16plex was missing.
  # This is fixed now, so code is commented.
  # sites <- c("K", "N-term", "N-term", "H", "S", "T")
  # sites2 <- c("K", "N-term", "P-N-term", "H", "S", "T")
  # where <- c("Anywhere", "Any N-term", "Protein N-term", "Anywhere", "Anywhere", "Anywhere")
  # tmtp <- data.frame(Id = paste0("TMTpro:", sites2),
  #                    UnimodId = 2016,
  #                    Name = "TMTpro",
  #                    Description = "TMTpro 16plex Tandem Mass Tag®",
  #                    Composition = "H(25) C(8) 13C(7) N 15N(2) O(3)",
  #                    AvgMass = 304.3127,
  #                    MonoMass = 304.207146,
  #                    Site = sites2,
  #                    Position = where,
  #                    Classification = "Isotopic label",
  #                    SpecGroup = 1,
  #                    NeutralLoss = FALSE,
  #                    LastModified = "2021-07-22 21:42:31",
  #                    Approved = FALSE,
  #                    Hidden = FALSE)
  # View(rbind(UniMod[grep("TMT", UniMod$Description),], tmtp))
  # UniMod <- rbind(UniMod, tmtp)
  # Remove:
  # - substitutions
  gS <- grepl("^[A-Z][a-z]{2}->[A-Z][a-z]{2} substitution$", UniMod$Description)
  UniMod_S <- UniMod[gS,]
  UniMod <- UniMod[!gS,]
  # - deltas
  gD <- grepl("^Delta:", UniMod$Name)
  UniMod_D <- UniMod[gD,]
  UniMod <- UniMod[!gD,]
  # For now these are removed as they are normally not expected.
  # However, we can keep them aside in case there is no match to main UniMod (not implemented currently).
  #
  parseMods <- \(mods, fixed = FALSE) { #mods <- varMods
    fixInd <- fixed + 1L
    x1 <- vapply(mods, \(x) { as.numeric(x[[1L]]) }, 1)
    x2 <- vapply(mods, \(x) { paste(x[2L:(length(x) - 2L)], collapse = "") }, "")
    x3 <- vapply(mods, \(x) { as.logical(toupper(x[[length(x) - 1L]])) }, TRUE)
    x4 <- vapply(mods, \(x) { as.integer(x[[length(x)]]) }, 1L)
    mods <- data.frame("Mass delta" = x1,
                       "Site" = x2,
                       "Enabled" = x3,
                       "Max occurences" = x4,
                       "Type" = c("Variable", "Fixed")[fixInd],
                       check.names = FALSE)
    mods <- mods[mods$"Mass delta" != 0,]
    mods <- mods[mods$Enabled,]
    mods$Enabled <- NULL
    mods$Site <- gsub(" ?\\([^)]+\\)", "", mods$Site) # Why do I need this again?
    mods$ID <- paste0(c("var", "fix")[fixInd], 1L:nrow(mods))
    tmp <- lapply(strsplit(mods$Site, ""), \(x) { #x <- strsplit(mods$Site, "")[4L]
      x <- unlist(x)
      w <- which(x == "n")
      if (length(w)) {
        x[w + 1L] <- paste0("n", x[w + 1L])
        x <- x[-w]
      }
      w <- which(x == "c")
      if (length(w)) {
        x[w + 1L] <- paste0("c", x[w + 1L])
        x <- x[-w]
      }
      w <- which(x == "[")
      if (length(w)) {
        x[w + 1L] <- paste0("[", x[w + 1L])
        x <- x[-w]
      }
      w <- which(x == "]")
      if (length(w)) {
        x[w + 1L] <- paste0("]", x[w + 1L])
        x <- x[-w]
      }
      return(x)
    })
    tmp <- listMelt(tmp, 1L:nrow(mods))
    mods <- mods[tmp$L1,]
    mods$Site <- tmp$value
    #
    # This bit I think was just for varmods...
    mods$Site[mods$Site == "N-Term Peptide"] <- "n^"
    mods$Site[mods$Site == "C-Term Peptide"] <- "^c"
    mods$Site[mods$Site == "N-Term Protein"] <- "[^"
    mods$Site[mods$Site == "C-Term Protein"] <- "]^"
    #
    mods <- mods[grep("^site_[0-9]+$", mods$Site, invert = TRUE),]
    return(mods)
  }
  allPTMs <- list()
  # - Fixed modifications
  pat <- topattern("msfragger.table.fix-mods=")
  g <- grep(pat, FP_Wrkflw, value = TRUE)
  fixMods <- unlist(strsplit(sub(pat, "", g), ";"))
  fixMods <- strsplit(fixMods, ",")
  allPTMs[["Fixed"]] <- parseMods(fixMods, TRUE)
  # - Variable modifications
  pat <- topattern("msfragger.table.var-mods=")
  g <- grep(pat, FP_Wrkflw, value = TRUE)
  varMods <- unlist(strsplit(sub(pat, "", g), ";"))
  varMods <- strsplit(varMods, ",")
  allPTMs[["Variable"]] <- parseMods(varMods, FALSE)
  # - Combine
  allPTMs <- plyr::rbind.fill(allPTMs)
  tst <- vapply(strsplit(allPTMs$Site, ""), \(x) {
    x <- x[x != " "]
    sum(!x %in% c("[", "n", "^", "c", "]", AA)) == 0L
  }, TRUE)
  if (sum(!tst)) {
    warning("Unexpected result when parsing FragPipe's modification sites, investigate!")
  }
  allPTMs <- allPTMs[tst,]
  allPTMs$Position <- "Anywhere"
  allPTMs$Position[grep("^\\[", allPTMs$Site)] <- "Protein N-term"
  allPTMs$Position[grep("^n", allPTMs$Site)] <- "Any N-term"
  allPTMs$Position[grep("^\\]", allPTMs$Site)] <- "Protein C-term"
  allPTMs$Position[grep("^c", allPTMs$Site)] <- "Any C-term"
  allPTMs$Site_long <- allPTMs$Site
  allPTMs$Site_long[allPTMs$Site_long == "[^"] <- "P-N-term"
  allPTMs$Site_long[allPTMs$Site_long == "]^"] <- "P-C-term"
  allPTMs$Site_long[grep("^n", allPTMs$Site_long)] <- "N-term"
  allPTMs$Site_long[grep("^c", allPTMs$Site_long)] <- "C-term"
  allPTMs$AA <- gsub("\\[|\\]|\\^|n|c", "", allPTMs$Site)
  allPTMs$AA[allPTMs$AA == "*"] <- ""
  w <- which(allPTMs$AA != "")
  allPTMs$Site_long[w] <- allPTMs$AA[w] # Restore AA specificity
  allPTMs$"Mass precision" <- nchar(gsub(".*\\.|0+$", "", allPTMs$"Mass delta"))
  allPTMs <- allPTMs[order(allPTMs$`Mass delta`, decreasing = FALSE),]
  tmp <- lapply(1L:nrow(allPTMs), \(i) { #i <- 1L #i <- 2L  #i <- 3L #i <- 4L
    prec <- allPTMs$`Mass precision`[i]
    dMass <- allPTMs$`Mass delta`[i]
    sites <- allPTMs$Site_long[i]
    # if (sites %in% c("N-term", "P-N-term")) { sites <- c("N-term", "P-N-term") } else {
    #   if (sites %in%  c("C-term", "P-C-term")) { sites <- c("C-term", "P-C-term") }
    # }
    # Sometimes UniMod is more generic:
    if (sites == "P-N-term") { sites <- c("N-term", "P-N-term") } else {
      if (sites == "P-C-term") { sites <- c("C-term", "P-C-term") }
    }
    pos <- allPTMs$Position[i]
    w <- which((abs(round(UniMod$MonoMass, prec) - round(as.numeric(dMass), prec)) <= 10^(1 - prec))
               &(UniMod$Site %in% sites)
               &(UniMod$Position == pos))
    # if (!length(w)) {
    #   w <- which((abs(round(UniMod$MonoMass+1.007276, prec) - round(as.numeric(dMass), prec)) <= 10L^(1-prec))
    #              &(UniMod$Site %in% sites)
    #              &(UniMod$Position == pos))
    # }
    if (length(w)) {
      #if (length(unique(UniMod$Description[w])) != 1L) {
      #  print(x)
      #  View(UniMod[w,])
      #  res <- NA
      #}
      res <- UniMod[w, c("UnimodId", "Site", "Position")]
      res$Position <- as.character(res$Position)
      res$Position[res$Position == "Anywhere"] <- ""
      res$Site[res$Site %in% paste0(c("N", "C"), "-term")] <- ""
      w <- which((res$Position != "") & (res$Site != ""))
      if (length(w)) {
        res$Site[w] <- sub("^ ", "", paste0(apply(res[w, c("Site", "Position")], 1L, paste, collapse = " ("), ")"))
      }
      w <- which((res$Position != "") & (res$Site == ""))
      if (length(w)) { res$Site[w] <- res$Position[w] }
      res <- aggregate(res$Site, list(res$UnimodId), list)
      colnames(res) <- c("UniMod", "Site_long")
      res$Row <- i
    } else {
      res <- data.frame(Row = i,
                        UniMod = NA)
      res$Site_long <- sites
    }
    return(res)
  })
  tmp <- plyr::rbind.fill(tmp)
  tmp2 <- listMelt(tmp$Site_long, 1L:nrow(tmp), c("Site_long", "tmpRow"))
  tmp2[, c("UniMod", "Row")] <- tmp[tmp2$tmpRow, c("UniMod", "Row")]
  kol <- c("ID", "Mass delta", "Site", "Max occurences", "Type", "Position", "AA")
  tmp2[, kol] <- allPTMs[tmp2$Row, kol]
  tmp2$"Full name" <- ""
  w <- which(!is.na(tmp2$UniMod))
  tmp2$"Full name"[w] <- UniMod$Name[match(tmp2$UniMod[w], UniMod$UnimodId)]
  w <- which(is.na(tmp2$UniMod))
  if (length(w)) {
    tmp2$"Full name"[w] <- paste0("Unknown", 1L:length(w), "_", c("", "+")[(tmp2$`Mass delta`[w] > 0) + 1L],
                                  tmp2$`Mass delta`[w])
  }
  kol <- base::setdiff(colnames(tmp2), c("tmpRow", "Row", "Mass delta", "Full name"))
  tmp2 <- aggregate(tmp2[, kol], list(tmp2$"Mass delta", tmp2$"Full name"), list)
  colnames(tmp2)[1L:2L] <- c("Mass delta", "Full name")
  tmp2$"Max occurences" <- vapply(tmp2$"Max occurences", max, 1L)
  tmp2$Type <- sapply(tmp2$Type, unique)
  stopifnot(!inherits(tmp2$Type, "list"))
  tmp2$UniMod <- sapply(tmp2$UniMod, unique)
  stopifnot(!inherits(tmp2$UniMod, "list"))
  allPTMs <- tmp2 
  # - Marks
  allPTMs$Mark <- tolower(substr(allPTMs$"Full name", 1L, 2L))
  ## Sometimes, some marks are duplicates, e.g. if you searched for "Acetyl (Protein N-term)" and "Acetyl (K)" together!
  ## We want to fix this so that each modification has a unique mark:
  #test <- aggregate(allPTMs$Mark, list(allPTMs$Mark), length)
  test <- aggregate(allPTMs$UniMod, list(allPTMs$Mark), \(x) { length(unique(x)) })
  W <- which(test$x > 1L)
  #allPTMs$Mark <- allPTMs$"Old mark"
  if (length(W)) {
    allPTMs$"Old mark" <- allPTMs$Mark
    for (i in W) {
      #i <- W[1L]
      w <- which(allPTMs$Mark == test$Group.1[i])
      m <- allPTMs[w,]
      # Simple case: multiple instances of a same isotopic label (e.g. 1, 2, 3 or 4 times 15N)
      labelTst <- tolower(gsub("^Label:[0-9]+|\\(|\\)", "", m$"Full name"))
      if (length(unique(labelTst)) == length(labelTst)) {
        allPTMs$Mark[w] <- labelTst
      } else {
        m$AA[lengths(m$AA) == 0L] <- "X"
        r <- { if ("Acetyl" %in% m$"Full name") { which(m$"Full name" == "Acetyl") } else { 1L } } # r is the one we will keep without changing it
        s <- base::setdiff(1L:nrow(m), r)
        tst <- lapply(s, \(x) {
          paste0(tolower(m$AA[[x]]), substr(m$Mark[[x]], 1L, 1L))
        })
        tst <- lapply(1L:length(tst), \(x) {
          rs <- tst[[x]]
          rs[!rs %in% allPTMs$Mark]
        })
        l <- length(tst)
        if (l > 1L) {
          for (i in 2L:l) {
            tst[[i]] <- tst[[i]][!tst[[i]] %in% unlist(tst[1L:(i - 1L)])]
          }
        }
        tst <- vapply(tst, \(x) {
          x <- unlist(x)
          x <- { if (length(x)) { x[1L] } else { "That didnae work, did it?" } }
          return(x)
        }, "")
        tst2 <- ((tst %in% allPTMs$Mark) | (tst == "That didnae work, did it?"))
        w2 <- which(!tst2)
        m$Mark[s][w2] <- tst[w2]
        w1 <- which(tst2)
        if (length(w1)) {
          # not tested
          s1 <- s[w1]
          rs <- c()
          kount <- 1L
          char <- c(as.character(0L:9L), letters)
          taken <- unique(c(allPTMs$Mark, m$Mark))
          for (j in s1) {
            tst <- paste0(tolower(m$AA[s1]), char[kount])
            while ((tst %in% taken) && (kount < length(char))) {
              kount <- kount + 1L
              tst <- paste0(tolower(m$AA[s1]), char[kount])
            }
            if (kount == length(char)) {
              stop("I am really out of options here, never thought this would go this far! Check the code just in case, this should never happen.")
            } else { rs <- c(rs, tst) }
          }
          m$Mark[[s1]] <- rs
        }
        allPTMs[w,] <- m
      }
    }
  }
  #
  allPTMs$"Mass shift" <- UniMod$MonoMass[match(allPTMs$UniMod, UniMod$UnimodId)]
  if ("Mass delta" %in% colnames(allPTMs)) {
    w <- which(is.na(allPTMs$"Mass shift"))
    allPTMs$"Mass shift"[w] <- allPTMs$"Mass delta"[w]
    # "Mass delta" is the one which was entered in the search - it may be slightly inaccurate
    # "Mass shift" is the monoisotopic which we pulled from UniMod and should be better
    w <- which(!is.na(allPTMs$"Mass shift"))
    tst <- allPTMs$"Mass shift"[w] - allPTMs$"Mass delta"[w]
    if (max(tst) > 0.001) {
      warning("Discrepancies detected between the masses entered in the FragPipe search and the mapped UniMod PTMs.")
    } # Could probably be more precise
  }
  test <- aggregate(allPTMs$"Full name", list(allPTMs$Mark), \(x) { length(unique(x)) })
  if (max(test$x) > 1L) {
    warning("The algorithm was not able to assign unique 2-character old-MaxQuant style short marks to each modification, so some modified sequences may be ambiguous.")
  }
  # Remove parentheses in names: parentheses are also used to wrap around names in the full sequence,
  # nested parentheses are ugly and break the code further down.
  g <- grep("\\([0-9]+\\)$", allPTMs$"Full name")
  if (length(g)) {
    allPTMs$"Full name"[g] <- vapply(g, \(i) {
      paste0(gsub("\\([0-9]+\\)$", "", allPTMs$"Full name"[i]),
             "*",
             gsub(".*\\(|\\)$", "", allPTMs$"Full name"[i]))
    }, "")
  }
  g <- grep("\\(.+\\)$", allPTMs$"Full name")
  if (length(g)) {
    allPTMs$"Full name"[g] <- gsub("\\(", "{", gsub("\\)", "}", allPTMs$"Full name"[g]))
  }
  allPTMs$"Type of search" <- "Closed"
  #
  # Load PSM files
  if (isWellBehaved) {
    PSMs <- lapply(1L:length(Samples), \(x) { #x <- 1L
      res <- data.table::fread(FP_PSMs[x], integer64 = "numeric", check.names = FALSE,
                               data.table = FALSE, fill = TRUE, sep = "\t")
      if (nrow(res)) {
        res$Sample <- Samples[x]
        res$PSMs_file <- FP_PSMs[x]
        return(res)
      }
      return(NA)
    })
    w <- which(vapply(PSMs, is.data.frame, TRUE))
    PSMs <- plyr::rbind.fill(PSMs[w])
  } else {
    PSMs <- data.table::fread(FP_PSMs, integer64 = "numeric", check.names = FALSE,
                              data.table = FALSE, fill = Inf, sep = "\t")
    PSMs2Mnsfst <- match(gsub_Rep("\\\\", "/", PSMs$`Spectrum File`),
                         paste0(FP_Dir, "/interact-", FP_Mnfst$`File name`, ".pep.xml"))
    w <- which(is.na(PSMs2Mnsfst))
    l <- length(w)
    if (l) {
      if (l == nrow(PSMs)) { PSMs2Mnsfst <- rep(1L, nrow(PSMs)) } else {
        stop("I don't know what to do!!!")
      }
    }
    PSMs$Sample <- FP_Mnfst$Path[PSMs2Mnsfst]
    #PSMs$Sample <- FP_Mnfst$FP_PSMs[PSMs2Mnsfst] # Corrected bug... based on code but maybe this breaks something downstream -> check!!!
    PSMs$PSMs_file <- FP_Mnfst$FP_PSMs[PSMs2Mnsfst]
  }
  # Check formats
  kols <- base::intersect(c("Charge", "Number of Enzymatic Termini", "Number of Missed Cleavages",
                            "Protein Start", "Protein End"),
                          colnames(PSMs))
  for (kol in kols) { PSMs[[kol]] <- as.integer(PSMs[[kol]]) }
  kols <- base::intersect(c("Retention", "Observed Mass", "Calibrated Observed Mass", 
                            "Observed M/Z", "Calibrated Observed M/Z", "Calculated Peptide Mass", 
                            "Calculated M/Z", "Delta Mass", "Expectation", "Hyperscore", 
                            "Nextscore", "PeptideProphet Probability", "Intensity", 
                            "Best Score with Delta Mass", "Best Score without Delta Mass"),
                          colnames(PSMs))
  if (FailIfNoQuant &&
      ((!"Intensity" %in% kols) ||
       (!sum(is.finite(PSMs$Intensity))) ||
       (!sum(PSMs$Intensity[is.finite(PSMs$Intensity)])))) {
    if (isTMT) {
      warning("No MS1 quantitative information identified in PSMs table! Using sum of reporter intensities!")
      PSMs$Intensity <- rowSums(PSMs[, tmtKols, drop = FALSE], na.rm = TRUE)
    } else {
      stop("No MS1 quantitative information identified in PSMs table!\nEither rerun FragPipe with IonQuant (or FreeQuant) turned on, or set FailIfNoQuant to FALSE if your workflow only requires identifications.")
    }
  }
  for (kol in kols) { PSMs[[kol]] <- as.numeric(PSMs[[kol]]) }
  # Let's now deal with those pesky un-localised open-search delta masses!!!
  if (OpenSearch) {
    OpenSearch <- ("Observed Modifications" %in% colnames(PSMs))
    if (OpenSearch) {
      PSMs$"Observed Modifications"[is.na(PSMs$"Observed Modifications")] <- ""
      OpenSearch <- sum(PSMs$`Observed Modifications` == "") > 0L
    }
  }
  if (OpenSearch) {
    # We need to map open search modification names to fold changes
    #allPTMs <- allPTMs[allPTMs$"Type of search" == "Closed",]
    wMod2 <- which(PSMs$`Observed Modifications` != "")
    row2orig <- gsub("; Mod", "; _;_Mod", gsub(", Mod", ", _;_Mod", PSMs$`Observed Modifications`[wMod2]))
    row2orig <- listMelt(strsplit(row2orig, "[,;] _;_"), wMod2, c("Full", "row"))
    row2orig$Full <- sub("^Mod[0-9]+: ", "", grep("^Mod[0-9]+: ", row2orig$Full, value = TRUE))
    u <- unique(row2orig$Full)
    tmp1 <- data.frame(Full = u,
                       Full_name = sub(" \\((Theoretical|PeakApex): .+", "", u))
    tmp1 <- aggregate(tmp1$Full, list(tmp1$Full_name), unique)
    colnames(tmp1) <- c("Full name", "Full")
    tmp1 <- tmp1[!tmp1$`Full name` %in% c("First isotopic peak", "Second isotopic peak", "Third isotopic peak",
                                          "Isotopic peak error"),]
    massKol <- c("Mass delta (UniMod)", "Mass delta (us)", "Mass delta")
    tmp1[, massKol] <- NA_real_
    tmp2 <- listMelt(strsplit(tmp1$`Full name`, "/"))
    tmp3 <- data.frame(value = unique(tmp2$value))
    tmp3[, c("UniMod", massKol[1L])] <- UniMod[match(tmp3$value, UniMod$Description), c("UnimodId", "MonoMass")]
    tmp3[[massKol[2L]]] <- NA
    g <- grep("^Unannotated mass-shift ", tmp3$value)
    if (length(g)) {
      tmp3[g, massKol[2L]] <- as.numeric(sub(" \\(.*", "", sub("^Unannotated mass-shift ", "", tmp3$value[g])))
    }
    wNA <- is.na(tmp3[[massKol[2L]]])
    g1 <- which(grepl("^Deletion of [A-Z]$", tmp3$value) & wNA)
    g2 <- which(grepl("^Addition of [A-Z]$", tmp3$value) & wNA)
    if (length(g1)) {
      tmp3[g1, massKol[2L]] <- - AA_table$`Residue monoisotopic mass`[match(sub("^Deletion of ", "", tmp3$value[g1]), AA_table$AA)]
    }
    if (length(g2)) {
      tmp3[g2, massKol[2L]] <- AA_table$`Residue monoisotopic mass`[match(sub("^Addition of ", "", tmp3$value[g2]), AA_table$AA)]
    }
    tmp3[[massKol[3L]]] <- vapply(1L:nrow(tmp3), \(i) { #i <- 1L #i <- 47L
      x <- as.numeric(tmp3[i, massKol[1L:2L]])
      x[!is.na(x)][1L]
    }, 1)
    tmp2[, c("UniMod", massKol)] <- tmp3[match(tmp2$value, tmp3$value),
                                         c("UniMod", massKol)]
    tmp2 <- data.table::as.data.table(tmp2)
    massFun <- \(x) {
      if (!is.na(x[1L])) { return(x[1L]) }
      return(x[2L])
    }
    tmp2 <- tmp2[,
                 .(UniMod = unique(UniMod),
                   `Mass delta (UniMod)` = massFun(`Mass delta (UniMod)`),
                   `Mass delta (us)` = massFun(`Mass delta (us)`),
                   `Mass delta` = massFun(`Mass delta`)),
                 by = .(row = L1)]
    tmp2 <- as.data.frame(tmp2)
    tmp1[, c("UniMod", massKol)] <- tmp2[match(1L:nrow(tmp1), tmp2$row),
                                         c("UniMod", massKol)]
    tmp1$Theor <- vapply(tmp1$Full, \(x) {
      x <- grep("Theoretical: ", x, value = TRUE)
      if (!length(x)) { return(NA_real_) }
      unique(as.numeric(sub("(, .+)?\\)$", "", sub(".+Theoretical: ", "", x))))
    }, 1)
    wNA <- which(!is.na(tmp1$Theor))
    if (length(wNA)) { tmp1$"Mass delta (FP)"[wNA] <- tmp1$Theor[wNA] }
    wNA <- which(is.na(tmp1$"Mass delta"))
    if (length(wNA)) { tmp1$"Mass delta"[wNA] <- tmp1$"Mass delta (FP)"[wNA] }
    tmp1 <- tmp1[, c("Full name", "Full", "Mass delta", "UniMod")]
    tmp1$Position <- "Unknown"
    #tmp1$Enabled <- NA
    #mp4$"Max occurences" <- NA_real_
    tmp1$AA <- lapply(1L:nrow(tmp1), \(x) { c() })
    tmp1$Type <- "Delta mass"
    tmp1$"Type of search" <- "Open"
    #
    ###############################################################################
    # Make 2-character PTM marks                                                  #
    ###############################################################################
    # Now, for historical reasons, I have stuck with the "2 characters PTM mark" format which older MaxQuant versions used to use.
    # It is debatable whether that still makes sense. Probably not.
    # Yeah...
    # But, I have no time to review all of my code base to change that just yet.
    # So... this means we need to generate marks for all of these PTMs...
    # Oh dear...
    # Well, let's try:
    markOpt <- c(letters, 0L:9L, "!", "?", ".", "*", "#", "^", "<", ">", "&", "$", "@", "-", "+", ";", ",", ":", "~", "=", "%")
    # Of course, we cannot ever hope to get 2-letter marks for all possible PTMs, because, e.g.:
    # nrow(UniMod)  # = 3704
    # length(markOpt)^2L # = 3025 -> we need more characters than 2, or more options per character (use different encoding?)
    # 
    #
    tmp1$Mark <- ""
    # - Marks from the same PTM from the closed search
    tmp1 <- tmp1[order(tmp1$UniMod),]
    w <- which(tmp1$UniMod %in% allPTMs$UniMod)
    if (length(w)) {
      tmp1$UniMod[w] <- allPTMs$Mark[match(tmp1$UniMod[w], allPTMs$UniMod)]
    }
    # - Marks based on amino-acid deletion or addition
    #   - Deletion
    g1 <- grep("^(.+/)?Deletion of [A-Z](/.+)?$", tmp1$`Full name`)
    #tmp1$`Full name`[g1]
    if (length(g1)) {
      tmp <- listMelt(strsplit(tmp1$`Full name`[g1], "/"), g1)
      tmp <- tmp[grep("^Deletion of [A-Z]$", tmp$value),]
      tmp$value <- sub("^Deletion of ", "", tmp$value)
      tmp$Mark <- paste0("-", tolower(tmp$value))
      m <- match(g1, tmp$L1)
      tmp1$Mark[g1] <- tmp$Mark[m]
      tmp1$AA[g1] <- tmp$value[m]
      tmp1$Position[g1] <- "Anywhere"
    }
    #   - Addition
    g2 <- grep("^(.+/)?Addition of [A-Z](/.+)?$", tmp1$`Full name`)
    if (length(g2)) {
      tmp <- listMelt(strsplit(tmp1$`Full name`[g2], "/"), g2)
      tmp <- tmp[grep("^Addition of [A-Z]$", tmp$value),]
      tmp$value <- sub("^Addition of ", "", tmp$value)
      tmp$Mark <- paste0("+", tolower(tmp$value))
      m <- match(g2, tmp$L1)
      tmp1$Mark[g2] <- tmp$Mark[m]
      tmp1$AA[g2] <- lapply(g2, \(x) { NULL })
      tmp1$Position[g2] <- "Anywhere"
    }
    # - Marks from the mass shift
    w <- which(tmp1$Mark == "")
    if (length(w)) {
      tmp1$Mark_opt2 <- ""
      w1a <- w[(tmp1$`Mass delta`[w] > 0) & (tmp1$`Mass delta`[w] < 10)]
      #w1b <- w[(tmp1$`Mass delta`[w] > 0) & (tmp1$`Mass delta`[w] < 100) & (tmp1$`Mass delta`[w] >= 10)] # No! We don't want confusion with UniMod IDs
      w2 <- w[(tmp1$`Mass delta`[w] < 0) & (tmp1$`Mass delta`[w] > -10)]
      if (length(w1a)) { tmp1$Mark_opt2[w1a] <- paste0("+", as.character(round(tmp1$`Mass delta`[w1a]))) }
      #if (length(w1b)) { tmp1$Mark_opt2[w1b] <- as.character(round(tmp1$`Mass delta`[w1b])) }
      if (length(w2)) { tmp1$Mark_opt2[w2] <- as.character(round(tmp1$`Mass delta`[w2])) }
      wNR <- w[!vapply(w, \(i) {
        if (i == 1L) { return(TRUE) }
        tmp1$Mark_opt2[i] %in% c(allPTMs$Mark, tmp1$Mark[1L:(i-1L)], tmp1$Mark_opt2[1L:(i-1L)], "!m")
      }, TRUE)]
      tmp1$Mark[wNR] <- tmp1$Mark_opt2[wNR]
    }
    # - Marks from the shortened PTM name
    w <- which(tmp1$Mark == "")
    if (length(w)) {
      tmp1$Mark_opt1 <- ""
      tmp1$Mark_opt1[w] <- tolower(substr(tmp1$`Full name`[w], 1L, 2L))
      wNR <- w[!vapply(w, \(i) {
        if (i == 1L) { return(TRUE) }
        tmp1$Mark_opt1[i] %in% c(allPTMs$Mark, tmp1$Mark[1L:(i-1L)], tmp1$Mark_opt1[1L:(i-1L)], "!m")
      }, TRUE)]
      tmp1$Mark[wNR] <- tmp1$Mark_opt1[wNR]
    }
    # - Automatically generated marks
    w <- which(tmp1$Mark == "")
    if (length(w)) {
      markOpt2 <- do.call(paste, c(as.data.frame(gtools::permutations(length(markOpt), 2L, markOpt)), sep = ""))
      markOpt2 <- base::setdiff(markOpt2, c(allPTMs$Mark, tmp1$Mark, "!m"))
      tmp1$Mark[w] <- markOpt2[1L:length(w)]
    }
    # - Generic marks for the rest... this will need a fix
    w <- which(tmp1$Mark == "")
    if (length(w)) {
      tmp1$Mark[w] <- "!m"
    }
    #
    tmp1$"Full name" <- gsub(" ", "_", gsub(",", ".", gsub("\\(", "{", gsub("\\)", "}", tmp1$"Full name"))))
    tmp1$"Type of search" <- "Open"
    tmp1$ID <- paste("open", as.character(1L:nrow(tmp1)))
    tmp1$"Mass shift" <- tmp1$"Mass delta"
    mrgCol <- c("Mass delta", "Full name", "UniMod", "ID", "Type", "Position", "AA", "Mark", "Mass shift", "Type of search")
    #View(tmp1[, mrgCol])
    allPTMs <- plyr::rbind.fill(allPTMs,
                                tmp1[, mrgCol])
    orig2name <- listMelt(tmp1$Full, tmp1$`Full name`, c("Orig", "Name"))
    orig2name$Mark <- tmp1$Mark[match(orig2name$Name, tmp1$`Full name`)]
    row2orig[, c("Name", "Mark")] <- orig2name[match(row2orig$Full, orig2name$Orig), c("Name", "Mark")]
    row2orig <- row2orig[!is.na(row2orig$Name),]
    wMod2 <- unique(row2orig$row)
  }
  # Filter peptides with ambiguous/non-classic amino acids 
  test <- gsub(paste(AA, collapse = "|"), "", PSMs$Peptide)
  w <- which(test == "")
  if (length(w) < nrow(PSMs)) {
    warning("Removing peptides with unknown amino acids!s")
    PSMs <- PSMs[w,]
  }
  # Create MQ-like data.frame
  kol <- c("Spectrum File", "Peptide", "Modified Peptide", "Charge", "Hyperscore",
           "Spectrum", "Intensity", "Retention", "Protein")
  #print(base::setdiff(kol, colnames(PSMs)))
  stopifnot(sum(!kol %in% colnames(PSMs)) == 0L)
  EV <- data.frame(Sequence = PSMs$Peptide,
                   Length = nchar(PSMs$Peptide), 
                   Charge = PSMs$Charge,
                   Intensity = PSMs$Intensity)
  #
  pKol <- base::intersect(c("PeptideProphet Probability", "Probability"),
                          colnames(PSMs))[1L]
  if (!length(pKol)) { stop("No PSMs quality column detected!") } else {
    EV$PEP <- 1 - PSMs[[pKol]]
  }
  a1 <- strsplit(PSMs$Protein, "[,;] ?")
  a2 <- strsplit(PSMs$"Mapped Proteins", "[,;] ?")
  a1 <- listMelt(a1, ColNames = c("prot", "row"))
  a2 <- listMelt(a2, ColNames = c("prot", "row"))
  tmp <- rbind(a1, a2)
  if (UniProtIDs) {
    tmp$prot <- sub("\\|[^|]+$", "",
                    sub("^[a-z]{2}\\|", "", tmp$prot))
  }
  tmp <- tmp[order(tmp$row, tmp$prot),]
  tmp <- data.table::as.data.table(tmp)
  tmp <- tmp[, .(prot2 = list(unique(prot))),
             by = .(row = row, prot = prot)]
  tmp <- tmp[, .(prots = paste(prot2, collapse = ";")),
             by = .(row = row)]
  tmp <- as.data.frame(tmp)
  EV$Proteins <- tmp$prots[match(1L:nrow(EV), tmp$row)]
  if (isWellBehaved) {
    EV$Experiment <- PSMs$Sample
    #EV$"Raw file path" <- FP_Mnfst$Path[match(EV$Experiment, FP_Mnfst$Samples)] # Buggy line
    EV$"Raw file" <- gsub_Rep("^(.*/)?interact-", "",
                              gsub_Rep("\\\\", "/",
                                       gsub_Rep("(\\.mod)?\\.pep\\.xml$", ".d", PSMs$`Spectrum File`)))
    u <- unique(EV$"Raw file")
    if (identical(u,  "interact.d")) {
      tst <- setNames(lapply(FP_Mnfst$`File name`, \(fl) { #fl <- FP_Mnfst$`File name`[1L]
        pat <- topattern(fl)
        grep(pat, PSMs$Spectrum)
      }), FP_Mnfst$`File name`)
      rg <- 1L:nrow(EV)
      tst2 <- identical(setNames(sort(unlist(tst)), NULL), rg)
      if (tst2) {
        tst <- listMelt(tst, ColNames = c("row", "file"))
        EV$"Raw file" <- tst$file[match(rg, tst$row)]
      } else {
        stop("Houston, we have a problem...")
        # This indicates that we were not able to map raw files to PSMs!!!
      }
    }
    EV$"Raw file" <- gsub_Rep(".+/|\\.(raw|mz(X?ML|BIN)|mgf|d)$", "", EV$"Raw file")
    u <- unique(EV$"Raw file")
    w <- lapply(FP_Mnfst$`File name`, \(x) { which(u == x) })
    l <- lengths(w)
    stopifnot(max(l) == 1L) # Would indicate that we have non unique file names, which this cannot deal with!
    EV$"Raw file path" <- FP_Mnfst$Path[match(EV$"Raw file", FP_Mnfst$`File name`)]
  } else {
    EV$Experiment <- FP_Mnfst$Experiment[PSMs2Mnsfst]
    EV$"Raw file path" <- FP_Mnfst$Path[PSMs2Mnsfst]
    EV$"Raw file" <- gsub_Rep(".+/|\\.(raw|mz(X?ML|BIN)|mgf|d)$", "", EV$"Raw file path")
  }
  # Modified sequence
  # - rewrite this with data.table on the model of DIANN_to_MQ()!!!
  #View(PSMs[, c("Assigned Modifications", "Observed Modifications")])
  #bckpMods <- EV[, c("Modified sequence", "Modified sequence_verbose", "Modifications")]
  #bckpMods[, c("Assigned Modifications", "Observed Modifications")] <- PSMs[, c("Assigned Modifications", "Observed Modifications")]
  #View(bckpMods)
  wMod <- which(PSMs$`Assigned Modifications` != "")
  tmp <- listMelt(strsplit(PSMs$`Assigned Modifications`[wMod], ", "), wMod, c("PTMsite", "row"))
  tmp <- do.call(paste, c(tmp[, c("row", "PTMsite")], sep = "_"))
  tmp <- data.frame(where = sub("\\(.*", "", tmp),
                    shift = sub(".*\\(", "", sub("\\)", "", tmp)))
  tmp$shift <- as.numeric(tmp$shift)
  tmpModTbl <- allPTMs[allPTMs$`Type of search` == "Closed",]
  m <- match(tmp$shift, tmpModTbl$`Mass delta`)
  prec <- 5L
  w <- which(is.na(m))
  while (length(w) || (prec > 1L)) {
    m[w] <- match(round(tmp$shift[w], prec), round(tmpModTbl$`Mass delta`, prec))
    prec <- prec-1L
    w <- which(is.na(m))
  }
  if (length(w)) { stop("Houston, we have a problem!") }
  tmp$Mark <- tmpModTbl$Mark[m]
  tmp$Name <- tmpModTbl$`Full name`[m]
  tmpMods <- data.table::as.data.table(tmp[, c("Mark", "Name", "where", "shift")])
  tmpVerb <- data.table::copy(tmpMods)
  tmpMods <- tmpMods[tmpMods$Mark %in% tmpModTbl$Mark[tmpModTbl$Type == "Variable"],]
  tmpMods <- tmpMods[, .(Mark = paste(Mark, collapse = ",")),
                     by = .(where = where)]
  tmpVerb <- tmpVerb[, .(shift = sum(shift, na.rm = TRUE),
                         Name = paste(Name, collapse = ",")),
                     by = .(where = where)]
  tmpVerb$Mark <- ""
  w <- which(tmpVerb$where %in% tmpMods$where)
  m <- match(tmpVerb$where[w], tmpMods$where)
  tmpVerb$Mark[w] <- tmpMods$Mark[m]
  # "z" sorting trick: I am lucky this is not an amino acid!
  tmp <- as.data.frame(tmpVerb)
  w <- which(tmp$Mark != "")
  tmp$Mark[w] <- paste0("z(", tmp$Mark[w], ")") 
  tmp$Name <- paste0("z(", tmp$Name, ")")
  tmp$row <- as.integer(sub("_.*", "", tmp$where))
  tmp$pos <- sub(".*_", "", tmp$where)
  tmp$pos <- vapply(sub("[A-Z]$", "", tmp$pos), \(x) {
    if (x == "N-term") { return(0) }
    if (x == "C-term") { return(Inf) }
    as.integer(x)
  }, 1)
  tmp2 <- listMelt(strsplit(paste0("<", EV$Sequence[wMod], ">"), ""), wMod, c("Mark", "row"))
  wN <- which(tmp2$Mark == "<")
  tmp2$Mark[wN] <- "N-term"
  wC <- which(tmp2$Mark == ">")
  tmp2$Mark[wC] <- "C-term"
  tmp2$pos <- base::sequence(base::tabulate(tmp2$row)) - 1L
  tmp2$pos2 <- as.character(tmp2$pos)
  tmp2$pos2[tmp2$Mark %in% c("N-term", "C-term")] <- ""
  tmp2$row2 <- paste0(as.character(tmp2$row), "_")
  tmp2$where <- do.call(paste, c(tmp2[, c("row2", "pos2", "Mark")], sep = ""))
  stopifnot(sum(!tmp$where %in% tmp2$where) == 0L)
  tmp2$shift <- 0
  tmp2$Name <- tmp2$Mark
  kol <- c("pos", "row", "Mark", "Name", "shift")
  tmpMods <- rbind(tmp2[, kol],
                   tmp[, kol])
  tmpMods <- data.table::as.data.table(tmpMods)
  tmpMods <- tmpMods[order(abs(tmpMods$shift)),]
  tmpMods <- tmpMods[order(tmpMods$pos),]
  tmpMods <- tmpMods[order(tmpMods$row),]
  w <- which(tmpMods$shift != 0)
  tmpMods$Mark[w] <- sub("^z", "", tmpMods$Mark[w])
  tmpMods$Name[w] <- sub("^z", "", tmpMods$Name[w])
  w <- which(tmpMods$Mark %in% c("N-term","C-term"))
  tmpMods$Mark[w] <- tmpMods$Name[w] <- "_"
  tmpVerb <- data.table::copy(tmpMods)
  tmpJust <- data.table::copy(tmpMods)
  tmpMods <- tmpMods[, .(modSeq = paste(Mark, collapse = "")),
                     by = .(row = row)]
  tmpVerb <- tmpVerb[, .(verbSeq = paste(Name, collapse = "")),
                     by = .(row = row)]
  tmpJust <- tmpJust[tmpJust$shift != 0,]
  tmpJust$Name <- gsub("[()]", "", tmpJust$Name)
  tmpJust$Name <- strsplit(tmpJust$Name, ",")
  tmpJust <- listMelt(tmpJust$Name, tmpJust$row)
  tmpJust <- data.table::as.data.table(tmpJust)
  tmpJust <- tmpJust[, .(n = .N),
                     by = .(row = L1, Mods = value)]
  tmpJust$n <- as.character(tmpJust$n)
  tmpJust$modTag <- tmpJust$Mods
  w2p <- which(tmpJust$n != "1")
  tmpJust$modTag[w2p] <- do.call(paste, c(tmpJust[w2p, c("n", "Mods")], sep = " "))
  tmpJust <- tmpJust[order(tmpJust$Mods),]
  tmpJust <- tmpJust[, .(Mods = paste(modTag, collapse = ",")), by = .(row = row),]
  tmpJust <- as.data.frame(tmpJust)
  tmpJust <- tmpJust[order(tmpJust$row),]
  EV$"Modified sequence_verbose" <- EV$"Modified sequence" <- paste0("_", EV$Sequence, "_")
  EV$"Modifications" <- ""
  m <- match(wMod, tmpMods$row)
  EV$"Modified sequence"[wMod] <- tmpMods$modSeq[m]
  m <- match(wMod, tmpVerb$row)
  EV$"Modified sequence_verbose"[wMod] <- tmpVerb$verbSeq[m]
  m <- match(wMod, tmpJust$row)
  EV$"Modifications"[wMod] <- tmpJust$Mods[m]
  #View(EV[wMod, c("Sequence", "Modified sequence", "Modified sequence_verbose", "Modifications")])
  if (OpenSearch) {
    EV$`Modified sequence`[wMod2] <- sub("_$", "(!m)_", EV$`Modified sequence`[wMod2])
    tmpVerb <- data.table::as.data.table(row2orig)
    tmpVerb <- tmpVerb[, .(n = .N),
                       by = .(row = row, Mods = Name)]
    tmpVerb$n <- as.character(tmpVerb$n)
    tmpVerb$modTag <- tmpVerb$Mods
    w2p <- which(tmpVerb$n != "1")
    tmpVerb$modTag[w2p] <- do.call(paste, c(tmpVerb[w2p, c("n", "Mods")], sep = " "))
    tmpVerb <- tmpVerb[order(tmpVerb$Mods),]
    tmpVerb <- tmpVerb[, .(Mods = paste(modTag, collapse = ",")), by = .(row = row),]
    tmpVerb <- as.data.frame(tmpVerb)
    tmpVerb <- tmpVerb[order(tmpVerb$row),]
    tmp <- data.frame(A = sub("_$", "", EV$`Modified sequence_verbose`[wMod2]),
                      B = paste0("(", tmpVerb$Mods[match(wMod2, tmpVerb$row)], ")_"))
    EV$`Modified sequence_verbose`[wMod2] <- do.call(paste, c(tmp, sep = ""))
    EV$Modifications[wMod2] <- sub("^,", "", paste0(EV$Modifications[wMod2], ",open search: ", tmp$B))
  }
  #
  #tst <- unique(sub("^[0-9]+ ", "", unlist(strsplit(EV$Modifications[wMdSq2], ","))))
  #sum(!tst %in% allPTMs$"Full name")
  # Mass and mass error columns
  EV$"m/z" <- PSMs$"Calibrated Observed M/Z"
  EV$"Theoretical m/z" <- PSMs$"Calculated M/Z"
  EV$Mass <- (EV$"m/z" - 1.007276466879) * EV$Charge # This used to be incorrect
  EV$"Uncalibrated - Calibrated m/z [Da]" <- PSMs$"Observed M/Z" - PSMs$"Calibrated Observed M/Z" # Somehow I get all 0s here?!
  EV$"Uncalibrated - Calibrated m/z [ppm]" <- 1000000*EV$"Uncalibrated - Calibrated m/z [Da]"/EV$"m/z"
  EV$"Mass error [Da]" <- PSMs$"Delta Mass"
  EV$"Mass error [ppm]" <- 1000000*EV$"Mass error [Da]"/EV$Mass
  EV$"Uncalibrated mass error [Da]" <- PSMs$"Observed Mass" - EV$Mass
  EV$"Uncalibrated mass error [ppm]" <- 1000000*EV$"Uncalibrated mass error [Da]"/EV$Mass
  if (OpenSearch) {
    EV$"Open search: Obs. modifications" <- PSMs$"Observed Modifications"
    EV$"Open search: Best score with Delta Mass" <- PSMs$"Best Score with Delta Mass"    
    EV$"Open search: Best score without Delta Mass" <- PSMs$"Best Score without Delta Mass"
  }
  # Retention time (in minutes, as in MaxQuant)
  EV$"Retention time" <- as.numeric(PSMs$Retention)/60
  # For now, I could not find ways to extract more retention time information
  # (it is not in the pepXML files)
  # So the below 3 columns cannot be created:
  #EV$"Retention time (start)" <- ...
  #EV$"Retention time (end)" <- ...
  #EV$"Retention length" <- ...
  #
  EV$Score <- PSMs$Hyperscore
  EV$"Delta score" <- PSMs$Hyperscore - PSMs$Nextscore
  EV$Expectation <- PSMs$Expectation
  #
  #tst <- as.data.frame(t(sapply(strsplit(PSMs$Spectrum, "\\."), unlist)))
  #sum(tst$V2 != tst$V3) == 0L
  #tst$V3 <- NULL
  #tst$V2 <- as.numeric(tst$V2)
  #tst$V4 <- as.numeric(tst$V4)
  #aggregate(tst$V2, list(tst$V1), min)
  #aggregate(tst$V2, list(tst$V1), max)
  #aggregate(PSMs$Retention, list(PSMs$`Spectrum File`), min)
  #aggregate(PSMs$Retention, list(PSMs$`Spectrum File`), max)
  #max(tst$V4)
  #tst2 <- aggregate(tst$V4, list(tst$V1, tst$V2), length) 
  #max(tst2$x)
  #
  #aggregate(tst$V4, list(tst$V4), length)
  # Conclusion:
  # - V2 is a scan number
  # - After checking, I can confirm that V4 is charge 
  EV$"MS/MS scan number" <- as.numeric(gsub("^[^.]+\\.|\\..+", "", PSMs$Spectrum))
  # Unfortunately, at least for DIA data these are pseudo-MS/MS spectrum numbers generated
  # by DIA-Umpire or DIA-Tracer.
  # They cannot be matched to those in the DiaNN report, if DiaNN was used for quantitation.
  #
  # Type:
  #   Even with files acquired in DIA mode, the MSFragger search is actually a DDA search on simulated DDA data
  #   which is created using DIA-Umpire/DIA Tracer.
  EV$Type <- "MULTI-MSMS"
  # Ion mobility" TO DO!
  EV$"Missed cleavages" <- PSMs$`Number of Missed Cleavages`
  EV$"Potential contaminant" <- "" # This can be fixed later
  EV$Reverse <- "" # This can be fixed later
  # Optional: filter by delta score
  EV <- EV[EV$`Delta score` > Min.Delta.Score,]
  #
  EV$id <- NULL
  EV <- cbind(data.frame(id = 1L:nrow(EV)), EV)
  #
  if (isTMT) {
    EV[, paste0("Reporter intensity ", TMTtbl$channel_code)] <- PSMs[, tmtKols]
  }
  EV$Search_ID <- FP_Workflow
  if (isActuallyDIANN) {
    DIANN <- DIANN_to_MQ(diannRep, cl = cl)
    Res <- list(Evidence = DIANN$Evidence,
                PTMs = DIANN$PTMs,
                FracMap = FP_Mnfst,
                WorkFlow = FP_Wrkflw,
                FP_Evidence = EV,
                FP_PTMs = allPTMs)
  } else {
    Res <- list(Evidence = EV,
                PTMs = allPTMs,
                FracMap = FP_Mnfst,
                WorkFlow = FP_Wrkflw)
  }
  if (isTMT) { Res[["TMT_annotations"]] <- TMTtbl }
  Res$Fasta <- sub("^ *database\\.db-path=", "", grep("fasta", FP_Wrkflw, value = TRUE))
  #
  if (stopCl) { parallel::stopCluster(cl) }
  return(Res)
}
