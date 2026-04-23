#' FP_to_MQ
#'
#' @description 
#' Reads FragPipe's workflow and manifest files
#' From them, detects and parses psm.tsv tables (one per Experiment) to a single MaxQuant evidence.txt-like table.
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
#' A list of (at least) 4 items:
#' - "Evidence": a PSMs table similar (nor identical!) to the evidence table MaxQuant creates
#' - "PTMs": a modifications table
#' - "FracMap": a map of the different MS raw files and their relationship to Experiments
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
  #FP_Workflow <- fpWorkflowFl_i; FP_Manifest <- fpManifestFl_i
  #FP_Workflow <- FP_WorkflowFl; FP_Manifest <- FP_ManifestFl
  #FP_Workflow <- wrkfl; FP_Manifest <- mnfst
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
  if ((is.null(cl))||(!inherits(cl, "cluster"))) {
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
  if ((!TESTING)&&((missing("FP_Workflow"))||(!file.exists(FP_Workflow)))) {
    stop("Invalid parameter \"FP_Workflow\": file not found!")
  }
  FP_Wrkflw <- readLines(FP_Workflow)
  pat <- topattern("diann.run-dia-nn=")
  isActuallyDIANN <- as.logical(toupper(gsub(pat, "", grep(pat, FP_Wrkflw, value = TRUE))))
  #
  FP_Dir <- gsub("\\\\", "",
                 gsub("\\\\\\\\", "/",
                      gsub("^workdir=", "",
                           grep("^workdir=", FP_Wrkflw, value = TRUE))))
  FixPaths <- FALSE
  if (!dir.exists(FP_Dir)) {
    FP_Dir_old <- FP_Dir
    FP_Dir <- dirname(FP_Workflow)
    if (FP_Dir == FP_Dir_old) {
      FP_Dir <- getwd()
      warning(paste0("FragPipe's output directory was not found: \"", FP_Dir, "\", assuming the local folder instead!"))
    } else {
      warning(paste0("FragPipe's output directory was not found: \"", FP_Dir, "\", assuming the parent folder of the FragPipe Workflow file instead!"))
    }
    FixPaths <- TRUE
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
  isTMT <- as.logical(gsub(pat, "", grep(pat, FP_Wrkflw, value = TRUE)))
  if (isTMT) {
    pat <- topattern("tmtintegrator.channel_num=TMT-")
    TMTplex <- as.integer(gsub(pat, "", grep(pat, FP_Wrkflw, value = TRUE)))
    TMTtblFl <- paste0(FP_Dir, "/experiment_annotation.tsv")
    if (!file.exists(TMTtblFl)) {
      stop(paste0("This is a TMT dataset but file \"", TMTtblFl, "\" could not be found in the results folder!"))
    }
    TMTtbl <- read.delim(TMTtblFl, check.names = FALSE)
    TMTtbl$plex[which(is.na(TMTtbl$plex))] <- TMTplex
    Channels <- c(126L, unlist(lapply(127L:134L, \(x) { paste0(as.character(x), c("N", "C")) })), "135N")
    TMTtbl$channel_code <- match(TMTtbl$channel, Channels) - 1L
    tmtKols <- as.character(TMTtbl$sample)
  }
  #
  # Parse FragPipe's samples manifest
  if ((!TESTING)&&((missing("FP_Manifest"))||(!file.exists(FP_Manifest)))) {
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
  FP_Mnfst$"File name" <- gsub(".+/|\\.((raw)|(mz(X?ML|BIN))|(mgf)|(d))$", "", FP_Mnfst$Path)
  FP_Mnfst$Experiment[which(is.na(FP_Mnfst$Experiment))] <- ""
  Exp <- Samples <- FP_Mnfst$Experiment
  Rep <- FP_Mnfst$Replicate
  if (length(Rep[which(!is.na(Rep))])) {
    Samples <- data.frame(Exp = Exp, Rep = Rep)
    # Rep seems to be duplicated somehow, but only sometimes... Why oh why!!!
    Samples1 <- apply(Samples, 1L, \(x) {
      x <- c(x[[1L]], x[[2L]])
      x <- x[which(!is.na(x))]
      return(paste(x, collapse = "_"))
    })
    Samples2 <- apply(Samples, 1L, \(x) {
      x <- c(x[[1L]], x[[2L]], x[[2L]])
      x <- x[which(!is.na(x))]
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
    FP_PSMs <- FP_PSMs[which(tst)]
    #FP_Fas <- FP_Fas[which(tst)]
    Samples <- Samples[which(tst)]
    FP_Mnfst2 <- FP_Mnfst2[which(tst),]
    Exp <- Exp[which(Exp %in% FP_Mnfst2$Experiment)]
    Rep <- Rep[which(Rep %in% FP_Mnfst2$Replicate)]
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
  #  SamplesNo <- unique(SamplesNo[w]); l <- length(SamplesNo)
  #  if (l > 1L) { SamplesNo <- paste0(paste(SamplesNo[1L:(l-1L)], collapse = ", "), " and ", SamplesNo[l]) }
  #  plur <- c("", "s")[(l > 1L)+1L]
  #  warning(paste0("The \"protein.fas\" file", plur, " for sample", plur, " ", SamplesNo,
  #                 " could not be found, skipping missing sample", plur, "!"))
  #  FP_PSMs <- FP_PSMs[which(tst)]
  #  FP_Fas <- FP_Fas[which(tst)]
  #  Samples <- Samples[which(tst)]
  #  FP_Mnfst <- FP_Mnfst[which(tst),]
  #  Exp <- Exp[which(Exp %in% FP_Mnfst$Experiment)]
  #  Rep <- Rep[which(Rep %in% FP_Mnfst$Replicate)]
  #}
  # Update values
  Exp <- unique(Exp)
  Samples <- unique(Samples)
  Rep <- unique(Rep)
  FP_PSMs <- gsub("//", "/", paste0(FP_Dir, "/", Samples, "/psm.tsv"))
  #
  # Process PTMs
  # - Load and post-process UniMod object
  data(modifications, package = "PTMods")
  UniMod <- modifications
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
  UniMod_S <- UniMod[which(gS),]
  UniMod <- UniMod[which(!gS),]
  # - deltas
  gD <- grepl("^Delta:", UniMod$Name)
  UniMod_D <- UniMod[which(gD),]
  UniMod <- UniMod[which(!gD),]
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
    mods <- mods[which(mods$"Mass delta" != 0),]
    mods <- mods[which(mods$Enabled),]
    mods$Enabled <- NULL
    mods$Site <- gsub(" ?\\([^\\)]+\\)", "", mods$Site) # Why do I need this again?
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
    mods$Site[which(mods$Site == "N-Term Peptide")] <- "n^"
    mods$Site[which(mods$Site == "C-Term Peptide")] <- "^c"
    mods$Site[which(mods$Site == "N-Term Protein")] <- "[^"
    mods$Site[which(mods$Site == "C-Term Protein")] <- "]^"
    #
    mods <- mods[grep("^site_[0-9]+$", mods$Site, invert = TRUE),]
    return(mods)
  }
  Modifs <- list()
  # - Fixed modifications
  g <- grep("^msfragger\\.table\\.fix-mods=", FP_Wrkflw, value = TRUE)
  fixMods <- unlist(strsplit(gsub("^msfragger\\.table\\.fix-mods=", "", g), ";"))
  fixMods <- strsplit(fixMods, ",")
  Modifs[["Fixed"]] <- parseMods(fixMods, TRUE)
  # - Variable modifications
  g <- grep("^msfragger\\.table\\.var-mods=", FP_Wrkflw, value = TRUE)
  varMods <- unlist(strsplit(gsub("^msfragger\\.table\\.var-mods=", "", g), ";"))
  varMods <- strsplit(varMods, ",")
  Modifs[["Variable"]] <- parseMods(varMods, FALSE)
  # - Combine
  Modifs <- plyr::rbind.fill(Modifs)
  tst <- vapply(strsplit(Modifs$Site, ""), \(x) {
    x <- x[which(x != " ")]
    sum(!x %in% c("[", "n", "^", "c", "]", AA)) == 0L
  }, TRUE)
  if (sum(!tst)) {
    warning("Unexpected result when parsing FragPipe's modification sites, investigate!")
  }
  Modifs <- Modifs[which(tst),]
  Modifs$Position <- "Anywhere"
  Modifs$Position[grep("^\\[", Modifs$Site)] <- "Protein N-term"
  Modifs$Position[grep("^n", Modifs$Site)] <- "Any N-term"
  Modifs$Position[grep("^\\]", Modifs$Site)] <- "Protein C-term"
  Modifs$Position[grep("^c", Modifs$Site)] <- "Any C-term"
  Modifs$Site_long <- Modifs$Site
  Modifs$Site_long[which(Modifs$Site_long == "[^")] <- "P-N-term"
  Modifs$Site_long[which(Modifs$Site_long == "]^")] <- "P-C-term"
  Modifs$Site_long[grep("^n", Modifs$Site_long)] <- "N-term"
  Modifs$Site_long[grep("^c", Modifs$Site_long)] <- "C-term"
  Modifs$AA <- gsub("\\[|\\]|\\^|n|c", "", Modifs$Site)
  Modifs$AA[which(Modifs$AA == "*")] <- ""
  w <- which(Modifs$AA != "")
  Modifs$Site_long[w] <- Modifs$AA[w] # Restore AA specificity
  Modifs$"Mass precision" <- nchar(gsub(".*\\.|0+$", "", Modifs$"Mass delta"))
  Modifs <- Modifs[order(Modifs$`Mass delta`, decreasing = FALSE),]
  tmp <- lapply(1L:nrow(Modifs), \(i) { #i <- 1L #i <- 2L  #i <- 3L #i <- 4L
    prec <- Modifs$`Mass precision`[i]
    dMass <- Modifs$`Mass delta`[i]
    sites <- Modifs$Site_long[i]
    # if (sites %in% c("N-term", "P-N-term")) { sites <- c("N-term", "P-N-term") } else {
    #   if (sites %in%  c("C-term", "P-C-term")) { sites <- c("C-term", "P-C-term") }
    # }
    # Sometimes UniMod is more generic:
    if (sites == "P-N-term") { sites <- c("N-term", "P-N-term") } else {
      if (sites == "P-C-term") { sites <- c("C-term", "P-C-term") }
    }
    pos <- Modifs$Position[i]
    w <- which((abs(round(UniMod$MonoMass, prec) - round(as.numeric(dMass), prec)) <= 10^(1 - prec))
               &(UniMod$Site %in% sites)
               &(UniMod$Position == pos))
    # if (!length(w)) {
    #   w <- which((abs(round(UniMod$MonoMass+1.007276, prec) - round(as.numeric(dMass), prec)) <= 10^(1-prec))
    #              &(UniMod$Site %in% sites)
    #              &(UniMod$Position == pos))
    # }
    if (length(w)) {
      #if (length(unique(UniMod$Description[w])) != 1) {
      #  print(x)
      #  View(UniMod[w,])
      #  res <- NA
      #}
      res <- UniMod[w, c("UnimodId", "Site", "Position")]
      res$Position <- as.character(res$Position)
      res$Position[which(res$Position == "Anywhere")] <- ""
      res$Site[which(res$Site %in% paste0(c("N", "C"), "-term"))] <- ""
      w <- which((res$Position != "")&(res$Site != ""))
      if (length(w)) {
        res$Site[w] <- gsub("^ ", "", paste0(apply(res[w, c("Site", "Position")], 1L, paste, collapse = " ("), ")"))
      }
      w <- which((res$Position != "")&(res$Site == ""))
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
  tmp2[, kol] <- Modifs[tmp2$Row, kol]
  tmp2$"Full name" <- ""
  w <- which(!is.na(tmp2$UniMod))
  tmp2$"Full name"[w] <- UniMod$Name[match(tmp2$UniMod[w], UniMod$UnimodId)]
  w <- which(is.na(tmp2$UniMod))
  if (length(w)) {
    tmp2$"Full name"[w] <- paste0("Unknown", 1L:length(w), "_", c("", "+")[(tmp2$`Mass delta`[w] > 0) + 1L],
                                  tmp2$`Mass delta`[w])
  }
  kol <- colnames(tmp2)[which(!colnames(tmp2) %in% c("tmpRow", "Row", "Mass delta", "Full name"))]
  tmp2 <- aggregate(tmp2[, kol], list(tmp2$"Mass delta", tmp2$"Full name"), list)
  colnames(tmp2)[1L:2L] <- c("Mass delta", "Full name")
  tmp2$"Max occurences" <- vapply(tmp2$"Max occurences", max, 1L)
  tmp2$Type <- sapply(tmp2$Type, unique)
  stopifnot(!"list" %in% class(tmp2$Type))
  tmp2$UniMod <- sapply(tmp2$UniMod, unique)
  stopifnot(!"list" %in% class(tmp2$UniMod))
  Modifs <- tmp2 
  # - Marks
  Modifs$Mark <- tolower(substr(Modifs$"Full name", 1L, 2L))
  ## Sometimes, some marks are duplicates, e.g. if you searched for "Acetyl (Protein N-term)" and "Acetyl (K)" together!
  ## We want to fix this so that each modification has a unique mark:
  #test <- aggregate(Modifs$Mark, list(Modifs$Mark), length)
  test <- aggregate(Modifs$UniMod, list(Modifs$Mark), \(x) { length(unique(x)) })
  W <- which(test$x > 1L)
  #Modifs$Mark <- Modifs$"Old mark"
  if (length(W)) {
    Modifs$"Old mark" <- Modifs$Mark
    for (i in W) {
      #i <- W[1L]
      w <- which(Modifs$Mark == test$Group.1[i])
      m <- Modifs[w,]
      # Simple case: multiple instances of a same isotopic label (e.g. 1, 2, 3 or 4 times 15N)
      labelTst <- tolower(gsub("^Label:[0-9]+|\\(|\\)", "", m$"Full name"))
      if (length(unique(labelTst)) == length(labelTst)) {
        Modifs$Mark[w] <- labelTst
      } else {
        m$AA[which(lengths(m$AA) == 0L)] <- "X"
        r <- { if ("Acetyl" %in% m$"Full name") { which(m$"Full name" == "Acetyl") } else { 1L } } # r is the one we will keep without changing it
        s <- setdiff(1L:nrow(m), r)
        tst <- lapply(s, \(x) {
          paste0(tolower(m$AA[[x]]), substr(m$Mark[[x]], 1L, 1L))
        })
        tst <- lapply(1L:length(tst), \(x) {
          rs <- tst[[x]]
          rs[which(!rs %in% Modifs$Mark)]
        })
        l <- length(tst)
        if (l > 1L) {
          for (i in 2L:l) {
            tst[[i]] <- tst[[i]][which(!tst[[i]] %in% unlist(tst[1L:(i - 1L)]))]
          }
        }
        tst <- vapply(tst, \(x) {
          x <- unlist(x)
          x <- { if (length(x)) { x[1L] } else { "That didnae work, did it?" } }
          return(x)
        }, "")
        tst2 <- ((tst %in% Modifs$Mark)|(tst == "That didnae work, did it?"))
        w2 <- which(!tst2)
        m$Mark[s][w2] <- tst[w2]
        w1 <- which(tst2)
        if (length(w1)) {
          # not tested
          s1 <- s[w1]
          rs <- c()
          kount <- 1L
          char <- c(as.character(0L:9L), letters)
          taken <- unique(c(Modifs$Mark, m$Mark))
          for (j in s1) {
            tst <- paste0(tolower(m$AA[s1]), char[kount])
            while (((tst) %in% taken)&&(kount < length(char))) {
              kount <- kount + 1L
              tst <- paste0(tolower(m$AA[s1]), char[kount])
            }
            if (kount == length(char)) {
              stop("I am really out of options here, never thought this would go this far! Check the code just in case, this should never happen.")
            } else { rs <- c(rs, tst) }
          }
          m$Mark[[s1]] <- rs
        }
        Modifs[w,] <- m
      }
    }
  }
  #
  Modifs$"Mass shift" <- UniMod$MonoMass[match(Modifs$UniMod, UniMod$UnimodId)]
  if ("Mass delta" %in% colnames(Modifs)) {
    w <- which(is.na(Modifs$"Mass shift"))
    Modifs$"Mass shift"[w] <- Modifs$"Mass delta"[w]
    # "Mass delta" is the one which was entered in the search - it may be slightly inaccurate
    # "Mass shift" is the monoisotopic which we pulled from UniMod and should be better
    w <- which(!is.na(Modifs$"Mass shift"))
    tst <- Modifs$"Mass shift"[w] - Modifs$"Mass delta"[w]
    if (max(tst) > 0.001) {
      warning("Discrepancies detected between the masses entered in the FragPipe search and the mapped UniMod PTMs.")
    } # Could probably be more precise
  }
  test <- aggregate(Modifs$"Full name", list(Modifs$Mark), \(x) { length(unique(x)) })
  if (max(test$x) > 1L) {
    warning("The algorithm was not able to assign unique 2-character old-MaxQuant style short marks to each modification, so some modified sequences may be ambiguous.")
  }
  # Remove parentheses in names: parentheses are also used to wrap around names in the full sequence,
  # nested parentheses are ugly and break the code further down.
  g <- grep("\\([0-9]+\\)$", Modifs$"Full name")
  if (length(g)) {
    Modifs$"Full name"[g] <- vapply(g, \(i) {
      paste0(gsub("\\([0-9]+\\)$", "", Modifs$"Full name"[i]),
             "*",
             gsub(".*\\(|\\)$", "", Modifs$"Full name"[i]))
    }, "")
  }
  g <- grep("\\(.+\\)$", Modifs$"Full name")
  if (length(g)) {
    Modifs$"Full name"[g] <- gsub("\\(", "{", gsub("\\)", "}", Modifs$"Full name"[g]))
  }
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
                              data.table = FALSE, fill = TRUE)
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
  kols <- c("Charge", "Number of Enzymatic Termini", "Number of Missed Cleavages", "Protein Start", "Protein End")
  kols <- kols[which(kols %in% colnames(PSMs))]
  for (kol in kols) { PSMs[[kol]] <- as.integer(PSMs[[kol]]) }
  kols <- c("Retention", "Observed Mass", "Calibrated Observed Mass", 
            "Observed M/Z", "Calibrated Observed M/Z", "Calculated Peptide Mass", 
            "Calculated M/Z", "Delta Mass", "Expectation", "Hyperscore", 
            "Nextscore", "PeptideProphet Probability", "Intensity", 
            "Best Score with Delta Mass", "Best Score without Delta Mass")
  kols <- kols[which(kols %in% colnames(PSMs))]
  if ((FailIfNoQuant)&&((!"Intensity" %in% kols)||(!length(is.all.good(PSMs$Intensity)))||(!sum(PSMs$Intensity, na.rm = TRUE)))) {
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
  }
  if (OpenSearch) {
    OpenSearch <- length(which(!is.na(PSMs$`Observed Modifications`))) > 0L
  }
  if (OpenSearch) {
    Modifs$"Type of search" <- "Closed"
    w <- which(PSMs$`Observed Modifications` != "")
    tmp <- gsub("; Mod", "; _;_Mod", gsub(", Mod", ", _;_Mod", PSMs$`Observed Modifications`[w]))
    tmp <- strsplit(tmp, "; _;_")
    tmp2 <- unique(unlist(tmp))
    tmp3 <- unique(unlist(strsplit(tmp2, ", _;_")))
    tmp3 <- gsub("^Mod[0-9]+: ", "", grep("^Mod[0-9]+: ", tmp3, value = TRUE))
    tmp4 <- data.frame(Full = tmp3,
                       `Mass delta` = NA,
                       PeakApex = NA,
                       check.names = FALSE)
    #tmp4 <- tmp4[grep("^Unannotated mass-shift |^Unidentified modification of ", tmp4$Full, invert = TRUE),]
    #tmp4 <- tmp4[grep("^((First)|(Second)|(Third)) isotopic peak|^Isotopic peak error$", tmp4$Full, invert = TRUE),]
    g <- grep("Theoretical: ", tmp4$Full)
    tmp4$"Mass delta"[g] <- as.numeric(gsub("(, .+)?\\)$", "", gsub(".+Theoretical: ", "", tmp4$Full[g])))
    g <- grep("PeakApex: ", tmp4$Full)
    tmp4$PeakApex[g] <- as.numeric(gsub("(, .+)?\\)$", "", gsub(".+PeakApex: ", "", tmp4$Full[g])))
    tmp4$"Full name" <- gsub(" \\(((Theoretical)|(PeakApex)): .+", "", tmp4$Full)
    w <- which(is.na(tmp4$"Mass delta"))
    tmp4$"Mass delta"[w] <- tmp4$PeakApex[w] 
    w <- which((is.na(tmp4$"Mass delta"))&(tmp4$"Full name" %in% UniMod$Description))
    tmp4$"Mass delta"[w] <- UniMod$MonoMass[match(tmp4$"Full name"[w], UniMod$Description)]
    tmp4$UniMod <- UniMod$UnimodId[match(tmp4$"Full name", UniMod$Description)]
    tmp4$Position <- "Unknown"
    tmp4$PeakApex <- NULL
    tmp4$Full <- NULL
    #tmp4$Enabled <- NA
    #mp4$"Max occurences" <- NA
    tmp4$AA <- lapply(1L:nrow(tmp4), \(x) { c() })
    tmp4$Type <- "Delta mass"
    tmp4$"Type of search" <- "Open"
    tmp4$Mark <- "!m" # We will mark all open-search delta masses as "!m" and add a column for delta masses 
    tmp4$"Full name" <- gsub(" ", "_", gsub(",", ".", gsub("\\(", "{", gsub("\\)", "}", tmp4$"Full name"))))
    Modifs <- plyr::rbind.fill(Modifs, tmp4)
  }
  # Filter peptides with ambiguous/non-classic amino acids 
  aaPat <- paste(AA, collapse = "|")
  parallel::clusterExport(cl, list("aaPat", "Modifs", "OpenSearch"), envir = environment())
  test <- gsub(aaPat, "", PSMs$Peptide)
  w <- which(test == "")
  if (length(w) < nrow(PSMs)) {
    warning("Removing peptides with unknown amino acids!s")
    PSMs <- PSMs[w,]
  }
  # Create MQ-like data.frame
  kol <- c("Spectrum File", "Peptide", "Modified Peptide", "Charge", "Hyperscore",
           "Spectrum", "Intensity", "Retention", "Protein")
  #print(kol[which(!kol %in% colnames(PSMs))])
  stopifnot(sum(!kol %in% colnames(PSMs)) == 0L)
  EV <- data.frame(Sequence = PSMs$Peptide,
                   Length = nchar(PSMs$Peptide), 
                   Charge = PSMs$Charge,
                   Intensity = PSMs$Intensity)
  #
  pKol <- c("PeptideProphet Probability", "Probability")
  pKol <- pKol[which(pKol %in% colnames(PSMs))[1L]]
  if (!length(pKol)) { stop("No PSMs quality column detected!") } else {
    EV$PEP <- 1 - PSMs[[pKol]]
  }
  a1 <- PSMs$Protein
  a2 <- PSMs$"Mapped Proteins"
  if (UniProtIDs) {
    a1 <- strsplit(gsub("^[a-z]{2}\\||\\|[^\\|]+$", "",
                        gsub("\\|[^\\|]+;[a-z]{2}\\|", ";", gsub(", ", ";", a1))), ";")
    a2 <- strsplit(gsub("^[a-z]{2}\\||\\|[^\\|]+$", "",
                        gsub("\\|[^\\|]+;[a-z]{2}\\|", ";", gsub(", ", ";", a2))), ";")
  }
  tmp <- cbind(a1, a2)
  f0 <- \(x) { paste(sort(unique(unlist(x))), collapse = ";") }
  f0 <- .bind_worker(f0,
                     list(tmp = tmp))
  tmp <- parallel::parApply(cl,
                            tmp,
                            1L,
                            f0)
  EV$Proteins <- tmp
  if (isWellBehaved) {
    EV$Experiment <- PSMs$Sample
    #EV$"Raw file path" <- FP_Mnfst$Path[match(EV$Experiment, FP_Mnfst$Samples)] # Buggy line
    EV$"Raw file" <- gsub_Rep("^(.*/)?interact-", "",
                              gsub_Rep("\\\\", "/",
                                       gsub_Rep("(\\.mod)?\\.pep\\.xml$", ".d", PSMs$`Spectrum File`)))
    EV$"Raw file" <- gsub_Rep(".+/|\\.((raw)|(mz(X?ML|BIN))|(mgf)|(d))$", "", EV$"Raw file")
    u <- unique(EV$"Raw file")
    w <- lapply(FP_Mnfst$`File name`, \(x) { which(u == x) })
    l <- lengths(w)
    stopifnot(max(l) == 1L) # Would indicate that we have non unique file names, which this cannot deal with!
    EV$"Raw file path" <- FP_Mnfst$Path[match(EV$"Raw file", FP_Mnfst$`File name`)]
  } else {
    EV$Experiment <- FP_Mnfst$Experiment[PSMs2Mnsfst]
    EV$"Raw file path" <- FP_Mnfst$Path[PSMs2Mnsfst]
    EV$"Raw file" <- gsub_Rep(".+/|\\.((raw)|(mz(X?ML|BIN))|(mgf)|(d))$", "", EV$"Raw file path")
  }
  # Modified sequence
  EV$"Modified sequence_verbose" <- EV$"Modified sequence" <- paste0("_", EV$Sequence, "_")
  wMdSq <- wMdSq2 <- which(PSMs$`Assigned Modifications` != "")
  if (OpenSearch) {
    wMdSq2 <- which((PSMs$`Assigned Modifications` != "")|(PSMs$`Observed Modifications` != ""))
  }
  if (length(wMdSq2)) {
    a1 <- strsplit(PSMs$Peptide[wMdSq2], "")
    a2 <- as.list(rep(NA, length(wMdSq2)))
    tmp <- strsplit(PSMs$`Assigned Modifications`[wMdSq], ", ?")
    f0 <- .bind_worker(.FP2MQ_modSeqWrkr1,
                       list(tmp = tmp,
                            Modifs = Modifs,
                            aaPat = aaPat))
    a2[which(wMdSq2 %in% wMdSq)] <- parallel::parLapply(cl, 
                                                        tmp,
                                                        f0,
                                                        mods = Modifs,
                                                        pat = aaPat)
    tmp <- cbind(a1, a2)
    if (OpenSearch) {
      a3 <- gsub("^Mod[0-9]+: ", "", gsub(", Mod[0-9]+: ", ";", gsub("; .+", "", PSMs$"Observed Modifications"[wMdSq2])))
      a3 <- gsub(" \\([^\\)]+\\)$", "", gsub(" \\([^\\)]+\\), ", ", ", a3)) 
      a3 <- gsub(" ", "_", gsub(",", ".", gsub("\\(", "{", gsub("\\)", "}", a3))))
      a3 <- strsplit(a3, ";")
      stopifnot(sum(!unique(unlist(a3)) %in% Modifs$"Full name") == 0L)
      tmp <- cbind(tmp, a3)
    }
    f0 <- .bind_worker(.FP2MQ_modSeqWrkr2,
                       list(tmp = tmp,
                            OpenSearch = OpenSearch))
    EV$"Modified sequence_verbose"[wMdSq2] <- parallel::parApply(cl,
                                                                 tmp,
                                                                 1L,
                                                                 f0,
                                                                 openSearch = OpenSearch)
    # Note: it seems that currently if several mods (different ones, or multiples of the same) occur twice
    # on the same location, MaxQuant writes it as:
    #   PEPTID(ModName1)(ModName2)R
    # We are not doing it like that right now.
    # Since we are not looking for perfect compatibility I guess this is ok?
    #
    #tst <- unique(unlist(strsplit(gsub("(_|\\))[A-Z]+(_|\\()", "_", EV$"Modified sequence_verbose"), "\\)_$|_|,")))
    #tst <- grep("^[A-Z]+$", tst, value = TRUE, invert = TRUE)
    #tst <- gsub("^[0-9]+ ", "", tst[which(tst != "")])
    #sum(!tst %in% Modifs$"Full name")
    tmp <- strsplit(EV$"Modified sequence_verbose"[wMdSq2], "\\(|\\)")
    f0 <- .bind_worker(.FP2MQ_modSeqWrkr3,
                       list(tmp = tmp,
                            Modifs = Modifs))
    EV$"Modified sequence"[wMdSq2] <- parallel::parSapply(cl,
                                                          tmp,
                                                          f0,
                                                          mods = Modifs)
  }
  EV$"Modified sequence" <- unlist(EV$"Modified sequence") # Because in some cases this has become a list and I don't know yet why
  # For now a corrective rather than a fix until I can identify the issue. Should trigger an error if the vector has the wrong length.
  EV$Modifications <- "Unmodified"
  a1 <- strsplit(gsub(paste0("_|\\)|", aaPat), "", EV$"Modified sequence"[wMdSq2]), "\\(")
  f0 <- \(x) { x[which(x != "")] }
  f0 <- .bind_worker(f0,
                     list(a1 = a1))
  tmp <- a1 <- parallel::parLapply(cl,
                                   a1,
                                   f0)
  if (OpenSearch) {
    tmp <- cbind(a1, a3)
    f0 <- .bind_worker(.FP2MQ_modSeqWrkr4,
                       list(tmp = tmp,
                            Modifs = Modifs))
    EV$Modifications[wMdSq2] <- parallel::parApply(cl,
                                                   tmp,
                                                   1L,
                                                   f0,
                                                   mods = Modifs)
  } else {
    f0 <- .bind_worker(.FP2MQ_modSeqWrkr5,
                       list(tmp = tmp,
                            Modifs = Modifs))
    EV$Modifications[wMdSq2] <- parallel::parSapply(cl,
                                                    tmp,
                                                    f0,
                                                    mods = Modifs)
  }
  #tst <- unique(gsub("^[0-9]+ ", "", unlist(strsplit(EV$Modifications[wMdSq2], ","))))
  #sum(!tst %in% Modifs$"Full name")
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
  #length(which(tst$V2 != tst$V3)) == 0L
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
  EV$"MS/MS scan number" <- as.numeric(gsub("^[^\\.]+\\.|\\..+", "", PSMs$Spectrum))
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
  EV <- EV[which(EV$`Delta score` > Min.Delta.Score),]
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
                FP_PTMs = Modifs)
  } else {
    Res <- list(Evidence = EV,
                PTMs = Modifs,
                FracMap = FP_Mnfst,
                WorkFlow = FP_Wrkflw)
  }
  if (isTMT) { Res[["TMT_annotations"]] <- TMTtbl }
  #
  if (stopCl) { parallel::stopCluster(cl) }
  return(Res)
}
