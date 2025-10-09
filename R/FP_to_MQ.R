#' FP_to_MQ
#'
#' @description 
#' Reads FragPipe's workflow and manifest files
#' From them, detects and parses psm.tsv tables (one per Experiment) to a single MaxQuant evidence.txt-like table.
#' As with the other functions of this type, the idea is not to get perfect conversion, but close enough that the data can be fed into this package's analysis scripts.
#' 
#' The output is a list of 4 items:
#' - Evidence: an evidence table similar to the one MaxQuant creates
#' - PTMs: a modifications table
#' - FracMap: a map of the different MS spectrum files and their relationship to Experiments
#' - WorkFlow: the FragPipe workflow
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
                     N.reserved = 1,
                     cl,
                     OpenSearch = FALSE) {
  TESTING <- FALSE
  #proteoCraft::DefArg(proteoCraft::FP_to_MQ); cl = parClust;TESTING <- TRUE
  #FP_Workflow <- fpWorkflowFl_i; FP_Manifest <- fpManifestFl_i
  #FP_Workflow <- FP_WorkflowFl; FP_Manifest <- FP_ManifestFl
  #FP_Workflow <- paste0(wd, "/fragpipe.workflow"); FP_Manifest <- paste0(wd, "/fragpipe-files.fp-manifest")
  #
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
  # Load FragPipe parameter file
  if ((!TESTING)&&((missing("FP_Workflow"))||(!file.exists(FP_Workflow)))) {
    stop("Invalid parameter \"FP_Workflow\": file not found!")
  }
  FP_Wrkflw <- readLines(FP_Workflow)
  pat <- proteoCraft::topattern("diann.run-dia-nn=")
  isActuallyDIANN <- as.logical(toupper(gsub(pat, "", grep(pat, FP_Wrkflw, value = TRUE))))
  #
  FP_Dir <- gsub("\\\\", "", gsub("\\\\\\\\", "/", gsub("^workdir=", "",
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
  pat <- proteoCraft::topattern("tmtintegrator.run-tmtintegrator=")
  isTMT <- as.logical(gsub(pat, "", grep(pat, FP_Wrkflw, value = TRUE)))
  if (isTMT) {
    pat <- proteoCraft::topattern("tmtintegrator.channel_num=TMT-")
    TMTplex <- as.integer(gsub(pat, "", grep(pat, FP_Wrkflw, value = TRUE)))
    TMTtblFl <- paste0(FP_Dir, "/experiment_annotation.tsv")
    if (!file.exists(TMTtblFl)) {
      stop(paste0("This is a TMT dataset but file \"", TMTtblFl, "\" could not found in the results folder!"))
    }
    TMTtbl <- read.delim(TMTtblFl, check.names = FALSE)
    TMTtbl$plex[which(is.na(TMTtbl$plex))] <- TMTplex
    Channels <- c(126, unlist(lapply(127:134, function(x) { paste0(as.character(x), c("N", "C")) })), "135N")
    TMTtbl$channel_code <- match(TMTtbl$channel, Channels)-1
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
    Samples1 <- apply(Samples, 1, function(x) {
      x <- c(x[[1]], x[[2]])
      x <- x[which(!is.na(x))]
      return(paste(x, collapse = "_"))
    })
    Samples2 <- apply(Samples, 1, function(x) {
      x <- c(x[[1]], x[[2]], x[[2]])
      x <- x[which(!is.na(x))]
      return(paste(x, collapse = "_"))
    })
    tst1 <- sum(dir.exists(paste0(FP_Dir, "/", Samples1))) == length(Samples1)
    tst2 <- sum(dir.exists(paste0(FP_Dir, "/", Samples2))) == length(Samples2)
    stopifnot(tst1+tst2 == 1)
    Samples <- get(paste0("Samples", 1:2)[which(c(tst1, tst2))])
  }
  FP_Mnfst$Samples <- Samples
  #
  FP_PSMs <- unique(gsub("//", "/", paste0(FP_Dir, "/", Samples, "/psm.tsv")))
  #FP_Fas <- paste0(FP_Dir, "/", Samples, "/protein.fas")
  tst <- file.exists(FP_PSMs)
  w <- which(!tst); l <- length(w)
  FP_Mnfst2 <- FP_Mnfst
  if (l) {
    #SamplesNo <- unique(SamplesNo[w]); l <- length(SamplesNo)
    #if (l > 1) { SamplesNo <- paste0(paste(SamplesNo[1:(l-1)], collapse = ", "), " and ", SamplesNo[l]) }
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
    stopifnot(length(Samples) == 0)
    stopifnot(length(Exp) == 0)
    stopifnot(dim(FP_Mnfst2)[1] == 0)
    FP_PSMs <- paste0(FP_Dir, "/psm.tsv")
    stopifnot(file.exists(FP_PSMs))
    isWellBehaved <- FALSE
  }
  #tst <- file.exists(FP_Fas)
  #w <- which(!tst); l <- length(w)
  #if (l) {
  #  SamplesNo <- unique(SamplesNo[w]); l <- length(SamplesNo)
  #  if (l > 1) { SamplesNo <- paste0(paste(SamplesNo[1:(l-1)], collapse = ", "), " and ", SamplesNo[l]) }
  #  plur <- c("", "s")[(l > 1)+1]
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
  UniMod <- unimod::modifications
  # Unfortunately, as of early 2023 this package did not have the whole of Unimod,
  # in particular TMT16plex was missing.
  # This should be fixed now but just in case.
  sites <- c("K", "N-term", "N-term", "H", "S", "T")
  sites2 <- c("K", "N-term", "P-N-term", "H", "S", "T")
  where <- c("Anywhere", "Any N-term", "Protein N-term", "Anywhere", "Anywhere", "Anywhere")
  tmtp <- data.frame(Id = paste0("TMTpro:", sites2),
                     UnimodId = 2016,
                     Name = "TMTpro",
                     Description = "TMTpro 16plex Tandem Mass Tag®",
                     Composition = "H(25) C(8) 13C(7) N 15N(2) O(3)",
                     AvgMass = 304.3127,
                     MonoMass = 304.207146,
                     Site = sites2,
                     Position = where,
                     Classification = "Isotopic label",
                     SpecGroup = 1,
                     NeutralLoss = FALSE,
                     LastModified = "2021-07-22 21:42:31",
                     Approved = FALSE,
                     Hidden = FALSE)
  UniMod <- rbind(UniMod, tmtp)
  # Remove:
  # - substitutions
  UniMod <- UniMod[grep("^[A-Z][a-z]{2}->[A-Z][a-z]{2} substitution$", UniMod$Description, invert = TRUE),]
  # - deltas
  g <- grepl("^Delta:", UniMod$Name)
  UniMod_D <- UniMod[which(g),]
  UniMod <- UniMod[which(!g),]
  #
  parseMods <- function(mods, fixed = FALSE) { #mods <- varMods
    x1 <- sapply(mods, function(x) { as.numeric(x[[1]]) })
    x2 <- sapply(mods, function(x) { paste(x[2:(length(x)-2)], collapse = "") })
    x3 <- sapply(mods, function(x) { as.logical(toupper(x[[length(x)-1]])) })
    x4 <- sapply(mods, function(x) { as.integer(x[[length(x)]]) })
    mods <- data.frame("Mass delta" = x1,
                       "Site" = x2,
                       "Enabled" = x3,
                       "Max occurences" = x4,
                       "Type" = c("Variable", "Fixed")[fixed + 1],
                       check.names = FALSE)
    mods <- mods[which(mods$"Mass delta" != 0),]
    mods <- mods[which(mods$Enabled),]
    mods$Enabled <- NULL
    mods$Site <- gsub(" ?\\([^\\)]+\\)", "", mods$Site) # Why do I need this again?
    mods$ID <- paste0(c("var", "fix")[fixed+1], 1:nrow(mods))
    tmp <- lapply(strsplit(mods$Site, ""), function(x) { #x <- strsplit(mods$Site, "")[4]
      x <- unlist(x)
      w <- which(x == "n")
      if (length(w)) {
        x[w+1] <- paste0("n", x[w+1])
        x <- x[-w]
      }
      w <- which(x == "c")
      if (length(w)) {
        x[w+1] <- paste0("c", x[w+1])
        x <- x[-w]
      }
      w <- which(x == "[")
      if (length(w)) {
        x[w+1] <- paste0("[", x[w+1])
        x <- x[-w]
      }
      w <- which(x == "]")
      if (length(w)) {
        x[w+1] <- paste0("]", x[w+1])
        x <- x[-w]
      }
      return(x)
    })
    tmp <- proteoCraft::listMelt(tmp, 1:nrow(mods))
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
  tst <- sapply(strsplit(Modifs$Site, ""), function(x) {
    x <- x[which(x != " ")]
    sum(!x %in% c("[", "n", "^", "c", "]", proteoCraft::AA)) == 0
  })
  if (sum(!tst)) { warning("Unexpected result when parsing FragPipe's modification sites, investigate!") }
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
  tmp <- lapply(1:nrow(Modifs), function(i) { #i <- 1 #i <- 2  #i <- 3 #i <- 4
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
    w <- which((abs(round(UniMod$MonoMass, prec) - round(as.numeric(dMass), prec)) <= 10^(1-prec))
               &(UniMod$Site %in% sites)
               &(UniMod$Position == pos))
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
        res$Site[w] <- gsub("^ ", "", paste0(apply(res[w, c("Site", "Position")], 1, paste, collapse = " ("), ")"))
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
  tmp2 <- proteoCraft::listMelt(tmp$Site_long, 1:nrow(tmp), c("Site_long", "tmpRow"))
  tmp2[, c("UniMod", "Row")] <- tmp[tmp2$tmpRow, c("UniMod", "Row")]
  kol <- c("ID", "Mass delta", "Site", "Max occurences", "Type", "Position", "AA")
  tmp2[, kol] <- Modifs[tmp2$Row, kol]
  tmp2$"Full name" <- ""
  w <- which(!is.na(tmp2$UniMod))
  tmp2$"Full name"[w] <- UniMod$Name[match(tmp2$UniMod[w], UniMod$UnimodId)]
  w <- which(is.na(tmp2$UniMod))
  if (length(w)) {
    tmp2$"Full name"[w] <- paste0("Unknown", 1:length(w), "_", c("", "+")[(tmp2$`Mass delta`[w] > 0)+1],
                                  tmp2$`Mass delta`[w])
    
  }
  kol <- colnames(tmp2)[which(!colnames(tmp2) %in% c("ID", "tmpRow", "Row", "Full name"))]
  tmp2 <- aggregate(tmp2[, kol], list(tmp2$ID, tmp2$"Full name"), list)
  colnames(tmp2)[2] <- "Full name"
  tmp2$Group.1 <- NULL
  tmp2$"Mass delta" <- sapply(tmp2$"Mass delta", unique)
  stopifnot(!"list" %in% class(tmp2$"Mass delta"))
  tmp2$"Max occurences" <- sapply(tmp2$"Max occurences", unique)
  stopifnot(!"list" %in% class(tmp2$"Max occurences"))
  tmp2$Type <- sapply(tmp2$Type, unique)
  stopifnot(!"list" %in% class(tmp2$Type))
  tmp2$UniMod <- sapply(tmp2$UniMod, unique)
  stopifnot(!"list" %in% class(tmp2$UniMod))
  Modifs <- tmp2 
  # - Marks
  Modifs$Mark <- tolower(substr(Modifs$"Full name", 1, 2))
  ## Sometimes, some marks are duplicates, e.g. if you searched for "Acetyl (Protein N-term)" and "Acetyl (K)" together!
  ## We want to fix this so that each modification has a unique mark:
  #test <- aggregate(Modifs$Mark, list(Modifs$Mark), length)
  test <- aggregate(Modifs$UniMod, list(Modifs$Mark), function(x) { length(unique(x)) })
  W <- which(test$x > 1)
  #Modifs$Mark <- Modifs$"Old mark"
  if (length(W)) {
    Modifs$"Old mark" <- Modifs$Mark
    for (i in W) {
      #i <- W[1]
      w <- which(Modifs$Mark == test$Group.1[i])
      m <- Modifs[w,]
      m$AA[which(sapply(m$AA, length) == 0)] <- "X"
      if ("Acetyl" %in% m$"Full name") { r <- which(m$"Full name" == "Acetyl") } else { r <- 1 } # r is the one we will keep without changing it
      s <- 1:nrow(m); s <- s[which(s != r)]
      tst <- lapply(s, function(x) {
        paste0(tolower(m$AA[[x]]), substr(m$Mark[[x]], 1, 1))
      })
      tst <- lapply(1:length(tst), function(x) {
        rs <- tst[[x]]
        rs[which(!rs %in% Modifs$Mark)]
      })
      l <- length(tst)
      if (l > 1) {
        for (i in 2:l) {
          tst[[i]] <- tst[[i]][which(!tst[[i]] %in% unlist(tst[1:(i-1)]))]
        }
      }
      tst <- sapply(tst, function(x) {
        x <- unlist(x)
        if (length(x)) { x <- x[1] } else { x <- "That didnae work, did it?" }
        return(x)
      })
      tst2 <- ((tst %in% Modifs$Mark)|(tst == "That didnae work, did it?"))
      w2 <- which(!tst2)
      m$Mark[s][w2] <- tst[w2]
      w1 <- which(tst2)
      if (length(w1)) {
        # not tested
        s1 <- s[w1]
        rs <- c()
        kount <- 1
        char <- c(0:9, letters)
        taken <- unique(c(Modifs$Mark, m$Mark))
        for (j in s1) {
          tst <- paste0(tolower(m$AA[s1]), char[kount])
          while (((tst) %in% taken)&&(kount < length(char))) {
            kount <- kount+1
            tst <- paste0(tolower(m$AA[s1]), char[kount])
          }
          if (kount == length(char)) {
            stop("I am really out of options here, never thought this would go this far! Check the code just in case, this should never happen.")
          } else {
            rs <- c(rs, tst)
          }
        }
        m$Mark[[s1]] <- rs
      }
      Modifs[w,] <- m
    }
  }
  #
  Modifs$"Mass shift" <- UniMod$MonoMass[match(Modifs$UniMod, UniMod$UnimodId)]
  if ("Mass delta" %in% colnames(Modifs)) {
    # "Mass delta" is the one which was entered in the search - it may be slightly inaccurate
    # "Mass shift" is the monoisotopic which we pulled from UniMod and should be better
    w <- which(!is.na(Modifs$"Mass shift"))
    tst <- Modifs$"Mass shift"[w] - Modifs$"Mass delta"[w]
    if (max(tst) > 0.001) {
      warning("Discrepancies detected between the masses entered in the FragPipe search and the mapped UniMod PTMs.")
    } # Could probably be more precise
  }
  test <- aggregate(Modifs$"Full name", list(Modifs$Mark), function(x) { length(unique(x)) })
  if (max(test$x) > 1) {
    warning("The algorithm was not able to assign unique re-character old-MaxQuant style short marks to each modification, so some modified sequences may be ambiguous.")
  }
  #
  # Load PSM files
  if (isWellBehaved) {
    PSMs <- lapply(1:length(Samples), function(x) { #x <- 1
      res <- data.table::fread(FP_PSMs[x], integer64 = "numeric", check.names = FALSE,
                               data.table = FALSE, fill = TRUE, sep = "\t")
      if (nrow(res)) {
        res$Sample <- Samples[x]
        res$PSMs_file <- FP_PSMs[x]
        return(res)
      } else { return(NA) }
    })  
    w <- which(sapply(PSMs, function(x) { "data.frame" %in% class(x) }))
    PSMs <- plyr::rbind.fill(PSMs[w])
  } else {
    PSMs <- data.table::fread(FP_PSMs, integer64 = "numeric", check.names = FALSE,
                              data.table = FALSE, fill = TRUE)
    PSMs2Mnsfst <- match(proteoCraft::gsub_Rep("\\\\", "/", PSMs$`Spectrum File`),
               paste0(FP_Dir, "/interact-", FP_Mnfst$`File name`, ".pep.xml"))
    w <- which(is.na(PSMs2Mnsfst))
    l <- length(w)
    if (l) {
      if (l == nrow(PSMs)) { PSMs2Mnsfst <- rep(1, nrow(PSMs)) } else {
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
  kols <- c("Retention", "Observed Mass", "Calibrated Observed Mass", "Observed M/Z", "Calibrated Observed M/Z",
            "Calculated Peptide Mass", "Calculated M/Z", "Delta Mass", "Expectation", "Hyperscore", "Nextscore",
            "PeptideProphet Probability", "Intensity", "Best Score with Delta Mass", "Best Score without Delta Mass")
  kols <- kols[which(kols %in% colnames(PSMs))]
  if ((FailIfNoQuant)&&((!"Intensity" %in% kols)||(!length(proteoCraft::is.all.good(PSMs$Intensity))))) {
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
    OpenSearch <- length(which(!is.na(PSMs$`Observed Modifications`))) > 0  
  }
  if (OpenSearch) {
    Modifs$"Type of search" <- "Closed"
    w <- which(PSMs$`Observed Modifications` != "")
    tmp <- gsub("; Mod", "; _;_Mod", gsub(", Mod", ", _;_Mod", PSMs$`Observed Modifications`[w]))
    tmp <- strsplit(tmp, "; _;_")
    tmp2 <- unique(unlist(tmp))
    tmp3 <- unique(unlist(strsplit(tmp2, ", _;_")))
    tmp3 <- gsub("^Mod[0-9]+: ", "", grep("^Mod[0-9]+: ", tmp3, value = TRUE))
    tmp4 <- data.frame(Full = tmp3, "Mass delta" = NA, PeakApex = NA, check.names = FALSE)
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
    tmp4$AA <- lapply(1:nrow(tmp4), function(x) { c() })
    tmp4$Type <- "Delta mass"
    tmp4$"Type of search" <- "Open"
    tmp4$Mark <- "!m" # We will mark all open-search delta masses as "!m" and add a column for delta masses 
    tmp4$`Full name` <- gsub(" ", "_", gsub(",", ".", gsub("\\(", "{", gsub("\\)", "}", tmp4$`Full name`))))
    Modifs <- plyr::rbind.fill(Modifs, tmp4)
  }
  # Filter peptides with ambiguous/non-classic amino acids 
  aaPat <- paste(proteoCraft::AA, collapse = "|")
  parallel::clusterExport(cl, list("aaPat", "Modifs", "OpenSearch"), envir = environment()) # Useful for later
  test <- gsub(aaPat, "", PSMs$Peptide)
  w <- which(test == "")
  if (length(w) < nrow(PSMs)) {
    warning("Removing peptides with unknown amino acids!s")
    PSMs <- PSMs[w,]
  }
  # Create MQ-like file
  kol <- c("Spectrum File", "Peptide", "Modified Peptide", "Charge", "Hyperscore",
           "Spectrum", "Intensity", "Retention", "Protein")
  #print(kol[which(!kol %in% colnames(PSMs))])
  stopifnot(sum(!kol %in% colnames(PSMs)) == 0)
  EV <- data.frame("Sequence" = PSMs$Peptide,
                   "Length" = nchar(PSMs$Peptide),
                   "Charge" = PSMs$Charge,
                   "Intensity" = PSMs$Intensity,
                   check.names = FALSE)
  #
  pKol <- c("PeptideProphet Probability", "Probability")
  pKol <- pKol[which(pKol %in% colnames(PSMs))[1]]
  if (!length(pKol)) { stop("No PSMs quality column detected!") } else { EV$PEP <- 1-PSMs[[pKol]] }
  a1 <- PSMs$Protein
  a2 <- PSMs$"Mapped Proteins"
  if (UniProtIDs) {
    a1 <- strsplit(gsub("^[a-z]{2}\\||\\|[^\\|]+$", "", gsub("\\|[^\\|]+;[a-z]{2}\\|", ";", gsub(", ", ";", a1))), ";")
    a2 <- strsplit(gsub("^[a-z]{2}\\||\\|[^\\|]+$", "", gsub("\\|[^\\|]+;[a-z]{2}\\|", ";", gsub(", ", ";", a2))), ";")
  }
  tmp <- cbind(a1, a2)
  f0 <- function(x) { paste(sort(unique(unlist(x))), collapse = ";") }
  environment(f0) <- .GlobalEnv
  tmp <- parallel::parApply(cl, tmp, 1, f0) # Good idea to check how reliable this is using ProtMatch
  EV$Proteins <- tmp
  if (isWellBehaved) {
    EV$Experiment <- PSMs$Sample
    #EV$"Raw file path" <- FP_Mnfst$Path[match(EV$Experiment, FP_Mnfst$Samples)] # Buggy line
    EV$"Raw file" <- proteoCraft::gsub_Rep("^(.*/)?interact-", "",
                                           proteoCraft::gsub_Rep("\\\\", "/",
                                                                 proteoCraft::gsub_Rep("(\\.mod)?\\.pep\\.xml$", ".d",
                                                                                       PSMs$`Spectrum File`)))
    EV$"Raw file" <- proteoCraft::gsub_Rep(".+/|\\.((raw)|(mz(X?ML|BIN))|(mgf)|(d))$", "", EV$"Raw file")
    u <- unique(EV$"Raw file")
    w <- lapply(FP_Mnfst$`File name`, function(x) { which(u == x) })
    l <- sapply(w, length)
    stopifnot(max(l) == 1) # Would indicate that we have non unique file names, which this cannot deal with!
    EV$"Raw file path" <- FP_Mnfst$Path[match(EV$"Raw file", FP_Mnfst$`File name`)]
  } else {
    EV$Experiment <- FP_Mnfst$Experiment[PSMs2Mnsfst]
    EV$"Raw file path" <- FP_Mnfst$Path[PSMs2Mnsfst]
    EV$"Raw file" <- proteoCraft::gsub_Rep(".+/|\\.((raw)|(mz(X?ML|BIN))|(mgf)|(d))$", "", EV$"Raw file path")
  }
  # Modified sequence
  EV$"Modified sequence_verbose" <- EV$"Modified sequence" <- paste0("_", EV$Sequence, "_")
  wMdSq <- wMdSq2 <- which(PSMs$`Assigned Modifications` != "")
  if (OpenSearch) { wMdSq2 <- which((PSMs$`Assigned Modifications` != "")|(PSMs$`Observed Modifications` != "")) }
  if (length(wMdSq2)) {
    a1 <- strsplit(PSMs$Peptide[wMdSq2], "")
    a2 <- as.list(rep(NA, length(wMdSq2)))
    f0 <- function(x) {
      #for (i in wMdSq) {
      #x <- strsplit(PSMs$`Assigned Modifications`[i], ", ?")
      #x <- strsplit(PSMs$`Assigned Modifications`[wMdSq[1]], ", ?")
      x <- gsub("\\)$", "", unlist(x))
      ptms <- as.data.frame(t(sapply(x, function(y) { unlist(strsplit(y, split = "\\(")) })))
      colnames(ptms) <- c("Site", "Mass shift")
      rownames(ptms) <- NULL
      ptms[, c("Site", "Position")] <- t(sapply(ptms$Site, function(aa) {
        if (aa == "N-term") { ps <- 0 }
        if (aa == "C-term") { ps <- l+1 }
        if (!aa %in% c("N-term", "C-term")) {
          ps <- gsub(paste0(aaPat, "$"), "", aa)
          aa <- gsub("^[0-9]+", "", aa)
        }
        return(c(aa, ps))
      }))
      ptms$"Mass shift" <- as.numeric(ptms$"Mass shift")
      ptms$Position <- as.numeric(ptms$Position)
      ptms$"Full name" <- NA
      w <- which(is.na(ptms$"Full name"))
      k <- 5
      while ((length(w))&&(k > 2)) {
        ptms$"Full name"[w] <- Modifs$"Full name"[match(round(ptms$`Mass shift`[w], k), round(Modifs$`Mass delta`, k))]
        w <- which(is.na(ptms$"Full name"))
        k <- k-1
      }
      if (length(w)) {
        stop("I am an error and I mean that some assigned PTM was not recognized... but actually you should probably convert me to a warning!")
      }
      #Modifs[match(round(ptms$`Mass shift`, 3), round(Modifs$`Mass delta`, 3)),]
      #}
      return(ptms)
    }
    environment(f0) <- .GlobalEnv
    a2[which(wMdSq2 %in% wMdSq)] <- parallel::parLapply(cl, strsplit(PSMs$`Assigned Modifications`[wMdSq], ", ?"), f0)
    tmp <- cbind(a1, a2)
    if (OpenSearch) {
      a3 <- gsub("^Mod[0-9]+: ", "", gsub(", Mod[0-9]+: ", ";", gsub("; .+", "", PSMs$"Observed Modifications"[wMdSq2])))
      a3 <- gsub(" \\([^\\)]+\\)$", "", gsub(" \\([^\\)]+\\), ", ", ", a3)) 
      a3 <- gsub(" ", "_", gsub(",", ".", gsub("\\(", "{", gsub("\\)", "}", a3))))
      a3 <- strsplit(a3, ";")
      stopifnot(sum(!unique(unlist(a3)) %in% Modifs$`Full name`) == 0)
      tmp <- cbind(tmp, a3)
    }
    f0 <- function(x) {
      #sapply(1:length(tmp), function(i) {
      #i <- 1
      #x <- tmp[i,]
      seq <- unlist(x[[1]])
      l <- length(seq)
      rs <- data.frame(Seq = c("N-term", seq, "C-term"), Mod = "")
      ptms <- x[[2]]
      if ("data.frame" %in% class(ptms)) {
        stopifnot(sum(rs$Seq[ptms$Position+1] != ptms$Site) == 0)
        rs$Mod[ptms$Position+1] <- ptms$"Full name"
        tmp2 <- c(rs$Mod[1], rs$Mod[2])
        tmp2 <- tmp2[which(tmp2 != "")]
        if (length(tmp2)) {
          tmp2 <- aggregate(tmp2, list(tmp2), length)
          tmp2 <- paste(gsub("^1 ", "", apply(tmp2[, c("x", "Group.1")], 1, paste, collapse = " ")), collapse = ",")
        } else { tmp2 <- "" }
        rs$Mod[2] <- tmp2; rs$Mod[1] <- ""
        rs$Mod[l+1] <- paste0(rs$Mod[l+1], rs$Mod[l+2]); rs$Mod[l+2] <- ""
      }
      rs$Seq[1] <- rs$Seq[l+2] <- "_"
      # NB: although MaxQuant does not label fixed modifications...
      # we will keep those labels because it's useful information and does not hurt
      #
      # Open search:
      # Add a simple catch-all mark for delta mass which will just indicate that one should look into the "Mass error [Da]" column
      if (OpenSearch) {
        obsptms <- unlist(x[[3]])
        if (length(obsptms)) {
          obsptms <- obsptms[which(obsptms != "")]
          rs$Mod[l+1] <- gsub("^,", "", paste(c(rs$Mod[l+1], obsptms), collapse = ","))
        }
      }
      w <- which(rs$Mod != "")
      rs$Mod[w] <- paste0("(", rs$Mod[w], ")")
      rs <- paste(apply(rs, 1, paste, collapse = ""), collapse = "")
      return(rs)
    }
    environment(f0) <- .GlobalEnv
    EV$"Modified sequence_verbose"[wMdSq2] <- parallel::parApply(cl, tmp, 1, f0)
    #
    #tst <- unique(unlist(strsplit(gsub("(_|\\))[A-Z]+(_|\\()", "_", EV$"Modified sequence_verbose"), "\\)_$|_|,")))
    #tst <- grep("^[A-Z]+$", tst, value = TRUE, invert = TRUE)
    #tst <- gsub("^[0-9]+ ", "", tst[which(tst != "")])
    #sum(!tst %in% Modifs$Full)
  }
  f0 <- function(x) {
    x <- unlist(x)
    l <- length(x)
    wmds <- (1:((l - 1)/2)) * 2
    mds <- strsplit(x[wmds], ",")
    w1 <- grep("^[0-9]+ ", mds, invert = TRUE)
    w2 <- grep("^[0-9]+ ", mds)
    mds2 <- list()
    if (length(w1)) { mds2[w1] <- mds[w1] }
    if (length(w2)) {
      mds2[w2] <- lapply(mds[w2], function(y) {
        y <- unlist(strsplit(y, " "))
        rep(y[2], as.integer(y[1]))
      })
    }
    mds <- paste0("(", sapply(mds2, function(y) {
      y <- Modifs$Mark[match(y, Modifs$`Full name`)]
      w <- which(y == "!m")
      if (length(w) > 1) { y <- y[c(which(y != "!m"), w[1])] }
      return(paste(y, collapse = ","))
    }), ")")
    x[wmds] <- mds
    return(paste(x, collapse = ""))
  }
  environment(f0) <- .GlobalEnv
  EV$"Modified sequence"[wMdSq2] <- parallel::parSapply(cl,
                                                        strsplit(EV$"Modified sequence_verbose"[wMdSq2], "\\(|\\)"),
                                                        f0)
  EV$"Modified sequence" <- unlist(EV$"Modified sequence") # Because in some cases this has become a list and I don't know yet why
  # For now a corrective rather than a fix until I can identify the issue. Should trigger an error if the vector has the wrong length.
  EV$Modifications <- "Unmodified"
  a1 <- strsplit(gsub(paste0("_|\\)|", aaPat), "", EV$"Modified sequence"[wMdSq2]), "\\(")
  f0 <- function(x) { x[which(x != "")] }
  environment(f0) <- .GlobalEnv
  tmp <- a1 <- parallel::parLapply(cl, a1, f0)
  if (OpenSearch) {
    tmp <- cbind(a1, a3)
    f0 <- function(x) {
      #x <- tmp[1,]
      x1 <- unlist(x[[1]])
      x2 <- unlist(x[[2]])
      x1 <- x1[which(x1 != "!m")]
      if (length(x1)) {
        mds <- aggregate(x1, list(x1), length)
        mds$Group.1 <- Modifs$`Full name`[match(mds$Group.1, Modifs$Mark)]
      }
      if (length(x2)) {
        x2 <- aggregate(x2, list(x2), length)
        if (length(x1)) { mds <- rbind(mds, x2) } else { mds <- x2 }
      }
      mds <- mds[order(mds$Group.1, decreasing = FALSE),]
      paste(apply(mds, 1, function(y) { gsub("^1 ", "", paste(rev(y), collapse = " ")) }), collapse = ",")
    }
    environment(f0) <- .GlobalEnv
    EV$Modifications[wMdSq2] <- parallel::parApply(cl, tmp, 1, f0)
  } else {
    f0 <- function(x) {
      #x <- tmp[1]
      x <- unlist(x)
      x <- aggregate(x, list(x), length)
      x <- x[order(x$Group.1, decreasing = FALSE),]
      x$Group.1 <- Modifs$`Full name`[match(x$Group.1, Modifs$Mark)]
      paste(apply(x, 1, function(y) { gsub("^1 ", "", paste(rev(y), collapse = " ")) }), collapse = ",")
    }
    environment(f0) <- .GlobalEnv
    EV$Modifications[wMdSq2] <- parallel::parSapply(cl, tmp, f0)
  }
  #tst <- unique(gsub("^[0-9]+ ", "", unlist(strsplit(EV$Modifications[wMdSq2], ","))))
  #sum(!tst %in% Modifs$`Full name`)
  # Mass and mass error columns
  EV$"m/z" <- PSMs$"Calibrated Observed M/Z"
  EV$"Theoretical m/z" <- PSMs$"Calculated M/Z"
  EV$Mass <- (EV$"m/z"-1.007276466879)*EV$Charge # This used to be incorrect
  EV$"Uncalibrated - Calibrated m/z [Da]" <- PSMs$"Observed M/Z" - PSMs$"Calibrated Observed M/Z" #Somehow I get all 0s here?!
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
  #length(which(tst$V2 != tst$V3)) == 0
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
  EV <- cbind(data.frame(id = 1:nrow(EV)),
              EV)
  #
  if (isTMT) {
    EV[, paste0("Reporter intensity ", TMTtbl$channel_code)] <- PSMs[, tmtKols]
  }
  EV$Search_ID <- FP_Workflow
  if (isActuallyDIANN) {
    DIANN <- proteoCraft::DIANN_to_MQ(diannRep,
                                      cl = cl)
    #DIANN <- DIANN_to_MQ(diannRep, cl = cl) # For testing
    Res <- list(Evidence = DIANN$Evidence,
                PTMs = DIANN$PTMs,
                FracMap = FP_Mnfst,
                WorkFlow = FP_Wrkflw,
                FP_Evidence = EV,
                FP_PTMs = Modifs,
                PTMs = diannRep)
  } else {
    Res <- list(Evidence = EV,
                PTMs = Modifs,
                FracMap = FP_Mnfst,
                WorkFlow = FP_Wrkflw,
                PTMs = FP_PSMs)
  }
  if (isTMT) { Res[["TMT_annotations"]] <- TMTtbl }
  #
  if (stopCl) { parallel::stopCluster(cl) }
  return(Res)
}
