#' Skyline_to_MQ
#'
#' @description 
#' Converts a Skyline tsv export table to a table with MaxQuant evidence.txt-like formating.
#' As with the other functions of this type, the idea is not to get perfect conversion, but close enough that the data can be fed into this package's analysis scripts.
#' 
#' The output is a list of 3 items:
#' - Evidence: a data.frame similar to the evidence.txt table MaxQuant creates
#' - PTMs: a modifications table
#' 
#' Surprisingly, the diaNN PSMs table does not include MZ values! Is there a way to do this by adding a command line argument?
#' For now they are either detected from the library file (same folder) or computed.
#' Cysteine alkylation is assumed to be Carbamidomethyl
#' 
#' @param Skyline_fl External Skyline .tsv export file, which can be a PSMs- or Transitions-level table. If transitions, it will aggregate information by file, modified sequence and charge into PSM-style table.
#' @param Fixed.mods A vector of those Fixed modifications which will not be reported in the results (reported by diaNN, but MQ ignores them). Default = "Carbamidomethyl"
#' @param N.clust A limit on the number of vCPUs to use. If left as NULL (default), uses the number of available clusters - 1, to a minimum of 1.
#' @param N.reserved Default = 1. Number of reserved vCPUs the function is not to use. Note that for obvious reasons it will always use at least one.
#' @param cl Already have a cluster handy? Why spend time making a new one, which on top of that may invalidate the old one. Just pass it along!
#' @param shinyOpt One of "popup" (default), "pane", "dialog", or "browser".
#' @param digPattern Default = "KR" for trypsin. Only used for the "Missed cleavages" column if that from Skyline is missing. Vector of amino acids after which the enzyme cleaves. 
#' 
#' @examples
#' temp <- Skyline_to_MQ(Skyline_fl = skyline_fl)
#' ev <- temp$Evidences
#' Modifs <- temp$PTMs
#' 
#' @import data.table
#' @export

Skyline_to_MQ <- function(Skyline_fl,
                          Fixed.mods = c("Carbamidomethyl"),
                          N.clust,
                          N.reserved = 1,
                          cl,
                          shinyOpt = "popup",
                          digPattern = "KR") {
  TESTING <- FALSE
  #proteoCraft::DefArg(proteoCraft::Skyline_to_MQ);TESTING <- TRUE
  #Skyline_fl <- skyline_fl
  #
  if (TESTING) {
    tm1 <<- Sys.time()
    # Note:
    # This is not a perfect alternative to missing but will work in most cases, unless x matches a function imported by a package 
    misFun <- function(x) { return(!exists(deparse(substitute(x)))) }
  } else { misFun <- missing }
  #
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
  shinyOpt <- tolower(shinyOpt)
  if (!shinyOpt %in% c("popup", "dialog", "pane", "browser")) { shinyOpt <- "popup" }
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
  Skyline <- data.table::fread(Skyline_fl, integer64 = "numeric", check.names = FALSE, data.table = FALSE)
  #
  # Filter peptides with ambiguous/non-classical amino acids 
  uSeq <- unique(Skyline$Peptide)
  test <- gsub(paste(proteoCraft::AA, collapse = "|"), "", uSeq)
  test <- test[match(Skyline$Peptide, uSeq)]
  w <- which(test == "")
  lX <- nrow(Skyline)-length(w) 
  if (lX) {
    u <- unique(unlist(strsplit(Skyline$Peptide, "")))
    u <- paste(u[which(!u %in% AA)], collapse = "-")
    warning(paste0("Removing ", lX, " peptide", c("", "s")[(lX > 1)+1], " with unknown amino acids ", u, "!"))
    Skyline <- Skyline[w,]
  }
  #
  if ("Missed Cleavages" %in% colnames(Skyline)) {
    Skyline$"Missed Cleavages" <- as.integer(Skyline$"Missed Cleavages")
  } else {
    digPat <- paste(c("[", digPattern, "]"), collapse = "")
    tmp <- gsub(paste0(digPat, "$"), "", Skyline$Peptide)
    Skyline$"Missed Cleavages" <- nchar(tmp) - nchar(gsub(digPat, "", tmp))
  }
  if ("Library Ion Mobility" %in% colnames(Skyline)) {
    Skyline$"Library Ion Mobility" <- as.numeric(gsub(" .*", "", Skyline$"Library Ion Mobility"))
  }
  #
  minCols <- c("Replicate Name", "Peptide", "Protein", "Missed Cleavages", "Precursor Charge", "Precursor Mz", 
               "Peptide Retention Time", "Min Start Time", "Max End Time", "Is Decoy")
  if (sum(!minCols %in% colnames(Skyline))) {
    stop(paste0("The following columns are required in the Skyline export:\n - ",
                paste(minCols, collapse = "\n - "), "\n"))
  }
  Skyline$"Is Decoy" <- as.logical(Skyline$"Is Decoy")
  #
  fileCols <- c("File Path", "File Name", "Replicate Name") # Order of decreasing priority
  w <- which(fileCols %in% colnames(Skyline))
  if (!length(w)) {
    stop(paste0("At least one of these file identifier columns is required in the Skyline export:\n - ",
                paste(fileCols, collapse = "\n - "), "\n"))
  }
  for (k in fileCols[w]) { Skyline[[k]] <- proteoCraft::gsub_Rep("\\\\", "/", Skyline[[k]]) }
  for (k in fileCols) {
    if (!k %in% colnames(Skyline)) {
      j <- fileCols[which((fileCols %in% colnames(Skyline))&(fileCols != k))][1]
      Skyline[[k]] <- Skyline[[j]]
    }
  }
  #
  modSeqCols <- paste0("Peptide Modified Sequence", c("", " Unimod Ids"))
  w <- which(modSeqCols %in% colnames(Skyline))
  if (!length(w)) {
    stop(paste0("At least one of these modified sequence columns is required in the Skyline export:\n - ",
                paste(modSeqCols, collapse = "\n - "), "\n"))
  }
  Skyline$"Modified sequence" <- Skyline[[modSeqCols[w][1]]]
  modSeqCol <- "Modified sequence"
  intCols <- c("Total Area MS1", "Area") # The order matters!
  intCols <- intCols[which(intCols %in% colnames(Skyline))]
  isTransTable <- (intCols[1] == "Total Area MS1")
  Skyline$Intensity <- Skyline[[intCols[1]]]
  intCol <- "Intensity"
  if (isTransTable) {
    warning("Column \"Total Area MS1\" is present -> we will assume that this is a Transitions table!")
    Skyline2 <- data.table::data.table(Skyline)
    aggrLst <- paste0("list(`Intensity` = sum(`Intensity`, na.rm = TRUE),
                    `Peptide Retention Time` = mean(`Peptide Retention Time`),
                    `Min Start Time` = mean(`Min Start Time`),
                    `Max End Time` = mean(`Max End Time`)",
                      c("", ",\n`Library Ion Mobility` = mean(`Library Ion Mobility`)")[("Library Ion Mobility" %in% colnames(Skyline))+1],
                      c("", ",\n`Missed Cleavages` = mean(`Missed Cleavages`)")[("Missed Cleavages" %in% colnames(Skyline))+1],
                      c("", ",\n`Explicit Ion Mobility` = mean(`Explicit Ion Mobility`)")[("Explicit Ion Mobility" %in% colnames(Skyline))+1],
                      c("", ",\n`Area` = paste(`Area`, collapse = \";\")")[("Area" %in% colnames(Skyline))+1],
                      c("", ",\n`Product Charge` = paste(`Product Charge`, collapse = \";\")")[("Product Charge" %in% colnames(Skyline))+1],
                      c("", ",\n`Detection Q Value` = mean(`Detection Q Value`)")[("Detection Q Value" %in% colnames(Skyline))+1],
                      c("", ",\n`Library Probability Score` = mean(`Library Probability Score`)")[("Library Probability Score" %in% colnames(Skyline))+1],
                      c("", ",\n`Library Score1` = mean(`Library Score1`)")[("Library Score1" %in% colnames(Skyline))+1],
                      ", `Is Decoy` = unique(`Is Decoy`))")
    byLst <- "list(`Replicate Name`, `File Name`, `File Path`, `Peptide`, `Modified sequence`, `Protein`, `Missed Cleavages`, `Precursor Charge`,
                  `Precursor Mz`, `Peptide Retention Time`, `Min Start Time`, `Max End Time`)"
    aggrCall <- paste0("Skyline2 <- Skyline2[, ", aggrLst, ", by = ", byLst, "]")
    aggrCall <- gsub("\n *", " ", aggrCall)
    #cat(aggrCall)         
    eval(parse(text = aggrCall))
    #k <- colnames(Skyline)[which(!colnames(Skyline) %in% colnames(Skyline2))]
    #print(k)
    Skyline <- Skyline2
  }
  optCols <- c("Missed Cleavages", "Explicit Ion Mobility")
  #colnames(Skyline)[which(!colnames(Skyline) %in% c(minCols, fileCols, modSeqCols, intCols, optCols, scoreCols))]
  #
  # Create MQ-like table
  EV <- data.frame("Experiment" = Skyline$"Replicate Name",
                   "Raw file" = Skyline$"File Name",
                   "Raw file path" = Skyline$"File Path",
                   "Sequence" = Skyline$Peptide,
                   "Modified sequence" = paste0("_",
                                                gsub("\\]", ")",
                                                     gsub("\\[", "(",
                                                          Skyline[[modSeqCol]])), "_"),
                   "Proteins"= gsub(" / ", ";", Skyline$Protein),
                   "Intensity" = as.numeric(Skyline[[intCol]]),
                   "Charge" = as.integer(Skyline$"Precursor Charge"),
                   "m/z" = as.numeric(Skyline$"Precursor Mz"),
                   "Missed cleavages" = Skyline$"Missed Cleavages",
                   "Retention time" = as.numeric(Skyline$"Peptide Retention Time"),
                   "Retention time (start)" = as.numeric(Skyline$"Min Start Time"),
                   "Retention time (end)" = as.numeric(Skyline$"Max End Time"),
                   check.names = FALSE)
  EV$"Retention length" <- EV$"Retention time (end)" - EV$"Retention time (start)"
  #
  isoProbs <- proteoCraft::IsotopeProbs
  #massH <- as.numeric(gsub("_.*", "", isoProbs$Monoisotopic[match("H", isoProbs$Atom)]))
  #massElectr <- 0.00054858
  #massProt <- massH + massElectr
  massProt <- 1.007276466879
  EV$Mass <- EV$Charge * (EV$"m/z" - massProt)
  # 
  # PEP
  if ("Detection Q Value" %in% colnames(Skyline)) {
    EV$PEP <- as.numeric(Skyline$"Detection Q Value")
  }
  #
  # Reverse
  EV$"Reverse" <- c("", "+")[Skyline$"Is Decoy" + 1]
  #
  # Score
  scoreCols <- c("Library Score1", "Library Probability Score") # The order matters slightly
  scoreCols <- scoreCols[which(scoreCols %in% colnames(Skyline))]
  if (length(scoreCols)) {
    EV$Score <- Skyline[[scoreCols[1]]]
  }
  #
  # Process modifications
  wMod <- grep("\\(", EV$`Modified sequence`)
  mWMod <- match(unique(EV$`Modified sequence`[wMod]), EV$`Modified sequence`)
  tmp <- data.frame(Seq = EV$Sequence[mWMod],
                    ModSeq = EV$`Modified sequence`[mWMod]#,
                    #Z = EV$Charge[mWMod]
                    )
  gU <- grep("UniMod:", tmp$ModSeq, ignore.case = TRUE)
  if (sum(gU)) {
    EV$"Modified sequence"[wMod] <- gsub("UniMod:", "UniMod:", EV$"Modified sequence"[wMod], ignore.case = TRUE)
    tmp$ModSeq[gU] <- gsub("UniMod:", "UniMod:", tmp$ModSeq[gU], ignore.case = TRUE)
  }
  tmp$Comp <- gsub("\\)[A-Z]*\\(", "___",
                   gsub("^_[A-Z]*\\(|\\)[A-Z]*_$|^_[A-Z]+_$", "",
                        gsub("\\[", "(",
                             gsub("\\]", ")", tmp$ModSeq))))
  # f0 <- function(x) { Peptides::mz(x[[1]], as.integer(x[[2]]), cysteins = 0) }
  # environment(f0) <- .GlobalEnv
  # tmp$Theoretical <- parallel::parApply(cl, tmp[, c("Seq", "Z")], 1, f0)
  tmp$Comp <- gsub("\\)[A-Z]*\\(", "___",
                   gsub("^_[A-Z]*\\(|\\)[A-Z]*_$|^_[A-Z]+_$", "",
                        gsub("\\[", "(",
                             gsub("\\]", ")", tmp$ModSeq))))
  tmp$Comp <- strsplit(tmp$Comp, "___")
  comps <- unique(unlist(tmp$Comp))
  gC <- sum(grepl("^UniMod:", comps, ignore.case = TRUE))
  l <- length(comps)
  modsType <- "No idea"
  if (gC == l) { modsType <- "UniMod" }
  if (modsType == "No idea") {
    if (shinyOpt == "popup") {
      #runApp1 <- "print(shiny::shinyApp(proteoCraft::DIANN_to_MQ_ui1(tstMap), proteoCraft::DIANN_to_MQ_server1, options = list(height = screenRes$height, width = screenRes$width)))"
      runApp2 <- "print(shiny::shinyApp(proteoCraft::DIANN_to_MQ_ui2, proteoCraft::DIANN_to_MQ_server2, options = list(height = screenRes$height, width = screenRes$width)))"
    }
    if (shinyOpt %in% c("dialog", "pane")) {
      if (shinyOpt == "pane") {
        myViewer <- shiny::paneViewer("Viewer", minHeight = "maximize")
      }
      if (shinyOpt == "dialog") {
        myViewer <- shiny::dialogViewer("Viewer", width = screenRes$width, height = screenRes$height)
      }
      #runApp1 <- "shiny::runGadget(shiny::shinyApp(proteoCraft::DIANN_to_MQ_ui1(tstMap), proteoCraft::DIANN_to_MQ_server1, options = list(height = screenRes$height, width = screenRes$width)), viewer = myViewer, stopOnCancel = FALSE)"
      runApp2 <- "shiny::runGadget(shiny::shinyApp(proteoCraft::DIANN_to_MQ_ui2, proteoCraft::DIANN_to_MQ_server2, options = list(height = screenRes$height, width = screenRes$width)), viewer = myViewer, stopOnCancel = FALSE)"
    }
    if (shinyOpt == "browser") {
      #runApp1 <- "print(shiny::shinyApp(proteoCraft::DIANN_to_MQ_ui1(tstMap), proteoCraft::DIANN_to_MQ_server1, options = list(height = \"100%\", width = \"100%\", launch.browser = TRUE)))"
      runApp2 <- "print(shiny::shinyApp(proteoCraft::DIANN_to_MQ_ui2, proteoCraft::DIANN_to_MQ_server2, options = list(height = \"100%\", width = \"100%\", launch.browser = TRUE)))"
    }
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
    # parallel::clusterExport(cl, "Mods", envir = environment())
    # f0 <- function(x) {
    #   x <- aggregate(x, list(x), length)
    #   res <- setNames(rep(0, length(Mods)), Mods)
    #   w <- which(Mods %in% x$Group.1)
    #   m <- match(Mods[w], x$Group.1)
    #   res[w] <- x$x[m]
    #   return(res)
    # }
    # environment(f0) <- .GlobalEnv
    # d.matr <- t(parallel::parSapply(cl, tmp$Mods, f0))
    # fit <- lm.fit(d.matr, tmp$DM)
    Mods <- data.frame(Name = Mods, #DM = fit$coefficients,
                       origMark = gsub(".*mod_", "", Mods))
    Mods$DM <- suppressWarnings(as.numeric(Mods$origMark))
    Mods$AA <- gsub("mod_.*", "", Mods$Name)
    precVal <- 1
    Mods$Options <- apply(Mods[, c("DM", "AA")], 1, function(x) {
      dm <- as.numeric(x[[1]])
      if (substr(x[[2]], 1, 1) == "_") {
        w <- which((abs(UniMod$MonoMass - dm) <= precVal)&(UniMod$Site %in% c(substr(x[[2]], 2, 2), "N-term")))
      } else {
        w <- which((abs(UniMod$MonoMass - dm) <= precVal)&(UniMod$Site == x[[2]]))
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
    #runApp2 <- "print(shiny::shinyApp(DIANN_to_MQ_ui2, DIANN_to_MQ_server2, options = list(height = screenRes$height, width = screenRes$width)))"
    eval(parse(text = runApp2))
    con <- file(fl, "r")
    Mods <- unserialize(con)
    suppressWarnings(suppressMessages(try(close(con), silent = TRUE))) # Only necessary in function mode
    unlink(fl)
    allPTMs <- aggregate(Mods$DM, list(Mods$Match), function(x) {
      rg <- c(min(x), max(x))
      stopifnot(rg[2]-rg[1] < 0.010)
      return(mean(x))
    })
    colnames(allPTMs) <- c("Full name", "Mass delta")
    tmp <- aggregate(Mods$origMark, list(Mods$Match), unique)
    allPTMs$origMark <- tmp$x[match(allPTMs$`Full name`, tmp$Group.1)]
    allPTMs$Position <- "Anywhere"
    allPTMs$Position[grep("C-term", allPTMs$"Full name", ignore.case = TRUE)] <- "Any C-term"
    allPTMs$Position[grep("N-term", allPTMs$"Full name", ignore.case = TRUE)] <- "Any N-term"
    allPTMs$Position[grep("protein C-term", allPTMs$"Full name", ignore.case = TRUE)] <- "Protein C-term"
    allPTMs$Position[grep("protein N-term", allPTMs$"Full name", ignore.case = TRUE)] <- "Protein N-term"
  }
  if (modsType == "UniMod") {
    allPTMs <- data.frame("UniMod" = comps,
                          origMark = comps,
                          check.names = FALSE)
    allPTMs$"Full name" <- UniMod$Name[match(gsub("UniMod:", "", allPTMs$UniMod),
                                             UniMod$UnimodId)]
  }
  if (modsType == "StrangerThings") {
    allPTMs <- data.frame("Full name" = comps,
                          origMark = comps,
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
    if (length(w)) { allPTMs$Site[w] <- lapply(allPTMs$Site[w], function(x) { paste0("n", x) }) }
    w <- which(allPTMs$Position == "Any C-term")
    if (length(w)) { allPTMs$Site[w] <- lapply(allPTMs$Site[w], function(x) { paste0("c", x) }) }
    w <- which(allPTMs$Position == "Protein N-term")
    if (length(w)) { allPTMs$Site[w] <- "[^" }
    w <- which(allPTMs$Position == "Protein C-term")
    if (length(w)) { allPTMs$Site[w] <- "^]" }
    allPTMs$Site_long <- lapply(strsplit(gsub("\\^", "", lapply(allPTMs$Site, paste, collapse = "")), ""), function(x) {
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
    for (i in W) { #i <- 7 #i <- W[1]
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
  if (modsType == "No idea") {
    allPTMs$"Old mark" <- allPTMs$origMark
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
  # #
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
  #
  EV <- cbind(data.frame(id = 1:nrow(EV)),
                         EV)
  Res <- list(Evidence = EV,
              PTMs = allPTMs)
  if (cleanUp) { parallel::stopCluster(cl) }
  return(Res)
}
