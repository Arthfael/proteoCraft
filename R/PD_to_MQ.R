#' PD_to_MQ
#'
#' @description 
#' Converts a Proteome Discoverer PSMs table to a MaxQuant evidence.txt-like table.
#' As with the other functions of this type, the idea is not to get perfect conversion, but close enough that the data can be fed into this package's analysis scripts.
#' 
#' Has not been updated for a long time, and is probably in need of a refresh!
#' 
#' The output is a list of two items:
#' - an evidence table similar to the one MaxQuant creates ("Evidences")
#' - a modifications table ("PTMs")
#' 
#' @param PD A PSMs table from Proteome Discover imported as a data frame.
#' @param Fixed.mods Fixed modifications which will not be reported in the results (PD reports isobaric tags and fixed modifications, MQ ignores them). Default covers "TMT0plex" to "TMT11plex", "iTRAQ0plex" to "iTRAQ8plex", "Carbamidomethyl" and "N-Ethylmaleimide".
#' @param Abs.quant.col Default = "Area". If not provided, and if reporter intensities are available, we will sum those instead. 
#' @param Isobaric Are there isobaric labelling reporter intensity columns in the table? Default = FALSE
#' @param Label This parameter is required if Isobaric is set to TRUE. Accepted values are: "TMT0plex", "TMTpro0plex", "TMT2plex", "TMT6to11plex", "TMTpro16plex", "iodoTMT6plex", "iodoTMT6plex", "aminoxyTMT0plex", "aminoxyTMT6plex", "iTRAQ4plex" or "iTRAQ8plex"
#' @param filter If set to TRUE, this will remove all PSMs flagged as "Rejected". Default = FALSE
#' @param MaxRank Will only accept search engine ranks lower or equal to the value of this parameter (default = 1)
#' 
#' @examples
#' temp <- PD_to_MQ(PD)
#' ev <- temp$Evidences
#' Modifs <- temp$PTMs
#' 
#' @export

PD_to_MQ <- function(PD,
                     Fixed.mods = c(paste0("TMT", c("", paste0(as.character(0L:11L), "plex"),
                                                    paste0("pro", as.character(0L:16L), "plex"))),
                                    paste0("iTRAQ", c("", paste0(as.character(0L:8L), "plex"))),
                                    "Carbamidomethyl", "N-Ethylmaleimide"),
                     Abs.quant.col = "Area",
                     Isobaric = FALSE,
                     Label,
                     filter = FALSE,
                     MaxRank = 1L) {
  if (Isobaric) {
    stopifnot(Label %in% c("TMT0plex", "TMTpro0plex", "TMT2plex", "TMT6to11plex", "TMTpro16plex", "iodoTMT0plex", "iodoTMT6plex", "aminoxyTMT0plex", "aminoxyTMT6plex", "iTRAQ4plex", "iTRAQ8plex"))
    LabelKlass <- c("TMT", "iodoTMT", "iTRAQ")
    LabelKlass <- LabelKlass[which(vapply(LabelKlass, \(x) { grepl(topattern(x), Label) }, TRUE))]
    # NB: the different iTRAQ channels are only isobaric at a lower resolution than for TMT:
    # In two given TMTs channels, including 10-11plex, while isotopologues are used for encoding channels,
    # the total number of isotopes of each atom is the same. Thus TMT is truly isobaric.
    # For iTRAQ the different channels have the same number of neutrons but not necessarily in the same atoms.
    # Thus one could say that iTRAQ is "isobaric at isotope but not isotopologue resolution".
    LabelMasses <- setNames(c(224.1524778966, 225.1558327344, 229.162932141, 324.2161407859, 329.2265950303, 296.2212261638, 301.2316804082, 144.1, 304.2),
                            c("TMT0plex", "TMT2plex", "TMT6to11plex", "iodoTMT0plex", "iodoTMT6plex", "aminoxyTMT0plex", "aminoxyTMT6plex", "iTRAQ4plex", "iTRAQ8plex"))
    
    if (! Label %in% names(LabelMasses)) { stop(paste0("Label ", Label, " is currently not supported... but it should be!")) }
    LabelMass <- LabelMasses[which(names(LabelMasses) == Label)]
  }
  if (filter) { PD <- PD[which(PD$"PSM Ambiguity" != "Rejected"),] }
  EV <- data.frame(Sequence = toupper(gsub("^\\[.+\\]\\.|\\.\\[.+\\]$", "", PD$"Annotated Sequence")))
  EV$Length <- nchar(EV$Sequence)
  test <- gsub(paste(AA, collapse = "|"), "", EV$Sequence)
  #w <- which(!grepl("X", EV$Sequence))
  w <- which(test == "")
  if (length(w) < nrow(EV)) {
    warning("Removing peptides with unknown amino acids")
    PD <- PD[w,]
    EV <- EV[w,]
  }
  PD$temp <- strsplit(PD$Modifications, "; ")
  temp <- gsub("\\(Prot\\)", "[Prot]", PD$Modifications)
  temp <- strsplit(gsub("\\)", "", gsub("\\(", "_", temp)), "; ")
  PD[,c("temp_mod", "temp_pos")] <- ""
  w <- which(lengths(temp) > 0L)
  PD[w, c("temp_mod", "temp_pos")] <- Isapply(temp[w], \(x) { #x <- temp[[1L]]
    x <- unlist(x)
    x <- as.data.frame(strsplit(x, "_"))
    x1 <- strsplit(as.character(x[2,]), "\\+")
    names(x1) <- x[1L,]
    x1 <- reshape2::melt(x1)
    x1 <- do.call(paste, c(x1, sep = ";"))
    return(x1)
  })
  Modifs <- data.frame(`Full name` = unique(unlist(strsplit(PD$temp_mod, ";"))), check.names = FALSE)
  Modifs$Type <- "Variable"
  Modifs$Type[which(Modifs$"Full name" %in% Fixed.mods)] <- "Fixed"
  Modifs$Mark <- tolower(substr(Modifs$"Full name", 1L, 2L))
  ## Identify affected Amino Acids
  temp <- data.frame(mod = unlist(strsplit(PD$temp_mod, ";")),
                     pos = unlist(strsplit(PD$temp_pos, ";")))
  w <- which(!temp$pos %in% c("N-Term", "N-Term[Prot]", "C-Term", "C-Term[Prot]"))
  temp$pos[w] <- substr(temp$pos[w], 1L, 1L)
  temp <- aggregate(temp$pos, list(temp$mod), unique)
  Modifs$AA <- temp$x[match(Modifs$"Full name", temp$Group.1)]
  ## Sometimes, some marks are duplicates, e.g. if you searched for "Acetyl (Protein N-term)" and "Acetyl (K)" together!
  ## We want to fix this so that each modification has a unique mark:
  test <- aggregate(Modifs$Mark, list(Modifs$Mark), length)
  W <- which(test$x > 1L)
  #Modifs$Mark <- Modifs$"Old mark"
  if (length(W)) {
    Modifs$"Old mark" <- Modifs$Mark
    for (i in W) {
      #i <- W[1L]
      w <- which(Modifs$Mark == test$Group.1[i])
      m <- Modifs[w,]
      m$AA[which(lengths(m$AA) == 0L)] <- "X"
      r <- if ("Acetyl" %in% m$"Full name") { which(m$"Full name" == "Acetyl") } else { 1L }
      s <- c(1L:nrow(m)); s <- s[which(s != r)]
      test <- apply(m[s, c("AA", "Mark")], 1L, \(x) { paste0(tolower(x[[1L]]), substr(x[[2L]], 1L, 1L)) })
      w1 <- which(!test %in% Modifs$Mark)
      m$Mark[s][w1] <- test
      w1 <- which(test %in% Modifs$Mark)
      if (length(w1) > 0L) {
        # not tested
        s <- s[w1]
        test <- sapply(s, list)
        kount <- 1L
        char <- c(as.character(0L:9L),
                  letters)
        taken <- unique(c(Modifs$Mark, m$Mark))
        for (j in s) {
          tst <- paste0(tolower(m$AA[s]), char[kount])
          while (((tst) %in% taken)&&(kount < length(char))) {
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
      Modifs[w,] <- m
    }
  }
  EV$Modifications <- "Unmodified"
  EV$"Modified sequence" <- paste0("_", EV$Sequence, "_")
  temp <- strsplit(PD$temp_mod, ";")
  w <- which(vapply(temp, \(x) { length(x[which(!x %in% Fixed.mods)]) }, 1L) > 0L)
  EV$Modifications[w] <- vapply(temp[w], \(x) {
    x <- x[which(! x %in% Fixed.mods)]
    x <- aggregate(x, list(x), length)
    x <- x[order(x$Group.1), 2L:1L]
    x <- apply(x, 1L, \(y) { gsub("^1 ", "", paste(y, collapse = " ")) } )
    return(paste(x, collapse = ","))
  }, "")
  temp <- PD[w,c("temp_mod", "temp_pos")]
  colnames(temp) <- c("mod", "pos")
  temp$mod <- strsplit(temp$mod, ";")
  temp$pos <- strsplit(temp$pos, ";")
  temp$seq <- strsplit(EV$Sequence[w], "")
  EV$"Modified sequence"[w] <- apply(temp, 1L, \(x) {
    sq <- data.frame(seq = c("_", unlist(x[3]), "_"))
    sq$pos <- 0L:(nrow(sq)-1L)
    sq$mod <- vector("list", nrow(sq))
    md <- unlist(x[1L]); wh <- unlist(x[2L])
    wh <- wh[which(!md %in% Fixed.mods)];md <- md[which(!md %in% Fixed.mods)]
    if (sum(c("N-Term", "N-Term[Prot]") %in% wh)) {
      sq$mod[1L] <- list(md[which(wh %in% c("N-Term", "N-Term[Prot]"))])
    }
    if (sum(c("C-Term", "C-Term[Prot]") %in% wh)) {
      list(sq$mod[nrow(sq)] <- md[which(wh %in% c("C-Term", "C-Term[Prot]"))])
    }
    md <- md[which(!wh %in% c("N-Term", "N-Term[Prot]", "C-Term", "C-Term[Prot]"))]
    wh <- wh[which(!wh %in% c("N-Term", "N-Term[Prot]", "C-Term", "C-Term[Prot]"))]
    if (length(md)) {
      wh <- Isapply(1L:length(wh), \(y) { c(substr(wh[y], 1L, 1L), substr(wh[y], 2L, nchar(wh[y]))) })
      wh[, 2L] <- as.numeric(wh[, 2L])
      test <- sum(apply(wh, 1L, \(y) { sq$seq[as.numeric(y[[2L]])+1L] != y[[1L]] })) == 0L
      if (test) {
        md <- aggregate(md, list(wh[, 2L]), sort)
        sq$mod[which(sq$pos %in% md[, 1L])] <- md[match(sq$pos[which(sq$pos %in% md[, 1L])], md[, 1L]), 2L]
      } else { stop("Amino acids and positions do not match!") }
    }
    wh2 <- which(lengths(sq$mod) > 0L)
    sq$mod[wh2] <- vapply(sq$mod[wh2], \(y) {
      return(paste0("(", paste(Modifs$Mark[match(unlist(y), Modifs$"Full name")], collapse = ","), ")"))
    }, "")
    sq$mod[which(vapply(sq$mod, is.null, TRUE))] <- ""
    return(paste(do.call(paste, c(sq[, c("seq", "mod")], sep = "")), collapse = ""))
  })
  EV$"Leading proteins" <- gsub("; ", ";", PD$"Master Protein Accessions")
  EV$Proteins <- gsub("; ", ";", PD$"Protein Accessions")
  EV$PEP <- PD$"Percolator PEP"
  EV$"Missed cleavages" <- PD$"# Missed Cleavages"
  EV$"Raw file" <- gsub("\\.raw$", "", PD$"Spectrum File", ignore.case = TRUE)
  EV$"Retention time" <- PD$"RT [min]"
  # For some peptides, we sometimes have a better estimate of Retention Time, as the apex of the peak could be detected:
  # This should not change the value by more than 1-2 min
  if ("Apex RT [min]" %in% colnames(PD)) { 
    w <- which(!is.na(PD$"Apex RT [min]"))
    EV$"Retention time"[w] <- PD$"Apex RT [min]"[w]
  }
  EV$Charge <- PD$Charge
  if (!Isobaric) { EV$"m/z" <- PD$"m/z [Da]" } else {
    ## Fix masses to the MQ format (isobaric tag mass removed)
    IsoTaest <- as.data.frame(sapply(unique(unlist(strsplit(PD$temp_mod, ";"))), \(m) {
      vapply(strsplit(PD$temp_mod, ";"), \(x) { length(which(unlist(x) == m)) }, 1L)
    }))
    KTaest <- nchar(EV$Sequence)-nchar(gsub("K", "", EV$Sequence))+1L
    taest <- sweep(IsoTaest, 1L, KTaest, "-")
    taest <- apply(taest, 2L, \(x) { sum(x == 0L) })
    if (max(taest) < nrow(PD)) { warning("There may be peptides with only partial isobaric labelling!") }
    w <- which(taest == max(taest))
    if (length(w) > 1L) { stop("No idea what's happening here!") } else {
      LabelDetails <- names(taest[which(taest == max(taest))])
      NIsobar <- IsoTaest[,w]
      EV$"m/z" <- PD$"m/z [Da]"-LabelMass*(NIsobar)/EV$Charge
    }
  }
  #massH <- as.numeric(gsub("_.*", "", isoProbs$Monoisotopic[match("H", isoProbs$Atom)]))
  #massElectr <- 0.00054858
  #massProt <- massH + massElectr
  massProt <- 1.007276466879
  EV$Mass <- EV$Charge * (EV$"m/z" - massProt)
  # (Using here the average of a hydrogen, not the mass of a proton:
  # these H+ come from the environment and are actually hydrons: not just protons but also some deuterons and tritons)
  EV$"Mass error [Da]" <- PD$"Deltam/z [Da]"
  EV$"Mass error [ppm]" <- PD$"DeltaM [ppm]"
  # Check whether the above is correct all the time!!!
  # Apparently, PD does include a MBR-like feature by default: where are those consensus features?
  if (!is.null(Abs.quant.col)) {
    if (Abs.quant.col %in% colnames(PD)) {
      EV$Intensity <- PD[[Abs.quant.col]]
    } else {
      warning(paste0("Column \"", Abs.quant.col, "\" could not be found in the Proteome Discover PSM table provided!"))
      Abs.quant.col <- NULL
    }
  }
  if (Isobaric) {
    g <- grep("^[0-9]+[N,C]?$", colnames(PD), value = TRUE)
    if (sum(c("131", "131N") %in% g) == 2L) {
      l <- length(which((!is.na(PD$"131"))&(!is.na(PD$"131N"))))
      if (l > 0) {
        warning("Channels named 131 and 131N are both present, and have non-NA values for the same PSMs, so will not be conflated.")
      } else {
        warning("Channels named 131 and 131N are both present and will be conflated.")
        w <- which((!is.na(PD$"131"))&(is.na(PD$"131N")))
        PD$"131N"[w] <- PD$"131"[w]
        PD$"131" <- NULL
        g <- g[which(g != "131")]
      }
    }
    g1 <- paste0("Reporter intensity ", (0L:(length(g)-1L)))
    EV[,g1] <- PD[, g]
    for (i in g1) {
      EV[[i]] <- as.numeric(EV[[i]])
      EV[which(is.na(EV[[i]])),i] <- 0L
    }
    if (is.null(Abs.quant.col)) {
      EV$Intensity <- rowSums(EV[, paste0("Reporter intensity ", (0L:(length(g)-1L)))])
    }
  }
  EV$"Isolation Interference [%]" <- PD$"Isolation Interference [%]"
  EV$"MS/MS scan number" <- PD$"First Scan"
  if ("Contaminant" %in% colnames(PD)) { EV$"Potential contaminant" <- c("", "+")[as.logical(toupper(PD$Contaminant)) + 1L] }
  EV$Reverse <- ""
  for (i in c("Activation Type", "Ion Inject Time [ms]", "Average Reporter S/N", "XCorr", "DeltaScore", "DeltaCn")) {
    if (i %in% colnames(PD)) { EV[[i]] <- PD[[i]] }
  }
  if ("Rank" %in% colnames(PD)) {
    EV$Rank <- PD$Rank
    t1 <- EV$Rank <= MaxRank
  } else { t1 <- rep(TRUE, nrow(EV)) }
  if ("Search Engine Rank" %in% colnames(PD)) {
    EV$"Search Engine Rank" <- PD$"Search Engine Rank"
    t2 <- EV$"Search Engine Rank" <= MaxRank
  } else { t2 <- rep(TRUE, nrow(EV)) }
  EV <- EV[which(t1 & t2),]
  PD <- PD[which(t1 & t2),]
  # Sometimes several search engines were used; are they creating duplicates?
  if (("Identifying Node" %in% colnames(PD))&&(length(unique(PD$"Identifying Node")) > 1L)) {
    EV$"Search Engine" <- PD$"Identifying Node"
    EV$test <- do.call(paste, c(EV[, c("Modified sequence", "Charge", "Retention time", "Raw file", "MS/MS scan number")], sep = "_-_"))
    evtemp <- data.frame(test = unique(EV$test))
    evtemp[, c("Modified sequence", "Charge", "Retention time", "Raw file", "MS/MS scan number")] <- Isapply(strsplit(unique(evtemp$test), "_-_"), unlist)
    evtemp$Spectrum <- do.call(paste, c(evtemp[, c("Raw file", "MS/MS scan number")], sep = "_-_"))
    evtemp$Identification <- do.call(paste, c(evtemp[, c("Modified sequence", "Charge")], sep = "_-_"))
    test <- aggregate(evtemp$Identification, list(evtemp$Spectrum), length)
    w <- which(test$x > 1L)
    if (length(w)) {
      l <- length(which(evtemp$Spectrum %in% test$Group.1[w]))
      warning(paste0("Removing ", l, " conflicting PSMs!"))
      evtemp <- evtemp[which(evtemp$Spectrum %in% test$Group.1[which(test$x == 1L)]),]
    }
    col <- colnames(EV)[which(!colnames(EV) %in% colnames(evtemp))]
    test <- aggregate(EV[,col], list(EV$test), \(x) { length(unique(x)) })
    test2 <- apply(test[,col], 2L, max, na.rm = TRUE)
    w <- which(test2 == 1L)
    test <- aggregate(EV[,col[w]], list(EV$test), unique)
    evtemp[, col[w]] <- test[match(evtemp$test, test$Group.1), col[w]]
    test <- aggregate(EV$"Search Engine", list(EV$test), \(x) { paste(sort(unique(x)), collapse = ";") })
    evtemp$"Search Engine" <- test$x[match(evtemp$test, test$Group.1)]
    w <- which((test2 != 1L)&(col != "Search Engine"))
    if ("PEP" %in% col[w]) {
      test <- aggregate(EV$PEP, list(EV$test), mean, na.rm = TRUE)
      evtemp$PEP <- test$x[match(evtemp$test, test$Group.1)]
    }
    if ("Proteins" %in% col[w]) {
      test <- setNames(strsplit(EV$Proteins, ";"), EV$test)
      test <- reshape2::melt(test)
      test <- aggregate(test$value, list(test$L1), \(x) { paste(sort(unique(x)), collapse = ";") })
      evtemp$Proteins <- test$x[match(evtemp$test, test$Group.1)]
    }
    col <- col[w][which(!col[w] %in% c("PEP", "Proteins"))]
    test3 <- sapply(col, \(x) { class(EV[[x]]) })
    w1 <- which(test3 %in% c("numeric", "double", "integer", "integer64"))
    w2 <- which(test3 == "character")
    if (length(w1)) {
      test <- aggregate(EV[, col[w1]], list(EV$test), mean, na.rm = TRUE)
      evtemp[, col[w1]] <- test[match(evtemp$test, test$Group.1), col[w1]]
    }
    if (length(w2)) {
      test <- aggregate(EV[, col[w2]], list(EV$test), paste, collapse = "___;___")
      evtemp[, col[w2]] <- test[match(evtemp$test, test$Group.1), col[w2]]
    }
    evtemp$test <- NULL
    evtemp$Charge <- as.integer(evtemp$Charge)
    evtemp$"Retention time" <- as.numeric(evtemp$"Retention time")
    evtemp$"MS/MS scan number" <- as.integer(evtemp$"MS/MS scan number")
    EV <- evtemp; rm(evtemp)
  }
  EV$Type <- "MULTI-MSMS"
  EV$id <- 1L:nrow(EV)
  return(list(Evidences = EV, PTMs = Modifs))
}
