# Check and process Fractions map
Exp.map$Use <- as.logical(Exp.map$Use)
MQ.Exp <- MQ.Exp[which(MQ.Exp %in% unique(unlist(Exp.map$MQ.Exp[which(Exp.map$Use)])))]
if (file.exists(FracMapPath)) {
  Frac.map %<o% read.csv(FracMapPath, check.names = FALSE)
  if (("Parent sample" %in% colnames(Frac.map))&&(!"MQ.Exp" %in% colnames(Frac.map))) {
    Frac.map$MQ.Exp <- Frac.map$"Parent sample"
  }
  Frac.map <- Frac.map[which(Frac.map$Use),]
  MQ.Exp <- MQ.Exp[which(MQ.Exp %in% unique(unlist(Frac.map$MQ.Exp[which(Frac.map$Use)])))]
  Exp.map <- Exp.map[which(vapply(Exp.map$MQ.Exp, function(x) { sum(x %in% MQ.Exp) }, 1) > 0),]
  Frac.map <- Frac.map[which(Frac.map$MQ.Exp %in% MQ.Exp),]
  #Frac.map$"Raw file" <- gsub("\\\\", "/", Frac.map$"Raw file") # Redundant
  #ev$"Raw file path" <- gsub_Rep("\\\\", "/", ev$"Raw file path")
  m <- match(ev$`Raw file`, rawFiles2)
  if (!sum(!is.na(m))) {
    stop()
    #ev$`Raw file` <- gsub_Rep("\\.[^\\.]+$", "", ev$`Raw file`)
    #m <- match(ev$`Raw file`, rawFiles2)
  }
  #rawFiles <- gsub_Rep("\\\\", "/", rawFiles)
  # if (!"Raw file path" %in% colnames(ev)) { ev$"Raw file path" <- rawFiles[match(ev$`Raw file`, rawFiles2)] } else {
  #   ev$"Raw file path" <- gsub_Rep("\\\\", "/", ev$"Raw file path")
  # }
  # if (SearchSoft != "DIANN") {
  #   Frac.map$"Raw file_Full" <- Frac.map$"Raw file"
  #   if (SearchSoft == "MAXQUANT") {
  #     Frac.map$"Raw file" <- gsub(".*/|\\.((raw)|(mzX?ML)|(d)|(dia))$", "", Frac.map$"Raw file", ignore.case = TRUE)
  #   } else {
  #     Frac.map$"Raw file" <- gsub("\\.((raw)|(mzX?ML)|(d)|(dia))$", "", Frac.map$"Raw file", ignore.case = TRUE)
  #   }
  # }
  if ("MQ.Exp" %in% colnames(ev)) {
    tst <- sum(is.na(ev$MQ.Exp)) == nrow(ev)
    if (tst) { ev$MQ.Exp <- NULL }
  }
  test <- sort(unique(Frac.map$MQ.Exp))
  if (sum(test != sort(MQ.Exp))) { stop("Column \"Experiment\" not defined properly in Fractions map!") }
  Frac.map <- Frac.map[which(Frac.map$MQ.Exp %in% MQ.Exp),]
  Frac.map$Experiment <- vapply(Frac.map$MQ.Exp, function(x) { #x <- Frac.map$MQ.Exp[1]
    x1 <- unlist(unique(Exp.map$Experiment[which(vapply(Exp.map$MQ.Exp, function(y) { x %in% unlist(y) }, TRUE))]))
    if (length(x1) == 1) { return(x1) } else {
      if (length(x1) > 1) { stop(paste0(x, " - each raw file must be mapped to exactly one Experiment!")) } else {
        return(NA)
      }
    }
  }, "")
  Frac.map <- Frac.map[which(!is.na(Frac.map$Experiment)),]
  if (("Replicate" %in% colnames(Frac.map))&&(sum(sort(unique(Frac.map$Replicate)) != sort(Rep)) != 0)) {
    stop("Replicates from Fractions map and Experiment map do not match!")
  }
  # if (SearchSoft == "DIANN") {
  #   if (length(which(!unique(ev$"Raw file path") %in% unique(Frac.map$"Raw file")))) {
  #     warning("Some evidences do not have corresponding raw files in Fractions map and shall thus be removed. Check that this is ok!")
  #     ev <- ev[which(ev$"Raw file path" %in% unique(Frac.map$"Raw file")),]
  #   }
  #   ev2fr %<o% match(ev$"Raw file path", Frac.map$"Raw file")
  #   if (length(which(!unique(Frac.map$"Raw file") %in% unique(ev$"Raw file path")))) {
  #     warning("There are some raw files for which not a single peptide evidence is present in the data!")
  #     Frac.map <- Frac.map[which(Frac.map$"Raw file" %in% unique(ev$"Raw file path")),]
  #   }
  # } else {
  #   ev$"Raw file" <- gsub_Rep("\\.((raw)|(mzX?ML)|(d)|(dia))$", "", ev$"Raw file", ignore.case = TRUE)
  #   Frac.map$"Raw file" <- gsub("\\.((raw)|(mzX?ML)|(d)|(dia))$", "", Frac.map$"Raw file", ignore.case = TRUE)
  #   if (length(which(!unique(ev$"Raw file") %in% unique(Frac.map$"Raw file")))) {
  #     warning("Some evidences do not have corresponding raw files in Fractions map and shall thus be removed. Check that this is ok!")
  #     ev <- ev[which(ev$"Raw file" %in% unique(Frac.map$"Raw file")),]
  #   }
  #   ev2fr %<o% match(ev$"Raw file", Frac.map$"Raw file")
  #   if (length(which(!unique(Frac.map$"Raw file") %in% unique(ev$"Raw file")))) {
  #     warning("There are some raw files for which not a single peptide evidence is present in the data!")
  #     Frac.map <- Frac.map[which(Frac.map$"Raw file" %in% unique(ev$"Raw file")),]
  #   }
  # }
  stopifnot("Raw file" %in% colnames(Frac.map),
            "Raw file path" %in% colnames(ev))
  ev2fr %<o% match(ev$"Raw file path", Frac.map$"Raw file")
  tst <- is.na(ev2fr)
  if (sum(tst)) {
    w <- which(tst)
    for (k in c("Raw file", "Raw files name")) { #k <- "Raw file" #k <-  "Raw files name"
      w2 <- which(Frac.map[[k]] %in% ev$"Raw file"[w])
      if (length(w2)) {
        m <- match(Frac.map[w2, k], ev$"Raw file")
        tst2 <- gsub(".*/|\\.[^\\.]+$", "", Frac.map$"Raw file"[w2]) == gsub(".*/|\\.[^\\.]+$", "", ev$"Raw file path"[m])
        stopifnot(sum(is.na(tst2)) == 0, sum(!tst2) == 0)
        Frac.map$"Raw file"[w2] <- ev$"Raw file path"[m]
        ev2fr %<o% match(ev$"Raw file path", Frac.map$"Raw file")
      }
    }
  }
  tst <- sum(is.na(ev2fr))
  if (tst) {
    warning(paste0("Removing ", tst, " PSMs not matching selected raw files..."))
    ev <- ev[which(!is.na(ev2fr)),]
    ev2fr <- ev2fr[which(!is.na(ev2fr))]
  }
  if (!"MQ.Exp" %in% colnames(ev)) {
    ev$MQ.Exp <- Frac.map$MQ.Exp[ev2fr]
    w <- which(is.na(ev$MQ.Exp))
    if (length(w)) { warning("Some PSMs do not have a corresponding Samples, is this expected?\n(This is normal if you decided to not use all samples...)") }
    if (LabelType == "Isobaric") { Iso <- sort(unique(Exp.map$Isobaric.set)) }
  }
  # if (exists("MQ.Exp2discrd")) {
  #   ev <- ev[which(!ev$MQ.Exp %in% MQ.Exp2discrd),]
  #   ev2fr %<o% match(ev$"Raw file path", Frac.map$"Raw file")
  #   Frac.map <- Frac.map[which(!Frac.map$MQ.Exp %in% MQ.Exp2discrd),]
  # }
  if ("Replicate" %in% colnames(Frac.map)) {
    Frac.map$Unique.Frac.ID <- apply(Frac.map[, c("Replicate", "MQ.Exp", "Fraction")], 1, function(x) {
      paste0("Rep", paste(x, collapse = "_"))
    })
    ev$Replicate <- Frac.map$Replicate[ev2fr]
  } else { Frac.map$Unique.Frac.ID <- do.call(paste, c(Frac.map[, c("MQ.Exp", "Fraction")], sep = "_")) }
  if (LabelType == "Isobaric") {
    if (!"Isobaric.set" %in% colnames(Frac.map)) {
      if (test.iso.set) { stop("The fractions map does not include an isobaric set column!") } else {
        Frac.map$Isobaric.set <- 1
        ev$Isobaric.set <- 1
      }
    } else { ev$Isobaric.set <- Frac.map$Isobaric.set[ev2fr] }
  }
  # Filter for ones to keep
  MQ.Exp <- MQ.Exp[which(MQ.Exp %in% Frac.map$MQ.Exp)]
  Exp.map <- Exp.map[which(vapply(Exp.map$MQ.Exp, function(x) { sum(x %in% MQ.Exp) }, 1) > 0),]
  Frac.map$Fraction <- as.numeric(gsub("^ | $", "", Frac.map$Fraction))
  ev <- ev[which(ev$MQ.Exp %in% MQ.Exp),]
  ev2fr %<o% match(ev$"Raw file path", Frac.map$"Raw file")
  # Final test
  test <- data.frame(Ev = sort(unique(ev$MQ.Exp)),
                     Fraction.map = sort(unique(Frac.map$MQ.Exp)),
                     Exp.map = sort(unique(unlist(Exp.map$MQ.Exp))),
                     MQ.Exp = sort(MQ.Exp))
  if (max(apply(test, 1, function(x) { length(unique(x)) })) > 1) {
    stop(paste0("The MQ.Exp object, the Fractions and Experiment map, and the ", evNm, " file do not match!!!"))
    #apply(test, 1, unique)
  }
  if ("Norma.groups" %in% colnames(Frac.map)) {
    warning("Column \"Norma.groups\" in the Fractions Map is deprecated, use \"PTM-enriched\" instead!")
  }
  Unique.Frac %<o% data.frame(Unique.Frac.ID = unique(Frac.map$Unique.Frac.ID))
  Unique.Frac$Raw.files <- lapply(Unique.Frac$Unique.Frac.ID, function(x) {
    Frac.map$"Raw file"[which(Frac.map$Unique.Frac.ID == x)]
  })
  ev$"Unique Frac" <- Frac.map$Unique.Frac.ID[ev2fr]
} else {
  kol <- "MQ.Exp"
  if (length(Exp) == 1) {
    ev$Experiment <- Exp
    kol <- c(kol, "Experiment")
  } else { stop("There are several Experiments in Exp.map, yet Fractions map is missing!") }
  if (LabelType == "Isobaric") {
    if (length(Iso) == 1) {
      ev$Isobaric.set <- Iso
      kol <- c(kol, Isobaric.set)
    } else { stop("There are several Isobaric sets in Exp.map, yet Fractions map is missing!") }
  }
  Frac.map %<o% set_colnames(aggregate(ev[, kol], list(ev$"Raw file path"), unique), c("Raw file", kol))
  # if (exists("MQ.Exp2discrd")) {
  #   ev <- ev[which(!ev$MQ.Exp %in% MQ.Exp2discrd),]
  #   Frac.map <- Frac.map[which(!Frac.map$MQ.Exp %in% MQ.Exp2discrd),]
  # }
}
# Filter
Exp.map$Use <- as.logical(Exp.map$Use)
Frac.map$Use <- as.logical(Frac.map$Use)
Exp.map <- Exp.map[which(Exp.map$Use),]
Frac.map <- Frac.map[which(Frac.map$Use),]
stopifnot(sum(sort(unique(Exp.map$Experiment)) != sort(unique(Frac.map$Experiment))) == 0)
stopifnot(sum(sort(unique(unlist(Exp.map$MQ.Exp))) != sort(unique(unlist(Frac.map$MQ.Exp)))) == 0)
ev2fr %<o% match(ev$"Raw file path", Frac.map$"Raw file") # Update it
tst <- listMelt(lapply(1:nrow(Exp.map), function(x) { Exp.map$MQ.Exp[x] }), 1:nrow(Exp.map))
tst$L1 <- as.integer(tst$L1)
tst <- tst[order(tst$L1, tst$value),]
MQ.Exp <- sort(unique(tst$value))
if (LabelType == "LFQ") { stopifnot(length(MQ.Exp) == length(unique(MQ.Exp)))}
ev <- ev[which(ev$MQ.Exp %in% MQ.Exp),]
ev <- ev[which(ev$"Raw file path" %in% Frac.map$"Raw file"),]
# Update values
ev2fr %<o% match(ev$"Raw file path", Frac.map$"Raw file")
rawFiles %<o% unique(ev$"Raw file path")
rawFiles2 %<o% unique(ev$`Raw file`)
#
if (!"MQ.Exp" %in% colnames(ev)) { ev$MQ.Exp <- Frac.map$MQ.Exp[ev2fr] }
if (!"Parent sample" %in% colnames(ev)) { ev$"Parent sample" <- ev$MQ.Exp } # Synonym for now
if (!"Experiment" %in% colnames(ev)) { ev$Experiment <- Frac.map$Experiment[ev2fr] }
if (!"Fraction" %in% colnames(ev)) { ev$Fraction <- Frac.map$Fraction[ev2fr] }
if (!"Replicate" %in% colnames(ev)) {
  tmp <- listMelt(Exp.map$MQ.Exp, Exp.map$Replicate)
  ev$Replicate <- tmp$L1[match(ev$MQ.Exp, tmp$value)]
}
tmp <- listMelt(Exp.map$MQ.Exp, Exp.map$Experiment)
ev$Experiment <- tmp$L1[match(ev$MQ.Exp, tmp$value)]

# Intensity columns at PSM level
ev.col %<o% c(Original = "Intensity")
if (LabelType == "Isobaric") {
  ev.ref %<o% c(Original = "Reporter intensity ")
  korkol <- grep(topattern(paste0(ev.ref["Original"], "corrected ")), colnames(ev), value = TRUE)
  ev <- ev[, which(!colnames(ev) %in% korkol)]
}

# Update Factors
FactorsLevels <- setNames(lapply(Factors, function(Fact) {
  unique(Exp.map[[Fact]])
}), Factors)
w <- which(vapply(FactorsLevels, length, 1) > 0)
Factors <- Factors[w]
FactorsLevels <- FactorsLevels[Factors]
