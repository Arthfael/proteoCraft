#' OpMS_to_MQ
#' 
#' @description
#' Very old/deprecated, but could be the basis of an updated function should we start working with OpenMS again.
#' Process OpenMS csv output file into a valid PSM file
#' Not yet tested for:
#' - Labelled samples.
#' - Cases where a single amino-acid has been assigned multiple PTMs.
#'
#' @param OpMS The name of a results csv file from OpenMS.
#' @param Fixed.mods Fixed modifications. Will not be reported in the results.
#' @param Contaminants Database of common contaminant proteins used for the search. Set either to "MQ" (default; in this case the next argument is used), or directly to a valid fasta database.

#' @param FixUniProtIDs If set to TRUE, UniProtKB IDs of the form "xy|ACCESS|NAME_SPECIE" will be changed into the shorter "ACCESS" (for accession) format.
#' @param Isobaric Are there isobaric labelling reporter intensity columns in the table? Default = FALSE
#' @param Label This parameter is required if Isobaric is set to TRUE. Accepted values are: "TMT0plex", "TMTpro0plex", "TMT2plex", "TMT6to11plex", "TMTpro16plex", "iodoTMT6plex", "iodoTMT6plex", "aminoxyTMT0plex", "aminoxyTMT6plex", "iTRAQ4plex" or "iTRAQ8plex"
#' @param MaxRank Will only accept search engine ranks lower or equal to the value of this parameter (default = 1)
#' 
#' @examples
#' temp <- OpMS_to_MQ("Results.csv")
#' ev <- temp$Evidences
#' Modifs <- temp$PTMs
#' 

OpMS_to_MQ <- function(OpMS,
                       Fixed.mods = c(paste0("tmt", c("", paste0(0:11, "plex"), paste0("pro", 0:16, "plex"))), #Not sure about that one, I need to test it.
                                      paste0("itraq", c("", paste0(0:8, "plex"))), #Same as above.
                                      "carbamidomethyl", "nethylmaleimide"),
                       Contaminants = "MQ",
                       MQ_loc = "C:/MaxQuant",
                       FixUniProtIDs = TRUE,
                       Isobaric = FALSE,
                       Label,
                       MaxRank = 0) {
  #
  # Useful functions:
  colconvert <- function(x, tag) {
    tag <- paste0("^#", tag)
    rs <- unlist(strsplit(gsub(tag, "", x[rev(grep(tag, x))[1]]), "\t"))
    return(rs[2:length(rs)])
  }
  dataextr <- function(x, tag) {
    rs <- as.data.frame(t(as.data.frame(sapply(x[grep(paste0("^", tag, "\t"), x)],
                                               strsplit, split = "\t"))))
    rs <- rs[, 2:ncol(rs)]
    rownames(rs) <- NULL
    colnames(rs) <- colconvert(x, tag)
    return(rs)
  }
  #
  # Column names mapping:
  ColMappings <- readLines("https://raw.githubusercontent.com/HUPO-PSI/mzIdentML/master/validator/trunk/src/main/java/psidev/psi/pi/validator/objectrules/AObjectRule.java")
  ColMappings <- grep(" +put\\(", ColMappings, value = TRUE)
  ColMappings <- gsub("^ +put\\(\"|\"\\);$", "", ColMappings)
  ColMappings <- as.data.frame(t(sapply(strsplit(ColMappings, "\", \""), unlist)))
  colnames(ColMappings) <- c("Code", "Name")
  #
  # Load and process "Results.csv"
  ResCSV <- readLines(OpMS)
  #ResCSV <- readLines("Results.csv")
  ResCSV <- ResCSV[which(!grepl("^#?UNASSIGNEDPEPTIDE\t|^#?PROTEIN\t|^#?RUN\t", ResCSV))]
  #write(ResCSV, "test.txt")
  for (i in c("MAP", "CONSENSUS", "PEPTIDE")) {
    assign(paste0(substr(i, 1, 1), tolower(substr(i, 2, 4))), dataextr(ResCSV, i))
  }
  for (i in as.character(sapply(c("rt_", "mz_", "intensity_", "charge_", "width_", "quality_"), function(x) {
    paste0(x, c("cf", Map$id))
  }))) { if (i %in% colnames(Cons)) { Cons[[i]] <- as.numeric(Cons[[i]]) } }
  for (i in c("rt", "mz", "score", "rank", "charge", "map_index", "E-Value", "b_score", "delta", "hyperscore", "mass",
              "nextscore", "q-value_score", "y_score")) {
    if (i %in% colnames(Pept)) { Pept[[i]] <- as.numeric(Pept[[i]]) }
  }
  #
  # Fix some column names:
  for (i in c("charge", "sequence", "score", "rank")) {
    colnames(Pept)[which(colnames(Pept) == i)] <- paste0(toupper(substr(i, 1, 1)),
                                                         substr(i, 2, nchar(i)))
  }
  colnames(Pept)[match("rt", colnames(Pept))] <- "MS.MS.retention.time"
  Pept$MS.MS.scan.number <- gsub(".+scan=", "", Pept$spectrum_reference)
  Pept$MS.MS.retention.time <- Pept$MS.MS.retention.time/60
  for (i in c("rt_cf", paste0("rt_", Map$id))) { if (i %in% colnames(Cons)) { Cons[[i]] <- Cons[[i]]/60 } }
  ColMappings <- ColMappings[which(ColMappings$Code %in% colnames(Pept)),]
  if (nrow(ColMappings) > 0) {
    for (i in ColMappings$Code) { #i <- ColMappings$Code[1]
      col <- ColMappings$Name[which(ColMappings$Code == i)]
      colnames(Pept)[which(colnames(Pept) == i)] <- col
      tst <- tryCatch(as.numeric(Pept[[col]]), error = function(e) { e }, warning = function(w) { w })
      if (class(tst) == "numeric") { Pept[[col]] <- tst }
    }
  }
  #
  # Map peptides to consensus:
  Cons$id <- 1:nrow(Cons)
  Pept$id <- 1:nrow(Pept)
  colnames(Cons)[which(colnames(Cons) == "mz_cf")] <- "m.z"
  colnames(Cons)[which(colnames(Cons) == "charge_cf")] <- "Charge"
  g_con <- grep("^CONSENSUS\t", ResCSV)
  g_pep <- grep("^PEPTIDE\t", ResCSV)
  Pept$Consensus.ID <- sapply(g_pep, function(x) {
    Cons$id[match(max(g_con[which(g_con < x)]), g_con)]
  })
  temp <- aggregate(Pept$id, list(Pept$Consensus.ID), paste, collapse = ";")
  Cons$Peptide.IDs <- NA
  w <- which(Cons$id %in% temp$Group.1)
  Cons$Peptide.IDs[w] <- temp$x[match(Cons$id[w], temp$Group.1)]
  Cons <- Cons[which(!is.na(Cons$Peptide.IDs)),]
  Pept[, c("Intensity", "Retention.time")] <- as.data.frame(t(apply(Pept[,c("Consensus.ID", "map_index")], 1, function(x) {
    unlist(Cons[match(x[[1]], Cons$id), paste0(c("intensity_", "rt_"), x[[2]])])
  })))
  #
  Pept$m.z <- Cons$m.z[match(Pept$Consensus.ID, Cons$id)]
  Pept$Mass <- Pept$Charge*(Pept$m.z-1.007276466879) # (monoisotopic)
  #
  # Correct modified sequence:
  Pept$Modified.sequence <- Pept$Sequence
  Modifs <- data.frame(Old.full.name = unique(unlist(sapply(strsplit(Pept$Modified.sequence, "\\(|\\)"), function(x) { #x <- strsplit(Pept$Modified.sequence, "\\(|\\)")[1]
    l <- length(x)
    return(x[which(!as.logical((1:l %% 2)))])
  }))))
  Modifs$Full.name <- tolower(Modifs$Old.full.name)
  for (i in 1:nrow(Modifs)) { #i <- 1
    Pept$Modified.sequence <- gsub(proteoCraft::topattern(Modifs$Old.full.name[i], start = FALSE),
                                   Modifs$Full.name[i], Pept$Modified.sequence)
  }
  Modifs$Type <- "Variable"
  Modifs$Type[which(Modifs$Full.name %in% tolower(Fixed.mods))] <- "Fixed"
  Modifs$Mark <- tolower(substr(Modifs$Full.name, 1, 2))
  g1 <- grep("^mod:[0-9]+$", Modifs$Full.name)
  g2 <- which(!grepl("^mod:[0-9]+$", Modifs$Full.name))
  if (length(g1) > 0) {
    require(unimod)
    Mds <- unimod::modifications
    for (i in g1) { #i <- g1[1]
      id <- gsub("^mod:0*", "", Modifs$Full.name[i])
      md <- Mds[which(Mds$UnimodId == id),]
      nm <- gsub("\\(|\\)", "", tolower(unique(md$Name))) # I checked, should always be length = 1
      Pept$Modified.sequence <- gsub(paste0("\\(", Modifs$Full.name[i], "\\)"),
                                     paste0("(", nm, ")"), Pept$Modified.sequence)
      Modifs$Full.name[i] <- nm
      ma <- substr(tolower(Modifs$Full.name[i]), 1, 2)
      if (ma %in% Modifs$Mark) {
        ma <- c(substr(ma, 1, 1), 0)
        while ((paste(ma, collapse = "") %in% Modifs$Mark)&&(ma[2] < 10)) {
          ma[2] <- ma[2]+1
        }
        ma <- paste(ma, collapse = "")
        if (ma %in% Modifs$Mark) { stop("Well that is unexpected!") }
      }
      Modifs$Mark[i] <- ma
    }
  }
  if (length(g2) > 0) {
    for (i in g2) { #i <- g2[1]
      Modifs$Full.name[i] <- tolower(Modifs$Full.name[i])
      Pept$Modified.sequence <- gsub(paste0("\\(", Modifs$Full.name[i], "\\)"),
                                     paste0("(", Modifs$Full.name[i], ")"),
                                     Pept$Modified.sequence, ignore.case = TRUE)
    }
  }
  ## Identify affected AAs
  Modifs$AA <- rep(list(""), nrow(Modifs))
  for (i in 1:nrow(Modifs)) { #i <- 1
    g <- grep(proteoCraft::topattern(paste0("(", Modifs$Full.name[i], ")"), start = FALSE), Pept$Modified.sequence)
    g1 <- grep(Modifs$Full.name[i], Pept$Modified.sequence)
    if (sum(g1 != g) > 0) { stop("Unexpected behaviour!")}
    if (length(g) > 0) {
      tmp <- Pept$Modified.sequence[g]
      for (j in c(".", AA, ")")) { tmp <- gsub(proteoCraft::topattern(j, start = FALSE), paste0(j, "_"), tmp) }
      tmp <- strsplit(gsub("\\(|\\)", "", tmp), "_")
      aa <- sort(unique(unlist(sapply(tmp, function(x) { x[which(x == Modifs$Full.name[i])-1]}))))
      if ("." %in% aa) {
        g2 <- grep(proteoCraft::topattern(paste0(".(", Modifs$Full.name[i], ")"), start = FALSE), Pept$Modified.sequence)
        tmp2 <- Pept$Modified.sequence[g2]
        for (j in c(".", AA, ")")) { tmp2 <- gsub(proteoCraft::topattern(j, start = FALSE), paste0(j, "_"), tmp2) }
        tmp2 <- strsplit(gsub("\\(|\\)", "", tmp2), "_")
        aa2 <- sort(unique(unlist(sapply(tmp2, function(x) { x[which(x == Modifs$Full.name[i])+1]}))))
        if (length(aa2) <= 4) { # (Arbitrary threshold... may change)
          aa2 <- paste0("_", aa2)
        } else { aa2 <- "_" }
        aa <- c(aa[which(aa != ".")], aa2)
      }
      Modifs$AA[i] <- list(aa)
    } else { stop("We seem to have corrupted the Sequence column!") }
  }
  # Modifications column
  Pept$Modifications <- gsub(paste(c(AA, "^\\."), collapse = "|"), "", Pept$Modified.sequence)
  g <- grepl("\\(|\\)", Pept$Modifications)
  Pept$Modifications[which(!g)] <- "Unmodified"
  Pept$Modifications[which(g)] <- sapply(strsplit(gsub("^\\(|\\)$", "", Pept$Modifications[which(g)]), "\\)\\("), function(x) { #x <- strsplit(gsub("^\\(|\\)$", "", temp), "\\)\\(")[1]
    x <- unlist(x)
    x <- aggregate(x, list(x), length)
    x$x[which(x$x == 1)] <- ""
    return(paste(gsub("^ ", "", apply(x, 1, function(x) {paste(x[[2]], x[[1]])})), collapse = ", "))
  })
  # Replace modification names by marks
  for (i in 1:nrow(Modifs)) {
    Pept$Modified.sequence <- gsub(proteoCraft::topattern(Modifs$Full.name[i], start = FALSE), Modifs$Mark[i], Pept$Modified.sequence)
  }
  Pept$Modified.sequence <- gsub("^__", "_", paste0("_", gsub("^\\.", "_", Pept$Modified.sequence), "_"))
  Pept$Sequence <- gsub("[^A-Z]", "", Pept$Sequence)
  Pept$Length <- nchar(Pept$Sequence)
  colnames(Pept)[which(colnames(Pept) == "accessions")] <- "Proteins"
  if (FixUniProtIDs) {
    temp <- reshape2::melt(setNames(strsplit(Pept$Proteins, ";"), Pept$id))
    temp$value <- gsub("^[^\\|]{2}\\||\\|[^\\|]+$", "", temp$value)
    temp <- aggregate(temp$value, list(temp$L1), function(x) {paste(x, collapse = ";")})
    Pept$Proteins <- temp$x[match(Pept$id, temp$Group.1)]
  }
  Map$Raw.file <- gsub("^.+/", "", Map$filename)
  raw.files <- unique(Map$Raw.file)
  Pept$Raw.file <- Map$Raw.file[match(Pept$map_index, Map$id)]
  Pept$Type <- "MULTI-MSMS"
  Pept$Reverse <- c("", "+")[grepl("decoy", Pept$target_decoy)+1]
  Pept$target_decoy <- NULL
  # Identify consensus features not derived from a peptide, which we will now be adding to our table:
  # (This is the Match-Between-Runs equivalent)
  kol <- c("Mass", "Sequence", "Modified.sequence", "Modifications", "aa_before",
           "aa_after", "Proteins", "Length", "Reverse")
  if (!is.null(Contaminants)) {
    Cont_ok <- FALSE
    if (Contaminants == "MQ") {
      if (dir.exists(MQ_loc)) {
        MQ_loc <- gsub("/+$", "", MQ_loc)
        g <- grep("^maxquant_", list.files(MQ_loc), value = TRUE, ignore.case = TRUE)
        if (length(g)) {
          if (length(g) > 1) {
            tst <- as.data.frame(t(sapply(strsplit(gsub("^MaxQuant_", "", g), "\\."), function(x) { as.numeric(unlist(x)) })))
            tst2 <- apply(tst, 2, function(x) { 
              x <- ceiling(log10(x+1)); x[which(x == 0)] <- 1
              return(x)
            })
            tst2 <- apply(tst2, 2, max)
            tst <- sapply(1:ncol(tst), function(x) {
              c1 <- as.character(tst[,x]) ; c2 <- tst2[x]-nchar(c1)
              c2 <- sapply(c2, function(y) {paste(rep(0, y), collapse = "")})
              return(sapply(1:length(c1), function(y) { paste(c2[y], c1[y], sep = "") }))
            })
            tst <- as.numeric(apply(tst, 1, paste, collapse = ""))
            g <- g[order(tst, decreasing = TRUE)[1]]
          }
          MQ_loc <- paste0(MQ_loc, "/", g)
        }
        MQ_loc <- paste0(MQ_loc, "/bin/conf/contaminants.fasta")
        if (file.exists(MQ_loc)) { 
          Cont_ok <- TRUE
          Contaminants <- proteoCraft::Format.DB(MQ_loc)
        } else { warning("The Contaminants fasta database was not found at the location specified in the \"Contaminants\" argument.\n-> no \"Potential.contaminant\" column was created!") }
      } else { warning("Invalid MaxQuant install folder (argument \"MQ_loc\"), no contaminants annotations database could be found.\n-> no \"Potential.contaminant\" column was created!") }
    } else {
      if (file.exists(Contaminants)) {
        Contaminants <- proteoCraft::Format.DB(Contaminants)
        Cont_ok <- TRUE
      } else { warning("The Contaminants fasta database was not found at the location specified in the \"Contaminants\" argument.\n-> no \"Potential.contaminant\" column was created!") }
    }
    if (Cont_ok) {
      temp <- reshape2::melt(setNames(strsplit(Pept$Proteins, ";"), Pept$id))
      temp$Contaminant <- temp$value %in% Contaminants$Protein.ID
      temp <- aggregate(temp$Contaminant, list(temp$L1), sum)
      Pept$Potential.contaminant <- temp$x[match(Pept$id, temp$Group.1)] > 0
      Pept$Potential.contaminant <- c("", "+")[Pept$Potential.contaminant + 1]
      kol <- c(kol, "Potential.contaminant")
    }
  }
  temp <- aggregate(Pept[, kol], list(Pept$Consensus.ID), unique)
  Cons[, kol] <- temp[match(Cons$id, temp$Group.1), kol]
  for (i in Map$id) { #i <- Map$id[1]
    w <- which(!Cons$id %in% Pept$Consensus.ID[which(Pept$map_index == i)])
    temp <- Cons[w, c("m.z", "Charge", paste0(c("rt_", "intensity_"), i), "id", kol)]
    temp <- temp[which(proteoCraft::is.all.good(temp$Intensity, 2)), , drop = FALSE]
    temp <- temp[which(temp$Intensity > 0), , drop = FALSE]
    if (nrow(temp) > 0) {
      colnames(temp)[3:5] <- c("Retention.time", "Intensity", "Consensus.ID")
      temp$Raw.file <- Map$Raw.file[which(Map$id == i)] 
      temp$map_index <- i
      temp$Type <- "MULTI-MATCH"
      temp$id <- 1:nrow(temp) + nrow(Pept)
      kol2 <- colnames(Pept)[which(!colnames(Pept) %in% colnames(temp))]
      for (k in kol2) { temp[[k]] <- NA }
      Pept <- rbind(Pept, temp)
    }
  }
  if ("Rank" %in% colnames(Pept)) { Pept <- Pept[which(Pept$Rank <= MaxRank),] }
  for (i in c("mass", "mz", "map_index", "protein_references", "search_identifier", "feature_id", "nextscore", "spectrum_reference",
              "Consensus.ID")) {
    Pept[[i]] <- NULL
  }
  return(list(Evidences = Pept, PTMs = Modifs))
}
