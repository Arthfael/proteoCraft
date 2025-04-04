#' Pepper_ProcessData
#'
#' @description
#' A function to take an ev data.frame and convert it to a Pepper-compatible table. This can then be one-hot encoded for the purpose
#' of training a model on it, or a pre-trained model can be applied to it.
#' 
#' 
#' Note:
#' An issue may be silent fixed modifications (TMT, cysteine alkylations) in MaxQuant. The final model should thus specify
#' what type of data it is suitable for.
#' 
#' @param Ev The evidences (PSMs) data frame.
#' @param Modifs Modifications table as output by DIANN_to_MQ, FP_to_MQ, PD_to_MQ... Should ideally contain a "Mass shift" column, but for MQ workflows this will be usually missing: two additional arguments are provided to deal with those cases: SearchSoft and MQFold.
#' @param experiments.map The experiments map.
#' @param Ev_Sample.col Name of the sample column in Ev
#' @param experiments.map_Sample.col Name of the sample column in experiments.map
#' @param path The path to which to save the file.
#' @param SearchSoft Search software used. Important because Modifs may miss the "Mass shift" column for data originating from some software (MaxQuant). If it is missing, and software is MaxQuant, then providing the next argument (MQFold) can allow parsing directly modifications and calculating the associated mass shift. 
#' @param MQFold Main MaxQuant folder, containing /bin/conf subfolder with modification xml(s). Only used if SearchSoft == "MaxQuant"
#' @param intCol Which intensity column to use? Default = "Intensity"
#' @param filter If TRUE, the data is filtered to create a training set. If FALSE (default), the data isn't filtered (used for preparing data to which a pre-trained model will be transferred).
#' 
#' 
#' 
#' @import data.table
#' @export

Pepper_ProcessData <- function(Ev,
                               Modifs,
                               experiments.map,
                               Ev_Sample.col,
                               experiments.map_Sample.col,
                               path = "",
                               SearchSoft = SearchSoft,
                               MQFold,
                               intCol = "Intensity",
                               filter = FALSE) {
  # Has received a pinch of data.table
  # Can probably be further accelerated by adding some more + maybe some parallelization.
  # But, not so bad for now.
  #
  # !!!This function currently assumes a tryptic digest!!!
  #
  #proteoCraft::DefArg(proteoCraft::Pepper_ProcessData)
  #Ev <- ev; experiments.map <- Exp.map; Ev_Sample.col = "Sample"; experiments.map_Sample.col = "Ref.Sample.Aggregate"; path = paste0(wd, "/ModPep4Pepper.tsv")
  #Ev <- ev; experiments.map <- SamplesMap; Ev_Sample.col = "Experiment"; experiments.map_Sample.col = "MQ.Exp"; path = paste0(wd, "/ModPep4Pepper.tsv")
  #Ev <- ev; experiments.map <- smplsMap; Ev_Sample.col = smplKol; experiments.map_Sample.col = smplKol2; path = paste0(pepDir, "/ModPep4Pepper.tsv") 
  #
  if (!"Mass shift" %in% colnames(Modifs)) {
    if ("UniMod" %in% colnames(Modifs)) {
      UniMod <- unimod::modifications
      Modifs$"Mass shift" <- UniMod$MonoMass[match(Modifs$UniMod, UniMod$UnimodId)]
    } else {
      if (SearchSoft == "MAXQUANT") {
        modFls <- paste0(MQFold, "/bin/conf/modifications", c("", ".local"), ".xml")
        modFls <- modFls[which(file.exists(modFls))]
        modFls <- lapply(modFls, function(modFl) { #modFl <- modFls[1]
          xml_lst <- as_list(xml2::read_xml(modFl))
          xml_lst <- xml_lst[[1]]
          xml_lst <- as.data.frame(t(sapply(xml_lst, function(x) {
            #x <- xml_lst[[1]]
            return(c(attr(x, "title"), attr(x, "composition")))
          })))
        })
        modFls <- plyr::rbind.fill(modFls)
        colnames(modFls) <- c("Name", "Composition")
        Modifs$Composition <- modFls$Composition[match(Modifs$`Full name`, modFls$Name)]
        Modifs$"Mass shift" <- sapply(strsplit(Modifs$Composition, " "), function(x) {
          #x <- strsplit(modifs$Formula, " ")[1]
          x <- unlist(x)
          x <- as.data.frame(t(sapply(strsplit(gsub("\\)$", "", x), "\\("), function(y) {
            if (length(y) == 1) { y <- c(y, 1) }
            return(y)
          })))
          m <- match(x[[1]], proteoCraft::IsotopeProbs$Atom)
          stopifnot(sum(is.na(m)) == 0)
          # For now the code above throws an error if an elements is missing from the table
          # If it ever does, I should expand the table to add isotopic probabilities for more elements!!!
          x <- sum(as.numeric(gsub("_.+", "", proteoCraft::IsotopeProbs$Monoisotopic[m]))*as.integer(x[[2]]))
          return(x)
        })
      } else { stop() }
    }
  }
  # Provide dummy "Gene names" if the columns is missing
  if (!"Gene names" %in% colnames(Ev)) { Ev$"Gene names" <- paste0("Gene", 1:nrow(Ev)) }
  # Filters
  # 1 - We only want unique peptides!
  kol <- c("Modified sequence", "Charge")
  if (filter) { w1 <- grep(";", Ev$Proteins, invert = TRUE) } else { w1 <- 1:nrow(ev) }
  tempEv <- Ev[w1,]
  rownames(tempEv) <- paste0("id_", tempEv$id)
  tempEv$UniqID <- tmpPpId <- do.call(paste, c(tempEv[, kol], sep = "_;_"))
  # Unique
  # 2 - We also only want primary sequences existing in 1 single modification*charge state
  if (filter) {
    tst2 <- data.table::data.table(tmpPpId = tmpPpId, Seq = tempEv$Sequence)
    tst2 <- tst2[, list(x = length(unique(tmpPpId))), by = list(Group.1 = Seq)]
    tst2 <- as.data.frame(tst2)
    w2 <- which(tempEv$Sequence %in% tst2$Group.1[which(tst2$x == 1)])
    tempEv <- tempEv[w2,]
  }
  # 3 - We do not want any missed cleavages
  if (filter) {
    g <- grep("[RK].+", tempEv$Sequence)
    tmp2 <- tmp <- Ev$Sequence[g] # Sequences with at least one missed tryptic cleavage
    for (aa in c("R", "K")) { tmp2 <- gsub(aa, paste0(aa, "_"), tmp2) }
    tmp2 <- unlist(strsplit(gsub("_$", "", tmp2), "_"))
    tmp <- unique(c(tmp, tmp2))
    w3 <- which(!tempEv$Sequence %in% tmp) # New ev filter - includes the previous ones
    tempEv <- tempEv[w3,]
  }
  unq <- unique(tempEv$UniqID)
  tempEv2 <- data.frame("DUMMY" = 1:length(unq), "ID" = unq, check.names = FALSE)
  colnames(tempEv2)[1] <- ""
  kolz <- setNames(c("ModSeq", "Peptide", "Protein", "Gene", "QueryCharge"),
                   c("Modified sequence", "Sequence", "Proteins", "Gene names", "Charge"))
  m <- match(tempEv2$ID, tempEv$UniqID)
  tempEv2[, kolz] <- tempEv[m, names(kolz)]
  rownames(tempEv2) <- rownames(tempEv)[m]
  # 4 - We want to exclude peptides with no siblings
  if (filter) {
    tst <- aggregate(tempEv2$ID, list(tempEv2$Protein), function(x) { length(unique(x)) })
    tempEv2 <- tempEv2[which(tempEv2$Protein %in% tst$Group.1[which(tst$x > 1)]),]
  }
  #
  #grep(";", tempEv2$Gene)
  tempEv2$Gene <- gsub(";.+", "", tempEv2$Gene) # Just in case
  # Intensities
  Samples <- unique(experiments.map[[experiments.map_Sample.col]])
  tempEv3 <- data.table::data.table(tempEv[, c(Ev_Sample.col, "UniqID", intCol)])
  colnames(tempEv3)[which(colnames(tempEv3) == Ev_Sample.col)] <- "Sample"
  tempEv3 <- tempEv3[which(!is.na(tempEv3$Sample)),]
  tempEv3 <- dcast(tempEv3, UniqID ~ Sample, value.var = intCol, fun.aggregate = sum, na.rm = TRUE)
  tempEv3 <- as.data.frame(tempEv3)
  tempEv2[, Samples] <- tempEv3[match(tempEv2$ID, tempEv3$UniqID), Samples]
  rm(tempEv, tempEv3)
  # Sequence modifications
  tempEv2$PeptideSequenceModifications <- tempEv2$Peptide
  g <- grep("\\(", tempEv2$ID)
  tmpMS <- Modifs$`Mass shift`
  tmpMdSq <- gsub("|_;_[0-9]+$", "", tempEv2$ID[g])
  for (aa in AA) { tmpMdSq <- gsub(aa, paste0("_", aa, "_"), tmpMdSq) }
  tmpMdSq <- gsub("\\(|\\)", "", gsub("_+", "_", tmpMdSq))
  tmpMdSq <- strsplit(paste0(tmpMdSq, "_"), "_")
  tmpMdSq <- sapply(tmpMdSq, function(x) {
    x <- unlist(x)
    w1 <- which(x %in% c(AA, ""))
    w2 <- which(!x %in% c(AA, ""))
    ms <- Modifs$`Mass shift`[match(x[w2], Modifs$Mark)]
    res <- data.frame(AA = x, mod = "")
    w3 <- sapply(w2, function(y) { max(w1[which(w1 < y)]) })
    tst <- aggregate(ms, list(w3), sum)
    tst <- tst[which(tst$x > 0),]
    w4 <- which(tst$x > 0)
    tst$x <- as.character(round(tst$x, 3))
    tst$x[w4] <- paste0("+", tst$x[w4])
    res$mod[tst$Group.1] <- tst$x
    res <- res[which(!res$AA %in% Modifs$Mark),]
    if (nrow(res)) { res <- do.call(paste, c(res, sep = "")) }
    res <- paste(res, collapse = "")
    return(res)
  })
  tempEv2$PeptideSequenceModifications[g] <- tmpMdSq
  tempEv2[, paste0("Charge ", as.character(1:6))] <- 0
  w <- which(tempEv2[, paste0("Charge ", as.character(1:6))] == 0, arr.ind = TRUE)
  w <- as.matrix(data.frame(row = 1:nrow(tempEv2), col = tempEv2$QueryCharge))
  tempEv2[, paste0("Charge ", as.character(1:6))][w] <- 1
  tempEv2 <- tempEv2[, c("Peptide", "PeptideSequenceModifications", "Protein", "Gene", "QueryCharge",
                       "Charge 1", "Charge 2", "Charge 3", "Charge 4", "Charge 5", "Charge 6",
                       Samples)]
  Res <- list(Data = tempEv2, PTMs = Modifs)
  if (nchar(path)) {
    kol <- colnames(tempEv2)
    tempEv <- as.matrix(tempEv2)
    tempEv <- cbind(1:nrow(tempEv), tempEv)
    colnames(tempEv) <- c("", kol)
    try(write.table(tempEv, path, na = "", quote = FALSE, sep = "\t"), silent = TRUE)
  }
  return(Res)
}

