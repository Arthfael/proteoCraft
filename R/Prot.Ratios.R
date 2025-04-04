#' Prot.Ratios
#'
#' @description
#' A function to calculate protein SILAC Ratios from individual peptides.
#' NB: Use with caution, this is an old function which can be considered deprecated.
#' Anyway for SILAC ratios the safest way is to accept MaxQuant ratios as more accurate than the ratio of intensity values.
#' Candidate for deletion!
#' 
#' @param Prot The protein/protein groups file.
#' @param Peptide.IDs The name of the Protein/Protein groups file's peptide IDs column. Default = "Peptide IDs"
#' @param Pep The peptides file.
#' @param id The name of the Peptides file IDs column. Default = "id"
#' @param Pep.Ratios.Nms Name of the peptides ratios column(s)
#' @param Pep.Ratios.to.log Should the input peptide ratios be converted to log scale (log2)? If they aren't, then they must, as the function will assume they are log scale. Default = TRUE
#' @param Prot.Ratios.to.log Should the output protein ratios be log-scale or not? Default = TRUE
#' @param Mods Which modifications should be included? If set to FALSE, will not filter any modifications.
#' @param Mod.Nms Default = "Modified sequence". The name of the column containing the Modified sequence in the peptides file. Can be set to a non-modified sequence if Mods = FALSE
#' @param Min.Pep.Nb How many peptides should at least be present? Should be at the very least 1, but really we don't advise going below 2 (set do default).
#' @param SD Should the Output also include the standard deviation? Default = TRUE
#' @param Pvalue Should the Output also include the P-value? Default = TRUE
#' @param In.Norm  Should the Input be re-normalised? Default = FALSE
#' @param Out.Norm Should the Output be re-normalised? Default = TRUE
#' 
#' @examples
#' test <- Prot.Ratios(Prot = prot, Peptide.IDs = "New.Peptide.IDs", Pep = pep, id = "New.Peptide.ID",
#'                     Pep.Ratios.Nms = paste0("log2.", sapply(Exp, function(x) {paste0("norm.Ratio.", x, ".", unique(unlist(strsplit(Ratios.map$Ratios, split = ";"))))})),
#'                     Pep.Ratios.to.log = FALSE, Prot.Ratios.to.log = TRUE,
#'                     Mods = Mod4Quant, Mod.Nms = "Modified.sequence", Min.Pep.Nb = 2, SD = TRUE, Pvalue = TRUE, In.Norm = TRUE, Out.Norm = TRUE)
#' colnames(test) <- gsub("Protein_", "", colnames(test))
#' prot[, colnames(test)] <- test
#' 
#' @export

Prot.Ratios <- function(Prot, Peptide.IDs = "Peptide IDs", Pep, id = "id",
                        Pep.Ratios.Nms, Pep.Ratios.to.log = TRUE, Prot.Ratios.to.log = TRUE,
                        Mods, Mod.Nms = "Modified sequence", Min.Pep.Nb = 2, SD = TRUE, Pvalue = TRUE, In.Norm = FALSE,
                        Out.Norm = TRUE) {
  test <- rep(0, nrow(Pep))
  if (!is.logical(Mods)) { # A clumsy way to do it, I know...
    a <- paste(c(paste(c(AA, "_"), collapse = "|"), paste0("\\(", Mods, "\\)")), collapse = "|")
    test <- sapply(gsub(a, "", Pep[[Mod.Nms]]), nchar)
  }
  pep1 <- Pep[which(test == 0),]
  options(warn = -1)
  if (Pep.Ratios.to.log) {
    for (i in Pep.Ratios.Nms) {
      pep1[[i]] <- log2(pep1[[i]])
    }
  }
  options(warn = 0)
  if (In.Norm) {
    for (i in Pep.Ratios.Nms) {
      a <- median(proteoCraft::is.all.good(pep1[[i]]))
      pep1[[i]] <- pep1[[i]] - a
    }
  }
  temp.ids <- strsplit(Prot[[Peptide.IDs]], split = ";")
  # NB: I have checked this code, and yes, if SD = TRUE, the first half of the output columns is the ratios and the second half is the SDs.
  # Pfew, I had a horrible doubt!
  res <- as.data.frame(t(sapply(temp.ids, function(x) {
    x <- which(pep1[[id]] %in% unlist(x))
    if (length(Pep.Ratios.Nms) > 1) {
      temp <- apply(pep1[x, Pep.Ratios.Nms], 2, function(y) {
        # Below the number of peptides has to be at least 2:
        if (length(which(!is.na(y))) >= Min.Pep.Nb) {
          y <- proteoCraft::is.all.good(y)
          if (length(y) > 0) {
            if (length(y) > 1) {
              if (SD) { y <- c(median(y), sd(y)) } else { y <- median(y) }
            }
          } else { y <- NA }
        } else { y <- NA }
        if ((SD)&&(length(y) == 1)) { y <- c(y, NA) }
        return(y)
      })
      if (SD) { temp <- c(temp[1,], temp[2,]) }
    } else {
      y <- pep1[x, Pep.Ratios.Nms]
      # Below the number of peptides has to be at least 2:
      if (length(which(!is.na(y))) >= Min.Pep.Nb) {
        y <- proteoCraft::is.all.good(y)
        if (length(y) > 0) {
          if (length(y) > 1) {
            if (SD) { y <- c(median(y), sd(y)) } else { y <- median(y) }
          }
        } else { y <- NA }
      } else { y <- NA }
      if ((SD)&&(length(y) == 1)) { y <- c(y, NA) }
      temp <- y
    }
    return(temp)
  })))
  if (SD) {
    colnames(res) <- c(paste0("Protein_", Pep.Ratios.Nms),
                       paste0("Protein_", Pep.Ratios.Nms, ": SD"))
  } else { colnames(res) <- paste0("Protein_", Pep.Ratios.Nms) }
  if (Out.Norm) {
    for (i in paste0("Protein_", Pep.Ratios.Nms)) {
      a <- proteoCraft::is.all.good(as.numeric(res[[i]]))
      if (length(a) > 0) {
        a <- median(a)
        res[[i]] <- res[[i]] - a
      }      
    }
  }
  if (!Prot.Ratios.to.log) {
    for (i in paste0("Protein_", Pep.Ratios.Nms)) {
      res[[i]] <- 2^res[[i]]
    }
  }
  if (Pvalue) {
    res2 <- sapply(temp.ids, function(x) {
      x <- which(pep1[[id]] %in% unlist(x))
      if (length(Pep.Ratios.Nms) > 1) {
        temp2 <- apply(pep1[x, Pep.Ratios.Nms], 2, function(y) {
          # Below the number of peptides has to be at least 2:
          if (length(which(!is.na(y))) >= Min.Pep.Nb) {
            y <- proteoCraft::is.all.good(y)
            if (length(y) > 1) {
              y <- -log10(t.test(x = y, y = NULL, alternative = "two.sided")$p.value)
            } else { y <- NA }
          } else { y <- NA }
          return(y)
        })
      } else {
        y <- pep1[x, Pep.Ratios.Nms]
        # Below the number of peptides has to be at least 2:
        if (length(which(!is.na(y))) >= Min.Pep.Nb) {
          y <- proteoCraft::is.all.good(y)
          if (length(y) > 1) {
            y <- -log10(t.test(x = y, y = NULL, alternative = "two.sided")$p.value)
          } else { y <- NA }
        } else { y <- NA }
        temp2 <- y
      }
      return(temp2)
    })
    if (length(Pep.Ratios.Nms) > 1) { res2 <- as.data.frame(t(res2)) }
    if (length(Pep.Ratios.Nms) > 1) {
      colnames(res2) <- paste0("Protein_", Pep.Ratios.Nms, ": -log10(peptides Pvalue)")
      res <- cbind(res, res2)
    } else { res[, paste0("Protein_", Pep.Ratios.Nms, ": SD")] <- res2 }
  }
  return(res)
}
