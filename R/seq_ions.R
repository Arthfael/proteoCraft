#' seq_ions
#'
#' @description 
#' This function is not used, more of an unfinished exercise.
#' It should calculate the selected sequence ions for the peptide sequence, at all allowed charges.
#' It assumes a monoisotopic precursor.
#' This function has the following shortcomings - at least:
#' - It cannot accept TMT/iTRAQ modifications, or SILAC mass shifts.
#' - It cannot yet calculate internal cleavage, neutral loss or immonium fragments.
#' - It cannot predict fragment relative intensities... but that's kinda expected: who can? EDIT: well ML can now.
#' 
#' @param Seq The peptide sequence. Modifications are now supported.
#' @param Z The parent peptide's charge.
#' @param types The type of sequence ions to compute. Available are "a", "b", "c", "x", "y", "z".
#' @param C_alkylation One of "NEM" (default), "IAA" or "custom". For the latter case, a mass shift should be provided.
#' @param C_mass_shift Only required if "C_alkylation" is "Custom", other wise will be ignored.
#' @param Mods A data.frame with two columns: "Code" should be the 2-letter code of the modifications as reported by MaxQuant in the "Modified.sequence" column. "Shift" should be the associated expected mass shift.
#' @param plot Should we print a plot of expected peak positions? FALSE by default.
#' 
#' @examples
#' Result <- seq_ions(Seq = "GLSDGEWQQVLNVWGK", Z = 3, types = c("b", "y"), plot = TRUE)
#'

seq_ions <- function(Seq, Z, types = c("b", "y"), C_alkylation = "NEM", C_mass_shift = NULL,
                     Mods, plot = FALSE) {
  types <- tolower(types)
  Seq <- gsub("^_|_$", "", Seq)
  test <- gsub(paste(proteoCraft::AA, collapse = "|"), "", Seq)
  if (!missing("Mods")) {
    colnames(Mods) <- tolower(colnames(Mods))
    stopifnot(sum(sort(colnames(Mods)) == c("code", "shift")) == 2)
  }
  if (nchar(test) > 0) {
    test2 <- gsub("\\([a-z,-]{2}\\)", "", test)
    if (nchar(test2) == 0) {
      if (missing("Mods")) {
        Mods <- data.frame(code = c("ox", "ac", "de", "gl", "ph"),
                           shift = c(15.9949146221,
                                     42.0105646863,
                                     0.9840155848,
                                     -17.0265491015,
                                     79.9663304084))
        temp <- paste(paste0("\\(", Mods$Code, "\\)"), collapse = "|")
        test3 <- gsub(temp, "", test)
        if (nchar(test3) == 0) {
          warning("The sequence looks modified but no Mods table was provided so default values will be used.")
        } else {
          stop("The sequence is modified with PTMs for which no default is available, and argument \"Mods\" is missing!")
        }
      }
    } else {
      test2 <- unique(unlist(strsplit(test2, "")))
      stop(paste("The sequence contains unallowed characters:", test2, collapse = ""))
    }
  }
  stopifnot(class(Seq) == "character", sum(!types %in% c("a", "b", "c", "x", "y", "z")) == 0)
  if (!toupper(C_alkylation) %in% c("IAA", "IODOACETAMIDE", "NEM", "NETHYLMALEIMIDE", "N-ETHYLMALEIMIDE", "CUSTOM")) {
    stop("I don't know this alkylator! Allowed values are \"NEM\" and \"IAA\".")
  }
  if (toupper(C_alkylation) %in% c("IAA", "IODOACETAMIDE")) { C_mass_shift <- 57.021464 }
  if (toupper(C_alkylation) %in% c("NEM", "NETHYLMALEIMIDE", "N-ETHYLMALEIMIDE")) { C_mass_shift <- 125.047678 }
  if (toupper(C_alkylation) == "CUSTOM") {
    if (is.null(C_mass_shift)) {
      stop("\"C_alkylation\" is set to \"custom\" but no value was provided for the \"C_mass_shift\" argument!")
    }
  }
  ########################################################################################################################
  # Update the code below once I have variable modifications done: the sequence test is not correct
  ########################################################################################################################
  test <- gsub(paste(proteoCraft::AA, collapse = "|"), "", gsub("\\([a-z,-]{2}\\)", "", Seq))
  stopifnot(nchar(test) == 0)
  seq <- unlist(strsplit(Seq, split = ""))
  w1 <- which(seq %in% proteoCraft::AA)
  w2 <- which(seq == "(")
  seq_tbl <- data.frame(AA = seq,
                        Mod = "")
  seq_tbl$AA[unlist(sapply(w2, function(x) { x:(x+3) }))] <- ""
  if (length(w2) > 0) {
    seq_tbl$Mod[w2-1] <- sapply(w2, function(x) { paste(seq[(1:2)+x], collapse = "") })
  }
  seq_tbl <- seq_tbl[which(seq_tbl$AA != ""),]
  L <- nrow(seq_tbl)
  seq_tbl$Raw.Mass <- proteoCraft::AA_table$`Residue monoisotopic mass`[match(seq_tbl$AA,
                                                                        proteoCraft::AA_table$AA == x)]
  seq_tbl$Mod.Mass <- seq_tbl$Raw.Mass
  w <- which(seq_tbl$AA == "C")
  seq_tbl$Mod.Mass[w] <- seq_tbl$Raw.Mass[w] + C_mass_shift
  if (length(w2) > 0) {
    for (i in 1:nrow(Mods)) {
      w <- which(seq_tbl$Mod == Mods$code[i])
      seq_tbl$Mod.Mass[w] <- seq_tbl$Mod.Mass[w] + Mods$shift[i]
    }
  }
  RES <- sapply(types, function(t) {
    # t <- "b"
    # First Z = 1
    if (t %in% c("a", "b", "c")) {
      r <- sapply(1:L, function(x) {
        sum(seq_tbl$Mod.Mass[1:x]) + OrgMassSpecR::MonoisotopicMass(formula = list(H=1))
      })
      r <- r + c(-OrgMassSpecR::MonoisotopicMass(formula = list(C=1, O=1)),
                 0,
                 OrgMassSpecR::MonoisotopicMass(formula = list(N=1, H=3)))[which(c("a", "b", "c") == t)]
    }
    if (t %in% c("x", "y", "z")) {
      r <- sapply(1:L, function(x) {
        sum(seq_tbl$Mod.Mass[x:L]) + OrgMassSpecR::MonoisotopicMass(formula = list(O=1, H=1))
      })
      r <- r + c(OrgMassSpecR::MonoisotopicMass(formula = list(C=1, O=1)),
                 OrgMassSpecR::MonoisotopicMass(formula = list(H=2)),
                 -OrgMassSpecR::MonoisotopicMass(formula = list(N=1, H=1)))[which(c("x", "y", "z") == t)]
      r <- rev(r)
    }
    res <- r
    # Now other charge states:
    if (Z > 1) {
      for (i in 2:Z) {
        s <- (r + OrgMassSpecR::MonoisotopicMass(formula = list(H=i-1)))/i
        res <- c(res, s)
      }
    }
    names(res) <- sapply(c(1:Z), function(x) { paste0(c(1:length(r)), "_Z=", x) })
    return(res)
  })
  RES <- as.data.frame(RES)
  if (plot) {
    scale <- max(unlist(RES))
    col <- rainbow(n = length(types))
    temp <- RES
    temp$Z <- gsub("^[0-9]+_Z=", "", row.names(temp))
    temp <- reshape2::melt(temp)
    #temp$Label <- apply(temp[, c("variable", "Z")], 1, paste, collapse = "_")
    temp$Label <- paste0(round(temp$value, 3), "...") 
    n <- length(unique(temp$variable))
    temp$Offset <- (c(0.5:(n-0.5))[sapply(temp$variable, function(x) {
      which(types == x)
    })])/(n)
    plot <- ggplot2::ggplot(temp) +
      ggplot2::geom_vline(ggplot2::aes(xintercept = value, colour = variable)) +
      ggplot2::geom_text(ggplot2::aes(label = Label, x = value, y = Offset, colour = variable),
                         cex = 2.5, angle = 45) +
      ggplot2::ylim(0, 1) + ggplot2::xlim(0, scale*1.05) + ggplot2::theme_bw() +
      ggplot2::facet_wrap(~Z)
    proteoCraft::poplot(plot)
  }
  return(RES)
}
