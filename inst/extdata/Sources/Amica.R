#### Amica input tables
# Write tables for Amica input:
## PG table
if (Param$Amica) {
  amicaTst <- try({
    cat("Also writing Amica table...\n")
    w <- which(as.logical(Exp.map$Use))
    AmicaDesign <- data.frame(groups = cleanNms(Exp.map[w, VPAL$column], rep = "."),
                              samples = cleanNms(Exp.map[w, RSA$column], rep = "."))
    AmicTbl <- data.frame(Majority.protein.IDs = PG$"Leading protein IDs",
                          Gene.names = PG$Genes,
                          razorUniqueCount = PG$"Razor + unique peptides",
                          Potential.contaminant = PG$"Potential contaminant")
    tmp <- data.frame(IDs = PG$"Peptide IDs", Razor = PG$"Peptide is razor")
    tmp$Razor <- lapply(strsplit(tmp$Razor, ";"), \(x) {
      as.logical(toupper(x))
    })
    tmp$IDs <- lapply(strsplit(tmp$IDs, ";"), as.numeric)
    tmp$RazorIDs <- apply(tmp[, c("IDs", "Razor")], 1L, \(x) {
      x[[1L]][which(x[[2L]])]
    })
    for (i in Exp.map$Ref.Sample.Aggregate[which(as.logical(Exp.map$Use))]) { #i <- Exp.map$Ref.Sample.Aggregate[which(as.logical(Exp.map$Use))][1L]
      i2 <- cleanNms(i, rep = ".")
      kol <- paste0("LFQIntensity_", i2)
      AmicTbl[[kol]] <- PG[[paste0(prtRfRoot, i)]]/log10(2L)
      AmicTbl[which(!is.finite(AmicTbl[[kol]])), kol] <- NaN
      kol <- paste0("razorUniqueCount_", i2)
      tmp$Tmp <- strsplit(PG[[paste0("Peptide IDs - ", i)]], ";")
      AmicTbl[[kol]] <- apply(tmp[, c("IDs", "Tmp")], 1L, \(x) {
        sum(x[[2L]] %in% x[[1L]])
      })
      AmicTbl[which(!is.finite(AmicTbl[[kol]])), kol] <- 0L
    }
    tmp <- Exp.map$Ref.Sample.Aggregate[which(as.logical(Exp.map$Use))]
    kol <- paste0("LFQIntensity_", cleanNms(tmp, rep = "."))
    temp <- Data_Impute2(AmicTbl[, kol],
                         Exp.map[match(tmp, Exp.map$Ref.Sample.Aggregate), VPAL$column])
    temp <- temp$Imputed_data
    colnames(temp) <- gsub("^LFQIntensity_", "ImputedIntensity_", colnames(temp))
    AmicTbl[, colnames(temp)] <- temp
    for (contr in myContrasts$Contrast) { #contr <- myContrasts$Contrast[1L]
      nm <- gsub(" - ", "__vs__", contr)
      kol <- paste0("P.Value_", nm)
      AmicTbl[[kol]] <- 10L^(-PG[[paste0(pvalue.col[which(pvalue.use)], contr)]])
      AmicTbl[which(!is.finite(AmicTbl[[kol]])), kol] <- NaN
      kol <- paste0("adj.P.Val_", nm)
      PVkol <- paste0(pvalue.col[which(pvalue.use)], contr)
      AmicTbl[[kol]] <- p.adjust(10L^(-PG[[PVkol]]), method = "BH")
      AmicTbl[which(!is.finite(AmicTbl[[kol]])), kol] <- NaN
      kol <- paste0("logFC_", nm)
      AmicTbl[[kol]] <- PG[[paste0("Mean ", Prot.Rat.Root, contr)]]
      AmicTbl[which(!is.finite(AmicTbl[[kol]])), kol] <- NaN
      kol <- paste0("AveExpr_", nm)
      m <- match(contr, myContrasts$Contrast)
      AmicTbl[[kol]] <- PG[[paste0("Mean ", prtRfRoot, myContrasts$A_full[m])]]/log10(2L)
      AmicTbl[which(!is.finite(AmicTbl[[kol]])), kol] <- NaN
    }
    tst <- apply(AmicTbl[, grep("^AveExpr_", colnames(AmicTbl), value = TRUE), drop = FALSE], 1L, \(x) {
      sum(is.finite(x))
    }) > 0L
    AmicTbl$quantified <- c("", "+")[((AmicTbl$razorUniqueCount >= 2L)&(tst))+1L]
    AmicTbl <- AmicTbl[which(AmicTbl$quantified == "+"),]
    dir <- paste0(wd, "/Amica")
    if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
    dirlist <- unique(c(dirlist, dir))
    data.table::fwrite(AmicTbl, paste0(wd, "/Amica/Amica_file.csv"), row.names = FALSE, na = "NaN", sep = "\t", quote = FALSE)
    data.table::fwrite(AmicaDesign, paste0(wd, "/Amica/Experimental_design.csv"), row.names = FALSE, sep = "\t", quote = FALSE, na = "NA")
    #data.table::fwrite(AmicTbl[1L:500L,], paste0(wd, "/Amica/Amica_file_short.csv"), row.names = FALSE, na = "NaN", quote = FALSE)
  })
  if (inherits(amicaTst, "try-error")) {
    warning("Update this chunk to allow writing Amica tables again!")
  }
}
