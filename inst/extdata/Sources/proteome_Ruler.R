if (protrul) {
  ProtRulRoot %<o% "log10(est. copies/cell) - "
  tempPG <- PG
  exprsRt <- paste0("Mean ", prtRfRoot)
  if (LocAnalysis) {
    tempPG <- tempPG[, grep(topattern(exprsRt), colnames(tempPG), invert = TRUE)]
    for (grp2 in SubCellFracAggr2$values) { #grp2 <- SubCellFracAggr2$values[1L]
      em2 <- Exp.map[which(Exp.map[[SubCellFracAggr2$column]] == grp2),]
      for (grp in unique(em2[[SubCellFracAggr$column]])) { #grp <- unique(em2[[SubCellFracAggr$column]])[1L]
        em <- em2[which(em2[[SubCellFracAggr$column]] == grp),]
        kol <- paste0(prtRfRoot, em$Ref.Sample.Aggregate)
        tempPG[[paste0(prtRfRoot, grp)]] <- apply(10^tempPG[, kol], 1L, \(x) {
          log10(sum(x[which(is.finite(x))]))
        })
      }
      tempPG[[paste0(exprsRt, grp2)]] <- apply(tempPG[, paste0(prtRfRoot, unique(em2[[SubCellFracAggr$column]]))], 1L, \(x) {
        mean(x[which(is.finite(x))])
      })
    }
  }
  temp <- try(Prot.Ruler(tempPG, db, exprsRt, NuclL = ProtRulNuclL), silent = TRUE)
  if (inherits(temp, "list")) {
    db <- temp$Database
    temp <- temp$Protein.groups
    kol <- c(grep(topattern(exprsRt), colnames(temp), value = TRUE), grep(topattern(ProtRulRoot), colnames(temp), value = TRUE))
    PG[, kol] <- temp[, kol]
    protrul <- TRUE
    l <- length(DatAnalysisTxt)
    DatAnalysisTxt[l] <- paste0(DatAnalysisTxt[l],
                                " Protein group copy numbers per cell were estimated using a variant of the proteome ruler logic, normalizing to scaled values of all identified histones.")
  } else {
    warning("Failed to run Prot.Ruler function; is the remote NCBI server available?")
    protrul <- FALSE
  }
  rm(temp)
}
