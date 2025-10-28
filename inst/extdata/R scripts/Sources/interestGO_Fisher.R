#################################################
#  Fisher extact test for GO terms of interest  #
#################################################
#
GO.tabs <- unlist(strsplit(Param$GO.tabs, ";"))
if (length(GO.tabs)) {
  # Reference
  dbGO <- listMelt(strsplit(db$`GO-ID`, ";"), db$`Protein ID`, c("GO", "Protein")) # Background proteome
  #
  # Columns used to determine whether a protein is found in a specific sample
  kol <- paste0(pep.ref["Original"], Exp.map[[RSA$column]])
  kol <- aggregate(kol, list(Exp.map[[VPAL$column]]), list)
  kol <- setNames(kol$x, kol$Group.1)
  #
  # Map pep rows to proteins and GO terms
  tmpPrt <- listMelt(strsplit(pep$Proteins, ";"), 1:nrow(pep), c("Protein", "Row"))
  tmpPrt$GO <- db$`GO-ID`[match(tmpPrt$Protein, db$`Protein ID`)]
  tmpGO <- listMelt(strsplit(tmpPrt$GO, ";"), tmpPrt$Row, c("GO", "Row"))
  #
  # Output list
  fishrTsts <- list()
  for (go in GO.tabs) { #go <- GO.tabs[1]
    fishrTsts[[go]] <- c()
    #
    goRws <- unique(tmpGO$Row[which(tmpGO$GO == go)]) # pep rows with annotations for the current term
    tstFish <- lapply(names(kol), function(grp) { #grp <- names(kol)[1]
      tmp <- pep[, kol[[grp]]] > 0
      w <- which(rowSums(tmp) >= 2)
      w1 <- w[which(w %in% goRws)] # Row of pep with non missing quant values for this group AND the current term
      w2 <- w[which(!w %in% goRws)] # Row of pep with non missing quant values for this group BUT not the current term
      p1 <- unique(tmpPrt$Protein[which(tmpPrt$Row %in% w1)])
      p2 <- unique(tmpPrt$Protein[which(tmpPrt$Row %in% w2)])
      n1 <- length(p1)
      n2 <- length(p2)
      p3 <- unique(dbGO$Protein[which((dbGO$GO == go)&(!dbGO$Protein %in% c(p1, p2)))])
      p4 <- unique(dbGO$Protein[which(!dbGO$Protein %in% c(p1, p2, p3))])
      n3 <- length(p3)
      n4 <- length(p4)
      dat <- data.frame("In group" = c(n1, n2),
                        "Not in group" = c(n3, n4),
                        row.names = c("with GO term", "without GO term"),
                        check.names = FALSE)
      fisher.test(dat)
    })
    fishrTsts[[go]] <- vapply(tstFish, function(x) { x$p.value }, 1)
  }
  fishrTsts <- do.call(rbind, fishrTsts)
  rownames(fishrTsts) <- GO_terms$Term[match( rownames(fishrTsts), GO_terms$ID)]
  colnames(fishrTsts) <- paste0("P-value - ", cleanNms(names(kol)))
  write.csv(fishrTsts, paste0(wd, "/Reg. analysis/GO enrich/Dataset/GO_terms_of_interest.csv"))
}




