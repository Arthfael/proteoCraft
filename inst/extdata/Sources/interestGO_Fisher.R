#################################################
#  Fisher exact test for GO terms of interest  #
#################################################
#
if (!exists("myFisherDir")) { myFisherDir <- paste0(wd, "/Reg. analysis/GO enrich/Dataset") }
if (!dir.exists(myFisherDir)) { dir.create(myFisherDir, recursive = TRUE) }
if (scrptType == "withReps") {
  GO.tabs <- unlist(strsplit(Param$GO.tabs, ";"))
}
if (scrptType == "noReps") {
  GO.tabs <- GO_filter
}
if (length(GO.tabs)) {
  # Reference from background proteome
  dbGO <- listMelt(strsplit(db$`GO-ID`, ";"),
                   db$`Protein ID`, c("GO", "Protein"))
  #
  # Columns used to determine whether a protein is found in a specific sample
  if (scrptType == "withReps") {
    kol <- paste0(pep.ref["Original"], Exp.map[[RSA$column]])
    kol <- aggregate(kol, list(Exp.map[[VPAL$column]]), list)
    kol <- setNames(kol$x, kol$Group.1)
  }
  if (scrptType == "noReps") {
    kol <- setNames(lapply(Exp, \(exp) { paste0(int.cols["Original"], " - ", exp) }), Exp)
  }
  #
  # Map pep rows to proteins and GO terms
  tmpPrt <- listMelt(strsplit(pep$Proteins, ";"), 1L:nrow(pep), c("Protein", "Row"))
  tmpPrt$GO <- db$`GO-ID`[match(tmpPrt$Protein, db$`Protein ID`)]
  tmpGO <- listMelt(strsplit(tmpPrt$GO, ";"), tmpPrt$Row, c("GO", "Row"))
  #
  #
  # Consider switching to using topGO - I have a similar 
  MinCount <- 10
  FishTst <- new("classicCount",
                 testStatistic = GOFisherTest,
                 name = "Fisher test")
  #
  # Output list
  fishrTsts <- list()
  for (go in GO.tabs) { #go <- GO.tabs[1L]
    fishrTsts[[go]] <- c()
    go_and_offspring <- union(go, GO_terms$Offspring[[match(go, GO_terms$ID)]])
    #
    goRws <- unique(tmpGO$Row[which(tmpGO$GO %in% go_and_offspring)]) # pep rows with annotations for the current term
    tstFish <- lapply(names(kol), \(grp) { #grp <- names(kol)[1L]
      tmp <- pep[, kol[[grp]], drop = FALSE] > 0
      w <- which(rowSums(tmp) >= 2L)
      w1 <- w[which(w %in% goRws)] # Row of pep with non missing quant values for this group AND the current term
      w2 <- w[which(!w %in% goRws)] # Row of pep with non missing quant values for this group BUT not the current term
      p1 <- unique(tmpPrt$Protein[which(tmpPrt$Row %in% w1)])
      p2 <- unique(tmpPrt$Protein[which(tmpPrt$Row %in% w2)])
      n1 <- length(p1)
      n2 <- length(p2)
      p3 <- unique(dbGO$Protein[which((dbGO$GO %in% go_and_offspring)&(!dbGO$Protein %in% c(p1, p2)))])
      p4 <- unique(dbGO$Protein[which(!dbGO$Protein %in% c(p1, p2, p3))])
      n3 <- length(p3)
      n4 <- length(p4)
      dat <- data.frame("In group" = c(n1, n2),
                        "Not in group" = c(n3, n4),
                        row.names = c("with GO term", "without GO term"),
                        check.names = FALSE)
      fisher.test(dat)
    })
    fishrTsts[[go]] <- vapply(tstFish, \(x) { x$p.value }, 1)
  }
  fishrTsts <- do.call(rbind, fishrTsts)
  w <- which(rownames(fishrTsts) %in% GO_terms$ID)
  rownames(fishrTsts)[w] <- GO_terms$Term[match(rownames(fishrTsts)[w], GO_terms$ID)]
  colnames(fishrTsts) <- paste0("P-value - ", cleanNms(names(kol)))
  fishrTsts <- fishrTsts[order(fishrTsts[[1L]], decreasing = FALSE), , drop = FALSE]
  #
  #
  Mappings <- data.table(dbGO)
  Mappings <- Mappings[, list(GO = list(unique(GO))), by = list(Protein = Protein)]
  Mappings <- as.data.frame(Mappings)
  Mappings <- setNames(Mappings$GO, Mappings$Protein)
  allProt <- setNames(seq_along(Mappings), names(Mappings))
  #
  IDs_flt <- setNames(lapply(names(kol), \(nm) { #nm <- names(kol)[1L]
    k <- kol[[nm]]
    w <- which(rowSums(pep[, k, drop = FALSE], na.rm = TRUE) > 0)
    unique(unlist(strsplit(pep$Proteins[w], ";")))
  }), names(kol))
  #
  #
  myGOdata <- setNames(lapply(Ontologies, \(ont) { #ont <- Ontologies[1L]
    my_topGO_tst <- setNames(lapply(names(kol), \(nm) { #nm <- names(kol)[1L]
      myFlt <- allProt[IDs_flt[[nm]]]
      require(AnnotationDbi)
      require(topGO)
      GOdata <- new("topGOdata",
          description = paste0(ont, " - ", nm),
          ontology = ont,
          allGenes = allProt,
          geneSel = \(x) { x %in% myFlt }, # Apply filter
          nodeSize = MinCount,
          annot = topGO::annFUN.gene2GO,
          gene2GO = Mappings)
      resultFisher <- topGO::getSigGroups(GOdata, FishTst)
      if (resultFisher@geneData[["Significant"]] <= 0) {
        return()
      }
      GO_tbl <- data.frame(ID = names(resultFisher@score),
                           Pvalue = resultFisher@score)
      return(GO_tbl)
    }), names(kol))
  }), Ontologies)
  myGOdata2 <- setNames(lapply(Ontologies, \(ont) {
    my_topGO_tst <- setNames(lapply(names(kol), \(nm) { #nm <- names(kol)[1L]
      x <- myGOdata[[ont]][[nm]]
      x <- x[which(x$ID %in% GO.tabs),]
      myGOdata[[ont]][[nm]] <- x
    }), names(kol))
  }), Ontologies)
  myGOdata2 <- setNames(lapply(Ontologies, \(ont) {
    pval <- lapply(names(kol), \(nm) { #nm <- names(kol)[1L]
      myGOdata[[ont]][[nm]]$Pvalue
    })
    pval <- do.call(cbind, pval)
    colnames(pval) <- names(kol)
    trm <- myGOdata[[ont]][[names(kol)[1L]]]$ID
    pval <- as.data.frame(pval)
    row.names(pval) <- trm
    return(pval)
  }), Ontologies)
  myGOdata2 <- do.call(rbind, myGOdata2)
  rownames(myGOdata2) <- sub("^[A-Z]{2}\\.", "", rownames(myGOdata2))
  w <- which(rownames(myGOdata2) %in% GO_terms$ID)
  rownames(myGOdata2)[w] <- GO_terms$Term[match(rownames(myGOdata2)[w], GO_terms$ID)]
  myGOdata2 <- myGOdata2[order(myGOdata2[[1L]], decreasing = FALSE), , drop = FALSE]
  #
  write.csv(fishrTsts, paste0(myFisherDir, "/GO_terms_of_interest.csv"))
  write.csv(myGOdata2, paste0(myFisherDir, "/GO_terms_of_interest_(topGO).csv"))
}
