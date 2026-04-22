#### Code chunk - Protein group profile plots and ranked abundance plots
#
# MAKE ME FASTER:
#   The last way I can make this faster (I think) is to delay saving the ggplots and plotlys,
#   and then parallelising the saving process.
#
#runRankAbundPlots %<o% TRUE
#runProfPlots %<o% TRUE
if (runRankAbundPlots||runProfPlots) {
  library(ggplot2)
  library(RColorBrewer)
  library(colorspace)
  library(ggplot2)
  library(plotly)
  library(AnnotationDbi)
  library(htmlwidgets)
  #
  # Change this hidden parameter to not trim subspecies tags
  #spLevel <- "subspecies"
  spLevel <- "species"
  #
  if (Annotate) {
    AllTerms %<o% unique(unlist(strsplit(db$`GO-ID`, ";")))
    AllTermNames %<o% unique(unlist(strsplit(db$GO, ";")))
  }
  #
  ggQuantLy %<o% list()
  ggProfLy %<o% list()
  #
  pepQuantTypes_ref %<o% setNames("PLACEHOLDER", "intensities")
  pepQuantTypes %<o% names(pepQuantTypes_ref)
  QuantTypes_ref %<o% setNames(c("PLACEHOLDER", "Sequence coverage [%] - ", "Spectral count - "),
                               c("LFQ", "Coverage", "Spectra"))
  if (!CreateMSMSKol) { QuantTypes_ref <- QuantTypes_ref[which(names(QuantTypes_ref) != "Spectra")] }
  QuantTypes %<o% names(QuantTypes_ref)
  lQ1 <- length(QuantTypes)
  lQ2 <- length(pepQuantTypes)
  #
  MainDir <- paste0(wd, "/Ranked abundance")
  MainDir2 <- paste0(wd, "/Profile plots")
  for (dr in c(MainDir, MainDir2)) {
    if (!dir.exists(dr)) { dir.create(dr, recursive = TRUE) }
  }
  if (scrptType == "withReps") { dirlist <- unique(c(dirlist, MainDir, MainDir2)) }
  for (quantType in QuantTypes) { #quantType <- QuantTypes[1L] #quantType <- "LFQ" #QuantType <- "Coverage"
    subDir <- paste0(MainDir, "/", quantType)
    if (scrptType == "withReps") { dirlist <- unique(c(dirlist, subDir)) }
    if (!dir.exists(subDir)) { dir.create(subDir, recursive = TRUE) }
  }
  #
  PG_varkol <- c("Leading protein IDs", "Protein IDs", "Common Name (short)", "id", "Label", "PEP")
  pep_varkol <- c("Modified sequence", "id", "Proteins", "PEP")
  if (prot.list.Cond) {
    PG_varkol <- union(PG_varkol, "In list")
    pep_varkol <- union(pep_varkol, "In list")
  }
  if ("GO-ID" %in% colnames(PG)) {
    PG_varkol <- union(PG_varkol, "GO-ID")
  }
  klstKol <- c()
  if ((exists("KlustKols"))&&(length(KlustKols))) {
    klstKol <- grep(" - Global$", KlustKols, value = TRUE)[1L]
    if (length(klstKol) == 1L) {
      PG_varkol <-  union(PG_varkol, klstKol)
    } 
  }
  tstOrg2 %<o% c()
  if (tstOrg) {
    tmp1 <- as.data.table(listMelt(setNames(strsplit(PG$`Protein IDs`, ";"), 1L:nrow(PG))))
    tmp2 <- as.data.table(listMelt(setNames(strsplit(pep$Proteins, ";"), 1L:nrow(pep))))
    tmp <- unique(c(tmp1$value, tmp2$value))
    tmp <- data.frame(ID = tmp,
                      Org = db[match(tmp, db$`Protein ID`),
                               dbOrgKol])
    tmp$Org <- gsub(" *\\([^\\)]*\\)", "", tmp$Org)
    if (spLevel == "species") { tmp$Org <- gsub(" subsp((\\.)|(ecies)) .*", "", tmp$Org) }
    tmp1$Org <- tmp$Org[match(tmp1$value, tmp$ID)]
    tmp2$Org <- tmp$Org[match(tmp2$value, tmp$ID)]
    f0 <- \(x) {
      x <- sort(unique(x))
      x <- c(x[which(x %in% Org$Organism)],
             x[which(!x %in% Org$Organism)])
      paste(x, collapse = ";")
    }
    tmp1 <- tmp1[, .(Org = f0(Org)), by = .(row = as.integer(L1))]
    tmp2 <- tmp2[, .(Org = f0(Org)), by = .(row = as.integer(L1))]
    PG$"Org. label" <- tmp1$Org[match(1L:nrow(PG), tmp1$row)]
    pep$"Org. label" <- tmp2$Org[match(1L:nrow(pep), tmp2$row)]
    PG$"Org. label"[which(PG$`Potential contaminant` == "+")] <- "Contaminant"
    pep$"Org. label"[which(pep$`Potential contaminant` == "+")] <- "Contaminant"
    if (prot.list.Cond) {
      PG$"Org. label"[which(PG$`In list` == "+")] <- "In list"
      pep$"Org. label"[which(pep$`In list` == "+")] <- "In list"
    }
    PG_varkol <- union(PG_varkol, "Org. label")
    pep_varkol <- union(pep_varkol, "Org. label")
    #
    tmp <- pep$"Org. label"
    tstOrg2 <- aggregate(tmp, list(tmp), length)
    tstOrg2 <- tstOrg2[order(tstOrg2$x, decreasing = TRUE),]
    tstOrg2 <- tstOrg2$Group.1[which(tstOrg2$x > 1L)]
    tstOrg2 <- tstOrg2[which(!tstOrg2 %in% c("Contaminant", "In list"))]
  }
  Org_Nms %<o% abbrOrg(tstOrg2)
  #
  myColors <- setNames("black", "-")
  myColors2 <- setNames(c("lightgrey", "brown"), c("-", "+"))
  lOrg <- length(tstOrg2)
  if (lOrg) {
    myColors <- c(myColors,
                  setNames(colorRampPalette(c("#47b7a4", "#eeee00"))(lOrg),
                           tstOrg2)) # We never expect more than a handful of organisms
    myColors[["Contaminant"]] <- "deepskyblue"
  }
  #myColorsB[[c("+", "In list")[(lOrg > 0)+1L]]] <- myColors[[c("+", "In list")[(lOrg > 0)+1L]]] <- "red"
  myColors["Specific"] <- "purple"
  myColors["In list"] <- "brown"
  colScale <- scale_colour_manual(name = "Category", values = myColors)
  #colScaleB <- scale_colour_manual(name = "Category", values = myColorsB)
  fillScale <- scale_fill_manual(name = "Category", values = myColors)
  colScale2 <- scale_colour_manual(name = "In list", values = myColors2)
  fillScale2 <- scale_fill_manual(name = "In list", values = myColors2)
  #
  if (scrptType == "withReps") {
    PG_ref <- paste0("Mean ", prtRfRoot)
    pep_ref <- paste0("Mean ", pep.ref[length(pep.ref)])
    Agg <- VPAL
    if ((WorkFlow == "PULLDOWN")||(sum(Exp.map$Use) < 12L)) {
      PG_ref <- prtRfRoot
      Agg <- RSA
    }
    mySamples <- Agg$values
    GO_PG_col %<o% unique(unlist(strsplit(Param$GO.tabs, ";")))
    GO_filt %<o% (length(GO_PG_col) > 0L)
    if (GO_filt) {
      if ((!exists("GO_terms"))&&(file.exists(paste0(wd, "/GO_terms.RData")))) { loadFun(paste0(wd, "/GO_terms.RData")) }
      GO_PG_col <- GO_PG_col[which(GO_PG_col %in% GO_terms$ID)]
      GO_filt <- length(GO_PG_col) > 0L
    }
    GO_filter <- GO_PG_col
  }
  if (scrptType == "noReps") {
    PG_ref <- rev(PG.int.cols[which(PG.int.cols != paste0("Imput. ", PG.int.cols["Original"]))])[1L]
    pep_ref <- paste0(int.cols["Original"], " - ")
    mySamples <- Exp
    if (GO_filt) {
      if ((!exists("GO_terms"))&&(file.exists(paste0(wd, "/GO_terms.RData")))) { loadFun(paste0(wd, "/GO_terms.RData")) }
    }
  }
  if (GO_filt) {
    names(GO_filter) <- GO_terms$Term[match(GO_filter, GO_terms$ID)]
    w <- which(is.na(names(GO_filter)))
    if (length(w)) {
      if (exists("CompGOTerms")) {
        w1 <- w[which(GO_filter[w] %in% CompGOTerms)]
        if (length(w1)) {
          names(GO_filter)[w1] <- names(CompGOTerms)[match(GO_filter[w1], CompGOTerms)]
          
        }
      }
      w <- which(is.na(names(GO_filter)))
      names(GO_filter)[w] <- GO_filter[w]
    }
    names(GO_filter) <- gsub(" *\\[.*", "", names(GO_filter))
  }
  QuantTypes_ref["LFQ"] <- PG_ref
  pepQuantTypes_ref["intensities"] <- pep_ref
  #
  runMark <- exists("CompGOTerms")
  #
  all_PG_kol <- unlist(lapply(QuantTypes_ref, \(ref) { paste0(ref, mySamples) }))
  all_pep_kol <- unlist(lapply(pepQuantTypes_ref, \(ref) { paste0(ref, mySamples) }))
  if (scrptType == "noReps") {
    # Currently left out for replicates script because there isn't any-more a 1:1 relationship between tests and sample groups!
    PG_regKol <- pep_regKol <- paste0("Regulated - ", mySamples)
    PG_regKol <- PG_regKol[which(PG_regKol %in% colnames(PG))]
    pep_regKol <- pep_regKol[which(pep_regKol %in% colnames(pep))]
  }
  PG_kntKol <- paste0("Peptides count - ", mySamples)
  PG_kntKol <- PG_kntKol[which(PG_kntKol %in% colnames(PG))]
  myPG <- PG[, c(all_PG_kol, PG_varkol, PG_kntKol)]
  if ((scrptType == "noReps")&&(length(PG_regKol))) {
    myPG[, PG_regKol] <- PG[, PG_regKol]
  }
  if (length(klstKol) == 1L) {
    colnames(myPG)[which(colnames(myPG) == klstKol)] <- "Cluster"
    PG_varkol[which(PG_varkol == klstKol)] <- "Cluster"
  }
  colnames(myPG)[which(colnames(myPG) == "Label")] <- "Protein Group"
  PG_varkol[which(PG_varkol == "Label")] <- "Protein Group"
  myPep <- pep[, c(all_pep_kol, pep_varkol)]
  if ((scrptType == "noReps")&&(length(pep_regKol))) {
    myPep[, pep_regKol] <- pep[, pep_regKol]
  }
  myPG$"Protein Group" <- factor(myPG$"Protein Group", levels = myPG$"Protein Group")
  myPep$"Modified sequence" <- gsub("^_|_$", "", myPep$"Modified sequence")
  myPep$"Modified sequence" <- factor(myPep$"Modified sequence", levels = myPep$"Modified sequence")
  if (prot.list.Cond) {
    myPG$`In list`[which(myPG$`In list` == "")] <- "-"
    myPep$`In list`[which(myPep$`In list` == "")] <- "-"
  }
  test <- aggregate(myPG$"Protein Group", list(myPG$"Protein Group"), length)
  w <- which(test$x > 1L)
  if (length(w)) {
    test <- test[w,]
    for (i in test$Group.1) {
      w <- which(myPG$"Protein Group" == i)
      myPG$"Protein Group"[w[2L:length(w)]] <- paste0(myPG$"Protein Group"[w[2L:length(w)]], "_", 2L:length(w))
    }
  }
  if (lOrg) {
    myPG$Category <- myPG$`Org. label`
    myPep$Category <- myPep$`Org. label`
  } else {
    myPG$Category <- "-"
    myPep$Category <- "-"
  }
  PG_varkol <- union(PG_varkol, "Category")
  pep_varkol <- union(pep_varkol, "Category")
  if (prot.list.Cond) {
    tst <- (lOrg > 0L)+1L
    myPG$Category[which(myPG$"In list" == "+")] <- c("+", "In list")[tst]
    myPep$Category[which(myPep$"In list" == "+")] <- c("+", "In list")[tst]
  }
  lev <- c("In list", "+", tstOrg2, "-", "Contaminant")
  lev1 <- lev[which(lev %in% unique(myPG$Category))]
  lev2 <- lev[which(lev %in% unique(myPep$Category))]
  myPG$Category <- factor(myPG$Category, levels = lev1)
  myPep$Category <- factor(myPep$Category, levels = lev2)
  #
  myFlt <- c()
  if (GO_filt||runMark) {
    if (GO_filt) { myFlt[names(GO_filter)] <- GO_filter }
    if (runMark) { myFlt[names(CompGOTerms)] <- CompGOTerms }
    tmp <- listMelt(strsplit(myPG$`GO-ID`, ";"), 1L:nrow(myPG), c("GO", "row"))
    for (goID in myFlt) { #goID <- GO_filter[1L]
      # Get children terms
      gofilter <- unique(unlist(c(goID,
                                  GOBPOFFSPRING[[goID]],
                                  GOCCOFFSPRING[[goID]],
                                  GOMFOFFSPRING[[goID]])))
      gofilter <- gofilter[which(!is.na(gofilter))]
      if (sum(gofilter %in% AllTerms)) {
        myPG[[goID]] <- "-"
        wtst <- unique(tmp$row[which(tmp$GO %in% gofilter)])
        myPG[wtst, goID] <- "+"
        PG_varkol <- union(PG_varkol, goID)
      } else {
        myFlt <- myFlt[which(myFlt != goID)] # Not setdiff, it strips names!!!
      }
    }
    if (length(myFlt)) {
      myGOcolors <- setNames(c("grey", rainbow(length(myFlt)-1L)), myFlt)
    }
  }
  # Custom colors
  myPG2 <- data.table(Category = myPG$Category,
                      PG = myPG$`Protein Group`)
  myPG2 <- myPG2[, list(x = unique(Category)),
                 by = list(Group.1 = PG)]
  myPG2 <- as.data.frame(myPG2)
  myPG2$Protein_group <- "darkgrey"
  w <- which(myPG2$x == "Contaminant")
  suppressWarnings(myPG2$Protein_group[w] <- brewer.pal(min(c(length(w), 12L)), "Blues"))
  w <- which(myPG2$x %in% tstOrg2)
  myPG2$Protein_group[w] <- rainbow_hcl(length(w))
  w <- which(myPG2$x %in% c("In list",  "+"))
  myPG2$Protein_group[w] <- "red"
  myPG$Protein_group <- myPG2$Protein_group[match(myPG$`Protein Group`, myPG2$Group.1)]
  PG_varkol <- union(PG_varkol, "Protein_group")
  #
  myPep2 <- data.table(Category = myPep$Category,
                       Pep = myPep$`Modified sequence`)
  myPep2 <- myPep2[, list(x = unique(Category)),
                   by = list(Group.1 = Pep)]
  myPep2 <- as.data.frame(myPep2)
  myPep2$Peptidoform <- "darkgrey"
  w <- which(myPep2$x == "Contaminant")
  suppressWarnings(myPep2$Peptidoform[w] <- brewer.pal(min(c(length(w), 12L)), "Blues"))
  w <- which(myPep2$x %in% tstOrg2)
  myPep2$Peptidoform[w] <- rainbow_hcl(length(w))
  w <- which(myPep2$x %in% c("In list",  "+"))
  myPep2$Peptidoform[w] <- "red"
  myPep$Peptidoform <- myPep2$Peptidoform[match(myPep$`Modified sequence`, myPep2$Group.1)]
  pep_varkol <- union(pep_varkol, "Peptidoform")
  #
  #
  # Define data frame of plots to draw in parallel:
  samplesList <- list(PG = sort(rep(mySamples, lQ1)),
                      pep = sort(rep(mySamples, lQ2)))
  samplesDF1 <- listMelt(samplesList, ColNames = c("values", "type"))
  samplesDF1$QuantType <- samplesDF1$ref <- ""
  wQ1 <- which(samplesDF1$type == "PG")
  samplesDF1$QuantType[wQ1] <- rep(QuantTypes, length(wQ1)/lQ1)
  samplesDF1$ref[wQ1] <- c(PG_ref,
                          "Sequence coverage [%] - ",
                          "Spectral count - ")[match(samplesDF1$QuantType[wQ1], QuantTypes)]
  wQ2 <- which(samplesDF1$type == "pep")
  samplesDF1$QuantType[wQ2] <- rep(pepQuantTypes, length(wQ2)/lQ2)
  samplesDF1$ref[wQ2] <- c(pep_ref)[match(samplesDF1$QuantType[wQ2], pepQuantTypes)]
  #
  PltTst <- setNames(c(TRUE, prot.list.Cond), c("All", "List"))
  # if (GO_filt) {
  #   tmp <- listMelt(strsplit(myPG$"GO-ID", ";"), 1L:nrow(myPG))
  #   for (goID in GO_filter) {
  #     w <- which(tmp$value == goID)
  #     PltTst[goID] <- length(w) > 0L
  #     if (PltTst[goID]) {
  #       myPG[[goID]] <- 1L:nrow(myPG) %in% tmp$L1[which(tmp$value == goID)]
  #       PG_varkol <- union(PG_varkol, goID)
  #     }
  #   }
  #   PG_varkol <- PG_varkol[which(PG_varkol != "GO-ID")]
  # }
  PltTst <- PltTst[which(PltTst)]
  samplesDF1 <- lapply(names(PltTst), \(x) {
    rs <- samplesDF1
    rs$subtype <- x
    return(rs)
  })
  samplesDF1 <- do.call(rbind, samplesDF1)
  #
  samplesDF_Lst <- list()
  if (runRankAbundPlots) {
    samplesDF_Lst$ranked <- samplesDF1
  }
  if ((runProfPlots)&&(moreThan1Exp)) {
    samplesDF2 <- aggregate(samplesDF1$values, list(samplesDF1$QuantType, samplesDF1$ref, samplesDF1$type), \(x){
      list(unique(x))
    })
    colnames(samplesDF2) <- c("QuantType", "ref", "type", "values")
    samplesDF2 <- samplesDF2[which(samplesDF2$type != "pep"),] # Those take just too bloody long!!!
    samplesDF2 <- samplesDF2[which(lengths(samplesDF2$values) > 1L),]
    #
    PltTst2 <- setNames(c(TRUE, prot.list.Cond, GO_filt, exists("CompGOTerms")),
                        c("All", "List", "GO", "Mark"))
    PltTst2 <- PltTst2[which(PltTst2)]
    samplesDF2 <- lapply(names(PltTst2), \(x) {
      rs <- samplesDF2
      rs$subtype <- x
      return(rs)
    })
    samplesDF2 <- do.call(rbind, samplesDF2)
    samplesDF_Lst$profiles <- samplesDF2
  }
  samplesDF <- do.call(rbind, samplesDF_Lst)
  #
  # Cluster operations
  source(parSrc)
  tmpFl1 <- tempfile(fileext = ".rds")
  tmpFl2 <- tempfile(fileext = ".rds")
  exports <- c("samplesDF", "tmpFl1", "tmpFl2", "wd", "MainDir", "MainDir2", "PG_varkol", "pep_varkol",
               "WorkFlow", "GO_filt", "runMark", "MakeRatios", "mySamples", "QuantTypes", "QuantTypes_ref", "pepQuantTypes", "pepQuantTypes_ref",
               "colScale", "fillScale", "colScale2", "fillScale2", "Exp", "scrptType", "tstOrg2")
  if (GO_filt) { exports <- union(exports, "GO_filter") }
  if (length(myFlt)) { exports <- union(exports, c("myFlt", "CompGOTerms", "myGOcolors")) }
  exports <- as.list(exports)
  clusterExport(parClust, exports, envir = environment())
  readr::write_rds(myPG, tmpFl1)
  readr::write_rds(myPep, tmpFl2)
  invisible(clusterCall(parClust, \() {
    library(ggplot2)
    library(RColorBrewer)
    library(colorspace)
    library(ggplot2)
    library(plotly)
    library(AnnotationDbi)
    library(htmlwidgets)
    library(proteoCraft)
    assign("myPG", readr::read_rds(tmpFl1), .GlobalEnv)
    assign("myPep", readr::read_rds(tmpFl2), .GlobalEnv)
    return()
  }))
  tmPlots <- try(parLapply(parClust, 1L:nrow(samplesDF), .plot_Rank_OR_Prof), silent = TRUE)
  #tmPlots <- try(parLapply(parClust, grep("^ranked\\.", rownames(samplesDF)), .plot_Rank_OR_Prof), silent = TRUE) # Test only ranked abundance plots
  #tmPlots <- try(parLapply(parClust, grep("^profiles\\.", rownames(samplesDF)), .plot_Rank_OR_Prof), silent = TRUE) # Test only profile plots
  # if (inherits(tmPlots, "try-error")) {
  #   warning("Going the slow way...")
  #   tmPlots <- try(lapply(1L:nrow(samplesDF), .plot_Rank_OR_Prof), silent = TRUE) 
  # }
  if (inherits(tmPlots, "try-error")) { stop(paste(c("Error: something's gone wrong:", tmPlots), collapse = "\n")) }
  unlink(tmpFl1)
  unlink(tmpFl2)
  # Ranked - something - plots
  tstL <- lengths(samplesDF$values)
  if (runRankAbundPlots) {
    for (quantType in QuantTypes) { #quantType <- QuantTypes[1L]
      w1 <- which((samplesDF$type == "PG")&(samplesDF$QuantType == quantType)&(samplesDF$subtype == "All")&(tstL == 1L))
      ggQuantLy[[quantType]] <- setNames(tmPlots[w1], samplesDF$values[w1])
    }
    for (quantType in pepQuantTypes) { #quantType <- pepQuantTypes[1L]
      w2 <- which((samplesDF$type == "pep")&(samplesDF$QuantType == quantType)&(samplesDF$subtype == "All")&(tstL == 1L))
      ggQuantLy[[paste0("peptides ", quantType)]] <- setNames(tmPlots[w2], samplesDF$values[w2])
    }
    saveFun(ggQuantLy, file = paste0(MainDir, "/quantPlots.RData"))
  }
  if (runProfPlots) {
    # Profile plots
    for (quantType in QuantTypes) {
      w1 <- which((samplesDF2$type == "PG")&(samplesDF2$QuantType == quantType)&(samplesDF2$subtype == "All")&(tstL > 1L))
      stopifnot(length(w1) == 1L)
      ggProfLy[[quantType]] <- tmPlots[[w1]]
    }
    for (quantType in pepQuantTypes) {
      w2 <- which((samplesDF2$type == "pep")&(samplesDF2$QuantType == quantType)&(samplesDF2$subtype == "All")&(tstL > 1L))
      lw2 <- length(w2)
      stopifnot(lw2 <= 1L)
      if (lw2) {
        ggProfLy[[paste0("peptides ", quantType)]] <- tmPlots[[w2]]        
      }
    }
    saveFun(ggProfLy, file = paste0(MainDir2, "/profilePlots.RData"))
  }
}
