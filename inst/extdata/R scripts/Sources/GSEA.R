# Gene-Set Enrichment Analysis (GSEA)
#    This does not need filtering of the data (e.g. based on significance in a given statistical test).
#    GSEA uses external annotations and correlates it with the average fold change per comparison,
#    to show trends as to whether a specific set in enriched or not.

usePar <- FALSE # For now I cannot make the parallel version work ---> this is off for the moment
if (usePar) { source(parSrc) }

keyType <- "UNIPROT"
idCol <- "Leading protein IDs"
if (dataType == "modPeptides") {
  #reNorm <- FALSE
  myData <- ptmpep
  ratRef <- pepRatRf
  idCol <- "Protein"
  namesRoot <- "Pep"
  ohDeer <- paste0(wd, "/Reg. analysis/", ptm, "/GSEA")
}
if (dataType == "PG") {
  #reNorm <- Norma.Prot.Ratio.classic
  myData <- PG
  ratRef <- Prot.Rat.Root
  idCol <- "Leading protein IDs"
  namesRoot <- "PG"
  ohDeer <- paste0(wd, "/Reg. analysis/GSEA")
}
if (!dir.exists(ohDeer)) { dir.create(ohDeer, recursive = TRUE) }
dirlist <- unique(c(dirlist, ohDeer))
log2Col <- paste0("Mean ", ratRef, VPAL$values)
log2Col <- log2Col[which(log2Col %in% colnames(myData))]
isOK <- length(log2Col) > 0
if (isOK) {
  #
  orgDBs <- data.frame(Full = c("Homo sapiens",
                                "Pan troglodytes",
                                "Macaca mulatta",
                                "Mus musculus",
                                "Rattus norvegicus",
                                "Canis familiaris",
                                "Sus scrofa",
                                "Bos taurus",
                                "Gallus gallus",
                                "Xenopus laevis",
                                "Danio rerio",
                                "Caenorhabditis elegans",
                                "Drosophila melanogaster",
                                "Anopheles egypti",
                                "Arabidopsis thaliana",
                                "Saccharomyces cerevisiae",
                                "Plasmodium falciparum",
                                "Escherichia coli strain K12",
                                "Escherichia coli strain Sakai",
                                "Myxococcus xanthus"),
                       db = c("org.Hs.eg.db",
                              "org.Pt.eg.db",
                              "org.Mmu.eg.db",
                              "org.Mm.eg.db",
                              "org.Rn.eg.db",
                              "org.Cf.eg.db",
                              "org.Ss.eg.db",
                              "org.Bt.eg.db",
                              "org.Gg.eg.db",
                              "org.Xl.eg.db",
                              "org.Dr.eg.db",
                              "org.Ce.eg.db",
                              "org.Dm.eg.db",
                              "org.Ag.eg.db",
                              "org.At.tair.db",
                              "org.Sc.sgd.db",
                              "org.Pf.plasmo.db",
                              "org.EcK12.eg.db",
                              "org.EcSakai.eg.db",
                              "org.Mxanthus.db"))
  if (Org$Organism %in% orgDBs$Full) { organism <- Org$Organism } else {
    organism <- dlg_list(c(orgDBs$Full, "none of these"),
                         orgDBs$Full[1], title = "Select organism")$res
  }
  isOK <- organism != "none of these"
}
if (isOK) {
  orgDBpkg <- orgDBs$db[match(organism, orgDBs$Full)]
  packs <- c("GO.db", "clusterProfiler", "pathview", "enrichplot", "DOSE", orgDBpkg)
  for (pck in packs) {
    if (!require(pck, character.only = TRUE)) { pak::pkg_install(pck) }
  }
  for (pck in packs) {
    if (usePar) {
      clusterCall(parClust, function() library(pck, character.only = TRUE))
    } else {
      library(pck, character.only = TRUE)
    }
  }
  eval(parse(text = paste0("myKeys <- keytypes(", orgDBpkg, ")")))
  if (!"UNIPROT" %in% myKeys) {
    if ((organism == "Arabidopsis thaliana")&&("TAIR" %in% colnames(db))) {
      if (!"TAIR" %in% colnames(myData)) {
        myData$TAIR <- db$TAIR[match(gsub(";.*", "", myData[[idCol]]), db$`Protein ID`)]
      }
      keyType <- idCol <- "TAIR"
    } else { isOK <- FALSE }
  }
}
if (isOK) {
  tmpDat <- myData[, c(idCol, log2Col)]
  f0 <- function(kol) { #kol <- log2Col[1]
    tmp <- setNames(tmpDat[[kol]],
                    gsub(";.*", "", tmpDat[[idCol]]))
    tmp <- na.omit(tmp)
    tmp <- sort(tmp, decreasing = TRUE)
    tmp <- tmp[which(nchar(names(tmp)) > 0)]
    gse <- gseGO(tmp, 
                 ont = "ALL", 
                 keyType = keyType, 
                 nPerm = 10000, 
                 minGSSize = 3, 
                 maxGSSize = 800, 
                 pvalueCutoff = 0.05, 
                 verbose = TRUE, 
                 OrgDb = orgDBpkg, 
                 pAdjustMethod = "none")
    return(gse)
  }
  if (usePar) {
    clusterExport(parClust, list("idCol", "log2Col", "keyType", "tmpDat"), envir = environment())
    environment(f0) <- .GlobalEnv
    gses <- parLapply(parClust, log2Col, f0)
  } else {
    gses <- lapply(log2Col, f0)
  }
  gses <- setNames(gses, proteoCraft::cleanNms(gsub(proteoCraft::topattern(ratRef), "", log2Col)))
  lapply(names(gses), function(grp) {
    gse <- gses[[grp]]
    plot <- dotplot(gse, showCategory = 10, split = ".sign") + facet_grid(.~.sign)
    #plot <- dotplot(gse, showCategory = 10, color = "pvalue", split = ".sign") + facet_grid(.~.sign)
    ggplot2::ggsave(paste0(ohDeer, "/", grp, " GSEA dotplot.jpeg"), plot, dpi = 300)
    ggplot2::ggsave(paste0(ohDeer, "/", grp, " GSEA dotplot.pdf"), plot, dpi = 300)
  })
}
