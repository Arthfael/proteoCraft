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
  if (scrptType == "withReps") { ratRef <- paste0("Mean ", pepRatRf) }
  if (scrptType == "noReps") { ratRef <- PTMs_ratRf[length(PTMs_ratRf)] }
  idCol <- "Protein"
  myData$Protein <- gsub(";.*", "", myData$Proteins)
  namesRoot <- "Pep"
  ohDeer <- paste0(wd, "/Reg. analysis/", ptm, "/GSEA")
}
if (dataType == "PG") {
  #reNorm <- Norma.Prot.Ratio.classic
  myData <- PG
  if (scrptType == "withReps") { ratRef <- paste0("Mean ", Prot.Rat.Root) }
  if (scrptType == "noReps") { ratRef <- PG.rat.cols }
  idCol <- "Leading protein IDs"
  namesRoot <- "PG"
  ohDeer <- paste0(wd, "/Reg. analysis/GSEA")
}
if (scrptType == "withReps") {
  log2Col <- paste0(ratRef, VPAL$values)
}
if (scrptType == "noReps") {
  log2Col <- paste0(ratRef, Exp)
}
log2Col <- log2Col[which(log2Col %in% colnames(myData))]
if (!dir.exists(ohDeer)) { dir.create(ohDeer, recursive = TRUE) }
if (exists("dirlist")) { dirlist <- unique(c(dirlist, ohDeer)) }
isOK <- length(log2Col) > 0
if (isOK) {
  if ((!exists("Org"))||(!"data.frame" %in% class(Org))||(nrow(Org) != 1)) {
    kol <- c("Organism_Full", "Organism")
    kol <- kol[which(kol %in% colnames(db))]
    tst <- sapply(kol, function(x) { length(unique(db[which(!as.character(db[[x]]) %in% c("", "NA")), x])) })
    kol <- kol[order(tst, decreasing = TRUE)][1]
    w <- which(db$`Potential contaminant` != "+")
    Org %<o% aggregate(w, list(db[w, kol]), length)
    colnames(Org) <- c("Organism", "Count")
    Org <- Org[which(Org$Count == max(Org$Count)[1]),]
    Org$Source <- aggregate(db$Source[which(db[[kol]] %in% Org$Organism)], list(db[which(db[[kol]] %in% Org$Organism), kol]), function(x) {
      unique(x[which(!is.na(x))])
    })$x
    Org$Source[which(is.na(Org$Source))] <- ""
  }
  # For now only the following 20 organisms are supported, because we need their annotations package:
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
  # I will need to make this more universal.
  # ...
  if (Org$Organism %in% orgDBs$Full) { organism <- Org$Organism } else {
    organism <- dlg_list(c(orgDBs$Full, "none of these"),
                         orgDBs$Full[1], title = "Select organism")$res
  }
  if (!length(organism)) { organism <- "none of these" }
  isOK <- organism != "none of these"
  # See https://guangchuangyu.github.io/2016/01/go-analysis-using-clusterprofiler/ for how to create annotations with the format clusterProfiler expects
  #BiocManager::install("AnnotationHub")
  #library(AnnotationHub)
  #hub <- AnnotationHub()
  #query(hub, organism)
  #myOrgAnnot <- hub[db$`Protein ID`]
  # ...
  # test this: can I use this and skip altogether the db package installation step?
}
if (isOK) {
  orgDBpkg <- orgDBs$db[match(organism, orgDBs$Full)]
  packs <- c("GO.db", "clusterProfiler", "pathview", "enrichplot", "DOSE", orgDBpkg)
  for (pck in packs) {
    if (!require(pck, character.only = TRUE)) {
      pak::pkg_install(pck, upgrade = FALSE, ask = FALSE)
    }
  }
  if (usePar) {
    invisible(clusterCall(parClust, function() {
      for (pck in packs) { library(pck, character.only = TRUE) }
      return()
    }))
  } else {
    for (pck in packs) {
      library(pck, character.only = TRUE)
    }
  }
  eval(parse(text = paste0("myKeys <- keytypes(", orgDBpkg, ")")))
  if (!"UNIPROT" %in% myKeys) {
    if ((organism == "Arabidopsis thaliana")&&("TAIR" %in% colnames(db))) {
      myData$TAIR <- gsub(";.*", "", db$TAIR[match(gsub(";.*", "", myData[[idCol]]), db$`Protein ID`)])
      keyType <- idCol <- "TAIR"
    } else { isOK <- FALSE }
  }
}
if (isOK) {
  tmpDat <- myData[, c(idCol, log2Col)]
  tmpDat <- tmpDat[which((nchar(tmpDat[[idCol]]) > 0)&(!is.na(tmpDat[[idCol]]))),]
  if (length(unique(tmpDat[[idCol]])) < nrow(tmpDat)) {
    tmpDat <- aggregate(tmpDat[, log2Col], list(tmpDat[[idCol]]), mean, na.rm = TRUE)
    colnames(tmpDat) <- c(idCol, log2Col)
  }
  f0 <- function(kol) { #kol <- log2Col[1]
    tmp <- setNames(tmpDat[[kol]],
                    gsub(";.*", "", tmpDat[[idCol]]))
    tmp <- tmp[which(!is.na(tmp))]
    tmp <- na.omit(tmp)
    tmp <- sort(tmp, decreasing = TRUE)
    tmp <- tmp[which(nchar(names(tmp)) > 0)]
    #View(tmp)
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
    return(list(GSE = gse,
                lFC = tmp))
  }
  if (usePar) {
    clusterExport(parClust, list("idCol", "log2Col", "keyType", "tmpDat"), envir = environment())
    environment(f0) <- .GlobalEnv
    gses <- parLapply(parClust, log2Col, f0)
  } else {
    gses <- lapply(log2Col, f0)
  }
  gses <- setNames(gses, proteoCraft::cleanNms(gsub(proteoCraft::topattern(ratRef), "", log2Col)))
  #
  d <- GOSemSim::godata(annoDb = orgDBpkg, ont = "BP") # It seems to make sense to use BP here since we are interested in which biological processes are reacting to the perturbation
  nCat <- 50
  #
  # GSEA dot plots
  nmRoot <- "GSEA dotplot"
  lapply(names(gses), function(grp) { #grp <- names(gses)[1]
    gse <- gses[[grp]]$GSE
    try({
      plot <- dotplot(gse, showCategory = nCat, split = ".sign", font.size = 4,
                      label_format = 500 # don't you dare wrap my labels!!!
      ) + facet_grid(.~.sign) +
        coord_fixed(0.025)
      #plot <- dotplot(gse, showCategory = nCat, color = "pvalue", split = ".sign") + facet_grid(.~.sign)
      #poplot(plot)
      ggplot2::ggsave(paste0(ohDeer, "/", grp, " ", nmRoot, ".jpeg"), plot, dpi = 300)
      ggplot2::ggsave(paste0(ohDeer, "/", grp, " ", nmRoot, ".pdf"), plot, dpi = 300)
    }, silent = TRUE)
  })
  #
  # GSEA enrichment map plots
  nmRoot <- "GSEA enrichment map"
  lapply(names(gses), function(grp) { #grp <- names(gses)[1]
    gse <- gses[[grp]]$GSE
    try({
      gse2 <- pairwise_termsim(gse, method = "Wang", semData = d)
      plot <- emapplot(gse2, showCategory = nCat)
      l <- length(plot$layers)
      w <- which(sapply(1:l, function(x) { "GeomTextRepel" %in% class(plot$layers[[x]]$geom) }))
      plot$layers[[w]]$aes_params$size <- 2
      #getMethod("emapplot", "gseaResult")
      #poplot(plot)
      ggplot2::ggsave(paste0(ohDeer, "/", grp, " ", nmRoot, ".jpeg"), plot, dpi = 300)
      ggplot2::ggsave(paste0(ohDeer, "/", grp, " ", nmRoot, ".pdf"), plot, dpi = 300)
    }, silent = TRUE)
  })
  #
  # GSEA category net plots
  if (!"Label" %in% colnames(db)) {
    db$Label <- do.call(paste, c(db[, c("Common Name", "Protein ID")], sep = "\n"))
  }
  nmRoot <- "GSEA category net plot"
  # For this plot it would be nice to be able to:
  # - Do it for GO terms of interest only
  # - Update labels to ones more informative
  # - Plot as interactive plotly
  #  (see commented discussion below about how to achieve this)
  lapply(names(gses), function(grp) { #grp <- names(gses)[1]
    gse <- gses[[grp]]$GSE
    try({
      lFC <- gses[[grp]]$lFC
      plot <- cnetplot(gse, categorySize = "pvalue", foldChange = lFC, showCategory = 10, colorEdge = TRUE,
                       #cex_label_category = 1.2, cex_label_gene = 0.8 # Those parameters do not work for me...
      )
      # ... so I used a hacky solution:
      l <- length(plot$layers)
      w <- which(sapply(1:l, function(x) { "GeomTextRepel" %in% class(plot$layers[[x]]$geom) }))
      plot$layers[[w]]$aes_params$size <- 1.6 # Downside: applies to both categories and proteins!
      # Edit labels
      plot$data$label <- plot$data$label
      w <- which(plot$data$label %in% db$`Protein ID`)
      plot$data$label[w] <- db$Label[match(plot$data$label[w], db$`Protein ID`)]
      #
      #poplot(plot)
      ggplot2::ggsave(paste0(ohDeer, "/", grp, " ", nmRoot, ".jpeg"), plot, dpi = 300)
      ggplot2::ggsave(paste0(ohDeer, "/", grp, " ", nmRoot, ".pdf"), plot, dpi = 300)
      #
      # plot2 <- plotly::ggplotly(plot, tooltip = "label")
      # htmlwidgets::saveWidget(plot2, paste0(ohDeer, "/", grp, " ", nmRoot, ".html"),
      #                         selfcontained = TRUE)
      # Doesn't work, even though the 
      # Maybe the solution would be to get the data from the ggplot created, including the segment layer, which uses its own data,
      # and rewrite my own ggplot2 call?
      # That way I could also easily filter for specific GO terms.
      # TBC...
    }, silent = TRUE)
  })
  #
  # GSEA ridge plots
  nmRoot <- "GSEA ridge plot"
  lapply(names(gses), function(grp) { #grp <- names(gses)[1]
    gse <- gses[[grp]]$GSE
    try({
      plot <- ridgeplot(gse, 100,
                        label_format = 500 # don't you dare wrap my labels!!!
      ) + labs(x = "enrichment distribution") +
        theme(axis.text.x = element_text(size = 4.5),
              axis.text.y = element_text(size = 4.5))
      #poplot(plot)
      ggplot2::ggsave(paste0(ohDeer, "/", grp, " ", nmRoot, ".jpeg"), plot, dpi = 300)
      ggplot2::ggsave(paste0(ohDeer, "/", grp, " ", nmRoot, ".pdf"), plot, dpi = 300)
    }, silent = TRUE)
  })
  #
  if (exists("DatAnalysisTxt")) {
    DatAnalysisTxt <- paste0(DatAnalysisTxt,
                             " Gene Set Enrichment Analysis was run using clusterProfiler.")
    
  }
  # See https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/ for more
}
