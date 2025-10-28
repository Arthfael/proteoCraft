# Gene-Set Enrichment Analysis (GSEA)
#    This does not need filtering of the data (e.g. based on significance in a given statistical test).
#    GSEA uses external annotations and correlates it with the average fold change per comparison,
#    to show trends as to whether a specific set in enriched or not.

source(parSrc)

keyType <- "UNIPROT"
idCol <- "Leading protein IDs"
if (!exists("GSEAmode")) { GSEAmode <- "standard" }
if (GSEAmode == "standard") {
  if (dataType == "modPeptides") {
    myData <- ptmpep
    if (scrptType == "withReps") { ratRef <- paste0("Mean ", pepRatRf) }
    if (scrptType == "noReps") { ratRef <- PTMs_ratRf[length(PTMs_ratRf)] }
    idCol <- "Protein"
    myData$Protein <- gsub(";.*", "", myData$Proteins)
    ohDeer <- paste0(wd, "/Reg. analysis/", ptm, "/GSEA")
  }
  if (dataType == "PG") {
    myData <- PG
    if (scrptType == "withReps") { ratRef <- paste0("Mean ", Prot.Rat.Root) }
    if (scrptType == "noReps") { ratRef <- PG.rat.cols }
    idCol <- "Leading protein IDs"
    ohDeer <- paste0(wd, "/Reg. analysis/GSEA")
  }
  if (scrptType == "withReps") {
    rankCol <- paste0(ratRef, VPAL$values)
  }
  if (scrptType == "noReps") {
    rankCol <- paste0(ratRef, Exp)
  }
  rankCol <- rankCol[which(rankCol %in% colnames(myData))]
  isOK <- length(rankCol) > 0
}
if (GSEAmode == "WGCNA") {
  if (dataType == "modPeptides") {
    warning("Parameters for modified peptides not written yet! Why are you even running this source?")
  }
  if (dataType == "PG") {
    myData <- PGmodMembership
    myData$id <- colnames(exprData)
    idCol <- "id"
    ohDeer <- wgcnaDirs[3]
    rankCol <- modNames
  }
  isOK <- TRUE
}
if (!dir.exists(ohDeer)) { dir.create(ohDeer, recursive = TRUE) }
if (exists("dirlist")) { dirlist <- unique(c(dirlist, ohDeer)) }
if (isOK) {
  if (Annotate) {
    if (!exists("GO_mappings")) { try(loadFun("GO_mappings.RData"), silent = TRUE) }
    if (!exists("GO_terms")) { try(loadFun("GO_terms.RData"), silent = TRUE) }
    if (sum(!c(exists("GO_mappings"), exists("GO_terms")))) {
      Src <- paste0(libPath, "/extdata/R scripts/Sources/GO_prepare.R") # Doing this earlier but also keep latter instance for now
      #rstudioapi::documentOpen(Src)
      source(Src, local = FALSE)
    }
    term2Prot <- GO_mappings$Protein
    term2Prot <- listMelt(strsplit(term2Prot$Protein, ";"), term2Prot$GO, c("protein", "term"))
    term2Prot <- term2Prot[, c("term", "protein")] # The order matters actually, not the name!
    term2name <- data.frame(term = GO_terms$ID,
                            name = GO_terms$Term)
  } else {
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
    if (isOK) {
      orgDBpkg <- orgDBs$db[match(organism, orgDBs$Full)]
      packs <- c("GO.db", "clusterProfiler", "BiocParallel", "pathview", "enrichplot", "DOSE", orgDBpkg)
      for (pck in packs) {
        if (!require(pck, character.only = TRUE)) {
          pak::pkg_install(pck, upgrade = FALSE, ask = FALSE)
        }
      }
      for (pck in packs) {
        library(pck, character.only = TRUE)
      }
      clusterExport(parClust, "packs", envir = environment())
      invisible(clusterCall(parClust, function() {
        for (pck in packs) { library(pck, character.only = TRUE) }
        return()
      }))
      eval(parse(text = paste0("myKeys <- keytypes(", orgDBpkg, ")")))
      if (!"UNIPROT" %in% myKeys) {
        if ((organism == "Arabidopsis thaliana")&&("TAIR" %in% colnames(db))) {
          myData$TAIR <- gsub(";.*", "", db$TAIR[match(gsub(";.*", "", myData[[idCol]]), db$`Protein ID`)])
          keyType <- idCol <- "TAIR"
        } else { isOK <- FALSE }
      }
    }
  }
}
if (isOK) {
  tmpDat <- myData[which((nchar(myData[[idCol]]) > 0)&(!is.na(myData[[idCol]]))),
                   c(idCol, rankCol)]
  if (length(unique(tmpDat[[idCol]])) < nrow(tmpDat)) {
    tmpDat <- aggregate(tmpDat[, rankCol], list(tmpDat[[idCol]]), mean, na.rm = TRUE)
    colnames(tmpDat) <- c(idCol, rankCol)
  }
  cpParam <- SerialParam()
  saveRDS(tmpDat, paste0(wd, "/tmpDat.RDS"))
  exports <- list("idCol", "rankCol", "keyType", "wd", "cpParam", "Annotate")
  if (Annotate) {
    saveRDS(term2Prot, paste0(wd, "/term2Prot.RDS"))
    saveRDS(term2name, paste0(wd, "/term2name.RDS"))
  } else {
    exports <- append(exports, "orgDBpkg")
  }
  clusterExport(parClust, exports, envir = environment())
  invisible(clusterCall(parClust, function(x) {
    require(clusterProfiler)
    tmpDat <<- readRDS(paste0(wd, "/tmpDat.RDS"))
    if (Annotate) {
      term2Prot <<- readRDS(paste0(wd, "/term2Prot.RDS"))
      term2name <<- readRDS(paste0(wd, "/term2name.RDS"))
    } else {
      require(orgDBpkg, character.only = TRUE)
    }
    return()
  }))
  #
  f0 <- function(kol, userAnnot = Annotate) { #kol <- rankCol[1]
    tmp <- setNames(tmpDat[[kol]],
                    gsub(";.*| - .*", "", tmpDat[[idCol]]))
    tmp <- tmp[which(!is.na(tmp))]
    tmp <- na.omit(tmp)
    tmp <- sort(tmp, decreasing = TRUE)
    tmp <- tmp[which(nchar(names(tmp)) > 0)]
    #View(tmp)
    suppressMessages({
      if (userAnnot) {
        gse <- GSEA(geneList = tmp,
                    TERM2GENE = term2Prot,
                    TERM2NAME = term2name,
                    nPerm = 10000,
                    minGSSize = 3,
                    maxGSSize = 800,
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "none",
                    verbose = TRUE,
                    BPPARAM = cpParam)
      } else {
        gse <- gseGO(tmp,
                     ont = "ALL", 
                     keyType = keyType, 
                     nPerm = 10000, 
                     minGSSize = 3, 
                     maxGSSize = 800, 
                     pvalueCutoff = 0.05, 
                     verbose = TRUE, 
                     OrgDb = orgDBpkg, 
                     pAdjustMethod = "none",
                     BPPARAM = cpParam)
      }
    })
    return(list(GSE = gse,
                lFC = tmp))
  }
  environment(f0) <- .GlobalEnv
  gses <- parLapply(parClust, rankCol, f0)
  #
  # Rename
  names(gses) <- rankCol
  if (GSEAmode == "standard") {
    names(gses) <- proteoCraft::cleanNms(gsub(proteoCraft::topattern(ratRef), "", names(gses)))
  }
  #
  unlink(paste0(wd, "/tmpDat.RDS"))
  if (Annotate) {
    unlink(paste0(wd, "/term2Prot.RDS"))
    unlink(paste0(wd, "/term2name.RDS"))
  }
  #
  #d <- GOSemSim::godata(annoDb = orgDBpkg, ont = "BP") # It seems to make sense to use BP here since we are interested in which biological processes are reacting to the perturbation
  nCat <- 50
  #
  # GSEA dot plots
  nmRoot <- "GSEA dotplot"
  invisible(lapply(names(gses), function(grp) { #grp <- names(gses)[1]
    gse <- gses[[grp]]$GSE
    try({
      plot <- clusterProfiler::dotplot(gse, showCategory = nCat, split = ".sign", font.size = 4,
                      label_format = 500 # don't you dare wrap my labels!!!
      ) + ggplot2::facet_grid(.~.sign) +
        ggplot2::coord_fixed(0.025)
      suppressMessages({
        plot <- plot + viridis::scale_fill_viridis()
        #plot <- dotplot(gse, showCategory = nCat, color = "pvalue", split = ".sign") + facet_grid(.~.sign)
        #poplot(plot)
        ggplot2::ggsave(paste0(ohDeer, "/", grp, " ", nmRoot, ".jpeg"), plot, dpi = 300)
        ggplot2::ggsave(paste0(ohDeer, "/", grp, " ", nmRoot, ".pdf"), plot, dpi = 300)
      })
    }, silent = TRUE)
  }))
  #
  # GSEA enrichment map plots
  nmRoot <- "GSEA enrichment map"
  invisible(lapply(names(gses), function(grp) { #grp <- names(gses)[1]
    gse <- gses[[grp]]$GSE
    try({
      gse2 <- pairwise_termsim(gse, method = "Wang", semData = d)
      plot <- clusterProfiler::emapplot(gse2, showCategory = nCat)
      suppressMessages({
        plot <- plot + viridis::scale_color_viridis(option = "cividis", direction = -1)
        l <- length(plot$layers)
        w <- which(sapply(1:l, function(x) { "GeomTextRepel" %in% class(plot$layers[[x]]$geom) }))
        plot$layers[[w]]$aes_params$size <- 2
        #getMethod("emapplot", "gseaResult")
        #poplot(plot)
        ggplot2::ggsave(paste0(ohDeer, "/", grp, " ", nmRoot, ".jpeg"), plot, dpi = 300)
        ggplot2::ggsave(paste0(ohDeer, "/", grp, " ", nmRoot, ".pdf"), plot, dpi = 300)
      })
    }, silent = TRUE)
  }))
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
  invisible(lapply(names(gses), function(grp) { #grp <- names(gses)[1]
    gse <- gses[[grp]]$GSE
    try({
      lFC <- gses[[grp]]$lFC
      plot <- clusterProfiler::cnetplot(gse, categorySize = "pvalue", foldChange = lFC, showCategory = 10,
                                        colorEdge = TRUE,
                                        #cex_label_category = 1.2, cex_label_gene = 0.8 # Those parameters do not work for me...
      )
      suppressMessages({
        plot <- plot + viridis::scale_color_viridis()
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
        # Doesn't work...
        # Maybe the solution would be to get the data from the ggplot created,
        # including the segment layer, which uses its own data,
        # and rewrite my own ggplot2 call?
        # That way I could also easily filter for specific GO terms.
        # TBC...
      })
    }, silent = TRUE)
  }))
  #
  # GSEA ridge plots
  nmRoot <- "GSEA ridge plot"
  invisible(lapply(names(gses), function(grp) { #grp <- names(gses)[1]
    gse <- gses[[grp]]$GSE
    try({
      plot <- clusterProfiler::ridgeplot(gse, 100,
                                         label_format = 500 # don't you dare wrap my labels!!!
      ) + ggplot2::labs(x = "enrichment distribution") +
        ggplot2::theme(axis.text.x = ggplot2::element_text(size = 4.5),
                       axis.text.y = ggplot2::element_text(size = 4.5))
      suppressMessages({
        plot <- plot + viridis::scale_fill_viridis()
        #poplot(plot)
        ggplot2::ggsave(paste0(ohDeer, "/", grp, " ", nmRoot, ".jpeg"), plot, dpi = 300)
        ggplot2::ggsave(paste0(ohDeer, "/", grp, " ", nmRoot, ".pdf"), plot, dpi = 300)
      })
    }, silent = TRUE)
  }))
  #
  if ((exists("DatAnalysisTxt"))&&(GSEAmode == "standard")) {
    DatAnalysisTxt <- paste0(DatAnalysisTxt,
                             " Gene Set Enrichment Analysis was run using clusterProfiler.")
  }
  # See https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/ for more
}
