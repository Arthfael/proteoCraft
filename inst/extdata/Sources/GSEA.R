# Gene-Set Enrichment Analysis (GSEA)
#    This does not need filtering of the data (e.g. based on significance in a given statistical test).
#    GSEA uses external annotations and correlates it with the average fold change per comparison,
#    to show trends as to whether a specific set in enriched or not.

if (GSEAmode == "standard") {
  cat("   GSEA\n   ----\n")
}
if (GSEAmode == "WGCNA") {
  cat("   GSEA on WGCNA results\n   ---------------------\n")
}

keyType <- "UNIPROT"
idCol <- "Leading protein IDs"
if (!exists("GSEAmode")) { GSEAmode <- "standard" }
stopifnot(GSEAmode %in% c("standard", "WGCNA"))
if (GSEAmode == "standard") {
  if (dataType == "modPeptides") {
    dataType2 <- Ptm
    myData <- ptmpep
    if (scrptType == "withReps") { ratRef <- ptms.ratios.ref }
    if (scrptType == "noReps") { ratRef <- PTMs_ratRf[length(PTMs_ratRf)] }
    idCol <- "Protein"
    myData$Protein <- gsub(";.*", "", myData$Proteins)
    ohDeer <- paste0(wd, "/Reg. analysis/", ptm, "/GSEA")
  }
  if (dataType == "PG") {
    dataType2 <- dataType
    myData <- PG
    if (scrptType == "withReps") { ratRef <- Prot.Rat.Root }
    if (scrptType == "noReps") { ratRef <- PG.rat.cols }
    idCol <- "Leading protein IDs"
    ohDeer <- paste0(wd, "/Reg. analysis/GSEA")
  }
  if (scrptType == "withReps") {
    rankCol <- paste0(ratRef, myContrasts$Contrast)
  }
  if (scrptType == "noReps") {
    rankCol <- paste0(ratRef, Exp)
  }
  rankCol <- intersect(rankCol, colnames(myData))
  isOK <- length(rankCol) > 0L
}
if (GSEAmode == "WGCNA") {
  if (dataType == "modPeptides") {
    stop("Parameters for modified peptides not written yet! Why are you even running this source with these parameters?")
    dataType2 <- Ptm
  }
  if (dataType == "PG") {
    dataType2 <- dataType
    myData <- PGmodMembership
    myData$id <- colnames(exprData)
    idCol <- "id"
    ohDeer <- wgcnaDirs[3L]
    rankCol <- modNames
  }
  isOK <- TRUE
}

GSEA_plotly_fl %<o% paste0(wd, "/Reg. analysis/GSEA/GSEA_plotly.RDS")
if ((!exists("GSEA_plotly")) && file.exists(GSEA_plotly_fl)) {
  loadFun(GSEA_plotly_fl)
}
if (!exists("GSEA_plotly")) {
  GSEA_plotly <- list()
}
if (!GSEAmode %in% names(GSEA_plotly)) { GSEA_plotly[[GSEAmode]] <- list() }
GSEA_plotly[[GSEAmode]][[dataType2]] <- list()

ohDeer <- paste0(ohDeer, c("", "/svg"))
for (dr in ohDeer) {
  if (!dir.exists(dr)) { dir.create(dr, recursive = TRUE) }
}
if (exists("dirlist")) { dirlist <- union(dirlist, ohDeer) }
packs <- c()
exports <- c("packs", "idCol", "rankCol", "keyType", "cpParam", "Annotate", "rdsFls")
if (isOK) {
  if (Annotate) {
    
    # Either we can use the annotations we already have
    if (!exists("GO_mappings_fl")) { GO_mappings_fl %<o% paste0(wd, "/GO_mappings.RDS") }
    if (!exists("GO_terms_fl")) { GO_terms_fl %<o% paste0(wd, "/GO_terms.RDS") }
    if (!exists("GO_mappings")) { try(loadFun(GO_mappings_fl), silent = TRUE) }
    if (!exists("GO_terms")) { try(loadFun(GO_terms_fl), silent = TRUE) }
    if (sum(!c(exists("GO_mappings"), exists("GO_terms")))) {
      Src <- paste0(libPath, "/extdata/Sources/GO_prepare.R") # Doing this earlier but also keep latter instance for now
      #rstudioapi::documentOpen(Src)
      source(Src)
    }
    term2Prot <- GO_mappings$Protein
    term2Prot <- listMelt(strsplit(term2Prot$Protein, ";"), term2Prot$GO, c("protein", "term"))
    term2Prot <- term2Prot[, c("term", "protein")] # The order matters actually, not the name!
    term2name <- data.frame(term = GO_terms$ID,
                            name = GO_terms$Term)
  } else {
    # Or we will have to get an annotations package (currently 20 organisms are supported)
    if ((!exists("Org"))||(!is.data.frame(Org))||(nrow(Org) != 1L)) {
      kol <- intersect(c("Organism_Full", "Organism"), colnames(db))
      tst <- sapply(kol, \(x) { length(unique(db[!as.character(db[[x]]) %in% c("", "NA"), x])) })
      kol <- kol[order(tst, decreasing = TRUE)][1L]
      w <- which(db$`Potential contaminant` != "+")
      Org %<o% aggregate(w, list(db[w, kol]), length)
      colnames(Org) <- c("Organism", "Count")
      isOK <- max(Org$Count) >= nrow(db) * 0.3
      if (isOK) {
        Org <- Org[which(Org$Count == max(Org$Count))[1L],]
        Org$Source <- aggregate(db$Source[db[[kol]] %in% Org$Organism], list(db[db[[kol]] %in% Org$Organism, kol]), \(x) {
          unique(x[!is.na(x)])
        })$x
        Org$Source[is.na(Org$Source)] <- ""
      }
    }
    if (isOK) {
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
      organism <- if (Org$Organism %in% orgDBs$Full) { Org$Organism } else {
        dlg_list(c(orgDBs$Full, "none of these"), orgDBs$Full[1L], title = "Organism")$res
      }
      if (!length(organism)) { organism <- "none of these" }
      isOK <- organism != "none of these" 
    }
    # See https://guangchuangyu.github.io/2016/01/go-analysis-using-clusterprofiler/ for how to create annotations with the format clusterProfiler expects
    #BiocManager::install("AnnotationHub")
    #library(AnnotationHub)
    #hub <- AnnotationHub()
    #query(hub, organism)
    #myOrgAnnot <- hub[db$`Protein ID`]
    if (isOK) {
      orgDBpkg <- orgDBs$db[match(organism, orgDBs$Full)]
      packs2 <- c("AnnotationDbi", orgDBpkg)
      for (pck in packs2) {
        if (!require(pck, character.only = TRUE)) {
          pak::pak(pck, upgrade = FALSE, ask = FALSE)
        }
      }
      for (pck in packs2) {
        library(pck, character.only = TRUE)
      }
      packs <- union(packs, packs2)
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
  packs3 <- c("GO.db", "clusterProfiler", "BiocParallel", "pathview", "enrichplot", "DOSE")
  for (pck in packs3) {
    if (!require(pck, character.only = TRUE)) {
      pak::pak(pck, upgrade = FALSE, ask = FALSE)
    }
  }
  for (pck in packs3) {
    library(pck, character.only = TRUE)
  }
  packs <- union(packs, packs3)
  tmpDat <- myData[(nchar(myData[[idCol]]) > 0L) & (!is.na(myData[[idCol]])),
                   c(idCol, rankCol)]
  if (length(unique(tmpDat[[idCol]])) < nrow(tmpDat)) {
    tmpDat <- aggregate(tmpDat[, rankCol], list(tmpDat[[idCol]]), mean, na.rm = TRUE)
    colnames(tmpDat) <- c(idCol, rankCol)
  }
  cpParam <- SerialParam()
  rdsFls <- paste0(wd, "/", c("tmpDat",
                              "term2Prot",
                              "term2name",
                              "gses",
                              "minDB"),
                   ".RDS")
  readr::write_rds(tmpDat, rdsFls[1L])
  if (Annotate) {
    readr::write_rds(term2Prot, rdsFls[2L])
    readr::write_rds(term2name, rdsFls[3L])
  } else {
    exports <- c(exports, "orgDBpkg")
  }
  source(parSrc)
  clusterExport(parClust, exports, envir = environment())
  invisible(clusterCall(parClust, \(x) {
    for (pck in packs) { library(pck, character.only = TRUE) }
    assign("tmpDat", readr::read_rds(rdsFls[1L]), envir = .GlobalEnv)
    if (Annotate) {
      assign("term2Prot", readr::read_rds(rdsFls[2L]), envir = .GlobalEnv)
      assign("term2name", readr::read_rds(rdsFls[3L]), envir = .GlobalEnv)
    }
    return()
  }))
  #
  cat("    - Running analysis...\n")
  f0 <- \(kol, userAnnot = Annotate) { #kol <- rankCol[1L]
    tmp <- setNames(tmpDat[[kol]],
                    gsub(";.*| - .*", "", tmpDat[[idCol]]))
    tmp <- tmp[!is.na(tmp)]
    tmp <- na.omit(tmp)
    tmp <- sort(tmp, decreasing = TRUE)
    tmp <- tmp[nchar(names(tmp)) > 0L]
    #View(tmp)
    gse <- suppressMessages({
      if (userAnnot) {
        GSEA(geneList = tmp,
             TERM2GENE = term2Prot,
             TERM2NAME = term2name,
             nPerm = 10000L,
             minGSSize = 3L,
             maxGSSize = 800L,
             pvalueCutoff = 0.05,
             pAdjustMethod = "none",
             verbose = TRUE,
             BPPARAM = cpParam)
      } else {
        gseGO(tmp,
              ont = "ALL", 
              keyType = keyType, 
              nPerm = 10000L, 
              minGSSize = 3L, 
              maxGSSize = 800L, 
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
  gses <- clusterApplyLB(parClust, rankCol, f0)
  #
  # Rename
  names(gses) <- rankCol
  if (GSEAmode == "standard") {
    names(gses) <- cleanNms(gsub(topattern(ratRef), "", names(gses)))
  }
  #
  unlink(rdsFls[1L])
  if (Annotate) {
    unlink(rdsFls[2L])
    unlink(rdsFls[3L])
  }
  #
  #d <- GOSemSim::godata(annoDb = orgDBpkg, ont = "BP") # It seems to make sense to use BP here since we are interested in which biological processes are reacting to the perturbation
  #
  cat("    - Drawing plots...\n")
  nmRoots <- paste0("GSEA ", c("dotplot",
                               "enrichment map",
                               "category net plot",
                               "ridge plot"))
  GSEA_plotFun <- \(grp, type) {
    switch(type,
           dotplot = GSEA_dotplotFun(grp),
           enrichment = GSEA_enrichFun(grp),
           net = GSEA_netFun(grp),
           ridge = GSEA_ridgeFun(grp))
  }
  GSEA_dotplotFun <- \(x,
                       nmRoot = nmRoots[1L],
                       nCat = 50L) { #x <- plotsDF$GSE[1L] #x <- plotsDF$GSE[9L]
    gse <- gses[[x]]$GSE
    if (!inherits(gse, "gseaResult")) { return() }
    try({
      ttl <- paste0(nmRoot, " _ ", x)
      svpth <- paste0(ohDeer, "/", ttl, ".", c("html", "svg"))
      plot <- clusterProfiler::dotplot(gse, showCategory = nCat, split = ".sign", font.size = 4L,
                                       label_format = 500L # don't you dare wrap my labels!!!
      ) + ggplot2::facet_grid(.~.sign) +
        ggplot2::coord_fixed(0.025) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                       plot.subtitle = ggplot2::element_text(hjust = 0.5),
                       axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(), 
                       axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(), 
                       panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                       plot.margin = ggplot2::margin(0, 0, 0, 0, "cm"))
      suppressMessages({
        plot <- plot + viridis::scale_fill_viridis()
        #plot <- dotplot(gse, showCategory = nCat, color = "pvalue", split = ".sign") + facet_grid(.~.sign)
        #poplot(plot)
        plotL <- plotly::ggplotly(plot)
        # Fix tooltip
        w <- which(vapply(1L:length(plotL$x$data), \(i) {
          x <- plotL$x$data[[i]]
          ("hoveron" %in% names(x)) && (x$hoveron == "points")
        }, TRUE))
        for (i in w) {
          plotL$x$data[[i]]$text <- gsub("I\\(enrichplot_point_shape\\): *21<br */>", "", plotL$x$data[[i]]$text)
        }
        plotL$x$layout$xaxis$autorange <- TRUE
        plotL$x$layout$yaxis$autorange <- TRUE
        plotL <- htmlwidgets::onRender(plotL, global_autorange)
        plotL <- plotly::config(plotL,
                                modeBarButtonsToRemove = c("select2d", "lasso2d"))
        plotL <- plotly::plotly_build(plotL)
        #plotL <- plotly::partial_bundle(plotL)
        plotL2 <- plotly::layout(plotL,
                                 title = list(text = x,
                                              automargin = TRUE,
                                              subtitle = list(text = nmRoot)))
        wd0 <- getwd()
        setwd(ohDeer[1L])
        htmlwidgets::saveWidget(plotly::partial_bundle(plotL2), svpth[1L], selfcontained = TRUE)
        #htmlwidgets::saveWidget(plotL2, svpth[1L], selfcontained = TRUE)
        setwd(wd0)
        #
        plot <- plot + ggtitle(x, subtitle = nmRoot)
        ggplot2::ggsave(svpth[2L], plot, dpi = 300L, width = 7L, height = 7L, unit = "in")
      })
      return(plotL)
    }, silent = TRUE)
  }
  GSEA_enrichFun <- \(x,
                      nmRoot = nmRoots[2L],
                      nCat = 50L) { #x <- plotsDF$GSE[1L] #x <- plotsDF$GSE[9L]
    gse <- gses[[x]]$GSE
    if (!inherits(gse, "gseaResult")) { return() }
    g <- grep("^NA(\\.[0-9]+)?$", rownames(gse@result), invert = TRUE)
    gse@result <- gse@result[g,]
    nCat <- min(c(nCat, nrow(gse@result)))
    try({
      gse2 <- pairwise_termsim(gse, method = "JC", semData = NULL)
      ttl <- paste0(nmRoot, " _ ", x)
      svpth <- paste0(ohDeer, "/", ttl, ".", c("html", "svg"))
      plot <- clusterProfiler::emapplot(gse2, showCategory = nCat) +
        ggplot2::theme_bw() +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                       plot.subtitle = ggplot2::element_text(hjust = 0.5),
                       axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(),
                       panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                       plot.margin = ggplot2::margin(0, 0, 0, 0, "cm"))
      suppressMessages({
        plot <- plot + viridis::scale_color_viridis(option = "cividis", direction = -1L)
        l <- length(plot$layers)
        w <- which(vapply(1L:l, \(x) { inherits(plot$layers[[x]]$geom, "GeomTextRepel") }, TRUE))
        plot$layers[[w]]$aes_params$size <- 4L
        #getMethod("emapplot", "gseaResult")
        #poplot(plot)
        plotL <- plotly::ggplotly(plot)
        # Fix tooltip
        w <- which(vapply(1L:length(plotL$x$data), \(i) {
          x <- plotL$x$data[[i]]
          ("hoveron" %in% names(x)) && (x$hoveron == "points")
        }, TRUE))
        for (i in w) {
          plotL$x$data[[i]]$text <- paste0(plot@data$name, "<br />",
                                           gsub("[xy]: *-?[0-9]+(\\.[0-9]+)?<br />", "", plotL$x$data[[i]]$text))
        }
        plotL$x$layout$xaxis$autorange <- TRUE
        plotL$x$layout$yaxis$autorange <- TRUE
        plotL <- htmlwidgets::onRender(plotL, global_autorange)
        plotL <- plotly::config(plotL,
                                modeBarButtonsToRemove = c("select2d", "lasso2d"))
        plotL <- plotly::plotly_build(plotL)
        #plotL <- plotly::partial_bundle(plotL)
        plotL2 <- plotly::layout(plotL,
                                 title = list(text = x,
                                              automargin = TRUE,
                                              subtitle = list(text = nmRoot)))
        wd0 <- getwd()
        setwd(ohDeer[1L])
        htmlwidgets::saveWidget(plotly::partial_bundle(plotL2), svpth[1L], selfcontained = TRUE)
        #htmlwidgets::saveWidget(plotL2, svpth[1L], selfcontained = TRUE)
        setwd(wd0)
        #
        plot <- plot + ggtitle(x, subtitle = nmRoot)
        ggplot2::ggsave(svpth[2L], plot, dpi = 300L, width = 7L, height = 7L, unit = "in")
      })
      return(plotL)
    }, silent = TRUE)
  }
  GSEA_netFun <- \(x,
                   nmRoot = nmRoots[3L]) { #x <- plotsDF$GSE[1L] #x <- plotsDF$GSE[9L]
    gse <- gses[[x]]$GSE
    if (!inherits(gse, "gseaResult")) { return() }
    try({
      lFC <- gses[[x]]$lFC
      ttl <- paste0(nmRoot, " _ ", x)
      svpth <- paste0(ohDeer, "/", ttl, ".", c("html", "svg"))
      plot <- clusterProfiler::cnetplot(gse, foldChange = lFC, showCategory = 10L,
                                        color_edge = "grey",
                                        #cex_label_category = 1.2, cex_label_gene = 0.8 # Those parameters do not work for me...
      ) +
        ggplot2::theme_bw() +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                       plot.subtitle = ggplot2::element_text(hjust = 0.5),
                       axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(), 
                       axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(), 
                       panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                       plot.margin = ggplot2::margin(0, 0, 0, 0, "cm"))
      suppressMessages({
        plot <- plot + viridis::scale_color_viridis()
        # ... so I used a hacky solution:
        l <- length(plot$layers)
        w <- which(vapply(1L:l, \(x) { inherits(plot$layers[[x]]$geom, "GeomTextRepel") }, TRUE))
        plot$layers[[w]]$aes_params$size <- 1.6 # Downside: applies to both categories and proteins!
        # Edit labels
        plot$data$label <- plot$data$label
        w <- which(plot$data$label %in% miniDB$`Protein ID`)
        plot$data$label[w] <- miniDB$Label[match(plot$data$label[w], miniDB$`Protein ID`)]
        #
        #poplot(plot)
        #
        plotL <- plotly::ggplotly(plot)
        # Fix tooltip
        w <- which(vapply(1L:length(plotL$x$data), \(i) {
          x <- plotL$x$data[[i]]
          ("hoveron" %in% names(x)) && (x$hoveron == "points")
        }, TRUE))
        for (i in w) {
          plotL$x$data[[i]]$text <- paste0(plot@data$label, "<br />",
                                           gsub("[xy]: *-?[0-9]+(\\.[0-9]+)?<br />", "",
                                                gsub("<br />I\\(\\.hilight\\):.*", "", plotL$x$data[[i]]$text)))
        }
        # Re-add lost traces
        # - Nodes
        plot2 <- ggplot2::ggplot_build(plot)
        w <- match("geom_point", names(plot@layers))
        pt_dat <- plot2$data[[w]]
        pt_dat$label <- paste0("GO term: ",
                               do.call(paste, c(plot@data[1L:nrow(pt_dat), c("label", "size")], sep = "\nsize: "))) # GO terms go first in the table
        plotL <- plotly::add_trace(plotL,
                                   x = pt_dat$x,
                                   y = pt_dat$y,
                                   type = "scatter",
                                   mode = "markers",
                                   text = pt_dat$label,
                                   marker = list(color = pt_dat$colour,
                                                 size = pt_dat$size*4),
                                   hoverinfo = "text",
                                   showlegend = FALSE,
                                   inherit = FALSE)
        # - segments
        w <- match("geom_segment", names(plot@layers))
        seg_dat <- plot2$data[[w]]
        length(plotL$x$data)
        plotL <- plotly::add_trace(plotL,
                                   x = as.vector(rbind(seg_dat$x, seg_dat$xend, NA)),
                                   y = as.vector(rbind(seg_dat$y, seg_dat$yend, NA)),
                                   type = "scatter",
                                   mode = "lines",
                                   line = list(color = "black",
                                               width = 1),
                                   hoverinfo = "skip",
                                   showlegend = FALSE,
                                   inherit = FALSE)
        plotL <- plotly::plotly_build(plotL)
        l <- length(plotL$x$data)
        plotL$x$data <- c(plotL$x$data[l],
                          plotL$x$data[-l])
        #
        plotL$x$layout$xaxis$autorange <- TRUE
        plotL$x$layout$yaxis$autorange <- TRUE
        plotL <- htmlwidgets::onRender(plotL, global_autorange)
        plotL <- plotly::config(plotL,
                                modeBarButtonsToRemove = c("select2d", "lasso2d"))
        plotL <- plotly::plotly_build(plotL)
        #plotL <- plotly::partial_bundle(plotL)
        plotL2 <- plotly::layout(plotL,
                                 title = list(text = x,
                                              automargin = TRUE,
                                              subtitle = list(text = nmRoot)))
        wd0 <- getwd()
        setwd(ohDeer[1L])
        htmlwidgets::saveWidget(plotly::partial_bundle(plotL2), svpth[1L], selfcontained = TRUE)
        #htmlwidgets::saveWidget(plotL2, svpth[1L], selfcontained = TRUE)
        setwd(wd0)
        #
        plot <- plot + ggtitle(x, subtitle = nmRoot)
        ggplot2::ggsave(svpth[2L], plot, dpi = 300L, width = 7L, height = 7L, unit = "in")
      })
      return(plotL)
    }, silent = TRUE)
  }
  GSEA_ridgeFun <- \(x,
                     nmRoot = nmRoots[4L]) { #x <- plotsDF$GSE[1L] #x <- plotsDF$GSE[9L]
    gse <- gses[[x]]$GSE
    if (!inherits(gse, "gseaResult")) { return() }
    try({
      lFC <- gses[[x]]$lFC
      ttl <- paste0(nmRoot, " _ ", x)
      svpth <- paste0(ohDeer, "/", ttl, ".", c("html", "svg"))
      plot <- enrichplot::ridgeplot(gse, fill = "pvalue", label_format = 500L # don't you dare wrap my labels!!!
      ) + ggplot2::labs(x = "enrichment distribution") +
        ggplot2::theme(axis.text.x = ggplot2::element_text(size = 5L),
                       axis.text.y = ggplot2::element_text(size = 5L),
                       plot.title = ggplot2::element_text(hjust = 0.5),
                       plot.subtitle = ggplot2::element_text(hjust = 0.5),
                       panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                       plot.margin = ggplot2::margin(0, 0, 0, 0, "cm"))
      suppressMessages({
        plot <- plot + viridis::scale_fill_viridis()
        #poplot(plot)
        plotL <- plotly::ggplotly(plot)
        #
        plotL$x$layout$xaxis$autorange <- TRUE
        plotL$x$layout$yaxis$autorange <- TRUE
        plotL <- htmlwidgets::onRender(plotL, global_autorange)
        plotL <- plotly::config(plotL,
                                modeBarButtonsToRemove = c("select2d", "lasso2d"))
        plotL <- plotly::plotly_build(plotL)
        #plotL <- plotly::partial_bundle(plotL)
        plotL2 <- plotly::layout(plotL,
                                 title = list(text = x,
                                              automargin = TRUE,
                                              subtitle = list(text = nmRoot)))
        wd0 <- getwd()
        setwd(ohDeer[1L])
        htmlwidgets::saveWidget(plotly::partial_bundle(plotL2), svpth[1L], selfcontained = TRUE)
        #htmlwidgets::saveWidget(plotL2, svpth[1L], selfcontained = TRUE)
        setwd(wd0)
        #
        plot <- plot + ggtitle(x, subtitle = nmRoot)
        ggplot2::ggsave(svpth[2L], plot, dpi = 300L, width = 7L, height = 7L, unit = "in")
      })
      return(plotL)
    }, silent = TRUE)
  }
  #
  plotsDF <- data.frame(GSE = rep(names(gses), 4L),
                        type = unlist(lapply(c("dotplot", "enrichment", "net", "ridge"), \(x) {  rep(x, length(gses)) })))
  clusterExport(parClust, list("GSEA_dotplotFun", "GSEA_enrichFun", "GSEA_netFun", "GSEA_ridgeFun",
                               "GSEA_plotFun", "plotsDF", "nmRoots", "ohDeer",
                               "global_autorange"), envir = environment())
  readr::write_rds(gses, rdsFls[4L])
  if (!"Label" %in% colnames(db)) {
    db$Label <- do.call(paste, c(db[, c("Common Name", "Protein ID")], sep = "\n"))
  }
  readr::write_rds(db[, c("Protein ID", "Label")], rdsFls[5L])
  invisible(clusterCall(parClust, \() {
    for (pck in packs) { library(pck, character.only = TRUE) }
    environment(GSEA_plotFun) <- .GlobalEnv
    environment(GSEA_dotplotFun) <- .GlobalEnv
    environment(GSEA_enrichFun) <- .GlobalEnv
    environment(GSEA_netFun) <- .GlobalEnv
    environment(GSEA_ridgeFun) <- .GlobalEnv
    assign("gses", readr::read_rds(rdsFls[4L]), envir = .GlobalEnv)
    assign("miniDB", readr::read_rds(rdsFls[5L]), envir = .GlobalEnv)
    return()
  }))
  unlink(rdsFls[4L])
  unlink(rdsFls[5L])
  # NB:
  # Do not use clusterApplyLB here, it fails if a "try-error" object is returned
  # (it does not distinguish between my try-errors and a failed node calculation)
  # Going for clusterApplyLB would thus require rewriting the output of the functions as
  #   list(outcome = ..., plot = plotL)
  # instead of directly returning plotL...
  temp <- parLapply(parClust, 1L:nrow(plotsDF), \(i) { #i <- 9L
    grp <- plotsDF$GSE[i]
    type <- plotsDF$type[i]
    GSEA_plotFun(grp, type)
  })
  #
  # Check results
  errorTst <- which(vapply(temp, \(x) { inherits(x, "try-error") }, TRUE))
  #
  # Assign to list object
  # - GSEA dot plots
  w <- which(plotsDF$type == "dotplot")
  w <- setdiff(w, errorTst)
  if (length(w)) {
    GSEA_plotly[[GSEAmode]][[dataType2]][[nmRoots[1L]]] <- setNames(temp[w], plotsDF$GSE[w])
  }
  # - GSEA enrichment map plots
  w <- which(plotsDF$type == "enrichment")
  w <- setdiff(w, errorTst)
  if (length(w)) {
    GSEA_plotly[[GSEAmode]][[dataType2]][[nmRoots[2L]]] <- setNames(temp[w], plotsDF$GSE[w])
  }
  # - GSEA category net plots
  w <- which(plotsDF$type == "net")
  w <- setdiff(w, errorTst)
  if (length(w)) {
    GSEA_plotly[[GSEAmode]][[dataType2]][[nmRoots[3L]]] <- setNames(temp[w], plotsDF$GSE[w])
  }
  # - GSEA ridge plots
  w <- which(plotsDF$type == "ridge")
  w <- setdiff(w, errorTst)
  if (length(w)) {
    GSEA_plotly[[GSEAmode]][[dataType2]][[nmRoots[4L]]] <- setNames(temp[w], plotsDF$GSE[w])
  }
  #
  if (exists("DatAnalysisTxt") && (GSEAmode == "standard")) {
    l <- length(DatAnalysisTxt)
    DatAnalysisTxt[l] <- paste0(DatAnalysisTxt[l],
                                " Gene Set Enrichment Analysis was run using clusterProfiler.")
  }
  # See https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/ for more
  #
}
# Final cleanup
# - global env
for (pck in rev(packs)) {
  try(detach(paste0("package:", pck), unload = TRUE), silent = TRUE)
}
# - the cluster (not applicable if mode is WGCNA)
if (GSEAmode == "standard") {
  invisible(clusterCall(parClust, \(x) {
    for (pck in rev(packs)) {
      try(detach(paste0("package:", pck), unload = TRUE), silent = TRUE)
    }
    rm(list = ls())
    gc()
    return()
  }))
}
cat("    - Saving results...\n")
saveFun(GSEA_plotly, GSEA_plotly_fl)
#loadFun(GSEA_plotly_fl)
cat("    - Done!\n\n")
