#### Code chunk - Protein group profile plots and ranked abundance plots
gc()
Script <- readLines(ScriptPath)
try({ stopCluster(parClust) }, silent = TRUE)
rm(list = ls()[which(!ls() %in% .obj)])
#runRankAbundPlots %<o% TRUE
#runProfPlots %<o% TRUE
if (runRankAbundPlots||runProfPlots) {
  if (Annotate) {
    AllTerms %<o% unique(unlist(strsplit(db$`GO-ID`, ";")))
    AllTermNames %<o% unique(unlist(strsplit(db$GO, ";")))
  }
  library(ggplot2)
  library(RColorBrewer)
  library(colorspace)
  library(ggplot2)
  library(plotly)
  library(AnnotationDbi)
  library(htmlwidgets)
  QuantTypes %<o% c("LFQ", "Coverage")
  if (CreateMSMSKol) { QuantTypes <- c(QuantTypes, "Spectra") }
  tstOrg2 %<o% c()
  if (tstOrg) {
    PG$temp <- PG[[pgOrgKol]]
    g <- grep(";", PG$temp)
    if (length(g)) {
      g2 <- g[grsep(Org$Organism[1], x = PG$temp[g])]
      if (length(g2)) {
        g <- g[which(!g %in% g2)]
        PG$temp[g2] <- Org$Organism[1]
      }
      if (length(g)) {
        PG$temp <- "Mixed"
      }
    }
    PG$temp[which(PG$`Potential contaminant` == "+")] <- "Contaminant"
    tstOrg2 <- aggregate(PG$temp, list(PG$temp), length)
    tstOrg2 <- tstOrg2[order(tstOrg2$x, decreasing = TRUE),]
    tstOrg2 <- tstOrg2$Group.1[which(tstOrg2$x > 1)]
    tstOrg2 <- tstOrg2[which(tstOrg2 != "Contaminant")]
  }
  abbrFun %<o% function(x) { #x <- tst[[1]]
    g1 <- grep("[A-Z]", x)
    g2 <- grep("[a-z]", x)
    w1 <- g2[which((g2-1) %in% g1)]
    x[w1] <- "."
    gsub("\\.[^ ]+", ".", paste(x, collapse = ""))
  }
  Org_Nm %<o% abbrFun(tstOrg2)
  MainDir <- paste0(wd, "/Ranked abundance")
  MainDir2 <- paste0(wd, "/Profile plots")
  if (scrptType == "withReps") { dirlist <- unique(c(dirlist, MainDir, MainDir2)) }
  ggQuant %<o% list()
  ggProf %<o% list()
  QuantLy %<o% list()
  ProfLy %<o% list()
  #
  myColors <- setNames("black", "-")
  myColors2 <- setNames(c("lightgrey", "brown"), c("-", "+"))
  # if (length(tstOrg2)) {
  #   myColorsB <- myColors <- setNames(colorRampPalette(c("blue", "green"))(length(tstOrg2)), tstOrg2)
  #   w <- which(names(myColorsB) != "In list")
  #   names(myColorsB)[w] <- gsub("[a-z]+ ", ". ", names(myColorsB)[w])
  #   names(myColorsB) <- gsub(";", "\n", names(myColorsB))
  #   myColorsB[["potential contaminant"]] <- myColors[["potential contaminant"]] <- "grey"
  # }
  if (length(tstOrg2)) {
    if (length(tstOrg2) >= 3) {
      if (length(tstOrg2) <= 12) {
        myColors <- c(myColors, setNames(brewer.pal(min(c(12, length(tstOrg2))), "Set3"), tstOrg2[min(c(12, length(tstOrg2)))]))
      }
    } else {
      myColors <- c(myColors, setNames(c("#8DD3C7", "#FFFFB3")[1:length(tstOrg2)], tstOrg2)) # We never expect more than a handful of organisms  
    }
    myColors[["Contaminant"]] <- "deepskyblue"
  }
  #myColorsB[[c("+", "In list")[(length(tstOrg2) > 0)+1]]] <- myColors[[c("+", "In list")[(length(tstOrg2) > 0)+1]]] <- "red"
  myColors["Specific"] <- "purple"
  myColors["In list"] <- "brown"
  colScale <- scale_colour_manual(name = "Category", values = myColors)
  #colScaleB <- scale_colour_manual(name = "Category", values = myColorsB)
  fillScale <- scale_fill_manual(name = "Category", values = myColors)
  colScale2 <- scale_colour_manual(name = "In list", values = myColors2)
  fillScale2 <- scale_fill_manual(name = "In list", values = myColors2)
  #
  if (scrptType == "withReps") {
    allAgg <- setNames(rep("VPAL" #"RSA"
                           , length(QuantTypes)), QuantTypes)
  }
  #
  if (scrptType == "withReps") {
    ref <- paste0("Mean ", prtRfRoot) #prtRfRoot
    GO_PG_col %<o% unique(unlist(strsplit(Param$GO.tabs, ";")))
    GO_filt %<o% (length(GO_PG_col) > 0)
    if (GO_filt) {
      if ((!exists("GO_terms"))&&(file.exists(paste0(wd, "/GO_terms.RData")))) { loadFun(paste0(wd, "/GO_terms.RData")) }
      GO_PG_col <- GO_PG_col[which(GO_PG_col %in% GO_terms$ID)]
      GO_filt <- length(GO_PG_col) > 0
    }
    GO_filter <- GO_PG_col
    #
    runMark <- exists("CompGOTerms")
  }
  if (scrptType == "noReps") {
    ref <- rev(PG.int.cols[which(PG.int.cols != paste0("Imput. ", PG.int.cols["Original"]))])[1]
  }
  #
  for (QuantType in QuantTypes) { #QuantType <- QuantTypes[1] #QuantType <- "LFQ" #QuantType <- "Coverage"
    source(parSrc, local = FALSE)
    invisible(clusterCall(parClust, function() {
      library(ggplot2)
      library(RColorBrewer)
      library(colorspace)
      library(ggplot2)
      library(plotly)
      library(AnnotationDbi)
      library(htmlwidgets)
      return()
    }))
    cat(" -> ", QuantType, "\n")
    QuantLy[[QuantType]] <- list()
    SubDir <- paste0(MainDir, "/", QuantType)
    if (scrptType == "withReps") { dirlist <- unique(c(dirlist, SubDir)) }
    if (!dir.exists(SubDir)) { dir.create(SubDir, recursive = TRUE) }
    kolnm <- c("log10 LFQ", "Coverage [%]", "Spectral count")[match(QuantType, QuantTypes)]
    #
    if (scrptType == "withReps") {
      Agg <- get(allAgg[QuantType])
      mySamples <- Agg$values
    }
    if (scrptType == "noReps") {
      mySamples <- Exp
    }
    #
    ref <- c(ref,
             "Sequence coverage [%] - ",
             "Spectral count - ")[match(QuantType, QuantTypes)]
    myKol <- paste0(ref, mySamples)
    Wh <- which(myKol %in% colnames(PG))
    myKol <- myKol[Wh]
    if (length(myKol)) {
      varkol <- c("Leading protein IDs", "Protein IDs", "Common Name (short)", "id", "Label")
      if (prot.list.Cond) { varkol <- c(varkol, "In list") }
      myData <- PG[, c(varkol, myKol)]
      if (prot.list.Cond) { myData$`In list`[which(myData$`In list` == "")] <- "-" }
      colnames(myData) <- gsub(topattern(ref), "", colnames(myData))
      test <- aggregate(myData$"Label", list(myData$"Label"), length)
      w <- which(test$x > 1)
      if (length(w)) {
        test <- test[w,]
        for (i in test$Group.1) {
          w <- which(myData$"Label" == i)
          myData$"Label"[w[2:length(w)]] <- paste0(myData$"Label"[w[2:length(w)]], "_", 2:length(w))
        }
      }
      myData <- set_colnames(reshape::melt.data.frame(myData, id.vars = varkol),
                             c(varkol[1:4], "Protein Group", c(c(), "In list")[prot.list.Cond], "Sample", "Y"))
      myData <- myData[which(is.all.good(myData$Y, 2)),]
      if (QuantType %in% c("Coverage", "Spectra")) { myData <- myData[which(myData$Y > 0),] }
      if (nrow(myData)) {
        mID <- match(myData$id, PG$id)
        if (length(tstOrg2)) { myData$Category <- PG$temp[mID] } else {
          myData$Category <- "-"
        }
        if (prot.list.Cond) {
          myData$Category[which(myData$"In list" == "+")] <- c("+", "In list")[(length(tstOrg2) > 0)+1]
        }
        lev <- c("In list", "+", tstOrg2, "-", "Contaminant")
        lev <- lev[which(lev %in% unique(myData$Category))]
        myData$Category <- factor(myData$Category, levels = lev)
        # For the "Band ID" workflow we want to delog LFQ to focus on the top proteins!
        if (scrptType == "withReps") {
          tst <- FALSE
        }
        if (scrptType == "noReps") {
          tst <- AnalysisParam$Type == "Band ID"
        }
        if ((QuantType == "LFQ")&&(tst)) {
          myData$Y <- 10^myData$Y
          kolnm <- "LFQ"
        }
        #
        myData$Value <- paste0(kolnm, ": ", myData$Y)
        if (GO_filt||runMark) {
          myFlt <- c()
          if (GO_filt) { myFlt <- unique(c(myFlt, GO_filter)) }
          if (runMark) { myFlt <- unique(c(myFlt, CompGOTerms)) }
          for (goID in myFlt) { #goID <- GO_filter[1]
            # Get children terms
            gofilter <- unique(unlist(c(goID,
                                        GOBPOFFSPRING[[goID]],
                                        GOCCOFFSPRING[[goID]],
                                        GOMFOFFSPRING[[goID]])))
            gofilter <- gofilter[which(!is.na(gofilter))]
            if (sum(gofilter %in% AllTerms)) {
              myData[[goID]] <- "-"
              wtst <- grsep2(gofilter, PG$`GO-ID`)
              #which(vapply(strsplit(PG$`GO-ID`[-wtst],";"), function(x) { sum(x %in% gofilter) }, 1) != 0)
              #which(vapply(strsplit(PG$`GO-ID`[wtst],";"), function(x) { sum(x %in% gofilter) }, 1) == 0)
              wtst2 <- which(myData$id %in% PG$id[wtst])
              myData[wtst2, goID] <- "+"
              #aggregate(myData$`GO:0043005`, list(myData$Sample), function(x) { setNames(aggregate(x, list(x), length)$x, c("-", "+")) })
              #View(PG[wtst, grep(topattern(PG.int.col), colnames(PG))])
            }
          }
        }
        myData2 <- data.table(Category = myData$Category,
                              PG = myData$`Protein Group`)
        myData2 <- myData2[, list(x = unique(Category)),
                           by = list(Group.1 = PG)]
        myData2 <- as.data.frame(myData2)
        myData2$Protein_group <- "darkgrey"
        w <- which(myData2$x == "Contaminant")
        suppressWarnings(myData2$Protein_group[w] <- brewer.pal(min(c(length(w), 12)), "Blues"))
        w <- which(myData2$x %in% tstOrg2)
        myData2$Protein_group[w] <- rainbow_hcl(length(w))
        w <- which(myData2$x %in% c("In list",  "+"))
        myData2$Protein_group[w] <- "red"
        myData$Protein_group <- myData2$Protein_group[match(myData$`Protein Group`, myData2$Group.1)]
        testReg <- QuantType == "LFQ"
        exports <- list("QuantType", "QuantTypes", "mySamples", "SubDir", "GO_filt", "wd",
                        "colScale", "fillScale", "colScale2", "fillScale2", "kolnm", "abbrFun",
                        "Exp", "scrptType", "testReg")
        if (GO_filt) {
          exports <- append(exports, "GO_filter")
        }
        clusterExport(parClust, exports, envir = environment())
        if (runRankAbundPlots) {
          cat("    - Drawing ranked abundance plots\n")
          readr::write_rds(myData, paste0(wd, "/tmp.RDS"))
          invisible(clusterCall(parClust, function() {
            myData <<- readr::read_rds(paste0(wd, "/tmp.RDS"))
            return()
          }))
          if (testReg) {
            rgKol <- paste0("Regulated - ", mySamples)
            myFilt <- which(rgKol %in% colnames(PG))
            rgKol <- rgKol[myFilt]
            tmpPG <- PG[, c("id", rgKol)]
            readr::write_rds(tmpPG, paste0(wd, "/tmp2.RDS"))
            invisible(clusterCall(parClust, function() {
              tmpPG <<- readr::read_rds(paste0(wd, "/tmp2.RDS"))
              return()
            }))
          } else { myFilt <- 1:length(mySamples) }
          #for (smpl in  mySamples[myFilt]) { #smpl <-  mySamples[myFilt][1]
          f0 <- function(smpl) { #smpl <- mySamples[myFilt][1]
            rgtst <- testReg 
            if (scrptType == "withReps") {
              smpl2 <- proteoCraft::cleanNms(smpl)
            }
            if (scrptType == "noReps") {
              smpl2 <- smpl
            }
            myData2 <- myData[which(myData$Sample == smpl),]
            myData2 <- myData2[order(myData2$Y, decreasing = TRUE),]
            myData2$"Protein Group" <- factor(myData2$"Protein Group", levels = myData2$"Protein Group")
            if (rgtst) {
              rg <- paste0("Regulated - ", smpl)
              myData2$Regulated <- tmpPG[match(myData2$id, tmpPG$id), rg]
              spc <- grep("^Specific", myData2$Regulated)
              if (length(spc)) {
                myData2sp <- myData2[spc,]
                myData2sp$Category <- "Specific"
              } else { rgtst <- FALSE }
            }
            intmin <- floor(min(myData2$Y))
            intmax <- ceiling(max(myData2$Y))
            intscale <- intmax-intmin
            #xmax <- max(c(max(c(0, which(myData2$Category == "+")))+round(nrow(myData2)/15), nrow(myData2)))
            xmax <- nrow(myData2)*18/15
            PltTst <- setNames(c(TRUE,
                                 "+" %in% myData2$"In list"), c("All", "List"))
            if (GO_filt) { for (goID in GO_filter) {
              PltTst[goID] <- goID %in% colnames(myData2)
            } }
            for (i in 1:length(PltTst)) { #i <- 1
              if (PltTst[i]) {
                myData3 <- myData2
                if (rgtst) {
                  myData3sp <- myData2sp
                }
                ttl <- paste0("Ranked abundance plots ",
                              c("LFQ", "coverage", "spectral counts")[match(QuantType, QuantTypes)], " - ", smpl2)
                if (i == 1) {
                  catnm <- "Category"
                  filt <- 1:nrow(myData3)
                }
                if (i == 2) {
                  ttl <- paste0(ttl, ", proteins of interest")
                  myData3$"In list" <- factor(myData3$"In list", levels = c("-", "+"))
                  catnm <- "In list"
                  filt <- which(myData3[[catnm]] == "+")
                }
                if (i > 2) {
                  goID <- names(PltTst)[i]
                  goID2 <- gsub("^GO:", "GO", goID)
                  nm <- AnnotationDbi::Term(goID)
                  myData3[[goID2]] <- factor(myData3[[goID]], levels = c("-", "+"))
                  ttl <- paste0(ttl, ", ", goID, " ", nm)
                  myColors3 <- setNames(c("lightgrey", "purple"), c("-", "+"))
                  colScale3 <- scale_colour_manual(name = goID, values = myColors3)
                  fillScale3 <- scale_fill_manual(name = goID, values = myColors3)
                  catnm <- goID
                  filt <- which(myData3[[catnm]] == "+")
                }
                ttl2 <- gsub("/|:|\\*|\\?|<|>|\\|", "-", ttl)
                myData3$y <- myData3$Y + intmax*0.01
                plot <- ggplot(myData3) +
                  geom_bar(stat = "identity", aes(x = `Protein Group`, y = Y, fill = .data[[catnm]], text = Value)) +
                  theme_bw() + ylab(kolnm) +
                  ggtitle(ttl, subtitle = paste0(length(unique(myData3$id)), " Protein Groups")) + 
                  coord_cartesian(xlim = c(1, xmax), ylim = c(intmin, intmax*1.2+intscale*0.05)) +
                  theme(panel.border = element_blank(),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        axis.line = element_line(colour = "black"),
                        axis.text.x = element_blank(), axis.ticks = element_blank(),
                        plot.margin = margin(r = 100))
                #poplot(plot)
                if (i == 1) { plot <- plot + fillScale }
                if (i == 2) { plot <- plot + fillScale2 }
                if (i > 2) { plot <- plot + fillScale3 }
                #poplot(plot)
                if (rgtst) {
                  plot <- plot +
                    ggplot2::geom_point(data = myData3sp,
                                        ggplot2::aes(`Protein Group`, Y+0.05, fill = Category),
                                        shape = 21)
                }
                #poplot(plot, 12, 22)
                plot_ly <- ggplotly(plot, tooltip = c("Protein Group", "text"))
                # a folder with external resources is created for each html plot!
                pth <- paste0(SubDir, "/", ttl2)
                plPath <- paste0(pth, ".html")
                setwd(SubDir) # For some reason, unless I do this the default selfcontained = TRUE argument gets ignored and
                tstPL <- try(saveWidget(partial_bundle(plot_ly), plPath), silent = TRUE)
                tstPL <- !("try-error" %in% class(tstPL))
                if ("try-error" %in% class(tstPL)) { tstPL <- try(saveWidget(plot_ly, plPath), silent = TRUE) }
                setwd(wd)
                #if ((i == 1)&&(!"try-error" %in% class(tstPL))) { system(paste0("open \"", plPath, "\"")) }
                plot <- plot +
                  geom_text(data = myData3[filt,], angle = 45, hjust = 0, cex = 3.5,
                            aes(`Protein Group`, y, colour = .data[[catnm]], label = `Protein Group`))
                #poplot(plot, 12, 22)
                evPlot <- proteoCraft::plotEval(plot)
              }
              if (i == 1) {
                Res <- list(plotly = plot_ly,
                            plotly_saved = tstPL,
                            plotly_path = plPath,
                            ggplot = list(title = ttl,
                                          path = pth,
                                          plot = evPlot,
                                          dpi = 150,
                                          width = 10,
                                          height = 10))
              }
            }
            return(Res)
          }
          tmp <- try(parLapply(parClust, mySamples[myFilt], f0), silent = TRUE)
          if ("try-error" %in% class(tmp)) {
            tmp <- try(lapply(mySamples[myFilt], f0), silent = TRUE)  
          }
          if ("try-error" %in% class(tmp)) {
            stop("Raaaaaahhhh!!!! This error is caused by how pandoc (called by htmlwidgets above) handles memory. I hoped that the non-parallel lapply above could avoid the issue encountered with parLapply, but this error proves me wrong. If you encounter this, you probably need more RAM.")
          } else {
            wN <- which(vapply(tmp, function(x) { x$plotly_saved }, TRUE))
            if (length(wN)){
              f1 <- function(x) {
                setwd(SubDir)
                htmlwidgets::saveWidget(x$plotly, x$plotly_path, selfcontained = TRUE)
                setwd(wd)
              }
              tst <- try(parLapply(parClust, tmp[wN], f1), silent = TRUE)
              if ("try-error" %in% class(tmp)) {
                tst <- try(lapply(tmp[wN], f1), silent = TRUE)  
              }
              if ("try-error" %in% class(tst)) {
                stop("Raaaaaahhhh!!!! This error is caused by how pandoc (called by htmlwidgets above) handles memory. I hoped that the non-parallel lapply above could avoid the issue encountered with parLapply, but this error proves me wrong. If you encounter this, you probably need more RAM.")
              }
            }
          }
          tmp <- setNames(tmp, mySamples[myFilt])
          QuantLy[[QuantType]] <- setNames(lapply(mySamples[myFilt], function(x) {
            tmp[[x]]$plotly
          }), mySamples[myFilt])
          ggQuant[[QuantType]] <- setNames(lapply(mySamples[myFilt], function(x) {
            tmp[[x]]$ggplot
          }), mySamples[myFilt])
          unlink(paste0(wd, "/tmp.RDS"))
          if (testReg) { unlink(paste0(wd, "/tmp2.RDS")) }
        }
        #}
        if ((runProfPlots)&&(length(myKol) > 1)) {
          SubDir <- paste0(MainDir2, "/", QuantType)
          if (!dir.exists(SubDir)) { dir.create(SubDir, recursive = TRUE) }
          baseTtl <- paste0("Protein group ", c("LFQ", "coverage", "spectral count")[match(QuantType, QuantTypes)],
                            " profiles")
          #
          cat("    - Drawing profile plot\n")
          #PltTst <- setNames(c(TRUE, "+" %in% myData$`In list`), c("All", "List"))
          PltTst <- setNames(c(TRUE, "+" %in% myData$`In list`, GO_filt, exists("CompGOTerms")),
                             c("All", "List", "GO", "Mark"))
          # if (GO_filt) {
          #   for (goID in GO_filter) {
          #     #PltTst[gsub(":", "", goID)] <- gsub(":", "", goID) %in% colnames(myData)
          #     PltTst$GO <- TRUE
          #   }
          # }
          wPlts <- which(PltTst)
          for (i in wPlts) { #i <- 1 #i <- 2 #i <- 3
            profData <- myData
            if (scrptType == "withReps") {
              profData$Sample <- proteoCraft::cleanNms(profData$Sample)
            }
            goOn <- TRUE
            if (names(PltTst)[i] == "All") {
              ttl <- baseTtl
              klstKol <- c()
              if ((exists("KlustKols"))&&(length(KlustKols))) {
                klstKol <- grep(" - Global$", KlustKols, value = TRUE)[1]
                if (length(klstKol)) {
                  profData$Cluster <- PG[match(profData$`Leading protein IDs`, PG$`Leading protein IDs`), klstKol]
                } 
              }
              profData$Category <- gsub(" *[\\(\\[].*", "", profData$Category)
              profData$Category[which(profData$Category == "Contaminant")] <- "Cont."
              lev1 <- unique(profData$Category)
              lev1 <- lev1[which(!lev1 %in% c("In list", "Cont."))]
              lev1 <- c(lev1, "In list", "Cont.")
              catnm <- "Category"
              profData$Category <- factor(profData$Category, levels = lev1)
              if (length(klstKol)) {
                uCl <- sort(unique(profData$Cluster))
                uCl <- uCl[which(!is.na(uCl))]
                w <- which(is.na(profData$Cluster))
                profData$Cluster[w] <- "Not clustered"
                profData$Cluster <- factor(profData$Cluster, levels = c(uCl, "Not clustered"))
                frm <- paste0("`", catnm, "`~Cluster")
              } else {
                frm <- paste0(".~`", catnm, "`")
              }
              frm <- as.formula(frm)
              myFacets <- facet_grid(frm)
              Ngl <- c(0, 90)[(length(levels(profData[[catnm]])) > 5) + 1]
            }
            if (names(PltTst)[i] == "List") {
              ttl <- paste0(baseTtl, ", proteins of interest")
              profData$"In list" <- factor(profData$"In list", levels = c("-", "+"))
              catnm <- "In list"
              frm <- as.formula(paste0(".~`", catnm, "`"))
              myFacets <- facet_grid(frm)
              Ngl <- c(0, 90)[(length(levels(profData[[catnm]])) > 5) + 1]
            }
            if (names(PltTst)[i] %in% c("GO", "Mark")) {
              if (names(PltTst)[i] == "GO") {
                ttl <- paste0(baseTtl, ", GO terms of interest")
                myFlt <- GO_filter
              }
              if (names(PltTst)[i] == "Mark") {
                ttl <- paste0(baseTtl, ", compartment markers")
                myFlt <- CompGOTerms
              }
              catnm <- "GO term"
              myFlt <- myFlt[which(myFlt %in% colnames(profData))]
              goOn <- length(myFlt)
              if (goOn) {
                profData <- lapply(myFlt, function(go) { #go <- myFlt[1]
                  dat <- profData[which(profData[[go]] == "+"),]
                  dat$"GO term" <- go
                  dat <- dat[, which(!colnames(dat) %in% myFlt)]
                  return(dat)
                })
                profData <- do.call(rbind, profData)
                names(myFlt) <- gsub(" \\[GO:.*", "", GO_terms$Term[match(myFlt, GO_terms$ID)])
                profData$"GO term" <- names(myFlt)[match(profData$"GO term", myFlt)]
                profData$"GO term" <- factor(profData$"GO term", levels = names(myFlt))
                frm <- as.formula(paste0("~`", catnm, "`"))
                myFacets <- facet_wrap(frm)
                Ngl <- 0
              }
            }
            if (goOn) {
              profData$Sample <- factor(profData$Sample)
              ttl2 <- gsub("/|:|\\*|\\?|<|>|\\|", "-", ttl)
              m <- match(profData$Category, c("-", "Contaminant", Org_Nm, "+", "In list"))
              profData$LineType <- c(rep("dotted", 2), rep("dashed", length(Org_Nm)), rep("solid", 2))[m]
              #profData$DotSize <- c(rep(0.2, 2), rep(0.3, length(Org_Nm)), rep(0.5, 2))[m]
              profData$DotSize <- 1
              #profData$Alpha <- c(rep(0.25, 2), rep(0.5, length(Org_Nm)), rep(1, 2))[m]
              profData$Alpha <- 1
              profData <- profData[which(!is.na(profData$Y)),]
              wTxt <- which(profData$Sample == rev(levels(profData$Sample))[1])
              yLim <- c(min(profData$Y, na.rm = TRUE),
                        max(profData$Y, na.rm = TRUE))
              # It's important to have fixed colours here (using scale_linetype_identity()) so we can consistently plot in the app
              # different protein subsets with always the same colors assigned to each.
              ggCall_txt <- "plot <- ggplot(profData, aes(x = Sample, y = Y, text1 = `Protein Group`, text2 = Value, colour = Protein_group)) +
                  geom_line(aes(group = id), alpha = 0.1) +
                  geom_point(aes(size = DotSize)) + ggtitle(ttl) + ylab(kolnm) + myFacets +
                  theme_bw() + scale_size_identity(guide = \"none\") + scale_linetype_identity(guide = \"none\") +
                  theme(legend.position = \"none\",
                        axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1),
                        strip.text = element_text(face = \"bold\", size = 8, lineheight = 0.8, angle = Ngl),
                        strip.background = element_rect(fill = \"lightblue\", colour = \"black\", linewidth = 1)) +
                  scale_alpha_identity(guide = \"none\") + scale_colour_identity()"
              suppressWarnings({
                eval(parse(text = ggCall_txt))
              })
              #poplot(plot, 12, 20)
              plot_txt <- plot +
                geom_text(data = profData[wTxt,], aes(label = `Protein Group`, x = Sample, y = Y, alpha = Alpha, color = id),
                          hjust = 0, cex = 1)
              #poplot(plot_txt, 12, 20)
              if ((QuantType == "LFQ")&&(names(PltTst)[i] %in% c("GO", "Mark"))) {
                profData3 <- profData
                # Normalize by row
                protMeans <- aggregate(profData3$Y, list(profData3$`Leading protein IDs`), mean, na.rm = TRUE)
                m <- match(profData3$`Leading protein IDs`, protMeans$Group.1)
                profData3$Y <- profData3$Y - protMeans$x[m]
                # Calculate enveloppe
                profData3 <- aggregate(profData3$Y, list(profData3$Sample, profData3$`GO term`), function(x) {
                  n <- length(x)
                  m <- mean(x)
                  ci <- qt(0.975, df = n-1)*sd(x)/sqrt(n)
                  c(m, m+ci, m-ci)
                })
                colnames(profData3)[1:2] <- c("Sample", "GO term")
                profData3[, c("Y", "Y + 95% CI", "Y - 95% CI")] <- do.call(as.data.frame, list(profData3$x))
                profData3$x <- NULL
                ttl3 <- paste0("Avg. norm. ", ttl)
                ggCall_txt3 <- gsub("\n +", "\n", gsub("^ +", "", unlist(strsplit(ggCall_txt, " +\\+ *\n?"))))
                g <- grep("geom_line\\(", ggCall_txt3)
                ggCall_txt3[g] <- "geom_ribbon(alpha = 0.1, linetype = \"dotted\", aes(ymin = `Y - 95% CI`, ymax = `Y + 95% CI`))"
                g <- grep("ylab\\(", ggCall_txt3)
                ggCall_txt3[g] <- "ylab(\"log10 avg. norm. profile\")"
                g <- grep("ggtitle\\(", ggCall_txt3)
                ggCall_txt3[g] <- "ggtitle(ttl3, subtitle = \"ribbon = 95% confidence interval\")"
                g <- grep("geom_point\\(|scale_[a-z]+_identity\\(", ggCall_txt3)
                ggCall_txt3 <- ggCall_txt3[-g]
                ggCall_txt3 <- c(ggCall_txt3, "geom_line()")
                ggCall_txt3 <- paste(ggCall_txt3, collapse = " + \n")
                ggCall_txt3 <- gsub(", text1 = `Protein Group`, text2 = Value, colour = Protein_group\\)",
                                    ", group = `GO term`, colour = `GO term`, fill = `GO term`)",
                                    ggCall_txt3)
                ggCall_txt3 <- paste0("rib_", gsub("profData", "profData3", ggCall_txt3))
                suppressWarnings({
                  eval(parse(text = gsub("  +", " ", ggCall_txt3)))
                })
                #poplot(rib_plot, 12, 20)
                pth3 <- paste0(SubDir, "/", ttl3)
                if (names(PltTst)[i] == "GO") {
                  nm <- "GO trends"
                }
                if (names(PltTst)[i] == "Mark") {
                  nm <- "Compartment trends"
                }
                ggProf[[nm]] <- list(title = ttl3,
                                     path = pth3,
                                     call = ggCall_txt3,
                                     plot = plotEval(rib_plot),
                                     dpi = 150,
                                     width = 10,
                                     height = 10)
                # This is a simple plot, no need to plotly it!
              }
              plotlyCall_txt <- "plotlyProfiles <- ggplotly(plot, tooltip = c(\"text1\", \"text2\"))"
              eval(parse(text = plotlyCall_txt))
              # Grey plot for shiny app
              # plotlyProfiles2 <- plot_ly(profData)
              # plotlyProfiles2 <- add_trace(plotlyProfiles2, x = ~Sample, y = ~Y,
              #                              #split = ~`Protein Group`, # This would be correct but makes it slow, so this small approximation seems ok for now
              #                              color = I("lightgrey"), type = "scatter",
              #                              mode = "lines+markers", text = ~`Protein Group`, connectgaps = FALSE,
              #                              name = "", showlegend = FALSE)
              if (names(PltTst)[i] == "All") {
                ProfLy[[QuantType]] <- list(Plot = plotlyProfiles,
                                            ggCall = ggCall_txt,
                                            plCall = plotlyCall_txt,
                                            data = profData,
                                            data2 = list(colName = kolnm,
                                                         facets = myFacets,
                                                         angle = Ngl),
                                            Ttl = ttl
                                            #,grey = plotlyProfiles2
                )
              }
              pth <- paste0(SubDir, "/", ttl2)
              plPath <- paste0(pth, ".html")
              setwd(SubDir)
              tst <- try(saveWidget(partial_bundle(plotlyProfiles), plPath), silent = TRUE)
              if ("try-error" %in% class(tst)) {
                tst <- try(saveWidget(plotlyProfiles, plPath), silent = TRUE)
              }
              setwd(wd)
              if ((names(PltTst)[i] == "All")&&(!"try-error" %in% class(tst))) { system(paste0("open \"", plPath, "\"")) }
              #system(paste0("open \"",pth, ".html"))
              evPlot <- proteoCraft::plotEval(plot)
              evPlot_txt <- proteoCraft::plotEval(plot_txt)
              rs <- list(title = ttl,
                         data = profData,
                         path = pth,
                         yRange = yLim,
                         plot = evPlot_txt,
                         plot_no_text = evPlot,
                         dpi = 150,
                         width = 10,
                         height = 10)
              nm <- QuantType
              if (names(PltTst)[i] != "All") { nm <- paste0(nm, " ", names(PltTst)[i]) }
              ggProf[[nm]] <- list(title = ttl,
                                   data = profData,
                                   path = pth,
                                   yRange = yLim,
                                   plot = evPlot_txt,
                                   plot_no_text = evPlot,
                                   dpi = 150,
                                   width = 10,
                                   height = 10)
              setwd(wd)
            }
          }
        }
      }
    }
    stopCluster(parClust)
  }
  PG$temp <- NULL
  saveFun(ggQuant, file = paste0(MainDir, "/ggQuantPlots.RData"))
  saveFun(QuantLy, file = paste0(MainDir, "/QuantPlots.RData"))
  if (runProfPlots){
    saveFun(ggProf, file = paste0(MainDir2, "/ggProfilePlots.RData"))
    saveFun(ProfLy, file = paste0(MainDir2, "/ProfilePlots.RData"))
  }
  source(parSrc, local = FALSE)
  for (QuantType in QuantTypes) { #QuantType <- QuantTypes[1]
    tst <- parSapply(parClust, ggQuant[[QuantType]], function(x) { #x <- ggQuant[[QuantType]][[1]]
      dr <- dirname(x$path)
      if (!dir.exists(dr)) { dir.create(dr, recursive = TRUE) }
      suppressMessages({
        ggplot2::ggsave(paste0(x$path, ".pdf"), x$plot, dpi = x$dpi, width = x$width, height = x$height)
        ggplot2::ggsave(paste0(x$path, ".jpeg"), x$plot, dpi = x$dpi, width = x$width, height = x$height)
      })
    })
  }
  if (runProfPlots) {
    tst <- parLapply(parClust, ggProf, function(x) {
      dr <- dirname(x$path)
      if (!dir.exists(dr)) { dir.create(dr, recursive = TRUE) }
      suppressMessages({
        ggplot2::ggsave(paste0(x$path, ".pdf"), x$plot, dpi = x$dpi, width = x$width, height = x$height)
        ggplot2::ggsave(paste0(x$path, ".jpeg"), x$plot, dpi = x$dpi, width = x$width, height = x$height)
      })
    })
    for (nm in names(ggProf)) {
      ReportCalls <- AddPlot2Report(Plot = ggProf[[nm]]$plot,
                                    Dir = dirname(ggProf[[nm]]$path),
                                    Title = ggProf[[nm]]$title)
    }
  }
  setwd(wd)
}
# To do!
# - Use saving save_Plotlys.R script to save time
# - Can I speed up the ggplotly part for profile plots?
# - Can this benefit from better parallelisation?
# - plotly::plotly_build()?
