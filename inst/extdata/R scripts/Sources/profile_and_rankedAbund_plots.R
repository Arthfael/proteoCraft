#### Code chunk - Protein group profile plots and ranked abundance plots
source(parSrc, local = FALSE)
rm(list = ls()[which(!ls() %in% .obj)])
invisible(clusterCall(parClust, function(x) { rm(list = ls());gc() }))
Script <- readLines(ScriptPath)
gc()
#runRankAbundPlots %<o% TRUE
#runProfPlots %<o% TRUE
if (runRankAbundPlots|runProfPlots) {
  library(ggplot2)
  library(RColorBrewer)
  library(colorspace)
  library(ggplot2)
  library(plotly)
  library(AnnotationDbi)
  library(htmlwidgets)
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
  QuantTypes %<o% c("LFQ", "Coverage")
  if (CreateMSMSKol) { QuantTypes <- c(QuantTypes, "Spectra") }
  tstOrg2 %<o% c()
  if (tstOrg) {
    PG$temp <- PG[[pgOrgKol]]
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
  tstOrg3 %<o% abbrFun(tstOrg2)
  MainDir <- paste0(wd, "/Ranked abundance")
  MainDir2 <- paste0(wd, "/Profile plots")
  if (scrptType == "withReps") { dirlist <- unique(c(dirlist, MainDir, MainDir2)) }
  ggQuant %<o% list()
  ggProf %<o% list()
  QuantLy %<o% list()
  ProfLy %<o% list()
  for (QuantType in QuantTypes) { #QuantType <- QuantTypes[1] #QuantType <- "LFQ" #QuantType <- "Coverage"
    cat(" + ", QuantType, "\n")
    QuantLy[[QuantType]] <- list()
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
    SubDir <- paste0(MainDir, "/", QuantType)
    if (scrptType == "withReps") { dirlist <- unique(c(dirlist, SubDir)) }
    if (!dir.exists(SubDir)) { dir.create(SubDir, recursive = TRUE) }
    kolnm <- c("log10 LFQ", "Coverage [%]", "Spectral count")[match(QuantType, QuantTypes)]
    if (scrptType == "withReps") {
      ref <- paste0("Mean ", prtRfRoot) #prtRfRoot
      GO_PG_col %<o% unique(unlist(strsplit(Param$GO.tabs, ";")))
      GO_filt %<o% length(GO_PG_col) > 0
      if (GO_filt) {
        if ((!exists("GO_terms"))&&(file.exists("GO_terms.RData"))) { loadFun("GO_terms.RData") }
        GO_PG_col <- GO_PG_col[which(GO_PG_col %in% GO_terms$ID)]
        GO_filt <- length(GO_PG_col) > 0
      }
      GO_filter <- GO_PG_col
    }
    if (scrptType == "noReps") {
      ref <- rev(PG.int.cols[which(PG.int.cols != paste0("Imput. ", PG.int.cols["Original"]))])[1]
    }
    ref <- c(ref,
             "Sequence coverage [%] - ",
             "Spectral count - ")[match(QuantType, QuantTypes)]
    if (scrptType == "withReps") {
      Agg <- get(c("VPAL",# "RSA",
                   "RSA",
                   "RSA")[match(QuantType, QuantTypes)])
      mySamples <- Agg$values
    }
    if (scrptType == "noReps") {
      mySamples <- Exp
    }
    myKol <- paste0(ref, mySamples)
    Wh <- which(myKol %in% colnames(PG))
    myKol <- myKol[Wh]
    if (length(myKol)) {
      varkol <- c("Leading protein IDs", "Protein IDs", "Common Name (short)", "id", "Label")
      if (prot.list.Cond) { varkol <- c(varkol, "In list") }
      temp <- PG[, c(varkol, myKol)]
      if (prot.list.Cond) { temp$`In list`[which(temp$`In list` == "")] <- "-" }
      colnames(temp) <- gsub(topattern(ref), "", colnames(temp))
      test <- aggregate(temp$"Label", list(temp$"Label"), length)
      w <- which(test$x > 1)
      if (length(w)) {
        test <- test[w,]
        for (i in test$Group.1) {
          w <- which(temp$"Label" == i)
          temp$"Label"[w[2:length(w)]] <- paste0(temp$"Label"[w[2:length(w)]], "_", 2:length(w))
        }
      }
      temp <- set_colnames(reshape::melt.data.frame(temp, id.vars = varkol),
                           c(varkol[1:4], "Protein Group", c(c(), "In list")[prot.list.Cond], "Sample", "Y"))
      temp <- temp[which(is.all.good(temp$Y, 2)),]
      if (QuantType %in% c("Coverage", "Spectra")) { temp <- temp[which(temp$Y > 0),] }
      if (nrow(temp)) {
        mID <- match(temp$id, PG$id)
        if (length(tstOrg2)) { temp$Category <- PG$temp[mID] } else {
          temp$Category <- "-"
        }
        if (prot.list.Cond) {
          temp$Category[which(temp$"In list" == "+")] <- c("+", "In list")[(length(tstOrg2) > 0)+1]
        }
        lev <- c("In list", "+", tstOrg2, "-", "Contaminant")
        lev <- lev[which(lev %in% unique(temp$Category))]
        temp$Category <- factor(temp$Category, levels = lev)
        # For the "Band ID" workflow we want to delog LFQ to focus on the top proteins!
        if (scrptType == "withReps") {
          tst <- FALSE
        }
        if (scrptType == "noReps") {
          tst <- AnalysisParam$Type == "Band ID"
        }
        if ((QuantType == "LFQ")&&(tst)) {
          temp$Y <- 10^temp$Y
          kolnm <- "LFQ"
        }
        #
        temp$Value <- paste0(kolnm, ": ", temp$Y)
        if (GO_filt) {
          for (goID in GO_filter) { #goID <- GO_filter[1]
            # Get children terms
            gofilter <- unique(unlist(c(goID,
                                        GOBPOFFSPRING[[goID]],
                                        GOCCOFFSPRING[[goID]],
                                        GOMFOFFSPRING[[goID]])))
            gofilter <- gofilter[which(!is.na(gofilter))]
            if (sum(gofilter %in% AllTerms)) {
              temp[[goID]] <- "-"
              wtst <- grsep2(gofilter, PG$`GO-ID`)
              #which(vapply(strsplit(PG$`GO-ID`[-wtst],";"), function(x) { sum(x %in% gofilter) }, 1) != 0)
              #which(vapply(strsplit(PG$`GO-ID`[wtst],";"), function(x) { sum(x %in% gofilter) }, 1) == 0)
              wtst2 <- which(temp$id %in% PG$id[wtst])
              temp[wtst2, goID] <- "+"
              #aggregate(temp$`GO:0043005`, list(temp$Sample), function(x) { setNames(aggregate(x, list(x), length)$x, c("-", "+")) })
              #View(PG[wtst, grep(topattern(PG.int.col), colnames(PG))])
            }
          }
        }
        temp2 <- data.table(Category = temp$Category,
                            PG = temp$`Protein Group`)
        temp2 <- temp2[, list(x = unique(Category)),
                       by = list(Group.1 = PG)]
        temp2 <- as.data.frame(temp2)
        temp2$Protein_group <- "darkgrey"
        w <- which(temp2$x == "Contaminant")
        suppressWarnings(temp2$Protein_group[w] <- brewer.pal(min(c(length(w), 12)), "Blues"))
        w <- which(temp2$x %in% tstOrg2)
        temp2$Protein_group[w] <- rainbow_hcl(length(w))
        w <- which(temp2$x %in% c("In list",  "+"))
        temp2$Protein_group[w] <- "red"
        temp$Protein_group <- temp2$Protein_group[match(temp$`Protein Group`, temp2$Group.1)]
        testReg <- QuantType == "LFQ"
        exports <- list("QuantType", "QuantTypes", "mySamples", "SubDir", "GO_filt", "GO_filter", "wd",
                        "colScale", "fillScale", "colScale2", "fillScale2", "kolnm", "abbrFun",
                        "Exp", "scrptType", "testReg")
        clusterExport(parClust, exports, envir = environment())
        if (runRankAbundPlots) {
          cat("   -> Drawing ranked abundance plots\n")
          saveRDS(temp, paste0(wd, "/tmp.RDS"))
          invisible(clusterCall(parClust, function() {
            temp <<- readRDS(paste0(wd, "/tmp.RDS"))
            return()
          }))
          if (testReg) {
            rgKol <- paste0("Regulated - ", mySamples)
            myFilt <- which(rgKol %in% colnames(PG))
            rgKol <- rgKol[myFilt]
            tmpPG <- PG[, c("id", rgKol)]
            saveRDS(tmpPG, paste0(wd, "/tmp2.RDS"))
            invisible(clusterCall(parClust, function() {
              tmpPG <<- readRDS(paste0(wd, "/tmp2.RDS"))
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
            temp2 <- temp[which(temp$Sample == smpl),]
            temp2 <- temp2[order(temp2$Y, decreasing = TRUE),]
            temp2$"Protein Group" <- factor(temp2$"Protein Group", levels = temp2$"Protein Group")
            if (rgtst) {
              rg <- paste0("Regulated - ", smpl)
              temp2$Regulated <- tmpPG[match(temp2$id, tmpPG$id), rg]
              spc <- grep("^Specific", temp2$Regulated)
              if (length(spc)) {
                temp2sp <- temp2[spc,]
                temp2sp$Category <- "Specific"
              } else { rgtst <- FALSE }
            }
            intmin <- floor(min(temp2$Y))
            intmax <- ceiling(max(temp2$Y))
            intscale <- intmax-intmin
            #xmax <- max(c(max(c(0, which(temp2$Category == "+")))+round(nrow(temp2)/15), nrow(temp2)))
            xmax <- nrow(temp2)*18/15
            PltTst <- setNames(c(TRUE,
                                 "+" %in% temp2$"In list"), c("All", "List"))
            if (GO_filt) { for (goID in GO_filter) {
              PltTst[goID] <- goID %in% colnames(temp2)
            } }
            for (i in 1:length(PltTst)) { #i <- 1
              if (PltTst[i]) {
                temp3 <- temp2
                if (rgtst) {
                  temp3sp <- temp2sp
                }
                ttl <- paste0("Ranked abundance plots ",
                                      c("LFQ", "coverage", "spectral counts")[match(QuantType, QuantTypes)], " - ", smpl2)
                if (i == 1) {
                  catnm <- "Category"
                  filt <- 1:nrow(temp3)
                }
                if (i == 2) {
                  ttl <- paste0(ttl, ", proteins of interest")
                  temp3$"In list" <- factor(temp3$"In list", levels = c("-", "+"))
                  catnm <- "In list"
                  filt <- which(temp3[[catnm]] == "+")
                }
                if (i > 2) {
                  goID <- names(PltTst)[i]
                  goID2 <- gsub("^GO:", "GO", goID)
                  nm <- AnnotationDbi::Term(goID)
                  temp3[[goID2]] <- factor(temp3[[goID]], levels = c("-", "+"))
                  ttl <- paste0(ttl, ", ", goID, " ", nm)
                  myColors3 <- setNames(c("lightgrey", "purple"), c("-", "+"))
                  colScale3 <- scale_colour_manual(name = goID, values = myColors3)
                  fillScale3 <- scale_fill_manual(name = goID, values = myColors3)
                  catnm <- goID
                  filt <- which(temp3[[catnm]] == "+")
                }
                ttl2 <- gsub("/|:|\\*|\\?|<|>|\\|", "-", ttl)
                temp3$y <- temp3$Y + intmax*0.01
                plot <- ggplot(temp3) +
                  geom_bar(stat = "identity", aes(x = `Protein Group`, y = Y, fill = .data[[catnm]], text = Value)) +
                  theme_bw() + ylab(kolnm) +
                  ggtitle(ttl, subtitle = paste0(length(unique(temp3$id)), " Protein Groups")) + 
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
                    ggplot2::geom_point(data = temp3sp,
                                        ggplot2::aes(`Protein Group`, Y+0.05, fill = Category),
                                        shape = 21)
                }
                #poplot(plot, 12, 22)
                plot_ly <- ggplotly(plot, tooltip = c("Protein Group", "text"))
                setwd(SubDir) # For some reason, unless I do this the default selfcontained = TRUE argument gets ignored and
                # a folder with external resources is created for each html plot!
                pth <- paste0(SubDir, "/", ttl2)
                plPath <- paste0(pth, ".html")
                tstPL <- try(saveWidget(partial_bundle(plot_ly), paste0(pth, ".html")), silent = TRUE)
                tstPL <- !("try-error" %in% class(tstPL))
                if ("try-error" %in% class(tstPL)) { tstPL <- try(saveWidget(plot_ly, paste0(ttl2, ".html")), silent = TRUE) }
                #if ((i == 1)&&(!"try-error" %in% class(tstPL))) { system(paste0("open \"", pth, ".html")) }
                plot <- plot +
                  geom_text(data = temp3[filt,], angle = 45, hjust = 0, cex = 3.5,
                            aes(`Protein Group`, y, colour = .data[[catnm]], label = `Protein Group`))
                #poplot(plot, 12, 22)
                evPlot <- plotEval(plot)
                #ggsave(filename = paste0(pth, ".jpeg"), plot, dpi = 600, width = 30, height = 10, units = "in")
                #ggsave(filename = paste0(pth, ".pdf"), plot, dpi = 600, width = 30, height = 10, units = "in")
                setwd(wd)
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
              # The hope here is maybe some succeeded and we can get away with 
              f1 <- function(x) {
                htmlwidgets::saveWidget(x$plotly, x$plotly_path, selfcontained = TRUE)
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
          cat("   -> Drawing profile plots\n")
          PltTst <- setNames(c(TRUE, "+" %in% temp$`In list`), c("All", "List"))
          if (GO_filt) { for (goID in GO_filter) {
            PltTst[gsub(":", "", goID)] <- gsub(":", "", goID) %in% colnames(temp)
          } }
          for (i in 1:length(PltTst)) { #i <- 1
            SubDir <- paste0(MainDir2, "/", QuantType)
            if (!dir.exists(SubDir)) { dir.create(SubDir, recursive = TRUE) }
            if (PltTst[i]) {
              temp2 <- temp
              if (scrptType == "withReps") {
                temp2$Sample <- proteoCraft::cleanNms(temp2$Sample)
              }
              temp2$Category <- gsub(" *[\\(\\[].*", "", temp2$Category)
              temp2$Category <- as.factor(temp2$Category)
              ttl <- paste0("Protein group ", c("LFQ", "coverage", "spectral count")[match(QuantType, QuantTypes)],
                            " profiles")
              if (i == 1) {
                lev <- levels(temp2$Category)
                w <- which(!lev %in% c("In list", "Contaminant"))
                lev <- c(c("In list", "Contaminant"), lev[w])
                w <- which(!lev %in% c("In list", "Contaminant"))
                lev2 <- lev
                lev2[w] <- sapply(strsplit(as.character(lev[w]), ""), abbrFun)
                lev2[which(lev2 == "Contaminant")] <- "Cont."
                lev2 <- gsub(";", "\n", lev2)
                temp2$Category <- as.character(temp2$Category)
                temp2$Category <- lev2[match(temp2$Category, lev)]
                temp2$Category <- factor(temp2$Category, levels = lev2)
                catnm <- "Category"
                filt <- 1:nrow(temp2)
              }
              if (i == 2) {
                ttl <- paste0(ttl, ", proteins of interest")
                temp2$"In list" <- factor(temp2$"In list", levels = c("-", "+"))
                catnm <- "In list"
                filt <- which(temp2[[catnm]] == "+")
              }
              if (i > 2) {
                goID <- names(PltTst)[i]
                goID2 <- gsub("^GO:", "GO", goID)
                nm <- AnnotationDbi::Term(goID)
                temp2[[goID]] <- factor(temp2[[goID]], levels = c("-", "+"))
                ttl <- paste0(ttl, ", ", goID, " ", nm)
                myColors3 <- setNames(c("lightgrey", "purple"), c("-", "+"))
                colScale3 <- scale_colour_manual(name = goID, values = myColors3)
                fillScale3 <- scale_fill_manual(name = goID, values = myColors3)
                catnm <- goID
                filt <- which(temp2[[catnm]] == "+")
              }
              ttl2 <- gsub("/|:|\\*|\\?|<|>|\\|", "-", ttl)
              m <- match(temp2$Category, c("-", "Contaminant", tstOrg3, "+", "In list"))
              temp2$LineType <- c(rep("dotted", 2), rep("dashed", length(tstOrg3)), rep("solid", 2))[m]
              #temp2$DotSize <- c(rep(0.2, 2), rep(0.3, length(tstOrg3)), rep(0.5, 2))[m]
              temp2$DotSize <- 1
              #temp2$Alpha <- c(rep(0.25, 2), rep(0.5, length(tstOrg3)), rep(1, 2))[m]
              temp2$Alpha <- 1
              Ngl <- c(0, 90)[(length(levels(temp2[[catnm]])) > 5) + 1]
              temp2 <- temp2[which(!is.na(temp2$Y)),]
              wTxt <- which(temp2$Sample == rev(levels(temp2$Sample))[1])
              frm <- as.formula(paste0("~`", catnm, "`"))
              suppressWarnings({
                plot <- ggplot(temp2, aes(text1 = `Protein Group`, text2 = Value, colour = Protein_group)) +
                  geom_line(aes(x = Sample, y = Y, group = id), alpha = 0.1) +
                  geom_point(aes(x = Sample, y = Y, size = DotSize)) + ggtitle(ttl) + ylab(kolnm) +
                  theme_bw() + facet_grid(frm) +
                  scale_size_identity(guide = "none") + scale_linetype_identity(guide = "none") +
                  theme(legend.position = "none",
                        axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1),
                        strip.text = element_text(face = "bold", size = 8, lineheight = 0.8, angle = Ngl),
                        strip.background = element_rect(fill = "lightblue", colour = "black", linewidth = 1)) +
                  scale_alpha_identity(guide = "none") + scale_colour_identity()
              })
              #poplot(plot, 12, 20)
              plotxt <- plot +
                geom_text(data = temp2[wTxt,], aes(label = `Protein Group`, x = Sample, y = Y, alpha = Alpha, color = id),
                          hjust = 0, cex = 2)
              #poplot(plotxt, 12, 20)
              setwd(SubDir)
              pth <- paste0(SubDir, "/", ttl2)
              plPath <- paste0(pth, ".html")
              plotlyProfiles <- ggplotly(plot, tooltip = c("text1", "text2"))
              # Grey plot for shiny app
              plotlyProfiles2 <- plot_ly(temp2)
              plotlyProfiles2 <- add_trace(plotlyProfiles2, x = ~Sample, y = ~Y,
                                           #split = ~`Protein Group`, # This would be correct but makes it slow, so this small approximation seems ok for now
                                           color = I("lightgrey"), type = "scatter",
                                           mode = "lines+markers", text = ~`Protein Group`, connectgaps = FALSE,
                                           name = "", showlegend = FALSE)
              ProfLy[[QuantType]] <- list(data = temp2,
                                          coloured = plotlyProfiles,
                                          grey = plotlyProfiles2)
              saveWidget(plotlyProfiles, paste0(pth, ".html"))
              tst <- try(saveWidget(partial_bundle(plotlyProfiles), paste0(pth, ".html")), silent = TRUE)
              if ((i == 1)&&(!"try-error" %in% class(tst))) { system(paste0("open \"", pth, ".html")) }
              #system(paste0("open \"",pth, ".html"))
              setwd(wd)
              evPlot <- plotEval(plotxt)
              #ggsave(paste0(pth, ".jpeg"), plotxt, dpi = 150)
              #ggsave(paste0(pth, ".pdf"), plotxt, dpi = 150)
              ggProf[[QuantType]] <- list(title = ttl,
                                          path = pth,
                                          plot = evPlot,
                                          dpi = 150,
                                          width = 10,
                                          height = 10)
              setwd(wd)
            }
          }
        }
      }
    }
  }
  PG$temp <- NULL
  invisible(clusterCall(parClust, function(x) { rm(list = ls());gc() }))
  saveFun(ggQuant, file = paste0(MainDir, "/ggQuantPlots.RData"))
  saveFun(ggProf, file = paste0(MainDir2, "/ggProfilePlots.RData"))
  saveFun(QuantLy, file = paste0(MainDir, "/QuantPlots.RData"))
  saveFun(ProfLy, file = paste0(MainDir2, "/ProfilePlots.RData"))
  for (QuantType in QuantTypes) { #QuantType <- QuantTypes[1]
    tst <- parSapply(parClust, ggQuant[[QuantType]], function(x) { #x <- ggQuant[[QuantType]][[1]]
      dr <- dirname(x$path)
      if (!dir.exists(dr)) { dir.create(dr, recursive = TRUE) }
      ggplot2::ggsave(paste0(x$path, ".pdf"), x$plot, dpi = x$dpi, width = x$width, height = x$height)
      ggplot2::ggsave(paste0(x$path, ".jpeg"), x$plot, dpi = x$dpi, width = x$width, height = x$height)
    })
  }
  tst <- parSapply(parClust, ggProf, function(x) {
    dr <- dirname(x$path)
    if (!dir.exists(dr)) { dir.create(dr, recursive = TRUE) }
    ggplot2::ggsave(paste0(x$path, ".pdf"), x$plot, dpi = x$dpi, width = x$width, height = x$height)
    ggplot2::ggsave(paste0(x$path, ".jpeg"), x$plot, dpi = x$dpi, width = x$width, height = x$height)
  })
  for (nm in names(ggProf)) {
    ReportCalls <- AddPlot2Report(Plot = ggProf[[nm]]$plot,
                                  Dir = dirname(ggProf[[nm]]$path),
                                  Title = ggProf[[nm]]$title)
  }
  setwd(wd)
}
