#' .plot_Rank_OR_Prof
#' 
#' Worker function used by the profile_and_rankedAbund_plots.R source to draw profile and ranked abundance plots.
#' Allows parallel drawing of all of these plots.
#'
#' @param ii Integer, index of the row of the samplesDF controlling what version of the plot is drawn.

.plot_Rank_OR_Prof <- function(ii) { #ii <- 1L #ii <- 2L #ii <- 5L #ii <- 7L #ii <- 13L #ii <- 18L #ii <- 23L
  quantType <- samplesDF$QuantType[ii]
  dataType <- samplesDF$type[ii]
  plotType <- samplesDF$subtype[ii]
  #tstReg <- ((quantType == "LFQ")&MakeRatios)
  ref <- samplesDF$ref[ii]
  kolnm <- gsub(" - $", "", ref)
  if (grepl("\\(Expr\\.\\)", kolnm)) { kolnm <- gsub("\\(Expr\\.\\)", " LFQ", kolnm) }
  if (grepl("^Sequence coverage ", kolnm)) { kolnm <- gsub("^Sequence coverage ", "Coverage ", kolnm) }
  smpls <- samplesDF$values[[ii]]
  toolTip <- paste0("text", as.character(1L:4L))
  tstPL <- FALSE
  myDPI <- 300L
  if (length(smpls) > 1L) { # Draw profile plot
    if (dataType == "PG") {
      subDir <- paste0(MainDir2, "/", quantType)
      varkol <- PG_varkol
      myData <- myPG
      yKol <- "Protein Group"
      colKol <- "Protein_group"
      txt4Kol <- "Peptides count"
      baseTtl <- paste0(quantType, " profiles")
    }
    if (dataType == "pep") {
      subDir <- paste0(MainDir2, "/Peptide ", quantType)
      varkol <- pep_varkol
      myData <- myPep
      yKol <- "Modified sequence"
      colKol <- "Peptidoform"
      txt4Kol <- "Proteins"
      baseTtl <- paste0("Pep. ", quantType, " profiles")
    }
    if (!dir.exists(subDir)) { dir.create(subDir, recursive = TRUE) }
    qKol <- grep(topattern(ref), colnames(myData), value = TRUE)
    if (dataType == "PG") {
      nKol <- gsub(topattern(ref), "Peptides count - ", qKol)
      myData2 <- reshape::melt(myData[, c(varkol, nKol)], id.vars = varkol)
    }
    #
    myData <- reshape::melt(myData[, c(varkol, qKol)], id.vars = varkol)
    colnames(myData)[which(colnames(myData) == "variable")] <- "Sample"
    myData$Sample <- gsub_Rep(topattern(ref), "", myData$Sample)
    if (scrptType == "withReps") {
      myData$Sample <- cleanNms(myData$Sample)
    }
    colnames(myData)[which(colnames(myData) == "value")] <- "Y"
    if (dataType == "PG") {
      myData[[txt4Kol]] <- myData2$value
    }
    # Process values
    if (dataType == "PG") {
      if ((quantType == "LFQ")&&(WorkFlow == "Band ID")) {
        kolnm <- "LFQ"
        myData$Y <- 10L^myData$Y
        myData <- myData[which(myData$Y > 0),]
      } else {
        if (quantType %in% c("Coverage", "Spectra")) {
          myData <- myData[which(myData$Y > 0),]
        }
      }
    }
    if (dataType == "pep") {
      if (quantType == "intensities") {
        #if (WorkFlow == "Band ID") { # Always log transform, otherwise it becomes difficult to read at peptides level!
        #  myData <- myData[which(myData$Y > 0),]
        #} else {
        kolnm <- "log10(intensity)"
        myData$Y <- log10(myData$Y)
        #}
      }
    }
    myData <- myData[which(is.all.good(myData$Y, 2L)),]
    myData[[kolnm]] <- myData$Y
    if (!nrow(myData)) { return(list(plotly_saved = FALSE,
                                     step = 1L)) }
    #
    profData <- myData
    goOn <- TRUE
    if (plotType == "All") {
      ttl <- baseTtl
      lev1 <- unique(profData$Category)
      lev1 <- lev1[which(!lev1 %in% c("In list", "Contaminant"))]
      lev1 <- union(lev1, c("In list", "Contaminant"))
      catnm <- "Category"
      profData$Category <- factor(profData$Category, levels = lev1)
      if ("Cluster" %in% colnames(profData)) {
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
      Ngl <- c(0, 90)[(length(levels(profData[[catnm]])) > 5L) + 1L]
    }
    if (plotType == "List") {
      ttl <- paste0(baseTtl, ", proteins of interest")
      profData$"In list" <- factor(profData$"In list", levels = c("-", "+"))
      catnm <- "In list"
      frm <- as.formula(paste0(".~`", catnm, "`"))
      myFacets <- facet_grid(frm)
      Ngl <- c(0, 90)[(length(levels(profData[[catnm]])) > 5L) + 1L]
    }
    if (plotType %in% c("GO", "Mark")) {
      if (plotType == "GO") {
        ttl <- paste0(baseTtl, ", GO terms of interest")
        myFlt2 <- GO_filter
      }
      if (plotType == "Mark") {
        ttl <- paste0(baseTtl, ", compartment markers")
        myFlt2 <- CompGOTerms
      }
      catnm <- "GO term"
      myFlt2 <- myFlt2[which(myFlt2 %in% colnames(profData))]
      goOn <- length(myFlt2)
      if (goOn) {
        profData <- lapply(myFlt2, \(go) { #go <- myFlt2[1L]
          dat <- profData[which(profData[[go]] == "+"),]
          if (!nrow(dat)) { return() }
          dat$"GO term" <- go
          dat <- dat[, which(!colnames(dat) %in% myFlt2)]
          return(dat)
        })
        profData <- do.call(rbind, profData)
        profData$"GO term" <- names(myFlt2)[match(profData$"GO term", myFlt2)]
        profData$"GO term" <- factor(profData$"GO term", levels = names(myFlt2))
        frm <- as.formula(paste0("~`", catnm, "`"))
        myFacets <- facet_wrap(frm)
        Ngl <- 0
      }
    }
    if (!goOn) {
      return(list(plotly_saved = FALSE,
                  step = 2L))
    }
    lvls <- mySamples
    if (scrptType == "withReps") {
      lvls <- cleanNms(mySamples, rep = " ")
    }
    profData$Sample <- factor(profData$Sample, levels = lvls)
    m <- match(profData$Category, c("-", "Contaminant", tstOrg2, "+", "In list"))
    profData$LineType <- c(rep("dotted", 2L), rep("dashed", length(tstOrg2)), rep("solid", 2L))[m]
    #profData$DotSize <- c(rep(0.2, 2L), rep(0.3, length(Org_Nms)), rep(0.5, 2L))[m]
    profData$DotSize <- 1L
    #profData$Alpha <- c(rep(0.25, 2L), rep(0.5, length(Org_Nms)), rep(1L, 2L))[m]
    profData$Alpha <- 1
    profData <- profData[which(!is.na(profData$Y)),]
    wTxt <- which(profData$Sample == rev(levels(profData$Sample))[1L])
    yLim <- c(min(profData$Y, na.rm = TRUE),
              max(profData$Y, na.rm = TRUE))
    # It's important to have fixed colours here (using scale_linetype_identity()) so we can consistently plot in the app
    # different protein subsets with always the same colors assigned to each.
    ggCall_txt <- "plot <- ggplot(profData, aes(x = Sample, y = Y, text1 = .data[[yKol]], text2 = .data[[kolnm]], text3 = PEP, text4 = .data[[txt4Kol]], colour = .data[[colKol]])) +
                  geom_line(aes(group = id), alpha = 0.1) +
                  geom_point(aes(size = DotSize)) + ggtitle(ttl) + ylab(kolnm) + myFacets +
                  theme_bw() + scale_size_identity(guide = \"none\") + scale_linetype_identity(guide = \"none\") +
                  theme(legend.position = \"none\",
                        axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1),
                        strip.text = element_text(face = \"bold\", size = 8L, lineheight = 0.8, angle = Ngl),
                        strip.background = element_rect(fill = \"lightblue\", colour = \"black\", linewidth = 1)) +
                  scale_alpha_identity(guide = \"none\") + scale_colour_identity()"
    suppressWarnings({
      eval(parse(text = ggCall_txt))
    })
    #poplot(plot, 12L, 22L)
    plot_txt <- plot +
      geom_text(data = profData[wTxt,], aes(label = .data[[yKol]], x = Sample, y = Y, alpha = Alpha, color = id),
                hjust = 0, cex = round(max(c(1, min(c(5, 100/length(wTxt))))), 1L))
    #poplot(plot_txt, 12L, 22L)
    if ((quantType == "LFQ")&&(plotType %in% c("GO", "Mark"))) {
      profData3 <- profData
      # Normalize by row
      protMeans <- aggregate(profData3$Y, list(profData3$`Leading protein IDs`), mean, na.rm = TRUE)
      m <- match(profData3$`Leading protein IDs`, protMeans$Group.1)
      profData3$Y <- profData3$Y - protMeans$x[m]
      # Calculate envelope
      profData3 <- aggregate(profData3$Y, list(profData3$Sample, profData3$`GO term`), \(x) {
        n <- length(x)
        m <- mean(x)
        ci <- qt(0.975, df = n-1L)*sd(x)/sqrt(n)
        c(m, m+ci, m-ci)
      })
      colnames(profData3)[1L:2L] <- c("Sample", "GO term")
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
      ggCall_txt3 <- gsub(", text1 = .+, colour = \\.data\\[\\[colKol\\]\\]\\)",
                          ", group = `GO term`, colour = `GO term`, fill = `GO term`)",
                          ggCall_txt3)
      ggCall_txt3 <- paste0("rib_", gsub("profData", "profData3", ggCall_txt3))
      suppressWarnings({
        eval(parse(text = gsub("  +", " ", ggCall_txt3)))
      })
      #poplot(rib_plot, 12L, 22L)
      pth3 <- paste0(subDir, "/", ttl3)
      if (plotType == "GO") {
        nm <- "GO trends"
      }
      if (plotType == "Mark") {
        nm <- "Compartment trends"
      }
      ggsave(paste0(pth3, ".jpeg"), rib_plot, dpi = 150L, width = 13L, height = 10L)
      ggsave(paste0(pth3, ".pdf"), rib_plot, dpi = 150L, width = 13L, height = 10L)
      myRes <- list(path = pth3,
                    title = ttl3,
                    ggCall = ggCall_txt3,
                    ggPlot = plotEval(rib_plot),
                    dpi = 150L,
                    width = 13L,
                    height = 10L)
    } else {
      plotlyCall_txt <- "plotlyProfiles <- ggplotly(plot, tooltip = toolTip)"
      eval(parse(text = plotlyCall_txt))
      pth <- paste0(subDir, "/", gsub("/|:|\\*|\\?|<|>|\\|", "-", ttl))
      plPath <- paste0(pth, ".html")
      currWD <- getwd()
      plPath2 <- paste0(currWD, "/", basename(plPath))
      if (file.exists(plPath)) { unlink(plPath) }
      if (file.exists(plPath2)) { unlink(plPath2) }
      tstPL <- try(saveWidget(partial_bundle(plotlyProfiles), plPath2), silent = TRUE)
      if (inherits(tstPL, "try-error")) { tstPL <- try(saveWidget(plotlyProfiles, plPath2), silent = TRUE) }
      if ((!inherits(tstPL, "try-error"))&&(file.exists(plPath2))) {
        tstPL <- file.rename(plPath2, plPath)
      } else { tstPL <- FALSE }
      ggsave(paste0(pth, ".jpeg"), plot, dpi = myDPI/2, width = 13L, height = 10L)
      ggsave(paste0(pth, ".pdf"), plot, width = 13L, height = 10L)
      ggsave(paste0(pth, "_lab.jpeg"), plot_txt, dpi = myDPI, width = 13L, height = 10L)
      ggsave(paste0(pth, "_lab.pdf"), plot_txt, width = 13L, height = 10L)
      evPlot <- plotEval(plot)
      evPlot_txt <- plotEval(plot_txt)
      myRes <- list(path = pth,
                    title = ttl,
                    ggCall = ggCall_txt,
                    ggPlot = evPlot_txt,
                    ggPlot_no_text = evPlot,
                    dpi = 150L,
                    width = 13L,
                    height = 10L,
                    data = profData,
                    data2 = list(colName = kolnm,
                                 yCol = yKol,
                                 color = colKol,
                                 labCol = txt4Kol,
                                 facets = myFacets,
                                 category = catnm,
                                 angle = Ngl),
                    yRange = yLim,
                    plotly_saved = tstPL,
                    plotly = plotlyProfiles,
                    plCall = plotlyCall_txt)
    }
  } else { # Draw ranked abundance plot
    smpl <- smpls; rm(smpls) # For clarity
    myKol <- paste0(ref, smpl)
    if (dataType == "PG") {
      subDir <- paste0(MainDir, "/", quantType)
      varkol <- PG_varkol
      xKol <- "Protein Group"
      myData <- myPG[, c(myKol, varkol)]
      txt4Kol <- "Peptides count"
      myData[[txt4Kol]] <- myPG[[paste0("Peptides count - ", smpl)]]
      varkol <- union(varkol, txt4Kol)
    }
    if (dataType == "pep") {
      subDir <- paste0(MainDir, "/Peptide ", quantType)
      varkol <- pep_varkol
      xKol <- "Modified sequence"
      txt4Kol <- "Proteins"  
      myData <- myPep[, c(myKol, varkol)]
    }
    if (!dir.exists(subDir)) { dir.create(subDir, recursive = TRUE) }
    if (!myKol %in% colnames(myData)) { return(list(plotly_saved = FALSE,
                                                    step = 1L)) }
    # if (tstReg) {
    #   rgKol <- paste0("Regulated - ", smpl)
    #   tstReg <- rgKol %in% colnames(myData)
    # }
    # if (tstReg) {
    #   myData2 <- myData[, c(rgKol, varkol)]
    #   colnames(myData2)[1L] <- gsub(topattern(ref), "", colnames(myData2)[1L])
    #   myData2 <- reshape::melt(myData2, id.vars = varkol)
    #   colnames(myData2) <- c(varkol, "Sample", "Reg")
    # }
    colnames(myData)[1L] <- gsub(topattern(ref), "", colnames(myData)[1L])
    myData <- reshape::melt(myData, id.vars = varkol)
    colnames(myData) <- c(varkol, "Sample", "Y")
    # if (tstReg) {
    #   myData$Regulated <- myData2$Reg
    #   rm(myData2)
    #   spcFlt <- grep("^Specific", myData$Regulated)
    #   if (length(spcFlt)) {
    #     myDataSp <- myData[spcFlt,]
    #     myDataSp$Category <- "Specific"
    #   }
    # }
    #
    # Process values
    if (dataType == "PG") {
      if ((quantType == "LFQ")&&(WorkFlow == "Band ID")) {
        kolnm <- "LFQ"
        myData$Y <- 10L^myData$Y
        myData <- myData[which(myData$Y > 0),]
      } else {
        if (quantType %in% c("Coverage", "Spectra")) {
          myData <- myData[which(myData$Y > 0),]
        }
      }
    }
    if (dataType == "pep") {
      if (quantType == "intensities") {
        #if (WorkFlow == "Band ID") { # Always log transform, otherwise it becomes difficult to read at peptides level!
        #  myData <- myData[which(myData$Y > 0),]
        #} else {
        kolnm <- "log10(intensity)"
        myData$Y <- log10(myData$Y)
        #}
      }
    }
    myData <- myData[which(is.all.good(myData$Y, 2L)),]
    if (!nrow(myData)) { return(list(plotly_saved = FALSE,
                                     step = 3L)) }
    #
    myData <- myData[order(myData$Y, decreasing = TRUE),]
    nrws <- nrow(myData)
    myData[[xKol]] <- factor(myData[[xKol]], levels = unique(myData[[xKol]]))
    myData[[kolnm]] <- myData$Y
    #
    if (scrptType == "withReps") {
      smpl2 <- cleanNms(smpl)
    }
    if (scrptType == "noReps") {
      smpl2 <- smpl
    }
    if (dataType == "PG") {
      ttl <- paste0("PG ranked by ", gsub(" *\\[.*", "", kolnm), " - ", smpl2)
    }
    if (dataType == "pep") {
      ttl <- paste0("Pep. ranked by ", gsub(" *\\[.*", "", kolnm), " - ", smpl2)
    }
    if (plotType == "All") {
      catnm <- "Category"
      myData <- myData[c(which(myData$Category == "Contaminant"),
                         which(!myData$Category %in% c("In list", "Contaminant")),
                         which(myData$Category == "In list")),]
      txtFilt <- 1L:nrws
    }
    if (plotType == "List") {
      ttl <- paste0(ttl, ", proteins of interest")
      myData$"In list" <- factor(myData$"In list", levels = c("-", "+"))
      catnm <- "In list"
      myData <- myData[order(myData$`In list`),]
      txtFilt <- which(myData[[catnm]] == "+")
    }
    # if (plotType %in% GO_filter) {
    #   goID <- plotType
    #   goID2 <- gsub("^GO:", "GO", goID)
    #   nm <- AnnotationDbi::Term(goID)
    #   myData[[goID2]] <- factor(myData[[goID]], levels = c("-", "+"))
    #   ttl <- paste0(ttl, ", ", goID, " ", nm)
    #   myColors3 <- setNames(c("lightgrey", "purple"), c("-", "+"))
    #   colScale3 <- scale_colour_manual(name = goID, values = myColors3)
    #   fillScale3 <- scale_fill_manual(name = goID, values = myColors3)
    #   catnm <- goID
    #   myData <- myData[order(myData[[catnm]]),]
    #   txtFilt <- which(myData[[catnm]] == "+")
    # }
    Sz <- max(c(1.5, min(c(0.01, round(600/nrow(myData), 2L)))))
    myData$xPos <- as.numeric(myData[[xKol]]) # not integer, we will add non integer shifts   
    wLst <- c()
    if ("Category" %in% names(myData)) {
      wLst <- which(myData$Category == "In list")
    }
    if ("In list" %in% names(myData)) {
      wLst <- which(myData$"In list" == "+")
    }
    wRst <- which(!1L:nrws %in% wLst)
    #
    pth <- paste0(subDir, "/", gsub("/|:|\\*|\\?|<|>|\\|", "-", ttl))
    #
    intmin <- min(myData$Y)
    intmax <- max(myData$Y)
    intscale <- intmax-intmin
    xmax <- nrws*5/4
    nSteps <- max(c(1, length(myFlt)))
    xStep <- xmax*0.033/nSteps
    yStep <- (intmax-intmin)*0.1/nSteps
    xSpan <- c(1-xStep*nSteps, xmax)
    ySpan <- c(intmin-yStep*nSteps, intmax*1.2+intscale*0.05)
    ySpan[1L] <- ySpan[1L]-(ySpan[2L]-ySpan[1L])*0.2
    #
    # if ((tstReg)&&(length(spcFlt))) {
    #   myDataSp$Y <- myDataSp$Y + intmax*0.01 # To offset the markers
    # }
    #
    plot <- ggplot() + theme_bw() + ylab(kolnm) + xlab(xKol) + scale_size_identity() +
      ggtitle(ttl,
              subtitle = paste0(length(unique(myData$id)), " ", tolower(xKol), "s"))  +
      theme(panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.line.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks = element_blank(),
            plot.margin = margin(r = 100L))
    #poplot(plot, 12L, 22L)
    if ((dataType == "PG")&&(length(myFlt))) {
      wGO <- setNames(lapply(myFlt, \(x) { which(myData[[x]] == "+") }), myFlt)
      GOdat <- lapply(seq_along(myFlt), \(i) { #i <- 1L
        x <- myData[wGO[[myFlt[i]]],]
        x$xPos <- x$xPos-i*xStep
        x$Y <- x$Y-i*yStep
        x$GO_term <- myGOcolors[myFlt[i]]
        x$Compartment <- names(myFlt)[i]
        return(x)
      })
      GOdat <- do.call(rbind, GOdat)
      myFltDF <- data.frame(GO_term = gsub(" *\\[.*", "", names(myFlt)),
                            xPos = nrws - xStep*(1L:nSteps+0.5),
                            Y = intmin - yStep*(1L:nSteps+0.5))
      plot <- plot +
        geom_segment(data = GOdat,
                     aes(x = xPos, xend = xPos-xStep, y = Y, yend = Y-yStep, text1 = .data[[xKol]], text2 = Compartment, color = GO_term)) +
        scale_color_identity()
      #poplot(plot, 12L, 22L)
    }
    #poplot(plot, 12L, 22L)
    # Here we could add a switch to choose between different geom:
    # geom_bar(stat = "identity", aes(x = xPos, y = Y, fill = .data[[catnm]], text1 = .data[[xKol]], text2 = .data[[kolnm]], text3 = PEP, text4 = .data[[txt4Kol]])) +
    # geom_tile(aes(x = xPos, y = Y, fill = .data[[catnm]], text1 = .data[[xKol]], text2 = .data[[kolnm]], text3 = PEP, text4 = .data[[txt4Kol]]), height = intscale*0.001) +
    plot <- plot +
      geom_point(data = myData[wRst,],
                 aes(x = xPos, y = Y, fill = .data[[catnm]], text1 = .data[[xKol]], text2 = .data[[kolnm]], text3 = PEP, text4 = .data[[txt4Kol]]),
                 size = Sz, shape = 21L, colour = "transparent")
    #poplot(plot, 12L, 22L)
    if (length(wLst)) {
      plot <- plot +
        geom_point(data = myData[wLst,],
                   aes(x = xPos, y = Y, fill = .data[[catnm]], text1 = .data[[xKol]], text2 = .data[[kolnm]], text3 = PEP, text4 = .data[[txt4Kol]]),
                   size = Sz*1.5, shape = 23L, colour = "transparent")
    }
    #poplot(plot, 12L, 22L)
    if (plotType == "All") { plot <- plot + fillScale }
    if (plotType == "List") { plot <- plot + fillScale2 }
    #if (plotType %in% GO_filt) { plot <- plot + fillScale3 }
    #poplot(plot, 12L, 22L)
    # if ((tstReg)&&(length(spcFlt))) {
    #   plot <- plot +
    #     geom_point(data = myDataSp, aes(xPos, Y, fill = .data[[catnm]]), shape = 21L)
    # }
    #poplot(plot, 12L, 22L)
    if (plotType == "All") {
      plot_ly <- ggplotly(plot, tooltip = toolTip)
      # a folder with external resources is created for each html plot!
      plPath <- paste0(pth, ".html")
      currWD <- getwd()
      plPath2 <- paste0(currWD, "/", basename(plPath))
      if (file.exists(plPath)) { unlink(plPath) }
      if (file.exists(plPath2)) { unlink(plPath2) }
      tstPL <- try(saveWidget(partial_bundle(plot_ly), plPath2), silent = TRUE)
      if (inherits(tstPL, "try-error")) { tstPL <- try(saveWidget(plot_ly, plPath2), silent = TRUE) }
      if ((!inherits(tstPL, "try-error"))&&(file.exists(plPath2))) {
        tstPL <- file.rename(plPath2, plPath)
      } else { tstPL <- FALSE }
    }
    if ((dataType == "PG")&&(length(myFlt))) {
      plot <- plot +
        geom_text(data = myFltDF, aes(x = xPos, y = Y, label = GO_term), size = 1.5, angle = -66, hjust = 0)
      #poplot(plot, 12L, 22L)
    }
    plot <- plot + ggnewscale::new_scale_color()
    if (plotType == "All") { plot <- plot + fillScale + colScale }
    if (plotType == "List") { plot <- plot + fillScale2 + colScale2 }
    #if (plotType %in% GO_filt) { plot <- plot + colScale3 }
    plot <- plot +
      coord_cartesian(xlim = xSpan, ylim = ySpan) +
      geom_text(data = myData[txtFilt,], angle = 45, hjust = 0, cex = 2.5,
                aes(xPos, Y + intscale*0.025, colour = .data[[catnm]], label = .data[[xKol]]))
    #poplot(plot, 12L, 22L)
    ggsave(paste0(pth, ".jpeg"), plot, dpi = myDPI, width = 25L, height = 10L)
    ggsave(paste0(pth, ".pdf"), plot, dpi = myDPI, width = 25L, height = 10L)
    evPlot <- plotEval(plot)
    myRes <- list(title = ttl,
                  path = pth,
                  ggPlot = evPlot,
                  dpi = myDPI,
                  width = 25L,
                  height = 10L,
                  plotly_saved = tstPL,
                  data2 = list(colName = kolnm,
                               xCol = xKol,
                               labCol = txt4Kol,
                               category = catnm))
    if (tstPL) {
      myRes$plotly_path <- plPath
      myRes$plotly <- plot_ly
    }
  }
  return(myRes)
}
