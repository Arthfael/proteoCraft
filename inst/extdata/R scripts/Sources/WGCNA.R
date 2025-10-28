##################
# WGCNA analysis #
##################
#
# WGCNA = Weighted Gene Correlation Network Analysis
#
# This implementation is based on the tutorial:
# https://fuzzyatelin.github.io/bioanth-stats/module-F21-Group1/module-F21-Group1.html
#
if (!require(WGCNA)) { pak::pak("WGCNA") }
library(WGCNA)
#
L <- length(RSA$values)
if (L <= 15) { warning("WGCNA is not recommended below 15 samples!") } else {
  if (L <= 30) { warning("Low number of samples, WGCNA modules reliability may be poor !") }
}
#
# Global parameters
hubThresh <- 0.8
sftThresh <- 0.8
repThresh <- 3
repProp <- 0.75
corThresh <- 0.5
pvalThresh <- 0.05
mergeThresh <- 0.25
#
wgcnaDir <- paste0(wd, "/WGCNA")
wgcnaDirs <- paste0(wgcnaDir, "/", c("Initial",
                                     "Traits",
                                     "Modules"))
if (scrptType == "withReps") { dirlist <- unique(c(dirlist, wgcnaDir, wgcnaDirs)) }
for (dr in wgcnaDirs) {
  if (!dir.exists(dr)) { dir.create(dr, recursive = TRUE) }
}
#
# Prepare data
# We will filter based on our own criteria:
# NB: really important to filter out all low quality data, otherwise the 
col <- paste0(Prot.Expr.Root, RSA$values)
minReps <- max(c(repThresh, round(length(Rep)*repProp))) # Let's be stringent and ask for 3 replicates
exprData <- as.data.frame(t(PG[, col]))
colnames(exprData) <- PG$Label
rownames(exprData) <- RSA$values
tmp <- Exp.map[match(rownames(exprData), Exp.map$Ref.Sample.Aggregate), VPAL$column]
m <- setNames(lapply(VPAL$value, function(x) { which(tmp == x) }), VPAL$value)
tmp2 <- !is.na(exprData)
tst <- setNames(lapply(m, function(x) {
  colSums(tmp2[x,])
}), VPAL$value)
tst <- do.call(cbind, tst)
tst <- as.data.frame(tst)
rownames(exprData) <- cleanNms(rownames(exprData))
whWGCNA <- which(apply(tst, 1, min) >= minReps)
exprData <- exprData[, whWGCNA]
exprMean <- rowMeans(PG[whWGCNA, col], na.rm = TRUE)
#
# Visualize data
sampleTree <- hclust(dist(exprData), method = "average") #Clustering samples based on distance 
#Plotting the cluster dendrogram
# Original code using base functions:
# par(cex = 0.6)
# par(mar = c(0,4,2,0))
# plot(sampleTree, main = "Sample clustering to detect outliers",
#      sub = "", xlab = "", cex.lab = 1.5,
#      cex.axis = 1.5, cex.main = 2)
# #draw on line to show cutoff height
# abline(h = max(sampleTree$height)*0.95, col = "red")
# dev.off()
# ggdendro option; problem: branch lengths are lost when it converts input from class "hclust" to "dendrogram"
# library(ggplot2)
# library(ggdendro)
# ggPlot <- ggdendrogram(sampleTree, rotate = FALSE, size = 1.5) +
#   ggtitle("Sample clustering to detect outliers") +
#   ylab("Height")
# poplot(ggPlot, 12, 22)
# Solution: fix branch lengths (they were originally arbitrary anyway)
L <- max(sampleTree$height)*0.025
sampleTree2 <- as.dendrogram(sampleTree)
dend_data <- ggdendro::dendro_data(sampleTree2, type = "rectangle")
w <- which(dend_data$segments$yend == 0)
dend_data$segments$yend[w] <- dend_data$segments$y[w] - L
dend_data$segments$label <- NA
dend_data$segments$label[w] <- dend_data$labels$label[match(dend_data$segments$x[w],
                                                            dend_data$labels$x)]
xMx <- max(dend_data$segments$x)
yMn <- min(c(dend_data$segments$y, dend_data$segments$yend))
yMx <- max(c(dend_data$segments$y, dend_data$segments$yend))
yMnRnd <-floor(yMn-L*0.1)
yMxRnd <- ceiling(yMx)
yLab <- data.frame(pos = yMnRnd:yMxRnd)
ttl <- "Samples clustering"
ggPlot <- ggplot(dend_data$segments) + theme_bw() +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_text(data = dend_data$segments[w,],
             aes(label = label, x = x, y = yend-L*0.1), angle = 90, hjust = 1, cex = 3) +
  geom_segment(y = yMx * 0.95, yend = yMx * 0.95, x = 1, xend = xMx, color = "red") +
  geom_segment(x = 0, xend = 0, y = yMnRnd, yend = yMxRnd) +
  geom_segment(data = yLab, x = -0.2, xend = 0, aes(y = pos, yend = pos)) +
  geom_text(data = yLab, x = -0.3, aes(y = pos, label = pos), hjust = 1, cex = 2) +
  labs(title = ttl, y = "Height", x = "") +
  xlim(-0.2, xMx + 1) + ylim(yMnRnd-0.5, yMxRnd+0.25) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
poplot(ggPlot, 6, 11)
suppressMessages({
  ggsave(paste0(wgcnaDirs[1], "/", ttl, ".jpeg"), ggPlot, dpi = 300)
  ggsave(paste0(wgcnaDirs[1], "/", ttl, ".pdf"), ggPlot, dpi = 300)
})
#
# Unlike in the tutorial, we will NOT remove outliers:
# We consider that these would already have been removed by then if necessary.
#
require(parallel)
try(stopCluster(parClust), silent = TRUE) # Free-up connections if necessary...
tst <- try(enableWGCNAThreads(N.clust), silent = TRUE) # Seems to be working based on plotting in tasks manager
if ("try-error" %in% class(tst)) {
  open_conns <- showConnections()
  sock_conns <- as.integer(rownames(open_conns[grep("sockconn", open_conns[,"class"]), , drop = FALSE]))
  for (i in sock_conns) try(close(getConnection(i)), silent = TRUE)
  tst <- try(enableWGCNAThreads(N.clust), silent = TRUE)
  if ("try-error" %in% class(tst)) {
    closeAllConnections()
    tst <- try(enableWGCNAThreads(N.clust), silent = TRUE)
  }
}
# Soft threshold
spt <- pickSoftThreshold(exprData) # Rather slow, although internally parallelized...
# Scale independence vs Mean connectivity
scaleFact <- max(spt$fitIndices$SFT.R.sq)/max(spt$fitIndices$mean.k.)
#mxSFT.R.sq <- max(spt$fitIndices$SFT.R.sq)
ttl <- "Scale independence vs Mean connectivity"
ggPlot <- ggplot(spt$fitIndices) +
  geom_text(aes(label = Power, x = Power, y = SFT.R.sq), color = "darkgreen") +
  geom_line(aes(x = Power, y = SFT.R.sq), color = "darkgreen", linetype = "dashed") +
  geom_text(aes(label = Power, x = Power, y = mean.k.*scaleFact), color = "purple") +
  geom_line(aes(x = Power, y = mean.k.*scaleFact), color = "purple", linetype = "dashed") +
  scale_y_continuous(name = "Scale independence",
                     sec.axis = sec_axis(transform = ~./scaleFact, name = "Mean connectivity")) +
  geom_hline(yintercept = #mxSFT.R.sq*
               sftThresh, color = "red") +
  ggtitle(ttl) + xlab("Soft Threshold (power)") +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.text.x = element_blank(),
                     axis.ticks.y = element_blank(),
                     axis.title.y = element_text(color = "darkgreen"),
                     axis.title.y.right = element_text(color = "purple"))
#poplot(ggPlot, 6, 11)
suppressMessages({
  ggsave(paste0(wgcnaDirs[1], "/", ttl, ".jpeg"), ggPlot, dpi = 300)
  ggsave(paste0(wgcnaDirs[1], "/", ttl, ".pdf"), ggPlot, dpi = 300)
})
pwrEst <- spt$powerEstimate
if (is.na(pwrEst)) { warning("Data is too low quality, skipping...") } else {
  # Shiny app to select threshold
  if (!require(shiny)) { pak::pak("shiny") }
  if (!require(shinyWidgets)) { pak::pak("shinyWidgets") }
  if (!require(shinyjs)) { pak::pak("shinyjs") }
  library(shiny)
  library(shinyWidgets)
  library(shinyjs)
  ggPlot2 <- ggPlot + geom_vline(xintercept = pwrEst, color = "red")
  appNm <- "Select WGCNA threshold"
  if (exists("appRunTest")) { rm(appRunTest) }
  ui <- fluidPage(
    useShinyjs(),
    setBackgroundColor( # Doesn't work
      color = c(#"#F8F8FF",
        "#E5EDE1"),
      gradient = "linear",
      direction = "bottom"
    ),
    extendShinyjs(text = jsToggleFS, functions = c("toggleFullScreen")),
    tags$head(tags$style(HTML("table {table-layout: fixed;"))), # So table widths can be properly adjusted!
    titlePanel(tag("u", appNm),
               appNm),
    br(),
    h5(tags$div(
      "Select the lowest power for which R^2 is above 0.8 - or stars to plateau, keeping mean connnectivity in a reasonable range (10-100),", tags$br(),
      "then click \"Save\".")),
    uiOutput("Current"),
    br(),
    actionBttn("saveBtn", "Save", icon = icon("save"), color = "success", style = "pill"),
    br(),
    sliderInput("Power", "Power", 1, max(spt$fitIndices$Power), pwrEst, 1, TRUE, width = "400px"),
    br(),
    plotOutput("IndepConnect", height = paste0(screenRes$width*0.4, "px")),
    br(),
  )
  server <- function(input, output, session) {
    POWER <- reactiveVal(pwrEst)
    updtPlot <- function(reactive = TRUE) {
      if (reactive) {
        pwr <- input$Power
        updtPlot <- POWER() != pwr
      } else {
        pwr <- pwrEst
        updtPlot <- TRUE
      }
      if (updtPlot) {
        POWER(pwr)
        ggPlot2 <<- ggPlot + geom_vline(xintercept = pwr, color = "red")
      }
      return(renderPlot({ ggPlot2 }))
    }
    updtVal <- function(reactive = TRUE) {
      if (reactive) {
        pwr <- input$Power
      } else {
        pwr <- pwrEst
      }
      m <- match(pwr, spt$fitIndices$Power)
      return(renderUI({ list(list(em("Current values:"), br(),
                                  em(paste0("- R^2 = ", signif(spt$fitIndices$SFT.R.sq[m], 3))), br(),
                                  em(paste0("- Mean connectivity = ", signif(spt$fitIndices$mean.k.[m], 3))), br())) }))
    }
    output$IndepConnect <- updtPlot(FALSE)
    output$Current <- updtVal(FALSE)
    spt$fitIndices
    observeEvent(input$Power, {
      output$IndepConnect <- updtPlot()
      output$Current <- updtVal()
    })
    observeEvent(input$saveBtn, {
      pwrEst <<- input$Power
      suppressWarnings({
        ggsave(paste0(wgcnaDirs[1], "/", ttl, ".jpeg"), ggPlot2, dpi = 300)
        ggsave(paste0(wgcnaDirs[1], "/", ttl, ".pdf"), ggPlot2, dpi = 300)
      })
      appRunTest <<- TRUE
      stopApp()
    })
    #observeEvent(input$cancel, { stopApp() })
    session$onSessionEnded(function() { stopApp() })
  }
  runKount <- 0
  while ((!runKount)||(!exists("appRunTest"))) {
    eval(parse(text = runApp), envir = .GlobalEnv)
    runKount <- runKount+1
  }
  #
  # Adjacency
  adjacency <- adjacency(exprData, power = pwrEst)
  rownames(adjacency) <- colnames(adjacency) <- colnames(exprData)
  #
  # Module construction
  TOM <- TOMsimilarity(adjacency) # TOM stands for Topological Overlap Matrix
  rownames(TOM) <- colnames(TOM) <- colnames(exprData)
  # (this is the rate limiting step!!!)
  TOMdissimil <- 1-TOM
  PGsTree <- hclust(as.dist(TOMdissimil), method = "average")
  # Plot
  ttl <- "Protein groups clustering on TOM-based dissimilarity"
  plot(PGsTree, xlab = "", sub = "",
       main = "Protein groups clustering on TOM-based dissimilarity", 
       labels = FALSE, hang = 0.04)
  # Identify modules
  Modules <- cutreeDynamic(dendro = PGsTree,
                           distM = TOMdissimil,
                           deepSplit = 2,
                           pamRespectsDendro = FALSE,
                           minClusterSize = 30)
  names(Modules) <- PGsTree$labels
  #table(Modules)
  modColors <- labels2colors(Modules) #assigns each module number a color
  #table(modColors) # returns the counts for each color (aka the number of PGs within each module)
  plotDendroAndColors(PGsTree, modColors,"Module",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = "Protein groups dendrogram and module colors")
  # Record plot
  p <- recordPlot()
  dev.off()
  # Replay plot to devices
  pdf(paste0(wgcnaDirs[1], "/", ttl, ".pdf"))
  replayPlot(p)
  dev.off()
  jpeg(paste0(wgcnaDirs[1], "/", ttl, ".jpeg"), quality = 100, width = 1600, height = 1200)
  replayPlot(p)
  dev.off()
  #
  # Identify each module's Eigengene ("Well, akshually, it's an Eigen-Protein-Group!")
  MElist <- moduleEigengenes(exprData, colors = modColors) 
  MEs <- MElist$eigengenes 
  #head(MEs)
  #
  MEdissimil <- 1-cor(MElist$eigengenes, use = "complete")
  METree <- hclust(as.dist(MEdissimil), method = "average") # Clustering Eigengenes 
  ttl <- "Cluster dendrogram"
  par(mar = c(0, 4, 2, 0),
      cex = 1)
  METree$labels <- substring(METree$labels, 3)
  plot(METree)
  abline(h = 0.25, col = "red") #a height of .25 corresponds to correlation of .75
  p <- recordPlot()
  dev.off()
  # Replay plot to devices
  pdf(paste0(wgcnaDirs[1], "/", ttl, ".pdf"))
  replayPlot(p)
  dev.off()
  jpeg(paste0(wgcnaDirs[1], "/", ttl, ".jpeg"), quality = 100, width = 1600, height = 1200)
  replayPlot(p)
  dev.off()
  #
  # Merge modules below default 0.25 line (75% similarity or more)
  merge <- mergeCloseModules(exprData, modColors, cutHeight = mergeThresh)
  mergedMEs <- merge$newMEs
  mergedModColors <- merge$colors # Update colors
  mergedModules <- match(mergedModColors, substring(colnames(mergedMEs), 3))
  mergedModHubs <- rep(FALSE, length(mergedModules))
  names(mergedModHubs) <- names(mergedModules) <- colnames(exprData)
  # From now on we only use these merged modules!
  #
  # New dendrogram
  ttl <- "Protein groups dendrogram and module colors for original and merged modules"
  plotDendroAndColors(PGsTree, cbind(modColors, mergedModColors), 
                      c("Original Module", "Merged Module"),
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = "Protein groups dendrogram and module colors for original and merged modules")
  p <- recordPlot()
  dev.off()
  # Replay plot to devices
  pdf(paste0(wgcnaDirs[1], "/", ttl, ".pdf"))
  replayPlot(p)
  dev.off()
  jpeg(paste0(wgcnaDirs[1], "/", ttl, ".jpeg"), quality = 100, width = 1600, height = 1200)
  replayPlot(p)
  dev.off()
  #
  # Relate samples to experimental factors
  w <- which(vapply(Factors, function(Fact) { length(FactorsLevels[[Fact]]) }, 1) > 1)
  myFact <- Factors[w]
  xpMap <- Exp.map[, c("Ref.Sample.Aggregate", myFact)]
  for (Fact in myFact) { # Values must be numeric but may be arbitrary (not quantitative)
    if (!is.numeric(xpMap[[Fact]])) {
      u <- unique(xpMap[[Fact]])
      if (length(u) == 2) {
        xpMap[[Fact]] <- as.integer(as.factor(xpMap[[Fact]]))
      } else {
        # One-hot encode if there are more than 2 levels
        for (i in u) {
          stopifnot(!i %in% colnames(xpMap))
          xpMap[[i]] <- c(0, 1)[(xpMap[[Fact]] == i) + 1]
        }
        xpMap[[Fact]] <- NULL
        myFact <- c(myFact[which(myFact != Fact)],
                    u)
      }
    }
  }
  myFact <- unname(myFact)
  xpMap$Sample <- cleanNms(Exp.map$Ref.Sample.Aggregate)
  xpMap <- xpMap[match(rownames(exprData), xpMap$Sample),]
  rownames(xpMap) <- xpMap$Sample
  xpMap$Ref.Sample.Aggregate <- NULL
  xpMap$Sample <- NULL
  # Define numbers of PGs and samples
  nPGs <- ncol(exprData)
  nSamples <- nrow(exprData)
  moduleTraitCor <- cor(mergedMEs, xpMap, use = "p") #p for pearson correlation coefficient 
  moduleTraitPVal <- corPvalueStudent(moduleTraitCor, nSamples) #calculate the p-value associated with the correlation
  # Identify for each trait the best modules
  traitModules <- setNames(lapply(myFact, function(Fact) {
    tst <- vapply(1:nrow(moduleTraitCor), function(x) {
      (abs(moduleTraitCor[x, Fact]) >= corThresh) & (moduleTraitPVal[x, Fact] <= pvalThresh)
    }, TRUE)
    substring(rownames(moduleTraitCor)[which(tst)], 3)
  }), myFact)
  traitModules <- traitModules[which(vapply(traitModules, length, 1) > 0)]
  #
  # Display correlations and their p-values
  textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                      signif(moduleTraitPVal, 1), ")", sep = "");
  dim(textMatrix) <- dim(moduleTraitCor)
  # Display the correlation values within a heatmap plot
  ttl <- "Module-exp. factor relationships"
  par(mar = c(12, 12, 4, 1))
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names(xpMap),
                 yLabels = substring(names(mergedMEs), 3),
                 ySymbols = names(mergedMEs),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.main = 1.6, cex.lab = 1.4, cex.text = 1.2,
                 zlim = c(-1,1),
                 main = ttl)
  p <- recordPlot()
  dev.off()
  # Replay plot to devices
  pdf(paste0(wgcnaDirs[1], "/", ttl, ".pdf"))
  replayPlot(p)
  dev.off()
  jpeg(paste0(wgcnaDirs[1], "/", ttl, ".jpeg"), quality = 100, width = 1600, height = 1200)
  replayPlot(p)
  dev.off()
  #
  # Let's check the relationship between our modules and our Traits
  MET <- MEs
  colnames(MET) <- substring(colnames(MET), 3)
  tmpXp <- xpMap[, myFact]
  colnames(tmpXp) <- paste0("[ ", colnames(tmpXp), " ]")
  MET <- orderMEs(cbind(MET, tmpXp))
  # Plot the relationships among the eigengenes and all traits
  ttl <- "Eigengenes vs exp. Factors"
  par(cex = 1.4, mar = c(1, 1, 1, 1))
  plotEigengeneNetworks(MET, "", marDendro = c(2, 6, 1, 2), marHeatmap = c(6, 6, 0, 2),
                        cex.lab = 0.8, xLabelsAngle = 90)
  p <- recordPlot()
  dev.off()
  # Replay plot to devices
  pdf(paste0(wgcnaDirs[1], "/", ttl, ".pdf"))
  replayPlot(p)
  dev.off()
  jpeg(paste0(wgcnaDirs[1], "/", ttl, ".jpeg"), quality = 100, width = 1600, height = 1200)
  replayPlot(p)
  dev.off()
  #
  shScale <- scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 24))
  if (length(traitModules)) {
    # Calculate modules membership and connectivity - sometimes called kME
    modNames <- substring(names(mergedMEs), 3)
    PGmodMembership <- signedKME(exprData, mergedMEs, outputColumnName = "")
    # Below, alternative, roughly equivalent way to create PGmodMembership:
    #PGmodMembership <- as.data.frame(cor(exprData, mergedMEs, use = "p"))
    # ... but signedKME() - according to chatGPT, applies the correct sign convention for signed networks,
    # optionally accounts for missing data, ensures column names match module colors,
    # and can handle networks with different sign types (signed, unsigned, signed hybrid).
    # At least in my examples I could verify that the outcome was identical.
    kIntraMod <- intramodularConnectivity(adjacency, mergedModColors)
    # Extract PGs which have high significance for each Trait (aka Factor)
    modMembPval <- as.data.frame(corPvalueStudent(as.matrix(PGmodMembership), nSamples))
    names(PGmodMembership) <- modNames
    names(modMembPval) <- paste0("Pval ", modNames)
    GSPvalues <- list()
    for (Fact in names(traitModules)) { #Fact <- names(traitModules)[2]
      dr <- paste0(wgcnaDirs[2], "/", Fact)
      if (scrptType == "withReps") { dirlist <- unique(c(dirlist, dr)) }
      if (!dir.exists(dr)) { dir.create(dr, recursive = TRUE) }
      # Isolate current factor from the others
      tmpFact <- as.data.frame(xpMap[[Fact]])
      names(tmpFact) <- Fact
      #Calculate the module membership and the associated p-values
      PGtraitSignificance <- as.data.frame(cor(exprData, tmpFact, use = "p"))
      GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(PGtraitSignificance), nSamples))
      # Calculate the significance and associated p-values
      names(PGtraitSignificance) <- paste0("GS.", names(tmpFact))
      names(GSPvalue) <- paste0("p.GS.", names(tmpFact))
      GSPvalues[[Fact]] <- GSPvalue
      #head(GSPvalue)
      #
      for (mod in traitModules[[Fact]]) { #mod <- traitModules[[Fact]][[1]]
        kl <- match(mod, modNames)
        modulePGs <- which(merge$colors == mod)
        # pick hubs
        hubs <- rownames(PGmodMembership)[which(PGmodMembership[, mod] > hubThresh)] # Value could be a parameter!
        mergedModHubs[match(hubs, names(mergedModHubs))] <- TRUE
        # Here we want to see a positive correlation (cor > 0.5) between module membership and significance
        tmp <- data.frame("X" = abs(PGmodMembership[modulePGs, kl]),
                          "Y" = abs(PGtraitSignificance[modulePGs, 1]),
                          "mean log10 expression" = exprMean[modulePGs],
                          "is a Hub?" = colnames(exprData)[modulePGs] %in% hubs,
                          check.names = FALSE)
        ttl <- paste0(Fact, " - ", mod, " module membership vs. significance")
        corXY <- WGCNA::cor(tmp$X, tmp$Y, use = "all.obs")
        displayAsZero = 1e-05
        if ((is.finite(corXY)) &&(abs(corXY) < displayAsZero)) { corXY = 0 }
        pvalXY <- signif(WGCNA::corPvalueStudent(corXY, sum(is.finite(tmp$X)&is.finite(tmp$Y))), 2)
        subttl <- paste0("cor. = ", signif(corXY, 2),
                         if (is.finite(corXY)) {
                           paste0(", p-val. = ", signif(pvalXY, 2))
                         } else {
                           ""
                         })
        tst <- sqrt(sum(hex2RGB(gplots::col2hex(mod))@coords^2)) < 0.2
        ggPlot <- ggplot(tmp) +
          geom_point(aes(x = X, y = Y, shape = `is a Hub?`, size = `mean log10 expression`),
                     color = c("black", "red")[tst + 1], fill = mod) +
          theme_bw() + ggtitle(ttl, subtitle = subttl) + shScale +
          xlab(paste("Module Membership in", mod, "module")) +
          ylab(paste0("Significance for factor \"", Fact, "\"")) +
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank())
        #poplot(ggPlot)
        suppressMessages({
          ggsave(paste0(dr, "/", ttl, ".jpeg"), ggPlot, dpi = 300)
          ggsave(paste0(dr, "/", ttl, ".pdf"), ggPlot, dpi = 300)          
        })
        #
        # save as list and do Preranked GSEA
        #
      }
      # Add all Factors to existing module eigengenes
      MET <- MEs
      colnames(MET) <- substring(colnames(MET), 3)
      tmpXp <- tmpFact
      colnames(tmpXp) <- paste0("[ ", colnames(tmpXp), " ]")
      MET <- orderMEs(cbind(MET, tmpXp))
      # Plot the relationships among the eigengenes and the trait to identify meta-modules:
      ttl <- paste0(Fact, " vs Eigengenes")
      par(cex = 0.9)
      plotEigengeneNetworks(MET, "", marDendro = c(0, 4, 1, 2), marHeatmap = c(5, 4, 1, 2),
                            cex.lab = 0.8, xLabelsAngle = 90)
      p <- recordPlot()
      dev.off()
      # Replay plot to devices
      pdf(paste0(dr, "/", ttl, ".pdf"))
      replayPlot(p)
      dev.off()
      jpeg(paste0(dr, "/", ttl, ".jpeg"), quality = 100, width = 1600, height = 1200)
      replayPlot(p)
      dev.off()
      #
      # -> Meta-modules:
      # = groups of modules with mutual correlations stronger than their correlation with the specified trait
      #
      ttl <- paste0(Fact, " - Eigengene dendrogram")
      par(cex = 1.0)
      plotEigengeneNetworks(MET, ttl, marDendro = c(0, 4, 2, 0),
                            plotHeatmaps = FALSE)
      p <- recordPlot()
      dev.off()
      # Replay plot to devices
      pdf(paste0(dr, "/", ttl, ".pdf"))
      replayPlot(p)
      dev.off()
      jpeg(paste0(dr, "/", ttl, ".jpeg"), quality = 100, width = 1600, height = 1200)
      replayPlot(p)
      dev.off()
      #
      # Plot the heatmap matrix
      ttl <- paste0(Fact, " - Eigengene adjacency heatmap")
      par(cex = 1.0, mar = c(1, 1, 1, 1))
      plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(5, 5, 2, 2),
                            plotDendrograms = FALSE, xLabelsAngle = 90)
      p <- recordPlot()
      dev.off()
      # Replay plot to devices
      pdf(paste0(dr, "/", ttl, ".pdf"))
      replayPlot(p)
      dev.off()
      jpeg(paste0(dr, "/", ttl, ".jpeg"), quality = 100, width = 1600, height = 1200)
      replayPlot(p)
      dev.off()
    }
    #
    if (exists("DatAnalysisTxt")) {
      if (GSEAmode == "standard") {
        DatAnalysisTxt <- paste0(DatAnalysisTxt,
                                 " Weighted Genes Correlation Networks Analysis was run using using package WGCNA.")
      }
    }
    # Perform GSEA analysis on each module of interest
    dataType <- "PG"
    GSEAmode <- "WGCNA"
    Src <- paste0(libPath, "/extdata/R scripts/Sources/GSEA.R")
    #rstudioapi::documentOpen(Src)
    source(Src, local = FALSE)
    #
    # Save results
    tmp <- data.frame("PG id" = PG$id[whWGCNA],
                      "Leading protein IDs" = PG$`Leading protein IDs`[whWGCNA],
                      "Common Name" = PG$`Common Name (short)`[whWGCNA],
                      "merged Module" = mergedModules,
                      "(original Module)" = Modules,
                      check.names = FALSE)
    tmp[, paste0("MM ", colnames(PGmodMembership))] <- PGmodMembership
    tmp[, paste0("Pval ", colnames(PGmodMembership))] <- PGmodMembership
    tmp[, colnames(kIntraMod)] <- kIntraMod
    tmp[, colnames(modMembPval)] <- modMembPval
    fwrite(tmp, paste0(wgcnaDir, "/WGCNA results.csv"), quote = FALSE, row.names = FALSE, eol = "\n", na = "NA")
  }
}
#
# Finally, since we killed our cluster, we recreate it for later
disableWGCNAThreads()
open_conns <- showConnections()
sock_conns <- as.integer(rownames(open_conns[grep("sockconn", open_conns[,"class"]), , drop = FALSE]))
for (i in sock_conns) try(close(getConnection(i)), silent = TRUE)
source(parSrc)
#
# Now, we have our modules, we may want to 