### Check P-values
#
print_gg %<o% function(p) { # Works for both ggplot and 
  grid::grid.newpage()
  grid::grid.draw(
    if (inherits(p, "ggplot")) { ggplot2::ggplot_gtable(ggplot2::ggplot_build(p)) } else { p } # already a gtable
  )
}
slim_plotly %<o% \(p) {
  p <- plotly::plotly_build(p)
  p$x$visdat <- NULL      # removes the source data closures
  p$x$cur_data <- NULL    # removes current data reference
  attr(p, "orig_gg") <- NULL   # drop embedded ggplot if present
  return(p)
}
source(parSrc)
if (dataType == "peptides") {
  pvalDir <- paste0(wd, "/Workflow control/Peptides/P-values")
  myData <- pep
  dataType <- "peptides"
  ratRef <- pep.rat.ref
  entityCol <- "peptide"
  descrTxt <- "Peptide:"
  rowsCol <- labCol <- "Modified sequence"
  appNm <- paste0(dtstNm, " - peptides P-values")
  appTtl <- "peptides - type of P-values"
}
if (dataType == "modPeptides") {
  pvalDir <- paste0(modDirs[1L], "/P-values")
  myData <- ptmpep
  entityCol <- "PG"
  ratRef <- ptms.ratios.ref
  descrTxt <- paste0(ptm, "modified peptide:")
  labCol <- "Name"
  rowsCol <- "Modified sequence"
  appNm <- paste0(dtstNm, " - ", ptm, "-mod. peptides P-values")
  appNm <- paste0(ptm, " - type of P-values")
}
if (dataType == "PG") {
  pvalDir <- paste0(wd, "/Workflow control/Protein groups/P-values")
  myData <- PG
  entityCol <- "PG"
  ratRef <- Prot.Rat.Root
  descrTxt <- "Protein:"
  labCol <- "Label"
  rowsCol <- "Leading protein IDs"
  appNm <- paste0(dtstNm, " - PG P-values")
  appTtl <- "PG - type of P-values"
}
if (!dir.exists(pvalDir)) { dir.create(pvalDir, recursive = TRUE) }
dirlist <- union(dirlist, pvalDir)
#
if ((dataType == "modPeptides") && (Ptm %in% names(PTMs_PVal_col))) {
  my_PVal_Col <- PTMs_PVal_col[[Ptm]]
} else {
  my_PVal_Col <- pvalue.col
}
my_PVal_Col <- my_PVal_Col[which(vapply(my_PVal_Col, \(type) { #type <- my_PVal_Col[1L]
  length(grep(topattern(type), colnames(myData))) > 0L
}, TRUE))]
stopifnot(length(my_PVal_Col) > 0L)
#
# Scatter plots:
pkol2 <- gsub(" - $", "", pval_Col)
whSingle <- which(!myContrasts$isDouble)
temp <- lapply(whSingle, \(i) { #i <- 1L
  nm <- myContrasts$Contrast[i]
  j <- paste0(pval_Col, nm)
  w <- which(j %in% colnames(myData))
  j <- j[w]
  if (!length(j)) { return(NA) }
  x <- myData[, j, drop = FALSE]
  colnames(x) <- pkol2[w]
  x[[entityCol]] <- myData[[labCol]]
  x$Contrast <- nm
  return(x)
})
temp <- temp[which(vapply(temp, is.data.frame, TRUE))]
temp <- plyr::rbind.fill(temp)
temp$Contr <- gsub_Rep(" - ", "\n- ", temp$Contrast)
tmpContr <- gsub(" - ", "\n- ", myContrasts$Contrast)
temp$Contr <- factor(temp$Contr, levels = tmpContr)
temp$Contrast <- factor(temp$Contrast, levels = myContrasts$Contrast)
kol <- setdiff(colnames(temp), c(entityCol, "Contrast", "Contr"))
temp <- temp[which(rowSums(temp[, kol], na.rm = TRUE) > 0),]
Comb <- gtools::combinations(length(kol), 2L, kol)
tmpFl <- tempfile(fileext = ".rds")
readr::write_rds(temp, tmpFl)
nr <- nrow(myContrasts)
#
clusterExport(parClust, list("Comb", "tmpFl", "pvalDir", "pval_Col", "gsub_Rep",
                             "plotEval", "slim_plotly", "entityCol", "my_PVal_Col", "nr"),
              envir = environment())
invisible(clusterCall(parClust, \() {
  require(ggplot2)
  require(viridis)
  require(scattermore)
  require(plotly)
  temp <- readr::read_rds(tmpFl)
  assign("temp", temp, envir = .GlobalEnv)
  return()
}))
unlink(tmpFl)
plotsList1 <- parLapply(parClust, 1L:nrow(Comb), \(i) { #i <- 1L
  X <- Comb[i, 1L]
  Y <- Comb[i, 2L]
  X2 <- names(my_PVal_Col)[match(paste0(X, " - "), my_PVal_Col)]
  Y2 <- names(my_PVal_Col)[match(paste0(Y, " - "), my_PVal_Col)]
  dat <- temp[, c(X, Y, entityCol, "Contrast", "Contr")]
  colnames(dat)[1L:2L] <- c("X", "Y")
  dat$"P-value, X axis" <- gsub_Rep(" -log10\\(Pvalue\\)$", "", Comb[i, 1L])
  dat$"P-value, Y axis" <- gsub_Rep(" -log10\\(Pvalue\\)$", "", Comb[i, 2L])
  ttl1 <- paste0("P-values scatter plot - ", X2, " VS ", Y2)
  ttl1a <- paste0(X2, " VS ", Y2, " (-log10)")
  Mx <- c(dat$X, dat$Y)
  Mx <- max(Mx[which(is.finite(Mx))])
  uX <- unique(dat$`P-value, X axis`)
  uY <- unique(dat$`P-value, Y axis`)
  myLim <- c(0, Mx)
  dat <- dat[which(is.finite(dat$X)&is.finite(dat$Y)),]
  plot1 <- ggplot(dat, aes(x = X, y = Y, colour = Contrast, text = PG)) +
    scale_color_viridis(begin = 0.4, end = 0.8, discrete = TRUE, option = "F") +
    #facet_grid(`P-value, Y axis`~Contrast+`P-value, X axis`) + # This was for when we saved all as one plot outside the parLapply
    #facet_grid(Contrast~`P-value, Y axis`) +
    facet_grid(Contr~., scales = "fixed") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.25) +
    coord_fixed(xlim = myLim,
                ylim = myLim) +
    ggtitle(ttl1a) + theme_bw() +
    scale_x_continuous(expand = c(0L, 0L)) +
    scale_y_continuous(expand = c(0L, 0L)) +
    xlab(uX) + ylab(uY)
  plot1ly <- ggplotly(plot1 + geom_point(size = 0.33) +
                        theme(strip.background = element_blank(),
                              strip.text = element_blank()),
                      tooltip = c("text", "X", "Y", "Contrast"),
                      width = 700L,
                      height = 100L+400L*nr)
  plot1 <- plot1 +
    geom_scattermore(pixels = c(1024L, 1024L), pointsize = 3.2) +
    theme(strip.text.y = element_text(angle = 0))
  Img1 <- paste0(pvalDir, "/", gsub(":", " - ", ttl1))
  w1 <- 5
  h1 <- w1*(length(unique(dat$Contrast)) + 0.5)/3
  suppressMessages({
    ggsave(paste0(Img1, ".jpeg"), plot1, dpi = 150L, width = w1,# height = h1,
           units = "in")
    ggsave(paste0(Img1, ".pdf"), plot1, dpi = 150L, width = w1,# height = h1,
           units = "in")
  })
  #system(paste0("open \"", Img1, ".jpeg\""))
  lst <- list(Plot = ggplot2::ggplot_gtable(ggplot2::ggplot_build(plot1)), # to print these use print_gg()
              Plotly = slim_plotly(plot1ly),
              Title = ttl1)
  return(lst)
})
invisible(clusterCall(parClust, \(x) { rm(list = ls());gc() }))
# i <- 1L
#print_gg(plotsList1[[i]]$Plot)
#vapply(plotsList, \(x) { x$Title }, "")
#
# P-values histogram:
nbin <- 40L
bd <- (0L:nbin)/nbin
pvalLvls <- gsub_Rep("('s)?( t-)?(test)? -log10\\(Pvalue\\) - $", "", pval_Col)
temp <- setNames(lapply(pval_Col, \(type) { #type <- pval_Col[1L]
  #
  # Estimate power at 0.05%
  kol <- grep(topattern(sub(" -log10\\(", " ", sub("\\) - $", " - ", type))), colnames(myData), value = TRUE)
  if (!length(kol)) { stop(type) }
  pvals <- unlist(myData[, kol, drop = FALSE])
  pvals <- pvals[which(is.finite(pvals))]
  power_est <- mean(pvals < 0.05)
  #
  # Data for plotting
  kol <- grep(topattern(type), colnames(myData), value = TRUE)
  if (!length(kol)) { stop(type) }
  temp <- myData[, kol, drop = FALSE]
  colnames(temp) <- cleanNms(gsub(topattern(type), "", colnames(temp)))
  temp <- dfMelt(temp)
  temp$value <- 10L^(-temp$value)
  temp <- temp[which(is.finite(temp$value)),]
  temp$Bin <- vapply(temp$value, \(x) { min(which(bd >= x))-1L }, 1L)
  res <- aggregate(temp$Bin, list(temp$variable, temp$Bin), length)
  colnames(res) <- c("Contrast", "Bin", "Count")
  grps <- unique(res$Contrast)
  res$Frequency <- NA_real_
  for (grp in grps) {
    w <- which(res$Contrast == grp)
    res$Frequency[w] <- res$Count[w]/sum(res$Count[w])
  }
  res$"P-value type" <- gsub_Rep("('s)?( t-)?(test)? -log10\\(Pvalue\\) - $", "", type)
  res$Low <- 0L
  return(list(Dat = res,
              pwrEst = power_est))
}), gsub("('s)?( t-)?(test)? -log10\\(Pvalue\\) - $", "", pval_Col))
pwrAnnot <- lapply(names(temp), \(x) { #x <- names(temp)[1L]
  data.frame("power estimate" = paste0(round(100L*temp[[x]]$pwrEst, 1L), "%"),
             "P-value type" = x,
             check.names = FALSE)
})
pwrAnnot <- do.call(rbind, pwrAnnot)
temp <- do.call(rbind, lapply(names(temp), \(x) { #x <- names(temp)[1L]
  temp[[x]]$Dat
}))
temp$"P-value type" <- factor(temp$"P-value type",
                              levels = pvalLvls)
pwrAnnot$`P-value type` <- factor(pwrAnnot$`P-value type`,
                                  pvalLvls)
temp$Contrast <- gsub_Rep(" - ", "\n- ", temp$Contrast)
temp$Contrast <- factor(temp$Contrast, levels = tmpContr)
pwrAnnot$Contrast <- factor(tmpContr[1L], levels = tmpContr)
ttl2 <- "P-values histogram"
plot2 <- ggplot(temp) +
  geom_rect(position = "identity",
            aes(xmin = (Bin-1L)/nbin, ymin = Low, xmax = Bin/nbin, ymax = Frequency,
                fill = `P-value type`),
            alpha = 0.5) +
  geom_text(data = pwrAnnot, label = "5% pVal pwr. est. =",
            x = 0.01, y = 0.01, hjust = 0, vjust = 0, color = "red", size = 4L) +
  geom_text(data = pwrAnnot, aes(label = `power estimate`),
            x = 0.9, y = 0.01, hjust = 1, vjust = 0, color = "red", size = 5L) +
  scale_fill_viridis(begin = 0.4, end = 0.8, discrete = TRUE, option = "G") +
  scale_x_continuous(expand = c(0L, 0L)) + scale_y_continuous(expand = c(0L, 0L)) +
  facet_grid(`P-value type`~Contrast) + ggtitle(ttl2) + theme_bw() +
  theme(strip.text.y = element_text(angle = -90, size = 12L),
        strip.text.x = element_text(size = 12L),
        axis.text = element_text(size = 10L))
#poplot(plot2, 12L, 22L)
Img2 <- paste0(pvalDir, "/", ttl2)
w2 <- ((length(whSingle)+1)*1.25)*2
h2 <- ((length(pval_Col)+0.2)*1.25)*2
suppressMessages({
  ggsave(paste0(Img2, ".jpeg"), plot2, dpi = 150L, width = w2, height = h2, units = "in")
  ggsave(paste0(Img2, ".pdf"), plot2, dpi = 150L, width = w2, height = h2, units = "in")
})
#system(paste0("open \"", Img2, ".jpeg\""))
#system(paste0("open \"", Img2, ".pdf\""))
ReportCalls <- AddPlot2Report(Title = ttl2, Dir = pvalDir)
#
plot2ly <- ggplotly(plot2)
plot2ly <- layout(plot2ly,
                  xaxis = list(range = list(0L, 1L)),
                  showlegend = FALSE)
plot2 <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(plot2))
plotsList2 <- list(list(Title = ttl2,
                        Plot = plot2,
                        Plotly = slim_plotly(plot2ly)))
Imgs1 <- list.files(pvalDir, "^P-values scatter plot - .*\\.jpeg$", full.names = TRUE)
Img2 <- gsub("(\\.jpeg)+$", ".jpeg", paste0(Img2, ".jpeg"))
Imgs1 <- Imgs1[which(Imgs1 != Img2)]
IMGS <- c(Img2, Imgs1)
plotsList <- append(plotsList2, plotsList1)
Imgs1Nms <- gsub(" t-test|'s", "",
                 gsub("Moderated", "Mod.",
                      gsub("Permutations", "Perm.",
                           gsub(".*/P-values scatter plot - |\\.jpeg", "", Imgs1))))
names(plotsList) <- c("Histogram", Imgs1Nms)
#
# Which type of P-values do we want to use?
msg <- "Confirm which type of P-values to use for t-test volcano plots"
if ((dataType == "modPeptides") && (Ptm %in% names(PTMs_PVal_use))) {
  my_PVal_Use <- PTMs_PVal_use[[Ptm]]
} else {
  my_PVal_Use <- setNames(names(my_PVal_Col) == Param$P.values.type,
                          names(my_PVal_Col))
  pval_Use <- my_PVal_Use[names(pval_Col)]
}
if (!sum(pval_Use)) { pval_Use["Moderated"] <- TRUE }
if (!sum(pval_Use)) { pval_Use[1L] <- TRUE }
#
source(parSrc, local = FALSE)
IMGsDims <- as.data.frame(t(parSapply(parClust, IMGS, \(x) { #x <- IMGs[1L]
  a <- jpeg::readJPEG(x)
  setNames(dim(a)[1L:2L], c("height", "width"))
})))
IMGsDims$height <- round(screenRes$width*IMGsDims$height/max(IMGsDims$height)*0.3)
IMGsDims$width <- round(screenRes$width*IMGsDims$width/max(IMGsDims$width)*0.3)
ui <- fluidPage(
  useShinyjs(),
  setBackgroundColor( # Doesn't work
    color = c(#"#F8F8FF",
      "#F1FDFF"),
    gradient = "linear",
    direction = "bottom"
  ),
  extendShinyjs(text = jsToggleFS, functions = c("toggleFullScreen")),
  titlePanel(tag("u", appTtl),
             appNm),
  br(),
  #
  fluidRow(column(2L,
                  selectInput("PVal", msg, names(pval_Col), names(pval_Col)[which(pval_Use)])),
           column(2L,
                  actionBttn("saveBtn", "Save", icon = icon("save"), color = "success", style = "pill"))),
  br(),
  fluidRow(column(5L,
                  selectInput("XY",
                              "Select comparison to display...",
                              Imgs1Nms,
                              grep(names(pval_Col)[which(pval_Use)], Imgs1Nms, value = TRUE)[1L]),
                  withSpinner(plotlyOutput("Img1", inline = TRUE))),
           column(7L,
                  withSpinner(plotOutput("Img2", inline = TRUE)))),
  br(),
  br()
)
#h0 <- paste0(round(screenRes$height*0.35), "px")
if (exists("appRunTst")) { rm(appRunTst) }
server <- \(input, output, session) {
  myIMG1 <- reactiveVal(1L)
  updtIMG1 <- \(reactive = TRUE) {
    i <- if (reactive) { myIMG1() } else { 1L }
    renderPlotly({ plotsList[[i+1L]]$Plotly })
  }
  output$Img1 <- updtIMG1(FALSE)
  output$Img2 <- renderPlot({ print_gg(plotsList$Histogram$Plot) },
                            height = 50L+200L*length(pval_Col),
                            width = 50L+200L*nrow(myContrasts))
  observeEvent(input$PVal, {
    assign("my_PVal_Use", setNames(names(my_PVal_Col) == input$PVal, names(my_PVal_Col)), envir = .GlobalEnv)
  }, ignoreInit = FALSE)
  observeEvent(input$XY, {
    i <- match(input$XY, Imgs1Nms)
    myIMG1(i)
    output$Img1 <- updtIMG1()
  }, ignoreInit = FALSE)
  observeEvent(input$saveBtn, {
    stopApp()
    assign("appRunTst", TRUE, envir = .GlobalEnv)
  })
  #observeEvent(input$cancel, { stopApp() })
  session$onSessionEnded(function() { stopApp() })
}
runKount <- 0L
while ((!runKount) || (!exists("appRunTst"))) {
  eval(parse(text = runApp), envir = .GlobalEnv)
  shinyCleanup()
  runKount <- runKount + 1L
}
#
#
if (dataType == "PG") {
  # Update fold changes
  # The more advanced statistical frameworks provide estimates of average logFC per contrast. Get those.
  # if (!exists("prExpr_roots")) {
  #   prExpr_roots <- setNames(Prot.Expr.Root, "Quantitation")
  # } else {
  #   Prot.Expr.Root <- prExpr_roots["Quantitation"]
  # }
  # prExpr_roots %<o% prExpr_roots
  # Get the log2FCs from all methods for comparison
  ratKol <- paste0(ratRef, myContrasts$Contrast)
  logFCs <- list(make_Rat2 = set_colnames(quantData_list$Data[, ratKol], myContrasts$Contrast))
  nmsConv <- data.frame(Name = c("limma", "DEqMS", "QFeatures", "ROTS", "MSstats"),
                        P.values.type = c("Moderated", "DEqMS", "MSqRob", "ROTS", "MSstats"))
  for (i in 1L:nrow(nmsConv)) { #i <- 5L
    nm <- nmsConv$Name[i]
    whContr <- 1L:nrow(myContrasts)
    if (nm %in% c("limma", "DEqMS")) { #nm <- "limma"
      repRat <- limmaFits[[dataType]][[nm]]$fit$coefficients[, myContrasts$Contrast[whContr], drop = FALSE] # Already log2!!!
      repRat <- as.data.frame(repRat)
    }
    if (nm == "QFeatures") {
      repRat <- MSqRob_infer[, grep(topattern("MSqRob logFC - "), colnames(MSqRob_infer), value = TRUE), drop = FALSE]
      # We produced log10 values for QFeatures, and we tested them as log10, so we need to change base here since we want log2FC!
      repRat <- repRat/log10(2L) # Convert to log2!!!
      colnames(repRat) <- sub("^MSqRob logFC - ", "", colnames(repRat))
    }
    if (nm == "MSstats") {
      repRat <- msstatsLFC # Already log2!!!
    }
    if (nm == "ROTS") {
      whContr <- which(!myContrasts$isDouble) # for ROTS we have only single contrasts
      tmp <- lapply(myContrasts$Contrast[whContr], \(contr) {
        rwnms <- rownames(ROTS_res[[contr]]$data)
        logfc <- as.data.frame(t(data.frame(ROTS_res[[contr]]$logfc)))
        colnames(logfc) <- rwnms
        return(logfc)
      })
      tmp <- plyr::rbind.fill(tmp) # Already log2!!!
      repRat <- as.data.frame(t(tmp))
      rownames(repRat) <- colnames(tmp)
      colnames(repRat) <- myContrasts$Contrast[whContr]
    }
    repRat <- as.data.frame(repRat)
    w <- which(!is.finite(as.matrix(repRat)), arr.ind = TRUE)
    repRat[w] <- NA_real_
    m <- match(myData[[rowsCol]], rownames(repRat))
    repRat <- repRat[m, myContrasts$Contrast[whContr]]
    rownames(repRat) <- myData[[rowsCol]]
    logFCs[[nm]] <- repRat
  }
  #vapply(logFCs, nrow, 1L)
  #vapply(logFCs, ncol, 1L)
  # To check the distribution of logFCs (they should be different but still roughly similar)
  tst <- lapply(names(logFCs), \(nm) {
    x <- melt(logFCs[[nm]])
    x$Type <- nm
    return(x)
  })
  tst <- do.call(rbind, tst)
  tst <- tst[which(is.finite(tst$value)),]
  colnames(tst)[1L:2L] <- c("Contrast", "logFC")
  tst$Contrast <- factor(tst$Contrast, levels = myContrasts$Contrast)
  tst$Type <- factor(tst$Type, levels = names(logFCs))
  ttl <- "Distribution of logFCs"
  plot <- ggplot(tst) +
    ggtitle(ttl) +
    geom_density(aes(x = logFC, color = Type)) +
    facet_grid(Type~Contrast) +
    theme_bw()
  #poplot(plot, 12L, 22L)
  ggsave(paste0(pvalDir, "/", ttl, ".jpeg"), plot, dpi = 100L)
  if (Param$P.values.type %in% nmsConv$P.values.type) {
    nm <- nmsConv$Name[match(Param$P.values.type, nmsConv$P.values.type)]
    myData <- myData[, which(!colnames(myData) %in% ratKol)]
    repRat <- logFCs[[nm]]
    myData[, paste0(ratRef, colnames(repRat))] <- repRat
  }
}


if (dataType == "peptides") {
  pep <- myData
  pvalue.col %<o% my_PVal_Col
  pvalue.use %<o% my_PVal_Use
  Param$P.values.type <- names(pvalue.col)[which(pvalue.use)]
}
if (dataType == "modPeptides") {
  ptmpep <- myData
  PTMs_PVal_col[[Ptm]] <- my_PVal_Col
  PTMs_PVal_use[[Ptm]] <- my_PVal_Use
  ptms.PVal <- my_PVal_Col[which(my_PVal_Use)]
}
if (dataType == "PG") {
  PG <- myData
  pvalue.col %<o% my_PVal_Col
  pvalue.use %<o% my_PVal_Use
  Param$P.values.type <- names(pvalue.col)[which(pvalue.use)]
}
