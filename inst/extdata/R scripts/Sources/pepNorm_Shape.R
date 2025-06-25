
# These normalisations are more like corrections: they will change the shape of the vector
#
# Initial values for outputs
tmpDat2 <- NA
wAG2 <- wAG1
Outcome <- TRUE
txt2 <- ""
#
normMeth <- normSequence[[nrmStp]]$Method #normMeth <- "LOESS" #normMeth <- "VSN"
pack <- "hexbin"
cran_req <- unique(c(cran_req, pack))
pack2 <- c("affy", "vsn")[match(normMeth, c("LOESS", "VSN"))]
bioc_req <- unique(c(bioc_req, pack2))
for (p in c(pack, pack2)) {
  if (!require(p, character.only = TRUE, quietly = TRUE)) { pak::pkg_install(p) }
}
for (p in c(pack, pack2)) {
  require(p, character.only = TRUE)
}
shpDr <- paste0(nrmDr, "/Step ", nrmStp, " - ", normMeth)
if (!dir.exists(shpDr)) { dir.create(shpDr, recursive = TRUE) }
dirlist <- unique(c(dirlist, shpDr))
#
currSamples <- allSamples[which(allSamples %in% colnames(tmpDat1))]
A <- parApply(parClust, tmpDat1[wAG1, currSamples], 1, function(x) { mean(proteoCraft::is.all.good(x)) })
#
# Impute
ImpGrps <- Exp.map[match(currSamples, Exp.map$Ref.Sample.Aggregate),
                   VPAL$column]
tmp <- proteoCraft::Data_Impute2(tmpDat1[, currSamples], ImpGrps)
tmpDat1a <- tmp$Imputed_data
wPos <- which(tmp$Positions_Imputed, arr.ind = TRUE)
tmpDat1a[, c("id", "Group")] <- tmpDat1[, c("id", "Group")] 

# Visualize - before
temp_plotA <- tmpDat1a[wAG1, c("id", "Group", currSamples)]
colnames(temp_plotA) <- proteoCraft::cleanNms(colnames(temp_plotA))
temp_plotA <- reshape::melt(temp_plotA, id.vars = c("id", "Group"))
temp_plotA$A <- rep(A, length(currSamples))
temp_plotA <- temp_plotA[which(proteoCraft::is.all.good(temp_plotA$value, 2)),]
temp_plotA$M <- temp_plotA$value - temp_plotA$A
annot <- data.frame(variable = unique(temp_plotA$variable))
annot$Median <- sapply(annot$variable, function(x) {
  paste0("Median: ", round(median(proteoCraft::is.all.good(temp_plotA$M[which(temp_plotA$variable == x)])), 3))
})
annot$IQR <- sapply(annot$variable, function(x) {
  paste0("IQR: ", round(IQR(proteoCraft::is.all.good(temp_plotA$M[which(temp_plotA$variable == x)])), 3))
})
annot$Amax <- max(proteoCraft::is.all.good(temp_plotA$A))*1.1
annot$Amin <- min(proteoCraft::is.all.good(temp_plotA$A))*1.1
annot$Mmax <- max(proteoCraft::is.all.good(temp_plotA$M))*1.1
annot$Mmin <- min(proteoCraft::is.all.good(temp_plotA$M))*1.1
annot1 <- annot[, c("variable", "Amax", "Mmin", "Mmax")] 
annot1 <- rbind(annot1, annot1)
annot1$Label <- c(annot$Median, annot$IQR)
MAttl <- "MA plot"
MAttlA <- paste0(MAttl, "_before")
ylim <- max(c(abs(c(annot$Mmax, annot$Mmin, (annot$Amax-annot$Amin)/4))))
annot1$Y <- ylim*0.9
w <- grep("^IQR: ", annot1$Label)
annot1$Y[w] <- -ylim*0.9
l1 <- length(unique(temp_plotA$variable))
l2 <- length(unique(temp_plotA$Group))
nkol <- max(c(1, round(sqrt(l1*l2))))
if ((l2 > 1)&&((nkol %% l2) != 0)) { nkol <- ceiling(nkol/l2)*l2 }
while (nkol > l1*l2) { nkol <- nkol-1 }
MAplotA <- ggplot(temp_plotA) +
  geom_scattermore(aes(x = A, y = M, colour = Group), size = 1, alpha = 1) +
  geom_hline(yintercept = 0, colour = "grey") + geom_smooth(aes(x = A, y = M), color = "red", linewidth = 0.8, linetype = "dashed") +
  geom_text(data = annot1, aes(x = Amax, y = Y, label = Label), hjust = 1, cex = 2) +
  scale_color_viridis_d(begin = 0.25) +
  facet_wrap(~variable+Group, ncol = nkol) + coord_fixed(log10(2)) + theme_bw() + ggtitle(MAttlA) +
  theme(legend.position = "bottom")
#poplot(MAplotA, 12, 22)
#
tmpDat2 <- tmpDat1[, currSamples]*NA
#
# Apply correction per peptides normalisation group
library(ggpubr)
for (lGrp in NormGrps$Group) { #lGrp <- NormGrps$Group[1]
  grpMtch <- match(NormGrps$IDs[[match(lGrp, NormGrps$Group)]],
                   tmpDat1$id[wAG1])
  grpMtch <- grpMtch[which(!is.na(grpMtch))]
  # Imputation to allow normalisation
  # We will fill gaps with normal distribution-based simulated data;
  # For each row, the mean of the distribution will be the mean of the row;
  # and the sd will be based on a LOESS regression estimate (nearest neighbours if missing):
  # Normalisation proper:
  if (normMeth == "LOESS") { tmpDat1b <- limma::normalizeCyclicLoess(as.matrix(tmpDat1a[grpMtch, currSamples])) }
  if (normMeth == "VSN") { tmpDat1b <- vsn::justvsn(as.matrix(10^tmpDat1a[grpMtch, currSamples]))/log2(10) }
  tmpDat2[grpMtch, currSamples] <- tmpDat1b[, currSamples]
  sd1 <- meanSdPlot(as.matrix(tmpDat1a[grpMtch, currSamples]), plot = FALSE)$gg
  sd2 <- meanSdPlot(as.matrix(tmpDat1b), plot = FALSE)$gg
  SDttl <- "mean SD plot"
  if (length(NormGrps$Group) > 1) { SDttl <- paste0(SDttl, " - ", lGrp) }
  SDttlB <- SDttlA <- SDttl
  SDttlA <- paste0(SDttlA, "_before")
  SDttlB <- paste0(SDttlB, "_after")
  SDplotA <- sd1 + theme_bw() + ggtitle(SDttl, subtitle = "Before")
  SDplotB <- sd2 + theme_bw() + ggtitle("", subtitle = "After")
  SDplot <- ggarrange(SDplotA, SDplotB, ncol = 2, nrow = 1)
  #poplot(SDplot, 6, 12)
  ggsave(paste0(shpDr, "/", SDttl, ".jpeg"), SDplot, dpi = 150, width = 12, height = 6, units = "in")
  ggsave(paste0(shpDr, "/", SDttl, ".pdf"), SDplot, dpi = 150, width = 12, height = 6, units = "in")
  ReportCalls <- AddPlot2Report(Plot = SDplot, Title = SDttl, Space = FALSE, Dir = shpDr)
}
#
wAG2 <- wAG1
# Remove imputed data:
if (nrow(wPos)) { tmpDat2[wPos] <- tmpDat1[, currSamples][wPos] }
# Visualize - after
temp_plotB <- tmpDat2[wAG2, currSamples]
colnames(temp_plotB) <- proteoCraft::cleanNms(colnames(temp_plotB))
temp_plotB[, c("id", "Group")] <- tmpDat1[, c("id", "Group")]
temp_plotB <- reshape::melt(temp_plotB, id.vars = c("id", "Group"))
temp_plotB$A <- rep(A, length(currSamples))
temp_plotB <- temp_plotB[which(proteoCraft::is.all.good(temp_plotB$value, 2)),]
temp_plotB$M <- temp_plotB$value - temp_plotB$A
annot <- data.frame(variable = unique(temp_plotB$variable))
annot$Median <- sapply(annot$variable, function(x) {
  paste0("Median: ", round(median(proteoCraft::is.all.good(temp_plotB$M[which(temp_plotB$variable == x)])), 3))
})
annot$IQR <- sapply(annot$variable, function(x) {
  paste0("IQR: ", round(IQR(proteoCraft::is.all.good(temp_plotB$M[which(temp_plotB$variable == x)])), 3))
})
annot$Amax <- max(proteoCraft::is.all.good(temp_plotB$A))*1.1
annot$Amin <- min(proteoCraft::is.all.good(temp_plotB$A))*1.1
annot$Mmax <- max(proteoCraft::is.all.good(temp_plotB$M))*1.1
annot$Mmin <- min(proteoCraft::is.all.good(temp_plotB$M))*1.1
annot2 <- annot[, c("variable", "Amax", "Mmin", "Mmax")] 
annot2 <- rbind(annot2, annot2)
annot2$Label <- c(annot$Median, annot$IQR)
MAttlB <- paste0(MAttl, "_after")
ylim <- max(c(abs(c(annot$Mmax, annot$Mmin, (annot$Amax-annot$Amin)/4))))
annot2$Y <- ylim*0.9
w <- grep("^IQR: ", annot2$Label)
annot2$Y[w] <- -ylim*0.9
l1 <- length(unique(temp_plotB$variable))
l2 <- length(unique(temp_plotB$Group))
nkol <- max(c(1, round(sqrt(l1*l2))))
if ((l2 > 1)&&((nkol %% l2) != 0)) { nkol <- ceiling(nkol/l2)*l2 }
while (nkol > l1*l2) { nkol <- nkol-1 }
MAplotB <- ggplot(temp_plotB) +
  geom_scattermore(aes(x = A, y = M, colour = Group), size = 1, alpha = 1) +
  geom_hline(yintercept = 0, colour = "grey") + geom_smooth(aes(x = A, y = M), color = "red", linewidth = 0.8, linetype = "dashed") +
  geom_text(data = annot2, aes(x = Amax, y = Y, label = Label), hjust = 1, cex = 2) +
  scale_color_viridis_d(begin = 0.25) +
  facet_wrap(~variable+Group, ncol = nkol) + coord_fixed(log10(2)) + theme_bw() + ggtitle(MAttlB) +
  theme(legend.position = "bottom")
tmp <- proteoCraft::is.all.good(c(temp_plotA$M, temp_plotB$M))
MAplotA <- MAplotA + ylim(min(tmp), max(tmp))
MAplotB <- MAplotB + ylim(min(tmp), max(tmp))
MAplot <- ggarrange(MAplotA, MAplotB, ncol = 2, nrow = 1)
#poplot(MAplot, 12, 22)
ggsave(paste0(shpDr, "/", MAttl, ".jpeg"), MAplot, dpi = 150, width = 24, units = "in")
ggsave(paste0(shpDr, "/", MAttl, ".pdf"), MAplot, dpi = 150, width = 24, units = "in")
ReportCalls <- AddPlot2Report(Plot = MAplot, Title = MAttl, Dir = shpDr)
#
appNm <- paste0(normMeth, " normalisation")
msg <- paste0("Keep results from ", normMeth, " normalisation? (untick to cancel correction)")
if ((!exists("KeepShapeCorrRes"))||(length(KeepShapeCorrRes) != 1)||(!is.logical(KeepShapeCorrRes))||(is.na(KeepShapeCorrRes))) {
  KeepShapeCorrRes <- TRUE
}
IMGs <- paste0(shpDr, "/", c(MAttl, SDttl), ".jpeg")
IMGsDims <- as.data.frame(t(parSapply(parClust, IMGs, function(x) { #x <- IMGs[1]
  a <- jpeg::readJPEG(x)
  setNames(dim(a)[1:2], c("height", "width"))
})))
IMGsDims$height <- screenRes$width*0.35*IMGsDims$height/max(IMGsDims$height)
IMGsDims$width <- screenRes$width*0.35*IMGsDims$width/max(IMGsDims$width)
if (exists("IHAVERUN")) { rm(IHAVERUN) }
ui <- fluidPage(
  useShinyjs(),
  setBackgroundColor( # Doesn't work
    color = c(#"#F8F8FF",
      "#E9F2F7"),
    gradient = "linear",
    direction = "bottom"
  ),
  extendShinyjs(text = jsToggleFS, functions = c("toggleFullScreen")),
  tags$head(tags$style(HTML("table {table-layout: fixed;"))), # So table widths can be properly adjusted!
  titlePanel(tag("u", appNm),
             appNm),
  br(),
  fluidRow(column(5,
                  checkboxInput("KeepResults", msg, KeepShapeCorrRes),
                  actionBttn("saveBtn", "Save", icon = icon("save"), color = "success", style = "pill"),
                  h4("Recommended criteria:"),
                  h5(HTML("&nbsp;Does the original MA plot look like it is skewed?")),
                  h5(HTML("&nbsp;&nbsp;-> If yes: accept the correction if the correction did remove the skew.")))),
  br(),
  fluidRow(withSpinner(imageOutput("MAplots", height = "505px"))),
  fluidRow(withSpinner(imageOutput("SDplots", height = "505px"))),
  br(),
  br()
)
server <- function(input, output, session) {
  output$MAplots <- renderPlot(MAplot, width = screenRes$width*0.8)
  output$SDplots <- renderPlot(SDplot, height = screenRes$width*0.8/2, width = screenRes$width*0.8)
  # output$MAplots <- renderImage({
  #   list(src = IMGs[1], height = 500, width = IMGsDims$width[1]*500/IMGsDims$height[1])
  # }, deleteFile = FALSE)
  # output$SDplots <- renderImage({
  #   list(src = IMGs[2], height = 500, width = IMGsDims$width[2]*500/IMGsDims$height[2])
  # }, deleteFile = FALSE)
  # output$Before <- renderPlotly(plot1)
  # output$After <- renderPlotly(plot2)
  observeEvent(input$saveBtn, {
    assign("KeepShapeCorrRes", as.logical(input[["KeepResults"]]), envir = .GlobalEnv)
    assign("IHAVERUN", TRUE, .GlobalEnv)
    stopApp()
  })
  #observeEvent(input$cancel, { stopApp() })
  session$onSessionEnded(function() { stopApp() })
}
while (!exists("IHAVERUN")) {
  eval(parse(text = runApp), envir = .GlobalEnv)
}
msg <- paste0(" -> ", normMeth, " correction for intensity range variance biases ", c("rejec", "accep")[KeepShapeCorrRes+1], "ted.\n")
if (KeepShapeCorrRes) {
  txt2 <- paste0("corrected for intensity range variance biases using ", c(paste0(normMeth, " regression"), "VSN")[match(normMeth, c("LOESS", "VSN"))])
}
Outcome <- KeepShapeCorrRes
cat(msg)
