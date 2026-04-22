### Check P-values
#
source(parSrc)
pvalDir <- paste0(wd, "/Workflow control/Protein groups/P-values")
if (!dir.exists(pvalDir)) { dir.create(pvalDir, recursive = TRUE) }
dirlist <- unique(c(dirlist, pvalDir))
#
w <- which(vapply(pvalue.col, \(type) { #type <- pvalue.col[1L]
  length(grep(topattern(type), colnames(PG))) > 0L
}, TRUE))
pvalueCol <- pvalue.col[w]
#
# Scatter plots:
pkol2 <- gsub(" - $", "", pvalueCol)
whSingle <- which(!myContrasts$isDouble)
temp <- lapply(whSingle, \(i) { #i <- 1L
  nm <- myContrasts$Contrast[i]
  j <- paste0(pvalueCol, nm)
  w <- which(j %in% colnames(PG))
  j <- j[w]
  if (!length(j)) { return(NA) }
  x <- PG[, j, drop = FALSE]
  colnames(x) <- pkol2[w]
  x$PG <- PG$Label
  x$Contrast <- nm
  return(x)
})
temp <- temp[which(vapply(temp, is.data.frame, TRUE))]
temp <- plyr::rbind.fill(temp)
temp$Contrast <- gsub_Rep(" - ", "\n- ", temp$Contrast)
tmpContr <- gsub(" - ", "\n- ", myContrasts$Contrast)
temp$Contrast <- factor(temp$Contrast, levels = tmpContr)
kol <- colnames(temp)
kol <- setdiff(kol, c("PG", "Contrast"))
temp <- temp[which(rowSums(temp[, kol]) > 0), ]
Comb <- gtools::combinations(length(kol), 2L, kol)
tmpFl <- tempfile(fileext = ".rds")
readr::write_rds(temp, tmpFl)
clusterExport(parClust, list("Comb", "tmpFl", "pvalDir", "pvalueCol", "gsub_Rep", "is.all.good", "plotEval"), envir = environment())
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
  X2 <- gsub(" -log10\\(Pvalue\\)$", "", X)
  Y2 <- gsub(" -log10\\(Pvalue\\)$", "", Y)
  dat <- temp[, c(X, Y, "PG", "Contrast")]
  colnames(dat)[1L:2L] <- c("X", "Y")
  dat$"P-value, X axis" <- gsub_Rep(" -log10\\(Pvalue\\)$", "", Comb[i, 1L])
  dat$"P-value, Y axis" <- gsub_Rep(" -log10\\(Pvalue\\)$", "", Comb[i, 2L])
  ttl1 <- paste0("P-values scatter plot - ", X2, " VS ", Y2)
  ttl1a <- paste0(ttl1, " (-log10)")
  Mx <- max(is.all.good(c(dat$X, dat$Y)))
  uX <- unique(dat$`P-value, X axis`)
  uY <- unique(dat$`P-value, Y axis`)
  plot1 <- ggplot(dat) +
    geom_scattermore(aes(x = X, y = Y, colour = Contrast),
                     pixels = c(1024L, 1024L), pointsize = 3.2) +
    scale_color_viridis(begin = 0.4, end = 0.8, discrete = TRUE, option = "F") +
    #facet_grid(`P-value, Y axis`~Contrast+`P-value, X axis`) + # This was for when we saved all as one plot outside the parLapply
    #facet_grid(Contrast~`P-value, Y axis`) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    coord_fixed() + ggtitle(ttl1a) + theme_bw() +
    theme(strip.text.y = element_text(angle = 0)) +
    xlim(0, Mx) + ylim(0L, Mx) +
    xlab(uX) + ylab(uY) +
    scale_x_continuous(expand = c(0L, 0L)) +
    scale_y_continuous(expand = c(0L, 0L))
  plot1ly <- plot_ly(data = dat, x = ~X, y = ~Y, color = ~Contrast, type = "scatter", mode = "markers",
                     text = ~paste("Protein:", PG,
                                   "<br>x = ", signif(X, 3L),
                                   "<br>y = ", signif(Y, 3L)),
                     hoverinfo = "text")
  plot1ly <- layout(plot1ly,
                    xaxis = list(range = list(0L, 5L), title = X),
                    yaxis = list(range = list(0L, 5L), title = Y),
                    showlegend = FALSE)
  #poplot(plot1, 12, 20)
  Img1 <- paste0(pvalDir, "/", gsub(":", " - ", ttl1))
  # This was for when we saved all as one plot outside the parLapply:
  #h1 <- (length(pvalueCol)-0.8)*2
  #w1 <- ((length(pvalueCol)-1L)*length(unique(dat$Contrast))+1L)*2L
  w1 <- 5
  h1 <- w1*(length(unique(dat$Contrast)) + 0.5)/3
  plot1 <- plotEval(plot1)
  suppressMessages({
    ggsave(paste0(Img1, ".jpeg"), plot1, dpi = 150L, width = w1,# height = h1,
           units = "in")
    ggsave(paste0(Img1, ".pdf"), plot1, dpi = 150L, width = w1,# height = h1,
           units = "in")
  })
  #system(paste0("open \"", Img1, ".jpeg\""))
  lst <- list(Plot = plot1,
              Plotly = plot1ly,
              Title = ttl1)
  return(lst)
})
invisible(clusterCall(parClust, \(x) { rm(list = ls());gc() }))
for (nm in 1L:length(plotsList1)) {
  ReportCalls <- AddPlot2Report(Plot = plotsList1[[nm]]$Plot,
                                Title = plotsList1[[nm]]$Title,
                                Dir = pvalDir)
}
#vapply(plotsList, \(x) { x$Title }, "")
#
# P-values histogram:
nbin <- 20L
bd <- (0L:nbin)/nbin

temp <- lapply(pvalueCol, \(type) { #type <- pvalueCol[1L]
  kol <- grep(topattern(type), colnames(PG), value = TRUE)
  if (!length(kol)) { stop(type) }
  temp <- PG[, kol, drop = FALSE]
  colnames(temp) <- cleanNms(gsub(topattern(type), "", colnames(temp)))
  temp <- dfMelt(temp)
  temp$value <- 10L^(-temp$value)
  temp <- temp[which(is.all.good(temp$value, 2L)),]
  temp$Bin <- vapply(temp$value, \(x) { min(which(bd >= x))-1L }, 1L)
  res <- aggregate(temp$Bin, list(temp$variable, temp$Bin), length)
  colnames(res) <- c("Contrast", "Bin", "Count")
  grps <- unique(res$Contrast)
  res$Frequency <- NA
  for (grp in grps) {
    w <- which(res$Contrast == grp)
    res$Frequency[w] <- res$Count[w]/sum(res$Count[w])
  }
  res$"P-value type" <- type
  res$Low <- 0L
  return(res)
})
temp <- plyr::rbind.fill(temp)
temp$"P-value type" <- gsub_Rep(" -log10\\(Pvalue\\) - $", "", temp$"P-value type")
dfltLvls <- c("Moderated t-test", "DEqMS mod. t-test", "Student's t-test", "Welch's t-test")
u <- unique(temp$"P-value type")
temp$"P-value type" <- factor(temp$"P-value type",
                              levels = c(dfltLvls, setdiff(u, dfltLvls)))
temp$Contrast <- gsub_Rep(" - ", "\n- ", temp$Contrast)
temp$Contrast <- factor(temp$Contrast, levels = tmpContr)
ttl2 <- "P-values histogram"
plot2 <- ggplot(temp) +
  geom_rect(position = "identity",
            aes(xmin = (Bin-1L)/nbin, ymin = Low, xmax = Bin/nbin, ymax = Frequency, fill = `P-value type`),
            colour = "black", alpha = 0.5) + ylim(0L, 1L) +
  scale_fill_viridis(begin = 0.4, end = 0.8, discrete = TRUE, option = "G") +
  facet_grid(`P-value type`~Contrast) + ggtitle(ttl2) + theme_bw() +
  theme(strip.text.y = element_text(angle = 0),
        strip.text.x = element_text(size = 7L),
        axis.text = element_text(size = 5L))
#poplot(plot2, 12, 20)
Img2 <- paste0(pvalDir, "/", ttl2)
w2 <- ((length(whSingle)+1)*1.25)*2
h2 <- ((length(pvalueCol)+0.2)*1.25)*2
suppressMessages({
  ggsave(paste0(Img2, ".jpeg"), plot2, dpi = 300L, width = w2, height = h2, units = "in")
  ggsave(paste0(Img2, ".pdf"), plot2, dpi = 300L, width = w2, height = h2, units = "in")
})
#system(paste0("open \"", Img2, ".jpeg\""))
ReportCalls <- AddPlot2Report(Title = ttl2, Dir = pvalDir)
#
# Which type of P-values do we want to use?
Imgs1 <- list.files(pvalDir, "^P-values scatter plot - .*\\.jpeg$", full.names = TRUE)
Img2 <- gsub("(\\.jpeg)+$", ".jpeg", paste0(Img2, ".jpeg"))
Imgs1 <- Imgs1[which(Imgs1 != Img2)]
IMGS <- c(Img2, Imgs1)
plot2ly <- ggplotly(plot2)
plot2ly <- layout(plot2ly,
                  xaxis = list(range = list(0L, 1L)),
                  showlegend = FALSE)
plot2 <- plot2 + coord_fixed(ratio = 0.75) + scale_x_continuous(limits = c(0L, 1L), expand = c(0L, 0L))
plot2 <- plotEval(plot2)
plotsList2 <- list(list(Title = ttl2,
                        Plot = plot2,
                        Plotly = plot2ly))
plotsList <- append(plotsList2, plotsList1)
msg <- "Confirm which type of P-values to use for t-test volcano plots"
pvalue.use %<o% setNames((names(pvalue.col) == Param$P.values.type), names(pvalue.col))
pvalueUse <- pvalue.use[names(pvalueCol)]
if (!sum(pvalueUse)) { pvalueUse["Moderated"] <- TRUE }
if (!sum(pvalueUse)) { pvalueUse[1L] <- TRUE }

appNm <- paste0(dtstNm, " - t-test P-values")
Imgs1Nms <- gsub(" t-test|'s", "", gsub("Moderated", "Mod.", gsub("Permutations", "Perm.", gsub(".*/P-values scatter plot - |\\.jpeg", "", Imgs1))))
#Mode <- "raster"
Mode <- "ggplot" # not working
Mode <- "plotly"
ui <- fluidPage(
  useShinyjs(),
  setBackgroundColor( # Doesn't work
    color = c(#"#F8F8FF",
      "#F1FDFF"),
    gradient = "linear",
    direction = "bottom"
  ),
  extendShinyjs(text = jsToggleFS, functions = c("toggleFullScreen")),
  titlePanel(tag("u", "Type of t-test P-values"),
             appNm),
  br(),
  #
  fluidRow(column(2L,
                  selectInput("PVal", msg, names(pvalueCol), names(pvalueCol)[which(pvalueUse)])),
           column(2L,
                  actionBttn("saveBtn", "Save", icon = icon("save"), color = "success", style = "pill"))),
  br(),
  br(),
  fluidRow(column(6L,
                  selectInput("XY", "Contrast to display...", Imgs1Nms, Imgs1Nms[1L]),
                  if (Mode == "raster") {
                    withSpinner(imageOutput("Img1", inline = TRUE))
                  },
                  if (Mode == "ggplot") {
                    withSpinner(plotOutput("Img1", inline = TRUE, height = "600px"))
                  },
                  if (Mode == "plotly") {
                    withSpinner(plotlyOutput("Img1", inline = TRUE, height = "600px"))
                  },
  ),
  column(6L,
         if (Mode == "raster") {
           withSpinner(imageOutput("Img2", inline = TRUE))
         },
         if (Mode == "ggplot") {
           withSpinner(plotOutput("Img2", inline = TRUE, height = "600px"))
         },
         if (Mode == "plotly") {
           withSpinner(plotlyOutput("Img2", inline = TRUE, height = "600px"))
         },
  )),
  br(),
  br()
)
#h0 <- paste0(round(screenRes$height*0.35), "px")
source(parSrc, local = FALSE)
IMGsDims <- as.data.frame(t(parSapply(parClust, IMGS, \(x) { #x <- IMGs[1L]
  a <- jpeg::readJPEG(x)
  setNames(dim(a)[1L:2L], c("height", "width"))
})))
IMGsDims$height <- round(screenRes$width*IMGsDims$height/max(IMGsDims$height)*0.3)
IMGsDims$width <- round(screenRes$width*IMGsDims$width/max(IMGsDims$width)*0.3)
fct <- 1L
server <- \(input, output, session) {
  myIMG1 <- reactiveVal(1L)
  updtIMG1 <- \(reactive = TRUE) {
    i <- if (reactive) { myIMG1() } else { 1L }
    if (Mode == "plotly") {
      rs <- renderPlotly({ plotsList[[i+1L]]$Plotly })
    }
    if (Mode == "ggplot") {
      rs <- renderPlot({ plotsList[[i+1L]]$Plot })
    }
    if (Mode == "raster") {
      rs <- renderImage({
        list(src = Imgs1[i],
             #height = IMGsDims$height[i+1L]*fct,
             width = IMGsDims$width[i+1L]*fct)
      }, deleteFile = FALSE)
    }
    return(rs)
  }
  output$Img1 <- updtIMG1(FALSE)
  if (Mode == "plotly") {
    output$Img2 <- renderPlotly({ plotsList[[1L]]$Plotly })
  }
  if (Mode == "ggplot") {
    output$Img2 <- renderPlot({ plotsList[[1L]]$Plot })
  }
  if (Mode == "raster") {
    output$Img2 <- renderImage({
      list(src = IMGS[1L], height = IMGsDims$height[1L], width = IMGsDims$width[1L])
    }, deleteFile = FALSE)
  }
  observeEvent(input$PVal, {
    assign("pvalue.use", setNames(names(pvalue.col) == input$PVal, names(pvalue.col)), envir = .GlobalEnv)
  }, ignoreInit = FALSE)
  observeEvent(input$XY, {
    i <- match(input$XY, Imgs1Nms)
    myIMG1(i)
    output$Img1 <- updtIMG1()
  }, ignoreInit = FALSE)
  observeEvent(input$saveBtn, { stopApp() })
  #observeEvent(input$cancel, { stopApp() })
  session$onSessionEnded(function() { stopApp() })
}
eval(parse(text = runApp), envir = .GlobalEnv)
shinyCleanup()
Param$P.values.type <- names(pvalue.col)[which(pvalue.use)]
