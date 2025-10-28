### Check P-values
#
pvalDir <- paste0(wd, "/Workflow control/Protein groups/P-values")
if (!dir.exists(pvalDir)) { dir.create(pvalDir, recursive = TRUE) }
dirlist <- unique(c(dirlist, pvalDir))
#
# Scatter plots:
pkol2 <- gsub(" - $", "", pvalue.col)
temp <- lapply(VPAL$values, function(vpal) { #vpal <- VPAL$values[1]
  j <- paste0(pvalue.col, vpal)
  w <- which(j %in% colnames(PG))
  j <- j[w]
  if (length(j)) {
    x <- set_colnames(PG[, j, drop = FALSE], pkol2[w])
    x$Group <- cleanNms(vpal)
  } else { x <- NA }
  return(x)
})
temp <- temp[which(vapply(temp, function(x) { "data.frame" %in% class(x) }, TRUE))]
temp <- plyr::rbind.fill(temp)
kol <- colnames(temp)
kol <- kol[which(kol != "Group")]
temp <- temp[which(rowSums(temp[, kol]) > 0), ]
Comb <- gtools::combinations(length(kol), 2, kol)
clusterExport(parClust, list("Comb", "temp", "pvalDir", "pvalue.col", "plotEval"), envir = environment())
plotsList1 <- parLapply(parClust, 1:nrow(Comb), function(i) { #i <- 1
  X <- Comb[i, 1]
  Y <- Comb[i, 2]
  X2 <- gsub(" -log10\\(Pvalue\\)$", "", X)
  Y2 <- gsub(" -log10\\(Pvalue\\)$", "", Y)
  dat <- temp[, c(X, Y, "Group")]
  colnames(dat) <- c("X", "Y", "Group")
  dat$"P-value, X axis" <- proteoCraft::gsub_Rep(" -log10\\(Pvalue\\)$", "", Comb[i, 1])
  dat$"P-value, Y axis" <- proteoCraft::gsub_Rep(" -log10\\(Pvalue\\)$", "", Comb[i, 2])
  ttl1 <- paste0("P-values scatter plot - ", X2, " VS ", Y2)
  ttl1a <- paste0(ttl1, " (-log10)")
  Mx <- max(proteoCraft::is.all.good(c(dat$X, dat$Y)))
  uX <- unique(dat$`P-value, X axis`)
  uY <- unique(dat$`P-value, Y axis`)
  plot1 <- ggplot2::ggplot(dat) +
    scattermore::geom_scattermore(ggplot2::aes(x = X, y = Y, colour = Group),
                                  pixels = c(1024, 1024), pointsize = 3.2) +
    viridis::scale_color_viridis(begin = 0.4, end = 0.8, discrete = TRUE, option = "F") +
    #ggplot2::facet_grid(`P-value, Y axis`~Group+`P-value, X axis`) + # This was for when we saved all as one plot outside the parLapply
    #ggplot2::facet_grid(Group~`P-value, Y axis`) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggplot2::coord_fixed() + ggplot2::ggtitle(ttl1a) + ggplot2::theme_bw() +
    ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 0)) +
    ggplot2::xlim(0, Mx) + ggplot2::ylim(0, Mx) +
    ggplot2::xlab(uX) + ggplot2::ylab(uY) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0))
  plot1ly <- plotly::plot_ly(data = dat, x = ~X, y = ~Y, color = ~Group, type = "scatter", mode = "markers")
  plot1ly <- plotly::layout(plot1ly,
                            xaxis = list(range = list(0, 5), title = X),
                            yaxis = list(range = list(0, 5), title = Y),
                            showlegend = FALSE)
  #proteoCraft::poplot(plot1, 12, 20)
  Img1 <- paste0(pvalDir, "/", gsub(":", " - ", ttl1))
  # This was for when we saved all as one plot outside the parLapply:
  #h1 <- (length(pvalue.col)-0.8)*2
  #w1 <- ((length(pvalue.col)-1)*length(unique(dat$Group))+1)*2
  w1 <- 5
  h1 <- w1*(length(unique(dat$Group)) + 0.5)/3
  plot1 <- plotEval(plot1)
  suppressMessages({
    ggplot2::ggsave(paste0(Img1, ".jpeg"), plot1, dpi = 150, width = w1,# height = h1,
                    units = "in")
    ggplot2::ggsave(paste0(Img1, ".pdf"), plot1, dpi = 150, width = w1,# height = h1,
                    units = "in")
  })
  #system(paste0("open \"", Img1, ".jpeg\""))
  lst <- list(Plot = plot1,
              Plotly = plot1ly,
              Title = ttl1)
  return(lst)
})
invisible(clusterCall(parClust, function(x) { rm(list = ls());gc() }))
for (nm in 1:length(plotsList1)) {
  ReportCalls <- AddPlot2Report(Plot = plotsList1[[nm]]$Plot,
                                Title = plotsList1[[nm]]$Title,
                                Dir = pvalDir)
}
#vapply(plotsList, function(x) { x$Title }, "")
#
# P-values histogram:
nbin <- 20
bd <- (0:nbin)/nbin
w <- which(vapply(pvalue.col, function(type) { #type <- pvalue.col[1]
  length(grep(topattern(type), colnames(PG))) > 0
}, TRUE))
temp <- lapply(pvalue.col[w], function(type) { #type <- pvalue.col[w][1]
  kol <- grep(topattern(type), colnames(PG), value = TRUE)
  if (!length(kol)) { stop(type) }
  temp <- PG[, kol, drop = FALSE]
  colnames(temp) <- cleanNms(gsub(topattern(type), "", colnames(temp)))
  temp <- reshape::melt(temp, measure.vars = colnames(temp))
  temp$value <- 10^(-temp$value)
  temp <- temp[which(is.all.good(temp$value, 2)),]
  temp$Bin <- vapply(temp$value, function(x) { min(which(bd >= x))-1 }, 1)
  res <- aggregate(temp$Bin, list(temp$variable, temp$Bin), length)
  colnames(res) <- c("Group", "Bin", "Count")
  grps <- unique(res$Group)
  res$Frequency <- NA
  for (grp in grps) {
    w <- which(res$Group == grp)
    res$Frequency[w] <- res$Count[w]/sum(res$Count[w])
  }
  res$"P-value type" <- type
  res[, RSA$names] <- Exp.map[match(res$Variable, cleanNms(Exp.map[[VPAL$column]])),
                              RSA$names]
  res$Low <- 0
  return(res)
})
temp <- plyr::rbind.fill(temp)
temp$"P-value type" <- gsub_Rep(" -log10\\(Pvalue\\) - $", "", temp$"P-value type")
ttl2 <- "P-values histogram"
plot2 <- ggplot(temp) +
  geom_rect(position = "identity",
            aes(xmin = (Bin-1)/nbin, ymin = Low, xmax = Bin/nbin, ymax = Frequency, fill = `P-value type`),
            colour = "black", alpha = 0.5) + ylim(0, 1) +
  scale_fill_viridis(begin = 0.4, end = 0.8, discrete = TRUE, option = "G") +
  facet_grid(`P-value type`~Group) + ggtitle(ttl2) + theme_bw() +
  theme(strip.text.y = element_text(angle = 0))
#poplot(plot2, 12, 20)
Img2 <- paste0(pvalDir, "/", ttl2)
w2 <- ((length(VPAL$values)+1)*1.25)*2
h2 <- ((length(pvalue.col)+0.2)*1.25)*2
suppressMessages({
  ggsave(paste0(Img2, ".jpeg"), plot2, dpi = 300, width = w2, height = h2, units = "in")
  ggsave(paste0(Img2, ".pdf"), plot2, dpi = 300, width = w2, height = h2, units = "in")
})
#system(paste0("open \"", Img2, ".jpeg\""))
ReportCalls <- AddPlot2Report(Title = ttl2, Dir = pvalDir)
#
# Which type of P-values do we want to use?
Imgs1 <- list.files(pvalDir, "^P-values scatter plot - .*\\.jpeg$", full.names = TRUE)
Img2 <- gsub("(\\.jpeg)+$", ".jpeg", paste0(Img2, ".jpeg"))
Imgs1 <- Imgs1[which(Imgs1 != Img2)]
IMGS <- c(Img2, Imgs1)
plot2ly <- plotly::ggplotly(plot2)
plot2ly <- plotly::layout(plot2ly,
                          xaxis = list(range = list(0, 1)),
                          showlegend = FALSE)
plot2 <- plot2 +  coord_fixed(ratio = 0.75) + scale_x_continuous(limits = c(0, 1), expand = c(0, 0))
plot2 <- plotEval(plot2)
plotsList2 <- list(list(Title = ttl2,
                        Plot = plot2,
                        Plotly = plot2ly))
plotsList <- append(plotsList2, plotsList1)
msg <- "Confirm which type of P-values to use for t-test volcano plots"
pvalue.use %<o% (names(pvalue.col) == Param$P.values.type)
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
  fluidRow(column(2, selectInput("PVal", msg, names(pvalue.col), Param$P.values.type)),
           column(2, actionBttn("saveBtn", "Save", icon = icon("save"), color = "success", style = "pill"))),
  br(),
  br(),
  fluidRow(column(6,
                  selectInput("XY", "Comparison to display...", Imgs1Nms, Imgs1Nms[1]),
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
  column(6,
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
IMGsDims <- as.data.frame(t(parSapply(parClust, IMGS, function(x) { #x <- IMGs[1]
  a <- jpeg::readJPEG(x)
  setNames(dim(a)[1:2], c("height", "width"))
})))
IMGsDims$height <- round(screenRes$width*IMGsDims$height/max(IMGsDims$height)*0.3)
IMGsDims$width <- round(screenRes$width*IMGsDims$width/max(IMGsDims$width)*0.3)
fct <- 1
server <- function(input, output, session) {
  myIMG1 <- reactiveVal(1)
  updtIMG1 <- function(reactive = TRUE) {
    if (reactive) { i <- myIMG1() } else { i <- 1 }
    if (Mode == "plotly") {
      rs <- renderPlotly({ plotsList[[i+1]]$Plotly })
    }
    if (Mode == "ggplot") {
      rs <- renderPlot({ plotsList[[i+1]]$Plot })
    }
    if (Mode == "raster") {
      rs <- renderImage({
        list(src = Imgs1[i],
             #height = IMGsDims$height[i+1]*fct,
             width = IMGsDims$width[i+1]*fct)
      }, deleteFile = FALSE)
    }
    return(rs)
  }
  output$Img1 <- updtIMG1(FALSE)
  if (Mode == "plotly") {
    output$Img2 <- renderPlotly({ plotsList[[1]]$Plotly })
  }
  if (Mode == "ggplot") {
    output$Img2 <- renderPlot({ plotsList[[1]]$Plot })
  }
  if (Mode == "raster") {
    output$Img2 <- renderImage({
      list(src = IMGS[1], height = IMGsDims$height[1], width = IMGsDims$width[1])
    }, deleteFile = FALSE)
  }
  observeEvent(input$PVal, {
    assign("pvalue.use", names(pvalue.col) == input$PVal, envir = .GlobalEnv)
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
