
# These normalisations are more like corrections: they will change the shape of the vector
#
if (!"Clean_name" %in% colnames(Exp.map)) {
  Exp.map$Clean_name <- cleanNms(Exp.map$Ref.Sample.Aggregate)
}
#
# Initial values for outputs
tmpDat2 <- NA
wAG2 <- wAG1
Outcome <- TRUE
txt2 <- ""
#
normMeth <- normSequence[[nrmStp]]$Method #normMeth <- "LOESS" #normMeth <- "VSN"
packs <- c()
if (normMeth == "VSN") {
  cran_req <- union(cran_req, "hexbin")
  bioc_req <- union(bioc_req, "vsn")
  packs <- union(packs, c("hexbin", "vsn"))
}
if (normMeth == "LOESS") {
  bioc_req <- union(bioc_req, "affy")
  packs <- union(packs, "affy")
}
if (normMeth == "GAM") {
  cran_req <- union(cran_req, "mgcv")
  packs <- union(packs, "mgcv")
}
for (p in packs) { if (!require(p, character.only = TRUE, quietly = TRUE)) { pak::pak(p) } }
for (p in packs) { require(p, character.only = TRUE) }
#
shpDr <- paste0(nrmDr, "/Step ", nrmStp, " - ", normMeth)
if (!dir.exists(shpDr)) { dir.create(shpDr, recursive = TRUE) }
dirlist <- unique(c(dirlist, shpDr))
#
currSamples <- intersect(allSamples, colnames(tmpDat1))
A <- parApply(parClust, tmpDat1[wAG1, currSamples], 1L, \(x) { mean(x[which(is.finite(x))]) })
#
tmpDat2 <- tmpDat1[, currSamples]*NA_real_
#
# Apply correction per peptides normalisation group
#library(ggpubr)
tstNorm <- try({
  for (lGrp in NormGrps$Group) { #lGrp <- NormGrps$Group[1]
    grpMtch <- match(NormGrps$IDs[[match(lGrp, NormGrps$Group)]],
                     tmpDat1$id[wAG1])
    grpMtch <- grpMtch[which(!is.na(grpMtch))]
    #
    dat <- as.matrix(tmpDat1[grpMtch, currSamples])
    if (normMeth == "LOESS") {
      tmpDat1b <- limma::normalizeCyclicLoess(dat, method = "fast")
    }
    if (normMeth == "VSN") {
      tmpDat1b <- vsn::justvsn(10L^dat)/log2(10L)
    }
    if (normMeth == "GAM") {
      avg <- rowMeans(replace(dat, !is.finite(dat), NA), na.rm = TRUE)
      A <- sweep(dat[, currSamples], 1L, avg, "+")/2
      M <- sweep(dat[, currSamples], 1L, avg, "-")
      xprtDat <- list(dat = as.data.frame(dat),
                      A = as.data.frame(A),
                      M = as.data.frame(M))
      tmpFl <- tempfile(fileext = ".RDS")
      readr::write_rds(xprtDat, tmpFl)
      clusterExport(parClust, list("tmpFl", "currSamples"), envir = environment())
      invisible(clusterCall(parClust, \() {
        library(mgcv)
        xprtDat <- readr::read_rds(tmpFl)
        assign("dat_", xprtDat$dat, envir = .GlobalEnv)
        assign("A_", xprtDat$A, envir = .GlobalEnv)
        assign("M_", xprtDat$M, envir = .GlobalEnv)
        return()
      }))
      unlink(tmpFl)
      tmpDat1b <- setNames(parLapply(parClust, currSamples, \(smpl) { #smpl <- currSamples[1L]
        A <- A_[[smpl]]
        M <- M_[[smpl]]
        res <- A * NA_real_ # Safety: start from pure NAs...
        w <- which(is.finite(A))
        A <- A[w]
        M <- M[w]
        n <- min(c(max(c(3L,
                         round(length(unique(A))/100))),
                   10L))
        fit <- mgcv::gam(M ~ s(A, bs = "cs", k = n),
                         method = "REML")
        # NB:
        # consider replacing with
        #   fit <- mgcv::bam(M ~ s(A, bs = "cs"), method = "fREML")
        # or
        #   fit <- mgcv::bam(M ~ s(A, bs = "cs"), method = "fREML", discrete = TRUE)
        # as alternatives if gam() takes too long.
        #
        trend <- predict(fit, newdata = data.frame(A = A))
        res[w] <- dat_[w, smpl] - trend # ... then add corrected data in the right place
        return(res)
      }), currSamples)
      tmpDat1b <- do.call(cbind, tmpDat1b)
    }
    tmpDat2[grpMtch, currSamples] <- tmpDat1b[, currSamples]
    # mean SD plots
    if (normMeth == "VSN") {
      sd1 <- vsn::meanSdPlot(as.matrix(tmpDat1[grpMtch, currSamples]), plot = FALSE)$gg
      sd2 <- vsn::meanSdPlot(as.matrix(tmpDat1b), plot = FALSE)$gg
      # - harmonize scales
      b1 <- ggplot_build(sd1)
      b2 <- ggplot_build(sd2)
      xLm <- c(b1$data[[1L]]$x, b2$data[[1L]]$x)
      yLm <- c(b1$data[[1L]]$y, b2$data[[1L]]$y)
      xLm <- xlim(min(xLm, na.rm = TRUE),
                  max(xLm, na.rm = TRUE))
      yLm <- ylim(min(yLm, na.rm = TRUE),
                  max(yLm, na.rm = TRUE))
      sd1 <- sd1 + xLm + yLm
      sd2 <- sd2 + xLm + yLm
      # - save plots
      SDttl <- "mean SD plot"
      if (length(NormGrps$Group) > 1L) { SDttl <- paste0(SDttl, " - ", lGrp) }
      SDttlB <- SDttlA <- SDttl
      SDttlA <- paste0(SDttlA, "_before")
      SDttlB <- paste0(SDttlB, "_after")
      SDplotA <- sd1 + theme_bw() + ggtitle(SDttl, subtitle = "Before")
      SDplotB <- sd2 + theme_bw() + ggtitle("", subtitle = "After")
      SDplot <- ggpubr::ggarrange(SDplotA, SDplotB, ncol = 2L, nrow = 1L)
      #poplot(SDplot, 6L, 12L)
      ggsave(paste0(shpDr, "/", SDttl, ".jpeg"), SDplot, dpi = 75L, width = 12L, height = 6L, units = "in")
      ggsave(paste0(shpDr, "/", SDttl, ".pdf"), SDplot, dpi = 75L, width = 12L, height = 6L, units = "in")
      ReportCalls <- AddPlot2Report(Plot = SDplot, Title = SDttl, Space = FALSE, Dir = shpDr)
    }
  }
}, silent = TRUE)
if (!inherits(tstNorm, "try-error")) {
  #
  wAG2 <- wAG1
  # Visualize
  stateLev <- c("before", "after")
  filtFun <- \(x) { x[which(is.finite(x))] }
  tmpAB <- lapply(1L:2L, \(i) { #i <- 1L
    if (i == 1L) { tmpDat <- tmpDat1[wAG1, currSamples] }
    if (i == 2L) { tmpDat <- tmpDat2[wAG2, currSamples] }
    colnames(tmpDat) <- smplNms <- cleanNms(currSamples)
    #w <- which(apply(tmpDat[, smplNms], 1L, \(x) { sum(is.finite(x)) })/length(smpls) >= 0.8) # If we wanted to filter for completeness
    w <- 1L:nrow(tmpDat)
    avg <- as.matrix(tmpDat[w, smplNms])
    avg <- rowMeans(replace(avg, !is.finite(avg), NA), na.rm = TRUE)
    tmpDat_M <- sweep(tmpDat[w, smplNms], 1L, avg, "-")/log10(2L)
    tmpDat_A <- sweep(tmpDat[w, smplNms], 1L,  avg, "+")/2
    tmpDat_M <- dfMelt(tmpDat_M, ColNames = c("Sample", "M (mean log2 FC)"))
    tmpDat_A <- dfMelt(tmpDat_A, ColNames = c("Sample", "A (mean log10 Intensity)"))
    tmpDat <- tmpDat_A
    tmpDat$"M (mean log2 FC)" <- tmpDat_M$"M (mean log2 FC)"
    rm(tmpDat_A, tmpDat_M)
    tmpDat <- tmpDat[which(is.finite(tmpDat$"M (mean log2 FC)")),]
    tmpDat <- tmpDat[which(is.finite(tmpDat$"A (mean log10 Intensity)")),]
    tmpDat[, c("Sample_group", "Replicate")] <- Exp.map[match(tmpDat$Sample, Exp.map$Clean_name),
                                                        c(VPAL$column, "Replicate")]
    tmpDat$state <- stateLev[i]
    tmpDat$state <- factor(tmpDat$state, levels = stateLev)
    annot_ <- data.frame(Sample = unique(tmpDat$Sample))
    annot_$Median <- vapply(annot_$Sample, \(x) {
      paste0("Median: ", round(median(filtFun(tmpDat$"M (mean log2 FC)"[which(tmpDat$Sample == x)])), 3L))
    }, "")
    annot_$IQR <- vapply(annot_$Sample, \(x) {
      paste0("IQR: ", round(IQR(filtFun(tmpDat$"M (mean log2 FC)"[which(tmpDat$Sample == x)])), 3L))
    }, "")
    annot_$Amax <- max(filtFun(tmpDat$"A (mean log10 Intensity)"))*1.1
    annot_$Amin <- min(filtFun(tmpDat$"A (mean log10 Intensity)"))*1.1
    annot_$Mmax <- max(filtFun(tmpDat$"M (mean log2 FC)"))*1.1
    annot_$Mmin <- min(filtFun(tmpDat$"M (mean log2 FC)"))*1.1
    annot <- annot_[, c("Sample", "Amax", "Mmin", "Mmax")] 
    annot <- rbind(annot, annot)
    annot$Label <- c(annot_$Median, annot_$IQR)
    annot$state <- stateLev[i]
    annot$state <- factor(annot$state, levels = stateLev)
    annot[, c("Sample_group", "Replicate")] <- Exp.map[match(annot$Sample, Exp.map$Clean_name),
                                                       c(VPAL$column, "Replicate")]
    return(list(data = tmpDat,
                annot = annot,
                annot_ = annot_))
  })
  datAB <- rbind(tmpAB[[1L]]$data,
                 tmpAB[[2L]]$data)
  annotA <- tmpAB[[1L]]$annot
  annotB <- tmpAB[[2L]]$annot
  annotA_ <- tmpAB[[1L]]$annot_
  annotB_ <- tmpAB[[2L]]$annot_
  ylim <- max(c(abs(c(annotA_$Mmax, annotA_$Mmin, annotB_$Mmax, annotB_$Mmin,
                      (annotA_$Amax-annotA_$Amin)/4,
                      (annotB$Amax-annotB$Amin)/4))))
  annotA$Y <- annotB$Y <- ylim*0.9
  annotA$Y[grep("^IQR: ", annotA$Label)] <- annotB$Y[grep("^IQR: ", annotA$Label)] <- -ylim*0.9
  annotAB <- rbind(annotA, annotB)
  # Consider drawing a binned arrows plot (up to 250-ish arrows, each corresponding to a local bin))
  # showing how the data is shifted by the transformation? 
  MAttl <- paste0(normMeth, " correction")
  MAfl <- paste0(shpDr, "/", MAttl)
  MAplot <- ggplot(datAB[1L:5L,]) +
    ggtitle(MAttl) +
    facet_grid(Replicate + state ~ Sample_group) + coord_fixed(0.05) + theme_bw() +
    geom_hex(aes(x = `A (mean log10 Intensity)`, y = `M (mean log2 FC)`, fill = state), linewidth = 1L, alpha = 1, bins = 100L) +
    geom_hline(yintercept = 0, colour = "grey") +
    scale_fill_viridis(begin = 0.25, end = 0.75, discrete = TRUE, option = "H") +
    theme(legend.position = "bottom",
          plot.margin = margin(0L, 0L, 0L, 0L)) +
    new_scale_color()
  #poplot(MAplot, 12L, 22L)
  MAplot <- if (normMeth == "LOESS") {
    MAplot +
      geom_smooth(aes(x = `A (mean log10 Intensity)`, y = `M (mean log2 FC)`, color = "LOESS", linetype = "LOESS"),
                  method = "loess", formula = y ~ x, se = FALSE)
  } else {
    MAplot +
      geom_smooth(aes(x = `A (mean log10 Intensity)`, y = `M (mean log2 FC)`, color = "GAM", linetype = "GAM"),
                  method = "gam", formula = y ~ s(x, bs = "cs"))
  }
  MAplot <- MAplot +
    scale_color_manual(name = "Method", values = c(LOESS = "green", GAM = "magenta")) +
    scale_linetype_manual(name = "Method", values = c(LOESS = "dotted", GAM = "dashed"))
  #
  tmpFl <- tempfile(fileext = ".RDS")
  readr::write_rds(datAB, tmpFl)
  clusterExport(parClust,
                list("MAplot", "tmpFl", "MAttl", "MAfl", "annotAB"),
                envir = environment())
  invisible(clusterCall(parClust, \() {
    assign("datAB", readr::read_rds(tmpFl), envir = .GlobalEnv)
    return()
  }))
  unlink(tmpFl)
  tmpMAplotFls <- setNames(parSapply(parClust, levels(datAB$Sample), \(smpl) { #smpl <- levels(datAB$Sample)[1L] 
    dat <- datAB[which(datAB$Sample == smpl),]
    ann <- annotAB[which(annotAB$Sample == smpl),]
    MAplot$data <- dat
    ttl <- paste0(MAttl, " - ", smpl)
    MAplot <- MAplot + ggtitle(ttl) +
      geom_text(data = ann, aes(x = Amax, y = Y, label = Label), hjust = 1, size = 3L)
    #poplot(MAplot, 12L, 22L)
    fl <- paste0(MAfl, " - ", smpl)
    ggsave(paste0(fl, ".jpeg"), MAplot, width = 10L, units = "in")
    ggsave(paste0(fl, ".pdf"), MAplot, width = 10L, units = "in")
    return(fl)
  }), levels(datAB$Sample))
  #
  appNm <- paste0(normMeth, " normalisation")
  msg <- paste0("Accept ", normMeth, " normalisation? (untick to cancel correction)")
  if ("Decision" %in% (normSequence[[nrmStp]])) {
    KeepShapeCorrRes <- normSequence[[nrmStp]]$Decision
  }
  if (!validLogicPar("KeepShapeCorrRes")) {
    KeepShapeCorrRes <- TRUE # Default: keep results only if data is very significantly different
  }
  if (exists("IHAVERUN")) { rm(IHAVERUN) }
  ui <- fluidPage(
    tags$head(
      tags$style(HTML("#MAplots img { min-width: 95%; height: auto; }")),
      tags$style(HTML("#SDplots img { min-width: 95%; height: auto; }"))
    ),
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
    fluidRow(column(5L,
                    checkboxInput("KeepResults", msg, KeepShapeCorrRes),
                    actionBttn("saveBtn", "Save", icon = icon("save"), color = "success", style = "pill"),
                    h4("Recommended criteria:"),
                    h5(HTML("&nbsp;Does the original MA plot look like it is skewed?")),
                    h5(HTML("&nbsp;&nbsp;-> If yes: accept the correction if the correction did remove the skew.")))),
    br(),
    fluidRow(column(12L,
                    selectInput("my_MA_plot", "", names(tmpMAplotFls), names(tmpMAplotFls)[1L]),
                    imageOutput("MAplot", inline = TRUE)),
             if (normMeth == "VSN") {
               fluidRow(column(12L,
                               imageOutput("SDplots", inline = TRUE)))
             },
    ),
    br(),
    br()
  )
  server <- function(input, output, session) {
    updtMA <- \(reactive = TRUE) {
      myMA <- { if (reactive) { input$my_MA_plot } else { names(tmpMAplotFls)[1L] } }
      renderImage(list(src = paste0(tmpMAplotFls[myMA], ".jpeg"), # For now any MA plots will do if there are several
                       contentType = "image/jpeg",
                       height = c("700px", "500px")[(normMeth == "VSN")+1L]),
                  deleteFile = FALSE)
    }
    output$MAplot <- updtMA(FALSE)
    #
    observeEvent(input$my_MA_plot, {
      output$MAplot <- updtMA()
    }, ignoreNULL = FALSE)
    #
    if (normMeth == "VSN") {
      output$SDplots <- renderPlot(SDplot, height = "500px")
    }
    observeEvent(input$saveBtn, {
      assign("KeepShapeCorrRes", as.logical(input[["KeepResults"]]), envir = .GlobalEnv)
      assign("IHAVERUN", TRUE, .GlobalEnv)
      stopApp()
    })
    session$onSessionEnded(function() { stopApp() })
  }
  while (!exists("IHAVERUN")) {
    eval(parse(text = runApp), envir = .GlobalEnv)
    shinyCleanup()
  }
  msg <- paste0(" -> ", normMeth, " correction for intensity range variance biases ",
                c("rejec", "accep")[KeepShapeCorrRes+1L], "ted.\n")
  if (KeepShapeCorrRes) {
    txt2 <- paste0("corrected for intensity range variance biases using ",
                   c(paste0(normMeth, " regression"), "VSN")[match(normMeth, c("LOESS", "VSN"))])
  }
  normSequence[[nrmStp]]$Decision <- Outcome <- KeepShapeCorrRes
  cat(msg)
} else {
  msg <- paste0(" -> ", normMeth, " correction for intensity range variance biases failed!\n")
  normSequence[[nrmStp]]$Decision <- Outcome <- FALSE
  cat(msg)
}

