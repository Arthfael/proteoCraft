#### Re-normalize protein group expression values
# In this revised version, we:
# a) Define nrmFlt, a normalization filter (i.e. which proteins are used to define the "center") for the relevant type of re-normalization:
#   - Norma.Prot.Ratio.classic -> all proteins,
#   - Norma.Prot.Ratio.to.Biot -> biotinylated proteins,
#   - Norma.Prot.Ratio.to.proteins -> user-defined list of proteins.
# b) Use this to normalize - within normalization groups - using:
#   - just the median or, ...
#   - ... if Norma.Prot.Ratio.Adv is TRUE, Levenberg-Marquardt.
# c) Add a small random error to singular rows (e.g. for when our filter only includes a single protein).
# d) Re-establish the original scale (check that the median of all values is unchanged after normalization).
# e) Re-calculate log ratios.
#
normPGs <- 0L
Norma.Prot.Ratio.classic %<o% FALSE
Norma.Prot.Ratio.to.Biot %<o% FALSE
Norma.Prot.Ratio.to.proteins %<o% FALSE
Norma.Prot.Ratio.Adv %<o% FALSE
if ((!exists("accept_PG_reNorm"))||(!is.logical(accept_PG_reNorm))||(length(accept_PG_reNorm) != 1L)) { accept_PG_reNorm <- TRUE }
accept_PG_reNorm %<o% accept_PG_reNorm
if (("Norma.Prot.Ratio" %in% colnames(Param))&&(Param$Norma.Prot.Ratio)) {
  dirPG <- paste0(wd, "/Workflow control/Protein groups/Renorm")
  if (!dir.exists(dirPG)) { dir.create(dirPG, recursive = TRUE) }
  dirlist <- unique(c(dirlist, dirPG))
  dirPep <- paste0(wd, "/Workflow control/Peptides/PG-based renorm")
  if (!dir.exists(dirPep)) { dir.create(dirPep, recursive = TRUE) }
  dirlist <- unique(c(dirlist, dirPep))
  #
  quantData <- quantData_list$Data
  quantData_norm <- quantData
  mySamples <- Exp.map$Ref.Sample.Aggregate
  ExpKol <- paste0(Prot.Expr.Root, mySamples)
  #
  if (!"Adv.Norma.Prot.Intens" %in% colnames(Param)) {
    ObjNm <- "Norma.Prot.Ratio.Adv"
    if ((ReUseAnsw)&&(ObjNm %in% AllAnsw$Parameter)) { ObjNm %<c% AllAnsw$Value[[match(ObjNm, AllAnsw$Parameter)]] } else {
      msg <- "Do you want to apply Levenberg-Marquardt (slow) re-normalisation at protein groups level?"
      tmp <- c(TRUE, FALSE)[match(dlg_message(msg, "yesno")$res, c("yes", "no"))]
      if (is.na(tmp)) { tmp <- FALSE }
      ObjNm %<c% tmp
      tmp <- AllAnsw[1L,]
      tmp[, c("Parameter", "Message")] <- c(ObjNm, msg)
      tmp$Value <- list(get(ObjNm))
      m <- match(ObjNm, AllAnsw$Parameter)
      if (is.na(m)) { AllAnsw <- rbind(AllAnsw, tmp) } else { AllAnsw[m,] <- tmp }
    }
  } else { Norma.Prot.Ratio.Adv <- Param$Adv.Norma.Prot.Intens }
  #
  # a) Define normalization filters
  Norma.Prot.Ratio.classic <- TRUE # The default
  if ((IsBioID2)&&("Norma.Prot.Ratio.to.Biot" %in% colnames(Param))&&
      (is.logical(Param$Norma.Prot.Ratio.to.Biot))&&
      (length(Param$Norma.Prot.Ratio.to.Biot) == 1)&&
      (Param$Norma.Prot.Ratio.to.Biot)) {
    nrmFlt <- which(PG$"Biot. peptides count" > 0L)
    if (length(w)) {
      Norma.Prot.Ratio.classic <- Norma.Prot.Ratio.to.proteins <- FALSE
      Norma.Prot.Ratio.to.Biot <- TRUE
      txTmp <- "biotinylated protein groups"
    }
  }
  if (("Norma.Prot.Ratio.to.proteins" %in% colnames(Param))&&
      (is.character(Param$Norma.Prot.Ratio.to.proteins))&&
      (!Param$Norma.Prot.Ratio.to.proteins %in% c("", "NA"))&&
      (!Norma.Prot.Ratio.to.Biot)) {
    Prot.Ratio.ref.Acc %<o% unique(unlist(strsplit(Param$Norma.Prot.Ratio.to.proteins, ";")))
    Prot.Ratio.ref.Acc <- Prot.Ratio.ref.Acc[which(Prot.Ratio.ref.Acc %in% db$"Protein ID")]
    if (length(Prot.Ratio.ref.Acc)) {
      Prot.Ratio.ref.Acc <- Prot.Ratio.ref.Acc[which(Prot.Ratio.ref.Acc %in% unique(unlist(strsplit(PG$"Protein IDs", ";"))))] # Having "Leading protein IDs" here was too stringent
      Norma.Prot.Ratio.classic <- FALSE
      l <- length(Prot.Ratio.ref.Acc)
      if (l) {
        temp <- strsplit(PG$"Protein IDs", ";")
        temp <- listMelt(temp, PG$id)
        temp$Norm <- temp$value %in% Prot.Ratio.ref.Acc
        temp <- temp$L1[which(temp$Norm)]
        nrmFlt <- which(PG$id %in% temp)
        if (length(nrmFlt)) {
          Norma.Prot.Ratio.classic <- Norma.Prot.Ratio.to.Biot <- FALSE
          Norma.Prot.Ratio.to.proteins <- TRUE
          txTmp <- paste0("the following proteins: ",
                          paste(db$`Full Name`[match(Prot.Ratio.ref.Acc, db$`Protein ID`)],
                                collapse = "-"))
        }
      }
    }
  }
  if (Norma.Prot.Ratio.classic) {
    nrmFlt <- 1L:nrow(quantData)
    Norma.Prot.Ratio.to.Biot <- Norma.Prot.Ratio.to.proteins <- FALSE
    txTmp <- "all protein groups"
  }
  #
  # nrmFlt is the rows filter we will use to ensure that we normalize to the correct proteins:
  # - Norma.Prot.Ratio.classic: the filter includes all rows
  # - Norma.Prot.Ratio.to.Biot: the filter includes biotinylated proteins
  # - Norma.Prot.Ratio.to.proteins: the filter includes specific user-defined protein(s)
  #
  normPGs <- sum(c(Norma.Prot.Ratio.to.Biot, Norma.Prot.Ratio.to.proteins, Norma.Prot.Ratio.classic))
  if (normPGs) {
    stopifnot(normPGs == 1L) # Only one method may run - otherwise it makes no sense!
    #
    # b) Normalize within groups
    grpKols <- setNames(lapply(RG$values, \(grp) { #grp <- RG$values[[1L]]
      em <- Exp.map[which(Exp.map[[RG$column]] == grp),]
      smpls <- em$Ref.Sample.Aggregate
      xpKol <- paste0(Prot.Expr.Root, smpls)
      w <- which(xpKol %in% colnames(quantData_norm))
      setNames(xpKol[w], smpls[w])
    }), RG$values)
    grpNorm_med <- setNames(lapply(RG$values, \(grp) {
      kol <- grpKols[[grp]]
      nms <- names(kol)
      dat <- quantData_norm[nrmFlt, kol]
      xpMed1 <- apply(dat[, kol], 2L, median, na.rm = TRUE)
      xpMed1 <- mean(xpMed1)-xpMed1
      return(setNames(xpMed1, nms))
    }), RG$values)
    mean_med <- mean(unlist(grpNorm_med))
    grpNorm_med <- setNames(lapply(RG$values, \(grp) { grpNorm_med[[grp]] - mean_med }), RG$values)
    if (Norma.Prot.Ratio.Adv) {
      tmpFl <- tempfile(fileext = ".rds")
      readr::write_rds(quantData_norm, tmpFl)
      ids <- PG$id[nrmFlt]
      exports <- list("tmpFl", "Exp.map", "RG", #"RRG",
                      "nrmFlt", "is.all.good", "Prot.Expr.Root", #"Prot.Rat.Root",
                      "ids", "AdvNorm.IL", "grpKols")
      clusterExport(parClust, exports, envir = environment())
      invisible(clusterCall(parClust, \() {
        quantData_norm <- readr::read_rds(tmpFl)
        assign("quantData_norm", quantData_norm, envir = .GlobalEnv)
        return()
      }))
      unlink(tmpFl)
      grpNorm_LM <- setNames(parLapply(parClust, RG$values, \(grp) {
        kol <- grpKols[[grp]]
        nms <- names(kol)
        dat <- quantData_norm[nrmFlt, kol]
        dat$id <- ids
        xpMed0 <- apply(dat[, kol], 2L, median, na.rm = TRUE)
        nrmDat <- AdvNorm.IL(dat, "id", kol, TRUE, 5L)
        nrmKol <- paste0("AdvNorm.", kol)
        xpMed2 <- apply(nrmDat[, nrmKol], 2L, median, na.rm = TRUE)
        xpMed2 <- xpMed2-xpMed0
        xpMed2 <- mean(xpMed2)-xpMed2
        return(setNames(xpMed2, nms))
      }), RG$values)
      mean_LM <- mean(unlist(grpNorm_LM))
      grpNorm_LM <- setNames(lapply(RG$values, \(grp) { grpNorm_LM[[grp]] - mean_LM }), RG$values)
      grpNorm <- grpNorm_LM
    } else {
      grpNorm <- grpNorm_med
    }
    # To check the correlation between both methods, and the effect of experimental factors on the discrepancy
    # (this can reveal insights into unexpected effects/trends in the data)
    #compareNorms <- FALSE
    #if (compareNorms) {
      tmp1 <- do.call(c, setNames(grpNorm_LM, NULL))
      tmp2 <- do.call(c, setNames(grpNorm_med, NULL))
      tst <- data.frame(sample = names(tmp1),
                        LM = setNames(tmp1, NULL))
      tst$median <- tmp2[tst$sample]
      tst[, RSA$names] <- Exp.map[match(tst$sample, Exp.map$Ref.Sample.Aggregate), RSA$names]
      nms <- VPAL$names
      if (length(Exp) == 1L) { nms <- nms[which(nms != "Experiment")] }
      tst$Samples_group <- do.call(paste, c(tst[, nms, drop = FALSE], sep = " "))
      l <- length(nms)
      if (l > 2L) {
        nms2 <- nms[3L:l]
        tst$Label <- do.call(paste, c(tst[, c(nms2, "Replicate")], sep = " "))
      }
      plot <- ggplot(tst)
      if (l == 1L) {
        plot <- plot + geom_point(aes(x = median, y = LM, shape = .data[[nms[1]]], color = Replicate))
      } else {
        plot <- plot + geom_point(aes(x = median, y = LM, shape = .data[[nms[1]]], color = .data[[nms[2]]]))
        plot <- if (l > 2L) {
          nms2 <- nms[3L:l]
          tst$Label <- do.call(paste, c(tst[, c(nms2, "Replicate")], sep = " "))
          plot + geom_text_repel(aes(x = median, y = LM, label = Label))
        } else {
          plot + geom_text_repel(aes(x = median, y = LM, label = Replicate))
        }
      }
      plot <- plot + ggtitle("Median normalisation (x) vs Levenberg-Marquardt (y)") +
        geom_abline(intercept = 0, slope = 1) + coord_fixed() +
        theme_bw()
      #poplot(plot)
      ggsave(paste0(dirPG, "/Effect of LM normalisation.jpeg"), plot, dpi = 100L)
      ggsave(paste0(dirPG, "/Effect of LM normalisation.pdf"), plot, dpi = 100L)
      #}
    #
    smplsNorm <- do.call(c, setNames(grpNorm, NULL))
    quantData_norm[, ExpKol] <- sweep(quantData_norm[, ExpKol],
                                      2L,
                                      smplsNorm[mySamples],
                                      "+")
    # If you are re-normalizing to a single protein, you get randomly singular values in one row.
    # This breaks stats so we need to fix that.
    #  - Expression:
    tst <- apply(quantData_norm[, ExpKol], 1L, \(x) {
      x <- x[which(!is.na(x))]
      length(unique(x))
    })
    if (1L %in% tst) {
      #aggregate(tst, list(tst), length)
      # Add a small value
      w <- which(tst > 0L)
      quantData_norm[w, ExpKol] <- quantData_norm[w, ExpKol] + runif(length(w)*length(ExpKol), min = 0, max = 1e-10)
    }
    #
    # Propagate effects to peptides
    pep.ref["Back-norm"] <- paste0("back-norm. ", pep.ref["Original"])
    prevRef <- match("Back-norm", names(pep.ref)) - 1L
    pepQuantData_norm <- sweep(pep[, paste0(pep.ref[prevRef], mySamples)],
                               2L,
                               10L^smplsNorm[mySamples],
                               "*")
    colnames(pepQuantData_norm) <- paste0(pep.ref["Back-norm"], mySamples)
    #
    ttl <- paste0("re-normalisation",
                  c(paste0(" (to ", c("biotinylated",
                                      "user-defined invariant"), " proteins)"),
                    "")[which(c(Norma.Prot.Ratio.to.Biot,
                                Norma.Prot.Ratio.to.proteins,
                                Norma.Prot.Ratio.classic))])
    l <- length(DatAnalysisTxt)
    insrt <- if (Norma.Prot.Ratio.Adv) { "using the Levenberg-marquart procedure to minimize the difference between " } else { "to the median of " }
    DatAnalysisTxt[l] <- gsub("\\.\\.\\.$",
                              paste0(" and re-normalized ", insrt, txTmp, "."),
                              DatAnalysisTxt[l])
    ReportCalls$Calls <- AddTxt2Report(paste0("Re-normalizing Protein groups-level expression values to the median value of ",
                                              txTmp))
  }
}
# Calculate protein ratios (now calculated here, not earlier as in earlier version)
RatKol <- paste0(Prot.Rat.Root, myContrasts$Contrast)
tmp1 <- make_Rat2(quantData,
                  experiment.map = Exp.map,
                  int.root = Prot.Expr.Root,
                  rat.root = Prot.Rat.Root)
quantData[, RatKol] <- tmp1[, RatKol]
#
if (normPGs) {
  tmp2 <- make_Rat2(quantData_norm,
                    experiment.map = Exp.map,
                    int.root = Prot.Expr.Root,
                    rat.root = Prot.Rat.Root)
  #summary(tmp1)
  #summary(tmp2)
  #summary(tmp1/log2(10))
  #summary(tmp2/log2(10))
  quantData_norm[, RatKol] <- tmp2[, RatKol]
  # Now check and visualize
  nrmPlots <- list()
  # - PGs
  #   - Expression values
  Samples <- cleanNms(mySamples)
  temp1 <- quantData[, ExpKol]
  colnames(temp1) <- Samples
  temp1 <- dfMelt(temp1)
  temp2 <- quantData_norm[, ExpKol]
  colnames(temp2) <- Samples
  temp2 <- dfMelt(temp2)
  temp1$Norm <- "Original"
  temp2$Norm <- "Normalised"
  temp <- rbind(temp1, temp2)
  rm(temp1, temp2)
  temp$Sample <- factor(temp$variable, levels = Samples)
  temp$variable <- NULL
  temp$Norm <- factor(temp$Norm, levels = c("Original", "Normalised"))
  temp <- temp[which(is.all.good(temp$value, 2L)),]
  ttlI1 <- paste0("Protein groups ", ttl, ", expression")
  intPlot1 <- ggplot(temp) +
    geom_violin(aes(x = Sample, y = value, colour = Sample, fill = Sample), alpha = 0.25) +
    geom_boxplot(aes(x = Sample, y = value, colour = Sample, fill = Sample), alpha = 0.5) +
    scale_color_viridis_d(begin = 0.25) +
    scale_fill_viridis_d(begin = 0.25) +
    facet_grid(Norm~.) + ggtitle(ttlI1) +
    theme_bw() + theme(axis.text.x = element_text(angle = 75, vjust = 1, hjust = 1))
  #print(intPlot1)
  nrmPlots[["PG_int"]] <- list(Path = paste0(dirPG, "/", ttlI1),
                               Plot = plotEval(intPlot1),
                               Ext = "jpeg")
  #suppressMessages({
  #  ggsave(paste0(dirPG, "/", ttlI1, ".jpeg"), intPlot1, dpi = 300L)
  #  ggsave(paste0(dirPG, "/", ttlI1, ".pdf"), intPlot1, dpi = 300L)
  #})
  #ReportCalls <- AddPlot2Report(Plot = intPlot1, Title = ttlI1)
  #
  #   - Ratios
  allContr <- gsub(topattern(Prot.Rat.Root), "", RatKol)
  temp1 <- data.table(quantData[, RatKol])
  colnames(temp1) <- allContr
  temp1 <- dfMelt(temp1)
  temp2 <- data.table(quantData_norm[, RatKol])
  colnames(temp2) <- allContr
  temp2 <- dfMelt(temp2)
  temp1$Norm <- "Original"
  temp2$Norm <- "Normalised"
  temp <- rbind(temp1, temp2)
  rm(temp1, temp2)
  allContr2 <- gsub(" - ", " -\n", allContr)
  temp$Contrast <- gsub_Rep(" - ", " -\n", temp$variable)
  temp$Contrast <- factor(temp$Contrast, levels = allContr2)
  temp$variable <- NULL
  temp$Norm <- factor(temp$Norm, levels = c("Original", "Normalised"))
  temp <- temp[which(is.all.good(temp$value, 2L)),]
  ttlR1 <- paste0("Protein groups ", ttl, ", ratios")
  ratPlot1 <- ggplot(temp) +
    geom_violin(aes(x = Contrast, y = value, colour = Contrast, fill = Contrast), alpha = 0.25) +
    geom_boxplot(aes(x = Contrast, y = value, colour = Contrast, fill = Contrast), alpha = 0.5) +
    scale_color_viridis_d(begin = 0.25) +
    scale_fill_viridis_d(begin = 0.25) +
    facet_grid(Norm~.) + ggtitle(ttlR1) +
    theme_bw() + theme(axis.text.x = element_text(angle = 75, vjust = 1, hjust = 1))
  #print(ratPlot1)
  nrmPlots[["PG_rat"]] <- list(Path = paste0(dirPG, "/", ttlR1),
                               Plot = plotEval(ratPlot1),
                               Ext = "jpeg")
  #suppressMessages({
  #  ggsave(paste0(dirPG, "/", ttlR1, ".jpeg"), ratPlot1, dpi = 300L)
  #  ggsave(paste0(dirPG, "/", ttlR1, ".pdf"), ratPlot1, dpi = 300L)
  #})
  #ReportCalls <- AddPlot2Report(Plot = ratPlot1, Title = ttlR1)
  #
  # - Peptides
  #   - Intensities
  temp1 <- data.table(pep[, paste0(pep.ref[prevRef], mySamples)])
  colnames(temp1) <- Samples
  temp1 <- dfMelt(temp1)
  temp2 <- data.table(pepQuantData_norm[, paste0(pep.ref["Back-norm"], mySamples)])
  colnames(temp2) <- Samples
  temp2 <- dfMelt(temp2)
  temp1$Norm <- "Original"
  temp2$Norm <- "Normalised"
  temp <- rbind(temp1, temp2)
  rm(temp1, temp2)
  temp$Sample <- factor(temp$variable, levels = Samples)
  temp$variable <- NULL
  temp$Norm <- factor(temp$Norm, levels = c("Original", "Normalised"))
  temp$value <- suppressWarnings(log10(temp$value)) # Peptide intensities are not log-transformed!
  temp <- temp[which(is.all.good(temp$value, 2L)),]
  ttlI2 <- paste0("Peptides ", ttl, " from PGs, intensity")
  intPlot2 <- ggplot(temp) +
    geom_violin(aes(x = Sample, y = value, colour = Sample, fill = Sample), alpha = 0.25) +
    geom_boxplot(aes(x = Sample, y = value, colour = Sample, fill = Sample), alpha = 0.5) +
    scale_color_viridis_d(begin = 0.25) +
    scale_fill_viridis_d(begin = 0.25) +
    facet_grid(Norm~.) + ggtitle(ttlI2) +
    theme_bw() + theme(axis.text.x = element_text(angle = 75, vjust = 1, hjust = 1))
  #print(intPlot2)
  nrmPlots[["Pep_int"]] <- list(Path = paste0(dirPep, "/", ttlI2),
                                Plot = plotEval(intPlot2),
                                Ext = "jpeg")
  #suppressMessages({
  #  ggsave(paste0(dirPep, "/", ttlI2, ".jpeg"), intPlot2, dpi = 300L)
  #  ggsave(paste0(dirPep, "/", ttlI2, ".pdf"), intPlot2, dpi = 300L)
  #})
  #ReportCalls <- AddPlot2Report(Plot = intPlot2, Title = ttlI2)
  #
  # - Peptides ratios
  temp1 <- make_Rat2(int.log = FALSE,
                     experiment.map = Exp.map,
                     int.root = pep.ref[prevRef],
                     rat.root = "")
  temp2 <- make_Rat2(pepQuantData_norm,
                     int.log = FALSE,
                     experiment.map = Exp.map,
                     int.root = pep.ref["Back-norm"],
                     rat.root = "")
  temp1 <- dfMelt(temp1)
  temp2 <- dfMelt(temp2)
  temp1$Norm <- "Original"
  temp2$Norm <- "Normalised"
  temp <- rbind(temp1, temp2)
  rm(temp1, temp2)
  temp$Contrast <- gsub_Rep(" - ", " -\n", temp$variable)
  temp$Contrast <- factor(temp$Contrast, levels = allContr2)
  temp$variable <- NULL
  temp$Norm <- factor(temp$Norm, levels = c("Original", "Normalised"))
  temp <- temp[which(is.all.good(temp$value, 2L)),]
  ttlR2 <- paste0("Peptides ", ttl, " from PGs, ratios")
  ratPlot2 <- ggplot(temp) +
    geom_violin(aes(x = Contrast, y = value, colour = Contrast, fill = Contrast), alpha = 0.25) +
    geom_boxplot(aes(x = Contrast, y = value, colour = Contrast, fill = Contrast), alpha = 0.5) +
    scale_color_viridis_d(begin = 0.25) +
    scale_fill_viridis_d(begin = 0.25) +
    facet_grid(Norm~.) + ggtitle(ttlR2) +
    theme_bw() + theme(axis.text.x = element_text(angle = 75, vjust = 1, hjust = 1))
  #print(ratPlot2)
  nrmPlots[["Pep_rat"]] <- list(Path = paste0(dirPep, "/", ttlR2),
                                Plot = plotEval(ratPlot2),
                                Ext = "jpeg")
  #suppressMessages({
  #  ggsave(paste0(dirPep, "/", ttlR2, ".jpeg"), ratPlot2, dpi = 300L)
  #  ggsave(paste0(dirPep, "/", ttlR2, ".pdf"), ratPlot2, dpi = 300L)
  #})
  #ReportCalls <- AddPlot2Report(Plot = ratPlot2, Title = ttlR2)
  #
  nrmPlots2 <- lapply(nrmPlots, \(x) {
    x$Ext <- "pdf"
    return(x)
  })
  nrmPlots <- c(nrmPlots, nrmPlots2)
  # Save plots
  tst <- parLapply(parClust, nrmPlots, \(x) {
    suppressMessages({
      ggplot2::ggsave(paste0(x$Path, ".", x$Ext), x$Plot, dpi = 300L)
    })
  })
  #
  msg <- "Accept re-normalisation results?"
  appNm <- paste0(dtstNm, " - PG re-normalisation")
  ui <- fluidPage(
    useShinyjs(),
    setBackgroundColor( # Doesn't work
      color = c(#"#F8F8FF",
        "#E1F7E7"),
      gradient = "linear",
      direction = "bottom"
    ),
    extendShinyjs(text = jsToggleFS, functions = c("toggleFullScreen")),
    titlePanel(tag("u", "Protein Groups (and peptides) re-normalisation"),
               appNm),
    br(),
    checkboxInput("accept_PG_reNorm", msg, accept_PG_reNorm),
    actionBttn("saveBtn", "Save", icon = icon("save"), color = "success", style = "pill"),
    # fluidRow(column(4L, withSpinner(plotOutput("pgIntensities"))),
    #          column(4L, withSpinner(plotOutput("pgRatios")))),
    #br(),
    # fluidRow(column(4L, withSpinner(plotOutput("pepIntensities"))),
    #          column(4L, withSpinner(plotOutput("pepRatios")))),
    fluidRow(column(6L, withSpinner(imageOutput("pgIntensities", width = "100%", height = "100%"))),
             column(6L, withSpinner(imageOutput("pgRatios", width = "100%", height = "100%")))),
    br(),
    br(),
    br(),
    fluidRow(column(6L, withSpinner(imageOutput("pepIntensities", width = "100%", height = "100%"))),
             column(6L, withSpinner(imageOutput("pepRatios", width = "100%", height = "100%")))),
    br(),
    br()
  )
  IMGs <- c(paste0(dirPG, "/", c(ttlI1, ttlR1), ".jpeg"),
            paste0(dirPep, "/", c(ttlI2, ttlR2), ".jpeg"))
  IMGsDims <- as.data.frame(t(parSapply(parClust, IMGs, \(x) { #x <- IMGs[1L]
    a <- jpeg::readJPEG(x)
    setNames(dim(a)[1L:2L], c("height", "width"))
  })))
  IMGsDims$height <- screenRes$width*0.35*IMGsDims$height/max(IMGsDims$height)
  IMGsDims$width <- screenRes$width*0.35*IMGsDims$width/max(IMGsDims$width)
  server <- \(input, output, session) {
    # output$pgIntensities <- renderPlot(intPlot1)
    # output$pgRatios <- renderPlot(ratPlot1)
    # output$pepIntensities <- renderPlot(intPlot2)
    # output$pepRatios <- renderPlot(ratPlot2)
    output$pgIntensities <- renderImage({
      list(src = IMGs[1L], height = IMGsDims$height[1L], width = IMGsDims$width[1L])
    }, deleteFile = FALSE)
    output$pgRatios <- renderImage({
      list(src = IMGs[2L], height = IMGsDims$height[2L], width = IMGsDims$width[2L])
    }, deleteFile = FALSE)
    output$pepIntensities <- renderImage({
      list(src = IMGs[3L], height = IMGsDims$height[3L], width = IMGsDims$width[3L])
    }, deleteFile = FALSE)
    output$pepRatios <- renderImage({
      list(src = IMGs[4L], height = IMGsDims$height[4L], width = IMGsDims$width[4L])
    }, deleteFile = FALSE)
    #
    observeEvent(input$accept_PG_reNorm, { assign("accept_PG_reNorm", as.logical(input[["accept_PG_reNorm"]]), envir = .GlobalEnv) })
    observeEvent(input$saveBtn, { stopApp() })
    #observeEvent(input$cancel, { stopApp() })
    session$onSessionEnded(\() { stopApp() })
  }
  eval(parse(text = runApp), envir = .GlobalEnv)
  shinyCleanup()
  #
  if (accept_PG_reNorm) {
    pep[, colnames(pepQuantData_norm)] <- pepQuantData_norm
    #summary(pep[, colnames(pepQuantData_norm)])
    # if (Param$Ratios.Thresholds == threshMsg) {
    #   pep.Ref.Ratios <- pep.Ref.Ratios.norm
    #   pep[, colnames(pep.Ref.Ratios)] <- pep.Ref.Ratios
    # }
    if (quantAlgo == "limpa") {
      # For limpa, re-normalisation requires re-running the quant from back-normalized peptides,
      # so the object is properly normalized "in-depth"!
      cat("Quant algorithm = limpa -> re-running quantitation from back-normalized peptides!\n")
      post_ReNorm_reRun <- TRUE
      #rstudioapi::documentOpen(quntSrc)
      source(quntSrc, local = FALSE)
      summary(quantData_list$Data)
    } else {
      quantData <- quantData_norm
    }
  } else {
    pep.ref <- pep.ref[setdiff(names(pep.ref), "Back-norm")]
  }
  ReportCalls <- AddSpace2Report()
}
