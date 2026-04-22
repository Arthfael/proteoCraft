### PCA plot for parameters app
# Create first PCA to check on sample relationships
if ((length(MQ.Exp) > 1L)||(LabelType == "Isobaric")) { # Should be always TRUE
  source(parSrc, local = FALSE)
  data <- ev
  colnames(data)[which(colnames(data) == "MQ.Exp")] <- "Parent sample"
  data <- data[which(data$Reverse != "+"),]
  data <- data[which((is.na(data$"Potential contaminant"))|(data$"Potential contaminant" != "+")),]
  kol <- if (LabelType == "Isobaric") {
    grep(paste0(topattern(ev.ref["Original"]), "[0-9]+$"), colnames(data), value = TRUE)
  } else { ev.col["Original"] }
  w <- which(rowSums(data[, kol, drop = FALSE], na.rm = TRUE) > 0)
  data <- data[w,]
  if (!"Fraction" %in% colnames(data)) { data$Fraction <- 1L }
  Fraction <- sort(unique(data$Fraction), decreasing = FALSE)
  Experiment <- Exp
  kols <- c("Parent sample", "Fraction", "Experiment")
  if (LabelType == "Isobaric") {
    X <- "Label"
    kols <- c("Fraction", "Parent sample", "Experiment") # The order matters!
    tst <- vapply(kols, \(x) { length(unique(data[[x]])) }, 1L)
    w1 <- which(tst > 1L)
    w2 <- which(tst >= 1L) 
    Y <- if (length(w1)) { kols[w1[1L]] } else { kols[w2[1L]] }
  }
  if (LabelType == "LFQ") {
    kols <- c("Parent sample", "Fraction", "Experiment") # The order matters!
    tst <- vapply(kols, \(x) { length(unique(data[[x]])) }, 1L)
    w1 <- which(tst > 1L)
    w2 <- which(tst >= 1L) 
    X <- kols[w1[1L]]
    Y <- if (length(w1) > 1L) { kols[w1[2]] } else { kols[w2[2]] }
  }
  kols <- setdiff(kols, c(X, Y))
  ReportCalls <- AddSpace2Report()
  ReportCalls$Calls <- AddTxt2Report("PSMs-level PCA plot:")
  ReportCalls$Calls <- append(ReportCalls$Calls, list())
  dir <- paste0(wd, "/Dimensionality red. plots/PCA")
  if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
  dirlist <- unique(c(dirlist, dir))
  LRepCalls <- length(ReportCalls$Calls)
  lsKl <- c("Modified sequence", Y)
  if (LabelType == "LFQ") { lsKl <- c(lsKl, X) }
  ls <- lapply(lsKl, \(kl) { data[[kl]] })
  tmp <- do.call(paste, c(data[, lsKl], sep = "---"))
  if (LabelType == "Isobaric") {
    kol2 <- gsub(topattern(ev.ref["Original"]), "", kol)
    data2 <- as.data.table(data[, kol])
    colnames(data2) <- kol2
    data2$Group <- tmp
    data2$MQ.Exp <- ev$MQ.Exp[w]
    data2 <- data2[, lapply(.SD, sum, na.rm = TRUE), keyby = list(Group = Group,
                                                                  MQ.Exp = MQ.Exp)]
    colnames(data2) <- c("Group", "MQ.Exp", kol2)
    data2 <- as.data.frame(data2)
    data2 <- melt(data2, id.vars = c("Group", "MQ.Exp"))
    data2$value <- log10(data2$value)
    data2 <- data2[which(is.all.good(data2$value, 2L)),]
    data2$IsoBarLab <- Exp.map$`Isobaric label details`[match(as.integer(data2$variable), Exp.map$`Isobaric label`)]
    data2$Parent_sample <- do.call(paste, c(data2[, c("MQ.Exp", "IsoBarLab")], sep = "_"))
    data2 <- data2[which(data2$Parent_sample %in% Exp.map$`Parent sample`),]
    data2 <- as.data.table(data2)
    data2 <- dcast(data2[, c("Group", "value", "Parent_sample")], Group~Parent_sample)
    data2 <- as.data.frame(data2)
    #data2[, lsKl] <- do.call(rbind, c(strsplit(data2$Group, "---")))
    #data2 <- data2[, which(!colnames(data2) %in% kol)]
    #data2[, lsKl] <- data[match(data2$Group, tmp), lsKl]
    data <- data2
    kol2 <- colnames(data)
    kol2 <- kol2[which(kol2 != "Group")]
    impGrps <- rep(1L, length(kol2))
  }
  if (LabelType == "LFQ") {
    kol2 <- if (X == "Parent sample") { get("MQ.Exp") } else { get(X) }
    data2 <- data.table(Intensity = data[[ev.col["Original"]]],
                        Group = tmp)
    data2 <- data2[, .(`log10(Intensity)` = sum(Intensity, na.rm = TRUE)),
                   keyby = .(Group)]
    data2$`log10(Intensity)` <- log10(data2$`log10(Intensity)`)
    data2 <- as.data.frame(data2)
    data2[, lsKl] <- data[match(data2$Group, tmp), lsKl]
    data2$Group <- NULL
    data <- spread(data2, X, "log10(Intensity)")
    impGrps <- rep(1L, length(MQ.Exp))
  }
  data <- data[, kol2]
  w <- which(is.na(data), arr.ind = TRUE)
  if (nrow(w)) {
    temp <- Data_Impute2(data, impGrps, is.log = FALSE)
    data <- temp$Imputed_data
  }
  pcA <- prcomp(t(data[, kol2]), scale. = TRUE)
  if (length(pcA$rotation)) {
    scoresA <- as.data.frame(pcA$x)
    if ("PC2" %in% colnames(scoresA)) {
      scoresA$Sample <- rownames(scoresA)
      rownames(scoresA) <- NULL
      pvA <- round(100*(pcA$sdev)^2 / sum(pcA$sdev^2), 0L)
      pvA <- pvA[which(pvA > 0)]
      pvA2 <- paste0("Original: ", paste(vapply(1L:length(pvA), \(x) {
        paste0("PC", x, ": ", pvA[x], "%")
      }, ""), collapse = ", "))
      scoresA$Label <- scoresA$Sample
      m <- if (LabelType == "Isobaric") {
        match(scoresA$Sample, Exp.map$`Parent sample`)
      } else {
        match(scoresA$Sample, Exp.map$MQ.Exp)
      }
      if (sum(is.na(m))) {
        warning("Mapping samples through MQ.Exp to sample groups failed, check code!\nMapping colors to samples instead of sample groups...")
        scoresA$Colour <- scoresA$Sample
        colKol <- "colKol"
      } else {
        tmp <- Exp.map[m, Factors[which(Factors != "Replicate")]]
        tmp <- tmp[, which(vapply(colnames(tmp), \(x) { length(unique(tmp[[x]])) > 1L }, TRUE)), drop = FALSE]
        scoresA$"Sample group" <- do.call(paste, c(tmp, sep = " "))
        colKol <- "Sample group"
        # tmp <- do.call(cbind, c(list(Exp.map[, "MQ.Exp", drop = FALSE],
        #                              tmp)))
        # scoresA$"Sample group" <- do.call(paste, c(tmp, sep = " "))
      }
      ttl <- "PCA plot - Samples (PSMs-level)"
      xLab <- paste0("PC1 = ", pvA[1L], "%")
      yLab <- paste0("PC2 = ", pvA[2L], "%")
      plot <- ggplot(scoresA, aes(x = PC1, y = PC2, colour = .data[[colKol]])) +
        geom_point() +
        ggpubr::stat_conf_ellipse(aes(fill = .data[[colKol]]),
                                  alpha = 0.1, geom = "polygon", show.legend = FALSE) +
        scale_color_viridis_d(begin = 0.25) +
        coord_fixed() + theme_bw() +
        xlab(xLab) + ylab(yLab) +
        geom_hline(yintercept = 0, colour = "black") + geom_vline(xintercept = 0, colour = "black") +
        ggtitle(ttl#, subtitle = pvA2
                ) +
        geom_text_repel(aes(x = PC1, y = PC2, label = Label, colour = .data[[colKol]]),
                        size = 2.5, show.legend = FALSE)
      #poplot(plot)
      suppressMessages({
        ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 300L, width = 20L, height = 20L, units = "in")
        ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300L, width = 20L, height = 20L, units = "in")
      })
      ReportCalls <- AddPlot2Report(Space = FALSE)
      Symb <- "circle"
      # Custom color scale
      scoresA$`Samples group` <- factor(scoresA[[colKol]])
      if ("PC3" %in% colnames(scoresA)) {
        plot_lyPSMsPCA <- plot_ly(scoresA, x = ~PC1, y = ~PC2, z = ~PC3,
                                  text = ~Label, type = "scatter3d", mode = "markers",
                                  color = ~`Samples group`, colors = "viridis",
                                  symbol = I(Symb))
      } else {
        plot_lyPSMsPCA <- plot_ly(scoresA, x = ~PC1, y = ~PC2,
                                  text = ~Label, type = "scatter", mode = "markers",
                                  color = ~`Samples group`, colors = "viridis",
                                  symbol = I(Symb))
      }
      plot_lyPSMsPCA %<o% layout(plot_lyPSMsPCA, title = ttl)
      setwd(dir)
      saveWidget(plot_lyPSMsPCA, paste0(dir, "/", ttl, ".html"), selfcontained = TRUE)
      #system(paste0("open \"", dir, "/", ttl, ".html"))
      setwd(wd)
    } else {
      msg <- "Not enough valid data to draw a PSM-level PCA plot!"
      ReportCalls <- AddMsg2Report(Offset = TRUE, Space = FALSE, Warning = TRUE)
    }
  }
  ReportCalls$Calls[[LRepCalls]] <- append(ReportCalls$Calls[[LRepCalls]],
                                           "body_add_par(Report, \"\", style = \"Normal\")")
} else {
  stop("Uh, I think you have the wrong analysis pipeline here...\nwhere are my sample groups and replicates!?!?!")
}
