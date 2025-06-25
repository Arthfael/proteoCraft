### PCA plot for parameters app
# Create first PCA to check on sample relationships
if ((length(MQ.Exp) > 1)||(LabelType == "Isobaric")) { # Should be always TRUE
  source(parSrc, local = FALSE)
  data <- ev
  colnames(data)[which(colnames(data) == "MQ.Exp")] <- "Parent sample"
  data <- data[which(data$Reverse != "+"),]
  data <- data[which((is.na(data$"Potential contaminant"))|(data$"Potential contaminant" != "+")),]
  if (LabelType == "Isobaric") {
    kol <- grep(paste0(topattern(ev.ref["Original"]), "[0-9]+$"), colnames(data), value = TRUE)
  } else { kol <- ev.col["Original"] }
  w <- which(rowSums(data[, kol, drop = FALSE], na.rm = TRUE) > 0)
  data <- data[w,]
  if (!"Fraction" %in% colnames(data)) { data$Fraction <- 1 }
  Fraction <- sort(unique(data$Fraction), decreasing = FALSE)
  Experiment <- Exp
  kols <- c("Parent sample", "Fraction", "Experiment")
  if (LabelType == "Isobaric") {
    X <- "Label"
    kols <- c("Fraction", "Parent sample", "Experiment") # The order matters!
    tst <- vapply(kols, function(x) { length(unique(data[[x]])) }, 1)
    w1 <- which(tst > 1)
    w2 <- which(tst >= 1) 
    if (length(w1)) { Y <- kols[w1[1]] } else { Y <- kols[w2[1]] }
  }
  if (LabelType == "LFQ") {
    kols <- c("Parent sample", "Fraction", "Experiment") # The order matters!
    tst <- vapply(kols, function(x) { length(unique(data[[x]])) }, 1)
    w1 <- which(tst > 1)
    w2 <- which(tst >= 1) 
    X <- kols[w1[1]]
    if (length(w1) > 1) { Y <- kols[w1[2]] } else { Y <- kols[w2[2]] }
  }
  kols <- kols[which(!kols %in% c(X, Y))]
  ReportCalls <- AddSpace2Report()
  ReportCalls$Calls <- append(ReportCalls$Calls,
                              "body_add_fpar(Report, fpar(ftext(\"PSMs-level PCA plot:\", prop = WrdFrmt$Section_title), fp_p = WrdFrmt$just))")
  ReportCalls$Calls <- append(ReportCalls$Calls, list())
  dir <- paste0(wd, "/Dimensionality red. plots/PCA")
  if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
  dirlist <- unique(c(dirlist, dir))
  LRepCalls <- length(ReportCalls$Calls)
  lsKl <- c("Modified sequence", Y)
  if (LabelType == "LFQ") { lsKl <- c(lsKl, X) }
  ls <- lapply(lsKl, function(kl) { data[[kl]] })
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
    data2 <- data2[which(is.all.good(data2$value, 2)),]
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
    impGrps <- rep(1, length(kol2))
  }
  if (LabelType == "LFQ") {
    if (X == "Parent sample") { kol2 <- get("MQ.Exp") } else { kol2 <- get(X) }
    data2 <- data.table(Intensity = data[[ev.col["Original"]]],
                        Group = tmp)
    data2 <- data2[, list(`log10(Intensity)` = sum(Intensity, na.rm = TRUE)),
                   keyby = Group]
    data2$`log10(Intensity)` <- log10(data2$`log10(Intensity)`)
    data2 <- as.data.frame(data2)
    data2[, lsKl] <- data[match(data2$Group, tmp), lsKl]
    data2$Group <- NULL
    data <- spread(data2, X, "log10(Intensity)")
    impGrps <- rep(1, length(MQ.Exp))
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
      pvA <- round(100*(pcA$sdev)^2 / sum(pcA$sdev^2), 0)
      pvA <- pvA[which(pvA > 0)]
      pvA <- paste0("Original: ", paste(vapply(1:length(pvA), function(x) {
        paste0("PC", x, ": ", pvA[x], "%")
      }, ""), collapse = ", "))
      scoresA$Label <- scoresA$Sample
      if (LabelType == "Isobaric") {
        m <- match(scoresA$Sample, Exp.map$`Parent sample`)
      } else {
        m <- match(scoresA$Sample, Exp.map$MQ.Exp)
      }
      if (sum(is.na(m))) {
        warning("Mapping samples through MQ.Exp to sample groups failed, check code!\nMapping colors to samples instead of sample groups...")
        scoresA$Colour <- scoresA$Sample
      } else {
        tmp <- Exp.map[m, Factors[which(Factors != "Replicate")]]
        tmp <- tmp[, which(vapply(colnames(tmp), function(x) { length(unique(tmp[[x]])) > 1 }, TRUE)), drop = FALSE]
        tmp <- do.call(cbind, c(list(Exp.map[, "MQ.Exp", drop = FALSE],
                                     tmp)))
        scoresA$Colour <- do.call(paste, c(tmp, sep = " "))
      }
      ttl <- "PCA plot - Samples (PSMs-level)"
      plot <- ggplot(scoresA) +
        geom_point(aes(x = PC1, y = PC2, colour = Colour)) +
        scale_color_viridis_d(begin = 0.25) +
        coord_fixed() + theme_bw() +
        geom_hline(yintercept = 0, colour = "black") + geom_vline(xintercept = 0, colour = "black") +
        ggtitle(ttl, subtitle = pvA) +
        geom_text_repel(aes(x = PC1, y = PC2, label = Label, colour = Colour),
                        size = 2.5, show.legend = FALSE)
      #poplot(plot)
      ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 300, width = 20, height = 20, units = "in")
      ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300, width = 20, height = 20, units = "in")
      ReportCalls <- AddPlot2Report(Space = FALSE)
      Symb <- "circle"
      # Custom color scale
      scoresA$`Samples group` <- factor(scoresA$Colour)
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
      saveWidget(plot_lyPSMsPCA, paste0(dir, "/", ttl, ".html"), selfcontained = TRUE)
      #system(paste0("open \"", dir, "/", ttl, ".html"))
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
