#### Peptidoforms-level, calculate quantitative values + draw PCA
## (future option, or when executing line-by-line: remove outliers/samples which will not be used)
## First visualize data: are there any clear outliers?
# Calculate single channel intensities and total intensity
#
# To do:
# - If outlier is the only reference sample in a reference group, unfortunately you will have to remove the whole group
# - Add Pearson correlation heatmap amongst those visualizations used to decide whether to remove any outliers, it is very good!
#
pep.ref %<o% setNames("Int. - ", "Original")
if (!"Use" %in% colnames(Exp.map)) { Exp.map$Use <- TRUE } else {
  if (is.character(Exp.map$Use)) {
    Exp.map$Use[which(Exp.map$Use == "T")] <- "TRUE"
    Exp.map$Use[which(Exp.map$Use == "F")] <- "FALSE"
    Exp.map$Use <- as.logical(Exp.map$Use)
    Exp.map$Use[which(is.na(Exp.map$Use))] <- TRUE
  }
}
source(parSrc, local = FALSE)
exports <- list("smpls", "Exp.map", "pep.ref", "LabelType", "wd", "is.all.good")
if (LabelType == "Isobaric") {
  tmp <- ev[, c("MQ.Exp", "Modified sequence",
                paste0(ev.ref[length(ev.ref)], as.character(sort(as.numeric(unique(Exp.map$"Isobaric label"))))))]
  exports <- append(exports, "ev.ref")
}
if (LabelType == "LFQ") {
  tmp <- ev[, c("MQ.Exp", "Modified sequence", ev.col[length(ev.col)])]
  exports <- append(exports, "ev.col")
}
readr::write_rds(tmp, paste0(wd, "/tmp.RDS"))
smpls <- unique(Exp.map$Ref.Sample.Aggregate[which(Exp.map$Use)])
clusterExport(parClust, exports, envir = environment())
invisible(clusterCall(parClust, \(x) {
  library(data.table)
  tmp <<- readr::read_rds(paste0(wd, "/tmp.RDS"))
  return()
}))
unlink(paste0(wd, "/tmp.RDS"))
tmp4 <- setNames(parLapply(parClust, smpls, \(smpl) { #smpl <- smpls[1L]
  m <- match(smpl, Exp.map$Ref.Sample.Aggregate)
  mqe <- unlist(Exp.map$MQ.Exp[m])
  w2 <- which(tmp$MQ.Exp %in% mqe)
  tmp2 <- data.frame(mod = NA, Intensity = NA)
  if (length(w2)) {
    if (LabelType == "Isobaric") {
      j <- as.character(sort(as.numeric(Exp.map$"Isobaric label"[m])))
      tmp3 <- tmp[w2, paste0(ev.ref[length(ev.ref)], j), drop = FALSE]
      for (k in j) {
        kk <- paste0(ev.ref[length(ev.ref)], j)
        tmp3[which(!is.all.good(tmp3[[kk]], 2L)), kk] <- NA
      }
      if (length(j) > 1) { tmp3 <- apply(tmp3, 1L, sum, na.rm = TRUE) } # Ultra-rare cases where the same parent sample is in different isobaric channels in different fractions
      tmp2 <- data.table(mod = tmp$"Modified sequence"[w2],
                         Intensity = unlist(tmp3))
    }
    if (LabelType == "LFQ") {
      tmp2 <- data.table(mod = tmp$"Modified sequence"[w2],
                         Intensity = tmp[w2, ev.col[length(ev.col)]])
      tmp2$Intensity[which(!is.all.good(tmp2$Intensity, 2L))] <- NA
    }
    tmp2 <- tmp2[, list(Intensity = sum(Intensity, na.rm = TRUE)), by = list(mod)]
    tmp2 <- as.data.frame(tmp2)
  }
  return(tmp2)
}), smpls)
for (smpl in smpls) { #smpl <- smpls[1L]
  tmp <- tmp4[[smpl]]
  pep[[paste0(pep.ref["Original"], smpl)]] <- 0
  w3 <- which(pep$"Modified sequence" %in% tmp$mod)
  pep[w3, paste0(pep.ref["Original"], smpl)] <- tmp$Intensity[match(pep$"Modified sequence"[w3], tmp$mod)]
}
kol <- paste0(pep.ref["Original"], RSA$values)
kol <- kol[which(kol %in% colnames(pep))]
data <- pep[, c("Modified sequence", kol)]
w <- which(rowSums(data[, kol], na.rm = TRUE) > 0)
data <- data[w,]
pc1 <- prcomp(t(data[, kol]), scale. = TRUE)
dir <- paste0(wd, "/Workflow control/Peptides/PCA plot")
if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
dirlist <- unique(c(dirlist, dir))
if (length(pc1$rotation)) {
  scores1 <- as.data.frame(pc1$x)
  if ("PC2" %in% colnames(scores1)) {
    rownames(scores1) <- gsub(topattern(pep.ref["Original"]), "", rownames(scores1))
    scores1[, RSA$names] <- Isapply(strsplit(rownames(scores1), "___"), unlist)
    scores1$Use <- Exp.map$Use[match(rownames(scores1), Exp.map$Ref.Sample.Aggregate)]
    rownames(scores1) <- NULL
    pv1 <- round(100*(pc1$sdev)^2 / sum(pc1$sdev^2), 0L)
    pv1 <- pv1[which(pv1 > 0)]
    pv1_ <- paste0("Original: ", paste(vapply(seq_along(pv1), \(x) {
      paste0("PC", x, ": ", pv1[x], "%")
    }, ""), collapse = ", "))
    w <- which(vapply(VPAL$names, \(x) { length(unique(scores1[[x]])) }, 1L) > 1L)
    w <- w[which(tolower(substr(names(w), 1L, 3L)) != "rep")]
    scores1$Samples_group <- do.call(paste, c(scores1[, VPAL$names[w], drop = FALSE], sep = " "))
    scores1$Label <- do.call(paste, c(scores1[, RSA$names, drop = FALSE], sep = " "))
    outlierAnnot_shape %<o% "Replicate"
    outlierAnnot_color %<o% "Samples_group"
    ttl <- "PCA plot - Preliminary - peptide level"
    xLab <- paste0("PC1 = ", pv1[1L], "%")
    yLab <- paste0("PC2 = ", pv1[2L], "%")
    plot <- ggplot(scores1, aes(x = PC1, y = PC2, colour = .data[[outlierAnnot_color]])) +
      geom_point(aes(shape = .data[[outlierAnnot_shape]])) +
      ggpubr::stat_conf_ellipse(aes(fill = .data[[outlierAnnot_color]]),
                                alpha = 0.2, geom = "polygon", show.legend = FALSE) +
      scale_color_viridis_d(begin = 0.25) +
      coord_fixed() + theme_bw() +
      xlab(xLab) + ylab(yLab) +
      geom_hline(yintercept = 0, colour = "black") + geom_vline(xintercept = 0, colour = "black") +
      ggtitle(ttl#, subtitle = pv1_
      ) +
      geom_text_repel(aes(label = Label), size = 2.5, show.legend = FALSE)
    #poplot(plot)
    suppressMessages({
      ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 300L, width = 10L, height = 10L, units = "in")
      ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300L, width = 10L, height = 10L, units = "in")
    })
    ReportCalls <- AddPlot2Report()
    nReps <- max(as.numeric(Rep))
    Symb <- rep(c("circle", "diamond", "square", "cross", "x"), nReps)[seq_len(nReps)]             
    Symb <- Symb[as.numeric(scores1[[outlierAnnot_shape]])]
    # Custom color scale
    scores1$"Samples group" <- factor(scores1$Samples_group)
    plot_lyPCA <- if ("PC3" %in% colnames(scores1)) {
      plot_ly(scores1, x = ~PC1, y = ~PC2, z = ~PC3,
              text = ~Label, type = "scatter3d", mode = "markers",
              color = ~get(outlierAnnot_color), colors = "viridis",
              symbol = I(Symb))
    } else {
      plot_ly(scores1, x = ~PC1, y = ~PC2,
              text = ~Label, type = "scatter", mode = "markers",
              color = ~`Samples group`, colors = "viridis",
              symbol = I(Symb))
    }
    plot_lyPCA %<o% layout(plot_lyPCA, title = ttl)
    renderPlotly({ plot_lyPCA <- plot_lyPCA })
    saveWidget(plot_lyPCA, paste0(wd, "/Workflow control/Peptides/PCA plot/", ttl, ".html"),
               selfcontained = TRUE)
    #system(paste0("open \"", wd, "/Workflow control/Peptides/PCA plot/", ttl, ".html"))
  } else {
    stop("There was only one component to the PCA, something must've gone wrong when generating the peptides table!") #(I think this will never happen, the previous check should be identical...?)
  }
} else { stop("There was only one component to the PCA, something must've gone wrong when generating the peptides table!") }
