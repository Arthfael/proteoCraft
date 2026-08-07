#### Peptidoforms-level, calculate quantitative values + draw PCA
#
# Calculate single channel intensities and total intensity
# For now only covers case with replicates (not inherently but because of some of the variables it uses),
# but could be redone to also cover the more general case.
#
# To do:
# - Add Pearson correlation heatmap amongst those visualizations used to decide whether to remove any outliers, it is very good!
#
require(parallel)
pep.ref %<o% setNames("Int. - ", "Original")
if (exists("scrptType")) {
  if (scrptType == "noReps") {
    stop("Write that bit!")
  }
  if (!"Use" %in% colnames(Exp.map)) { Exp.map$Use <- TRUE } else {
    if (is.character(Exp.map$Use)) {
      Exp.map$Use[which(Exp.map$Use == "T")] <- "TRUE"
      Exp.map$Use[which(Exp.map$Use == "F")] <- "FALSE"
      Exp.map$Use <- as.logical(Exp.map$Use)
      Exp.map$Use[which(is.na(Exp.map$Use))] <- TRUE
    }
    tmp_EM <- Exp.map
    refCol <- "Ref.Sample.Aggregate"
  }
} else {
  # You will need:
  #  - VPAL and RSA with at least the values and names fields
  #  - a FracMap data.frame with at least:
  #      - a "Samples" column, matching the values in ev$MQ.Exp,
  #      - a "Name" column, matching the values in ev$Experiment
  #      - optionally a "Use" logical column
  if (exists("LabelType")) {
    if (LabelType != "LFQ") { stop("This case isn't covered yet!") }
  } else {
    LabelType <- "LFQ"
  }
  stopifnot(exists("FracMap"),
            "Samples" %in%  colnames(FracMap))
  kol <- intersect(c("Group", "Sample", "Replicate"), colnames(FracMap))
  stopifnot(length(kol) > 0L)
  tmp_EM <- FracMap[, kol, drop = FALSE]
  if ((!"Use" %in% colnames(tmp_EM)) || is.logical(tmp_EM$Use) || sum(is.na(tmp_EM$Use)) || (!sum(tmp_EM$Use))) {
    tmp_EM$Use <- TRUE    
  }
  tmp_EM$MQ.Exp <- FracMap$Samples # (Not "Sample"!)
  refCol <- "Name"
  ev.col <- "Intensity"
}
source(parSrc, local = FALSE)
exports <- list("smpls", "tmp_EM", "pep.ref", "LabelType", "wd", "refCol")
if (LabelType == "Isobaric") {
  tmp <- ev[, c("MQ.Exp", "Modified sequence",
                paste0(ev.ref[length(ev.ref)], as.character(sort(as.numeric(unique(tmp_EM$"Isobaric label"))))))]
  exports <- append(exports, "ev.ref")
}
if (LabelType == "LFQ") {
  tmp <- ev[, c("MQ.Exp", "Modified sequence", ev.col[length(ev.col)])]
  exports <- append(exports, "ev.col")
}
readr::write_rds(tmp, paste0(wd, "/tmp.RDS"))
smpls <- unique(tmp_EM[which(tmp_EM$Use), refCol])
clusterExport(parClust, exports, envir = environment())
invisible(clusterCall(parClust, \(x) {
  library(data.table)
  tmp <<- readr::read_rds(paste0(wd, "/tmp.RDS"))
  return()
}))
unlink(paste0(wd, "/tmp.RDS"))
tmp4 <- setNames(parLapply(parClust, smpls, \(smpl) { #smpl <- smpls[1L]
  m <- match(smpl, tmp_EM[[refCol]])
  mqe <- unlist(tmp_EM$MQ.Exp[m])
  w2 <- which(tmp$MQ.Exp %in% mqe)
  tmp2 <- data.frame(mod = NA, Intensity = NA_real_)
  if (length(w2)) {
    if (LabelType == "Isobaric") {
      j <- as.character(sort(as.numeric(tmp_EM$"Isobaric label"[m])))
      tmp3 <- tmp[w2, paste0(ev.ref[length(ev.ref)], j), drop = FALSE]
      for (k in j) {
        kk <- paste0(ev.ref[length(ev.ref)], j)
        tmp3[which(!is.finite(tmp3[[kk]])), kk] <- NA_real_
      }
      if (length(j) > 1L) { tmp3 <- apply(tmp3, 1L, sum, na.rm = TRUE) } # Ultra-rare cases where the same parent sample is in different isobaric channels in different fractions
      tmp2 <- data.table(mod = tmp$"Modified sequence"[w2],
                         Intensity = unlist(tmp3))
    }
    if (LabelType == "LFQ") {
      tmp2 <- data.table(mod = tmp$"Modified sequence"[w2],
                         Intensity = tmp[w2, ev.col[length(ev.col)]])
      tmp2$Intensity[which(!is.finite(tmp2$Intensity))] <- NA_real_
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
pc1 <- stats::prcomp(t(data[, kol]), scale. = TRUE)
dir <- paste0(wd, "/Workflow control/Peptides/PCA plot")
if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
if (exists("dirlist")) { dirlist <- unique(c(dirlist, dir)) } 
if (length(pc1$rotation)) {
  scores1 <- as.data.frame(pc1$x)
  if ("PC2" %in% colnames(scores1)) {
    rownames(scores1) <- gsub(topattern(pep.ref["Original"]), "", rownames(scores1))
    scores1[, RSA$names] <- Isapply(strsplit(rownames(scores1), "___"), unlist)
    scores1$Use <- tmp_EM$Use[match(rownames(scores1), tmp_EM[[refCol]])]
    rownames(scores1) <- NULL
    pv1 <- round(100*(pc1$sdev)^2L / sum(pc1$sdev^2L), 0L)
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
    require(ggplot2)
    require(ggrepel)
    require(plotly)
    require(htmlwidgets)
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
    pcaDir <- paste0(wd, "/Workflow control/Peptides/PCA plot")
    if (!dir.exists(pcaDir)) { dir.create(pcaDir, recursive = TRUE) }
    setwd(pcaDir)
    saveWidget(plot_lyPCA, paste0(wd, "/Workflow control/Peptides/PCA plot/", ttl, ".html"),
               selfcontained = TRUE)
    setwd(wd)
    #system(paste0("open \"", wd, "/Workflow control/Peptides/PCA plot/", ttl, ".html"))
  } else {
    stop("There was only one component to the PCA, something must've gone wrong when generating the peptides table!") #(I think this will never happen, the previous check should be identical...?)
  }
} else { stop("There was only one component to the PCA, something must've gone wrong when generating the peptides table!") }
