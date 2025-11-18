#### Peptidoforms-level, calculate quantitative values + draw PCA
## (future option, or when executing line-by-line: remove outliers/samples which will not be used)
## First visualize data: are there any clear outliers?
# Calculate single channel intensities and total intensity
#
# To do:
# - If outlier is the only reference sample in a reference group, unfortunately you will have to remove the whole group
# - Add Pearson correlation heatmap amongst those visualizations used to decide whether to remove any outliers, it is very good!
#
pep.ref %<o% setNames("Evidence intensities - ", "Original")
if (!"Use" %in% colnames(Exp.map)) { Exp.map$Use <- TRUE } else {
  if (class(Exp.map$Use) == "character") {
    Exp.map$Use[which(Exp.map$Use == "T")] <- "TRUE"
    Exp.map$Use[which(Exp.map$Use == "F")] <- "FALSE"
    Exp.map$Use <- as.logical(Exp.map$Use)
    Exp.map$Use[which(is.na(Exp.map$Use))] <- TRUE
  }
}
source(parSrc, local = FALSE)
exports <- list("smpls", "Exp.map", "pep.ref", "LabelType", "wd")
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
invisible(clusterCall(parClust, function(x) {
  library(data.table)
  tmp <<- readr::read_rds(paste0(wd, "/tmp.RDS"))
  return()
}))
unlink(paste0(wd, "/tmp.RDS"))
tmp4 <- setNames(parLapply(parClust, smpls, function(smpl) { #smpl <- smpls[1]
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
        tmp3[which(!proteoCraft::is.all.good(tmp3[[kk]], 2)), kk] <- NA
      }
      if (length(j) > 1) { tmp3 <- apply(tmp3, 1, sum, na.rm = TRUE) } # Ultra-rare cases where the same parent sample is in different isobaric channels in different fractions
      tmp2 <- data.table(mod = tmp$"Modified sequence"[w2],
                         Intensity = unlist(tmp3))
    }
    if (LabelType == "LFQ") {
      tmp2 <- data.table(mod = tmp$"Modified sequence"[w2],
                         Intensity = tmp[w2, ev.col[length(ev.col)]])
      tmp2$Intensity[which(!proteoCraft::is.all.good(tmp2$Intensity, 2))] <- NA
    }
    tmp2 <- tmp2[, list(Intensity = sum(Intensity, na.rm = TRUE)), by = list(mod)]
    tmp2 <- as.data.frame(tmp2)
  }
  return(tmp2)
}), smpls)
for (smpl in smpls) { #smpl <- smpls[1]
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
    pv1 <- round(100*(pc1$sdev)^2 / sum(pc1$sdev^2), 0)
    pv1 <- pv1[which(pv1 > 0)]
    pv1_ <- paste0("Original: ", paste(vapply(seq_along(pv1), function(x) {
      paste0("PC", x, ": ", pv1[x], "%")
    }, ""), collapse = ", "))
    w <- which(vapply(VPAL$names, function(x) { length(unique(scores1[[x]])) }, 1) > 1)
    w <- w[which(tolower(substr(names(w), 1, 3)) != "rep")]
    scores1$Samples_group <- do.call(paste, c(scores1[, VPAL$names[w], drop = FALSE], sep = " "))
    scores1$Label <- do.call(paste, c(scores1[, RSA$names, drop = FALSE], sep = " "))
    outlierAnnot_shape %<o% "Replicate"
    outlierAnnot_color %<o% "Samples_group"
    ttl <- "PCA plot - Preliminary - peptide level"
    xLab <- paste0("PC1 = ", pv1[1], "%")
    yLab <- paste0("PC2 = ", pv1[2], "%")
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
      ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, dpi = 300, width = 10, height = 10, units = "in")
      ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300, width = 10, height = 10, units = "in")
    })
    ReportCalls <- AddPlot2Report()
    nReps <- max(as.numeric(Rep))
    Symb <- rep(c("circle", "diamond", "square", "cross", "x"), nReps)[seq_len(nReps)]             
    Symb <- Symb[as.numeric(scores1[[outlierAnnot_shape]])]
    # Custom color scale
    scores1$"Samples group" <- factor(scores1$Samples_group)
    if ("PC3" %in% colnames(scores1)) {
      plot_lyPCA <- plot_ly(scores1, x = ~PC1, y = ~PC2, z = ~PC3,
                            text = ~Label, type = "scatter3d", mode = "markers",
                            color = ~get(outlierAnnot_color), colors = "viridis",
                            symbol = I(Symb))
    } else {
      plot_lyPCA <- plot_ly(scores1, x = ~PC1, y = ~PC2,
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

# Peptides heatmap function - useful later and for checking normalisations
pepHtmp %<o% function(intProt = prot.list_pep,
                      Pep = pep,
                      ref = pep.ref[length(pep.ref)],
                      dstDir,
                      ttlRoot = "Peptides log2 heatmap",
                      DB = db,
                      Experiment = Exp,
                      RSA = Ref.Sample.Aggregate,
                      is.log = FALSE,
                      cl,
                      N.clust,
                      N.reserved = 1) {
  TESTING <- FALSE
  #TESTING <- TRUE
  if (TESTING) {
    # Note:
    # This is not a perfect alternative to missing but will work in most cases, unless x matches a function imported by a package 
    misFun <- function(x) { return(!exists(deparse(substitute(x)))) }
  } else { misFun <- missing }
  stopifnot(!is.null(intProt), length(intProt) > 0,
            !is.null(Pep), length(Pep) > 0)
  if ((is.null(ttlRoot))||(!"character" %in% class(ttlRoot))||(!nchar(ttlRoot))) {
    warning("Argument \"ttlRoot\" cannot be NULL!")
    ttlRoot <- "Peptides log2 heatmap"
  }
  if (!dir.exists(dstDir)) { dir.create(dstDir, recursive = TRUE) }
  if (exists("dirlist")) { assign("dirlist", unique(c(dirlist, dstDir)), envir = .GlobalEnv) }
  #
  g <- paste0(ref, RSA$values)
  g <- g[which(g %in% colnames(Pep))]
  g1 <- as.data.frame(t(as.data.frame(strsplit(gsub(topattern(ref), "", g), "___"))))
  colnames(g1) <- RSA$names
  g1$Col <- g
  test <- RSA$names[which(RSA$names != "Replicate")]
  for (i in rev(test)) { g1 <- g1[order(g1[[i]]),] }
  g <- g1$Col
  g1 <- gsub(topattern(ref), "", g)
  g1 <- cleanNms(g1, Experiment = Experiment)
  DB <- DB[match(intProt, DB$"Protein ID"), c("Common Name", "Protein ID", "Sequence")]
  tmpPep <- Pep[grsep2(intProt, Pep$Proteins), c("Proteins", "Sequence", "Modified sequence", g)]
  cleanUp <- FALSE
  if (nrow(tmpPep)) {
    # Create cluster (some steps are slow otherwise)
    if (misFun(cl)) {
      require(parallel)
      dc <- detectCores()
      if (misFun(N.reserved)) { N.reserved <- 1 }
      if (misFun(N.clust)) {
        N.clust <- max(c(dc-N.reserved, 1))
      } else {
        if (N.clust > max(c(dc-N.reserved, 1))) {
          warning("More cores specified than allowed, I will ignore the specified number! You should always leave at least one free for other processes, see the \"N.reserved\" argument.")
          N.clust <- max(c(dc-N.reserved, 1))
        }
      }
      cl <- makeCluster(N.clust, type = "SOCK")
      cleanUp <- TRUE
    }
    clusterExport(cl, list("DB", "tmpPep", "dir", "g", "g1", "ttlRoot", "is.log", "dstDir"), envir = environment())
    l <- length(intProt)
    tempDat <- data.frame(Protein = rep(intProt, 2),
                          Method = c(rep("Mean", l), rep("ZSc", l)))
    f0 <- function(x) { #x <- tempDat[1,]
      plp <- unlist(x[[1]])
      meth <- x[[2]]
      Plp <- paste(DB[which(DB$"Protein ID" == plp),
                      c("Common Name", "Protein ID")], collapse = " - ")
      grs <- proteoCraft::grsep2(plp, tmpPep$Proteins)
      if (length(grs)) {
        Seq <- DB$Sequence[match(plp, DB$"Protein ID")]
        temp <- tmpPep[grs, c("Sequence", "Modified sequence", g)]
        temp$tst1 <- vapply(temp$Sequence, function(x) { nchar(unlist(strsplit(Seq, x))[1]) }, 1)
        temp$tst2 <- vapply(temp$Sequence, nchar, 1)
        temp$tst3 <- vapply(temp$"Modified sequence", nchar, 1)
        temp <- temp[order(temp$tst1, temp$tst2, temp$tst3),]
        w <- which(rowSums(temp[, g], na.rm = TRUE) == 0)
        if (length(w)) {
          warning("Peptide(s) found with all invalid intensity values!")
          w <- which(rowSums(temp[, g], na.rm = TRUE) > 0)
          temp <- temp[w,]
        }
        if (nrow(temp)) {
          colnames(temp)[match(g, colnames(temp))] <- g1
          # Create heatmap
          temp <- temp[, c("Modified sequence", g1)]
          if (!is.log) {
            warning(" Converting data to log2...")
            temp[, g1] <- suppressWarnings(log2(temp[, g1]))  
          }
          w <- which(!is.finite(as.matrix(temp[, g1])), arr.ind = TRUE)
          temp[, g1][w] <- NA
          M <- rowMeans(temp[, g1], na.rm = TRUE)
          if (meth == "ZSc") { SD <- apply(temp[, g1], 1, sd, na.rm = TRUE) }
          temp[, g1] <- sweep(temp[, g1], 1, M, "-")
          if (meth == "ZSc") { temp[, g1] <- sweep(temp[, g1], 1, SD, "/") }
          temp2 <- proteoCraft::dfMelt(temp[, g1], c("Sample", "value"))
          temp2$"Modified sequence" <- temp$"Modified sequence"
          temp2$Sample <- as.character(temp2$Sample)
          temp2$Xmin <- match(temp2$Sample, colnames(temp))-1
          temp2$Xmax <- temp2$Xmin+1
          temp2$Ymax <- nrow(temp):1
          temp2$Ymin <- temp2$Ymax-1
          hlab <- aggregate(temp2[, c("Xmin", "Xmax")], list(temp2$Sample), unique)
          colnames(hlab)[1] <- "Sample"
          hlab$X <- rowMeans(hlab[, c("Xmin", "Xmax")])
          vlab <- aggregate(temp2[, c("Ymin", "Ymax")], list(temp2$"Modified sequence"), unique)
          colnames(vlab)[1] <- "Modified sequence"
          vlab$Y <- rowMeans(vlab[, c("Ymin", "Ymax")])
          Xscale <- max(temp2$Xmax)
          Yscale <- max(temp2$Ymax)
          # Create heatmap plot
          Xlim <- c(-1, Xscale+15)
          Ylim <- c(-6, Yscale+20)
          temp2a <- temp2[, c("Xmin", "Ymin", "value")]
          Splits <- 20
          XScale2 <- Xscale*0.1/Splits
          temp2b <- data.frame(Xmin = Xscale/2 - XScale2*((-Splits/2):(Splits/2)),
                               Ymin = -5)
          temp2b$Xmax <- temp2b$Xmin + XScale2
          temp2b$value <- min(temp2a$value, na.rm = TRUE) + (Splits:0)*(max(temp2a$value, na.rm = TRUE)-min(temp2a$value, na.rm = TRUE))/Splits
          temp2b$Label <- ""
          temp2b$Label[c(1, Splits/2+1, Splits+1)] <- round(temp2b$value[c(1, Splits/2+1, Splits+1)], 1)
          # Create graph
          ttl <- paste0(ttlRoot, " - ", gsub("/", "-", Plp), c("", " (Z-scored)")[(meth == "ZSc")+1])
          heatmap.plot <- ggplot2::ggplot(temp2a) +
            ggplot2::geom_rect(ggplot2::aes(xmin = Xmin, xmax = Xmin+1, ymin = Ymin, ymax = Ymin+1, fill = value)) +
            ggplot2::geom_rect(data = temp2b, ggplot2::aes(xmin = Xmin, xmax = Xmax, ymin = Ymin, ymax = Ymin+1, fill = value)) +
            ggplot2::geom_text(data = hlab, ggplot2::aes(x = X, label = Sample), y = Yscale+1, angle = 60, hjust = 0, size = 2.5) +
            ggplot2::geom_text(data = vlab, ggplot2::aes(y = Y, label = `Modified sequence`), x = Xscale+0.5, hjust = 0, size = 1.7) +
            ggplot2::geom_text(data = temp2b, ggplot2::aes(x = Xmin, label = Label), y = -3.5, hjust = 0.5, size = 2) +
            ggplot2::ggtitle(paste0(ttlRoot, ", ", c("mean-normalized", "Z-scored")[(meth == "ZSc")+1]),
                             subtitle = Plp) +
            ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(), 
                           axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(), 
                           panel.background = ggplot2::element_rect(fill = "transparent", color = NA),
                           plot.margin = ggplot2::margin(0, 0, 0, 0, "cm")) +
            ggplot2::coord_fixed(0.5) +
            ggplot2::scale_fill_gradient2(low = "red", mid = "black", high = "green", na.value = "lightblue") +
            ggplot2::xlab(NULL) + ggplot2::ylab(NULL) + ggplot2::theme(legend.position = "none") +
            ggplot2::xlim(Xlim[1], Xlim[2]) + ggplot2::ylim(Ylim[1], Ylim[2])
          #proteoCraft::poplot(heatmap.plot)
          suppressMessages({
            ggplot2::ggsave(paste0(dstDir, "/", ttl, ".jpeg"), heatmap.plot,
                            dpi = 600, width = 20, height = 12, units = "in")
            ggplot2::ggsave(paste0(dstDir, "/", ttl, ".pdf"), heatmap.plot,
                            dpi = 600, width = 20, height = 12, units = "in")
          })
          #system(paste0("open \"", dstDir, "/", ttl, ".jpeg", "\""))
          #system(paste0("open \"", dstDir, "/", ttl, ".pdf", "\""))
        }
      }
      return()
    }
    environment(f0) <- .GlobalEnv
    invisible(parApply(cl, tempDat, 1, f0))
  }
  if (cleanUp) { stopCluster(cl) }
}
# Another useful function for checking peptide normalisations
pepPlotFun %<o% function(df1,
                         df2,
                         ttl,
                         dstDir,
                         save = TRUE,
                         xpMap = Exp.map,
                         VPAL = Volcano.plots.Aggregate.Level) {
  tst1 <- df1
  tst2 <- df2
  colnames(tst1) <- gsub(".* - ", "", colnames(tst1))
  colnames(tst2) <- gsub(".* - ", "", colnames(tst2))
  tst1 <- proteoCraft::dfMelt(tst1)
  tst2 <- proteoCraft::dfMelt(tst2)
  tst1$Norm <- "Original"
  tst2$Norm <- "Re-normalized"
  tst1 <- rbind(tst1, tst2)
  rm(tst2)
  tst1 <- tst1[which(is.finite(tst1$value)),]
  g <- grepl("_REF\\.to\\.REF_", tst1$variable)
  tst1$Type <- c("Samples", "References")[g+1]
  tst1$Group <- xpMap[match(tst1$variable, xpMap$Ref.Sample.Aggregate),
                      VPAL$column]
  tst1$variable <- cleanNms(tst1$variable)
  w <- which(g)
  tst1$variable[w] <- gsub_Rep("_REF\\.to\\.REF_.*", "", tst1$variable[w])
  tst1$Group <- cleanNms(tst1$Group)
  tst1$Group[which(is.na(tst1$Group))] <- "References"
  tst1$Norm <- as.factor(tst1$Norm)
  tst1$Type <- as.factor(tst1$Type)
  plot <- ggplot(tst1) +
    geom_density(stat = "density", alpha = 0.1, aes(x = value, colour = variable, fill = variable)) +
    scale_y_continuous(expand = c(0, 0)) +
    facet_grid(Norm~Group) + ggtitle(ttl) + theme_bw()
  if (grepl("ratio", ttl, ignore.case = TRUE)) { ntrcpt <- 0 } else {
    ntrcpt <- median(tst1$value[which((tst1$Group != "References")&
                                        (tst1$Norm == "Original"))])
  }
  plot <- plot + geom_vline(xintercept = ntrcpt, linetype = "dashed")
  print(plot) # This type of QC plot does not need to pop up, the side panel is fine
  if (save) {
    suppressMessages({
      ggsave(paste0(dir[1], "/", ttl, ".jpeg"), plot, dpi = 150)
      ggsave(paste0(dir[1], "/", ttl, ".pdf"), plot, dpi = 150)
    })
  }
}
