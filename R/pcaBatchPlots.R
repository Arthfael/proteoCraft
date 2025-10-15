#' pcaBatchPlots
#'
#' @description
#' A function to make PCA (ggplot and plotly) plots for the purpose of monitoring the effects of batch corrections. 
#' 
#' @param dat Input data matrix, expected to be log-transformed!
#' @param root Just a name to stick to the 
#' @param batches Relevant batches, default = myBatches; should be valid column names of map!
#' @param map Map connecting samples to batches. Note that this may produce incorrect results if dat has repeated column names!
#' @param SamplesCol Name of the sample colum in the map, default = "Ref.Sample.Aggregate"
#' @param intRoot Root of the expression columns in the file, default = "log10 Intensity - "; if a character of length 0, all columns are assumed to be expression columns!
#' @param dir Directory where to save the plots, default = btchDir
#' @param ttl Plot title root, default = "PCA plot - SVA batch corr."
#' @param openMe Should we open the plotly plot? (default = FALSE)
#' @param isRef Whether a column in dat is an IRS reference column (mixed sample).
#' @param make_Avg In case isRef is missing, should we plot also the average of each batch? FALSE by default.
#'
#' @examples
#'   scoresLst <- PCAlyLst <- PCsLst <- list()
#'   tmp <- pcaBatchPlots(imputedExpdata,
#'                        "original",
#'                        map = smplsMap,
#'                        SamplesCol = "New")
#'   PCAlyLst[["original"]] <- tmp$PlotLy
#'   scoresLst[["original"]] <- tmp$Scores
#'   PCsLst[["original"]] <- tmp$PCs
#'   tmp <- pcaBatchPlots(imputedExpdata2,
#'                        "corrected",
#'                        map = smplsMap,
#'                        SamplesCol = "New")
#'   PCAlyLst[["corrected"]] <- tmp$PlotLy
#'   scoresLst[["corrected"]] <- tmp$Scores
#'   PCsLst[["corrected"]] <- tmp$PCs
#'
#' @export

pcaBatchPlots <- function(dat, # Expected to be log-transformed!
                          root,
                          batches = myBatches,
                          map = Exp.map,
                          SamplesCol = "Ref.Sample.Aggregate",
                          intRoot = "log10 Intensity - ",
                          dir = btchDir,
                          ttl = "PCA plot - SVA batch corr.",
                          openMe = FALSE,
                          isRef,
                          make_Avg = FALSE) {
  #proteoCraft::DefArg(pcaBatchPlots)
  suppressWarnings({
    make_Avg <- as.logical(make_Avg)
    if (is.na(make_Avg)) { make_Avg <- FALSE }
  })
  if (missing(map)) { stop() }
  if (nchar(intRoot)) {
    pat <- proteoCraft::topattern(intRoot)
    dat <- dat[, grep(proteoCraft::topattern(intRoot), colnames(dat), value = TRUE)]
    colnames(dat) <- gsub(pat, "", colnames(dat))
  }
  stopifnot(batches %in% colnames(map),
            ncol(dat) > 1,
            nrow(dat) > 0)
  map$Sample_name <- map[[SamplesCol]]
  colnames(dat) <- proteoCraft::cleanNms(colnames(dat))
  map$Sample_name <- proteoCraft::cleanNms(map[[SamplesCol]])
  m <- match(colnames(dat), map$Sample_name)
  stopifnot(sum(is.na(m)) == 0)
  map <- map[m,]
  if (length(batches) > 1) {
    whichBatch <- map$myBatch <- do.call(paste, c(map[, batches], sep = " "))
  } else {
    whichBatch <- map$myBatch <- map[, batches]
  }
  pc0 <- prcomp(t(dat), scale. = TRUE)
  scores0 <- as.data.frame(pc0$x)
  kol <- colnames(scores0)
  scores0$Sample <- rownames(scores0) <- proteoCraft::cleanNms(rownames(scores0))
  scores0[, paste(batches, collapse = " ")] <- whichBatch
  scores0[, batches] <- map[, batches]
  scores0$Batch <- whichBatch
  #
  pv0 <- round(100*(pc0$sdev)^2 / sum(pc0$sdev^2), 0)
  pv0 <- pv0[which(pv0 > 0)]
  pv0_ <- paste0(root, ": ", paste(sapply(1:length(pv0), function(x) {
    paste0("PC", x, ": ", pv0[x], "%")
  }), collapse = ", "))
  #print(pv0_)
  scores0$Corrected <- root
  scores0$Type <- "Individual sample"
  #scores0$Batch <- do.call(paste, c(scores0[, batches], sep = " / "))
  if (!missing(isRef)) {
    make_Avg <- TRUE
    refRoot <- "IntRef."
    refType <- "Internal reference"
    wR <- which(isRef)
    if (length(wR)) {
      grps <- aggregate(colnames(dat)[wR], list(whichBatch[wR]), function(x) { list(unique(x)) })
      l <- vapply(grps$x, length, 1)
    }
  }
  if (make_Avg) {
    if ((missing(isRef))||(!length(wR))||(min(l) < 1)) {
      refRoot <- "Avg."
      refType <- "Batch average"
      grps <- aggregate(colnames(dat), list(whichBatch), function(x) { list(unique(x)) })
    }
    datAv0 <- setNames(lapply(grps$x, function(x) {
      apply(scores0[match(x, scores0$Sample), kol, drop = FALSE],
            2, function(y) {
              mean(proteoCraft::is.all.good(y))
            })
    }), grps$Group.1)
    datAv0 <- do.call(rbind, c(datAv0))
    datAv0 <- as.data.frame(datAv0)
    datAv0$Batch <- grps$Group.1
    datAv0[[batches]] <- map[match(grps$Group.1, map$myBatch), batches]
    datAv0$Sample <- #row.names(datAv0) <-
      paste0(refRoot, datAv0[[batches]])
    datAv0$Corrected <- root
    datAv0$Type <- "Batch average"
    scores0 <- rbind(scores0,
                     datAv0)
  }
  plotlyPCA <- list()
  if ("PC2" %in% colnames(scores0)) {
    for (btch in batches) { #btch <- batches[1]
      nm1 <- paste0(ttl, " (", root, ", color = ", btch, ")")
      scores1 <- scores0
      scores1$Batch <- scores1[[btch]]
      scores1$Batch <- factor(scores1$Batch)
      xLab <- paste0("PC1 = ", pv0[1], "%")
      yLab <- paste0("PC2 = ", pv0[2], "%")
      plot <- ggplot2::ggplot(scores1, aes(x = PC1, y = PC2, color = Batch)) +
        ggplot2::geom_point(ggplot2::aes(shape = Type)) +
        ggpubr::stat_conf_ellipse(aes(fill = Batch),
                                  alpha = 0.2, geom = "polygon", show.legend = FALSE) +
        ggplot2::coord_fixed() + ggplot2::theme_bw() +
        ggplot2::xlab(xLab) + ggplot2::ylab(yLab) +
        ggplot2::geom_hline(yintercept = 0, colour = "black") +
        ggplot2::geom_vline(xintercept = 0, colour = "black") +
        ggplot2::ggtitle(nm1#, subtitle = pv0_
                         ) +
        ggplot2::scale_color_viridis_d(option = "D") +
        ggrepel::geom_text_repel(ggplot2::aes(label = Sample),
                                 size = 2.5, show.legend = FALSE)
      #poplot(plot)
      if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
      ggplot2::ggsave(paste0(dir, "/", nm1, ".jpeg"), plot, dpi = 300, width = 10, height = 10, units = "in")
      ggplot2::ggsave(paste0(dir, "/", nm1, ".pdf"), plot, dpi = 300, width = 10, height = 10, units = "in")
      #
      if ("PC3" %in% colnames(scores1)) {
        if (make_Avg) {
          plotlyPCA1 <- plotly::plot_ly(scores1, x = ~PC1, y = ~PC2, z = ~PC3,
                                        color = ~Batch, colors = "viridis",
                                        text = ~Sample, type = "scatter3d", mode = "markers", #symbol = ~Type,
                                        symbols = as.character(seq_along(unique(scores1$Type))))
        } else {
          plotlyPCA1 <- plotly::plot_ly(scores1, x = ~PC1, y = ~PC2, z = ~PC3,
                                        color = ~Batch, colors = "viridis",
                                        text = ~Sample, type = "scatter3d", mode = "markers")
        }
        plotlyPCA1 <- plotly::add_trace(plotlyPCA1, scores1, x = ~PC1, y = ~PC2, z = ~PC3,
                                        type = "scatter3d", mode = "text", showlegend = FALSE)
      } else {
        if (make_Avg) {
          plotlyPCA1 <- plotly::plot_ly(scores1, x = ~PC1, y = ~PC2,
                                        color = ~Batch, colors = "viridis",
                                        text = ~Sample, type = "scatter", mode = "markers", #symbol = ~Type,
                                        symbols = as.character(seq_along(unique(scores1$Type))))
        } else {
          plotlyPCA1 <- plotly::plot_ly(scores1, x = ~PC1, y = ~PC2,
                                        color = ~Batch, colors = "viridis",
                                        text = ~Sample, type = "scatter", mode = "markers")
        }
        plotlyPCA1 <- plotly::add_trace(plotlyPCA1, scores1, x = ~PC1, y = ~PC2,
                                        type = "scatter", mode = "text", showlegend = FALSE)
      }
      plotlyPCA1 <- plotly::layout(plotlyPCA1, title = nm1)
      plotlyPCA[[btch]] <- plotlyPCA1
      suppressWarnings(htmlwidgets::saveWidget(plotlyPCA1, paste0(dir, "/", nm1, ".html")))
      if (openMe) { system(paste0("open \"", dir, "/", nm1, ".html")) }
    }
  } else { stop("PCA failed, investigate!") }
  w <- which(colnames(scores0) %in% kol)
  colnames(scores0)[w] <- paste0(root, "_", colnames(scores0)[w])
  return(list(Scores = scores0,
              PlotLy = plotlyPCA,
              PCs = pc0))
}
