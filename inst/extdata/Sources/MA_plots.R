###########################################################################################################
# MA plots                                                                                                #
###########################################################################################################
# This will create MA plots to check whether any average logFC or logFC variance vs intensity biases exist.
# We will restrict ourselves to well observed peptides (80% non missing)!
#
#
# What are the options here?
# We can have:
# - label-free, one experiment => skip
# - label-free, multiple MQ experiments (+/- fractions)
# - isobaric (+/- MQ experiments/fractions)
#

# Update:
ev2fr %<o% match(ev$"Raw file path", Frac.map$"Raw file") # Update again
#
MAplotFls %<o% c()
if ((length(MQ.Exp) <= 1L) && (LabelType != "Isobaric")) { stop("Whut? Something's gone wrong!") } # Should be always TRUE
source(parSrc, local = FALSE)
invisible(clusterCall(parClust, \() {
  library(ggplot2)
  library(scattermore)
  library(ggnewscale)
  library(grid)
  library(mgcv)
  return()
}))
#
data <- ev
colnames(data)[which(colnames(data) == "MQ.Exp")] <- "Parent sample"
if (("PTM-enriched" %in% colnames(Frac.map))&&(sum(Modifs$"Full name" %in% Frac.map$"PTM-enriched"))) {
  data$"PTM-enrich." <- Frac.map$"PTM-enriched"[ev2fr]
  data$"PTM-enrich."[which(is.na(data$"PTM-enrich."))] <- ""
}
data <- data[which(data$Reverse != "+"),]
data <- data[which((is.na(data$"Potential contaminant"))|(data$"Potential contaminant" != "+")),]
quntKol <- if (LabelType == "Isobaric") {
  grep(paste0(topattern(ev.ref["Original"]), "[0-9]+$"), colnames(data), value = TRUE)
} else { ev.col["Original"] }
w <- which(rowSums(data[, quntKol, drop = FALSE], na.rm = TRUE) > 0)
data <- data[w,]
if (!"Fraction" %in% colnames(data)) { data$Fraction <- 1L }
Fraction <- sort(unique(data$Fraction), decreasing = FALSE) # actually useful!
Experiment <- Exp # This too!
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
  Y <- if (length(w1) > 1L) { kols[w1[2L]] } else { kols[w2[2L]] }
}
# X and Y each have length 1!
kols <- setdiff(kols, c(X, Y))
if (length(kols) > 1L) {
  data$Group <- do.call(paste, c(data[, kols], sep = " "))
  Grpkol <- "Group"
} else { Grpkol <- kols }
grps <- sort(unique(data[[Grpkol]]))
ReportCalls <- AddSpace2Report()
ReportCalls <- AddTxt2Report(paste0("MA plot", c("", "s")[(length(grps) > 1L)+1L], ":"))
ReportCalls$Objects$MA_groups <- c()
ReportCalls$Plots$MA_plots <- list()
ReportCalls$Calls <- append(ReportCalls$Calls, list())
LRepCalls <- length(ReportCalls$Calls)
dir <- paste0(wd, "/Workflow control/MA plots")
if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
dirlist <- unique(c(dirlist, dir))
#
for (grp in grps) { #grp <- grps[1L]
  data2 <- data[which(data[[Grpkol]] == grp),]
  lsKl <- c("Modified sequence", Y)
  if (LabelType == "LFQ") { lsKl <- c(lsKl, X) }
  if ("PTM-enrich." %in% colnames(data2)) { lsKl <- c(lsKl, "PTM-enrich.") }
  ls <- lapply(lsKl, \(kl) { data2[[kl]] })
  tmp <- do.call(paste, c(data2[, lsKl], sep = "---"))
  if (LabelType == "Isobaric") {
    smpls <- gsub(topattern(ev.ref["Original"]), "", quntKol)
    data2a <- as.data.table(data2[, quntKol])
    colnames(data2a) <- quntKol
    data2a$Group <- tmp
    data2a <- data2a[, lapply(.SD, sum, na.rm = TRUE), keyby = Group]
    colnames(data2a) <- c("Group", quntKol)
    data2a <- as.data.frame(data2a)
    data2a[, smpls] <- log10(data2a[, quntKol])
    data2a <- data2a[, which(!colnames(data2a) %in% quntKol)]
    data2a[, lsKl] <- data2[match(data2a$Group, tmp), lsKl]
    data2a$Group <- NULL
    data2 <- data2a
    smplKl <- "Label"
  }
  if (LabelType == "LFQ") {
    smpls <- if (X == "Parent sample") { get("MQ.Exp") } else { get(X) }
    data2a <- data.table(Intensity = data2[[ev.col["Original"]]], Group = tmp)
    data2a <- data2a[, .(`A (mean log10 Intensity)` = sum(Intensity, na.rm = TRUE)),
                     keyby = Group]
    data2a$`A (mean log10 Intensity)` <- log10(data2a$`A (mean log10 Intensity)`)
    data2a <- as.data.frame(data2a)
    data2a[, lsKl] <- data2[match(data2a$Group, tmp), lsKl]
    data2a$Group <- NULL
    data2 <- spread(data2a, X, "A (mean log10 Intensity)")
    smplKl <- "Parent sample"
  }
  kol_ <- c("Modified sequence", Y)
  if ("PTM-enrich." %in% colnames(data2)) { kol_ <- c(kol_, "PTM-enrich.") }
  
  #w <- which(apply(data2[, smpls], 1L, \(x) { sum(is.finite(x)) })/length(smpls) >= 0.8) # If we wanted to filter for completeness
  w <- 1L:nrow(tmpDat)
  if (length(w)) {
    avg <- as.matrix(data2[w, smpls])
    avg <- rowMeans(replace(avg, !is.finite(avg), NA), na.rm = TRUE)
    tmpDat_M <- (data2[w,] - avg)/log10(2L)
    tmpDat_A <- (data2[w,] + avg)/2
    tmpDat_M <- dfMelt(tmpDat_M, ColNames = c("Sample", "M (mean log2 FC)"))
    tmpDat_A <- dfMelt(tmpDat_A, ColNames = c("Sample", "A (mean log10 Intensity)"))
    tmpDat <- tmpDat_A
    tmpDat$"M (mean log2 FC)" <- tmpDat_M$"M (mean log2 FC)"
    rm(tmpDat_A, tmpDat_M)
    # removing redundant filters: we are already filtering which data we plot!
    #tmpDat <- tmpDat[which(is.finite(tmpDat$"M (mean log2 FC)")),]
    #tmpDat <- tmpDat[which(is.finite(tmpDat$"A (mean log10 Intensity)")),]
    tst <- vapply(c(Y, X), \(x) { length(unique(tmpDat[[x]])) }, 1L)
    wrpKl <- c(Y, X)
    wrpKl2 <- paste0("`", wrpKl, "`")
    if (1L %in% tst) {
      w <- which(tst > 1L)
      wrpKl <- wrpKl[w]
      wrpKl2 <- wrpKl2[w]
      #if ((length(wrpKl) == 1)&&(grepl(" ", wrpKl))) { wrpKl <- paste0("`", wrpKl, "`") }
    } else {
      if ((tst[1L] >= tst[2L]*3) || (tst[1L] <= tst[2L]/3)) {
        wrpKl2 <- paste0("`", paste(wrpKl, collapse = " | "), "`")
        wrpKl <- paste(wrpKl, collapse = " | ")
        tmpDat[[wrpKl]] <- as.factor(do.call(paste, c(tmpDat[, c(X, Y)], sep = " | ")))
      }
    }
    lst <- lapply(wrpKl, \(k) { tmpDat[[k]] })
    annot <- aggregate(tmpDat$"M (mean log2 FC)", lst, \(x) {
      x <- x[which(is.finite(x))]
      c(paste0("Median: ", signif(median(x), 3L)),
        paste0("IQR: ", signif(IQR(x), 3L)))
    })
    annot[, c("Median", "IQR")] <- do.call(as.data.frame, annot)
    annot$x <- NULL
    colnames(annot)[1L:length(wrpKl)] <- wrpKl
    filtFun <- \(x) { x[which(is.finite(x))] }
    annot$Amax <- max(filtFun(tmpDat$"A (mean log10 Intensity)"))*1.1
    annot$Amin <- min(filtFun(tmpDat$"A (mean log10 Intensity)"))*1.1
    annot$Mmax <- max(filtFun(tmpDat$"M (mean log2 FC)"))*1.1
    annot$Mmin <- min(filtFun(tmpDat$"M (mean log2 FC)"))*1.1
    annot2 <- annot[, c(wrpKl, "Amax", "Mmin", "Mmax")] 
    annot2 <- rbind(annot2, annot2)
    annot2$Tag <- c(annot$Median, annot$IQR)
    tmpDat[[Y]] <- factor(tmpDat[[Y]], levels = get(Y))
    MAttl <- "MA plot - PSMs"
    if (length(grps) > 1L) { MAttl <- paste0(MAttl, " - ", grp) }
    ylim <- max(c(abs(c(annot$Mmax, annot$Mmin, (annot$Amax-annot$Amin)/4))))
    annot2$Y <- ylim*0.9
    w <- grep("^IQR: ", annot2$Tag)
    annot2$Y[w] <- -ylim*0.9
    if ("PTM-enrich." %in% colnames(tmpDat)) {
      tmpDat <- tmpDat[order(tmpDat$"PTM-enrich."),]
      tmpDat$"PTM-enrich." <- factor(tmpDat$"PTM-enrich.", levels = c("", Modifs$`Full name`))
      tmpDat$Alpha <- c(0.1, 1)[(tmpDat$"PTM-enrich." != "")+1L]
      myColors <- setNames(c("darkblue", "red"), c("Not-enriched/Flow-through", "Enriched"))
      plot <- ggplot(tmpDat[1L:5L,]) +
        geom_scattermore(aes(x = `A (mean log10 Intensity)`, y = `M (mean log2 FC)`, color = `PTM-enrich.`),
                         shape = 16L, size = 1L)
      if (length(unique(data$`PTM-enrich.`)) == 1L) { plot <- plot + guides(color = "none") }
    } else {
      plot <- ggplot(tmpDat[1L:5L,]) +
        geom_scattermore(aes(x = `A (mean log10 Intensity)`, y = `M (mean log2 FC)`, color = .data[[smplKl]]),
                         size = 1L)
    }
    plot <- plot +
      geom_hline(yintercept = 0, color = "grey", linewidth = 0.3, linetype = "dashed") +
      scale_fill_viridis(begin = 0.25, end = 0.75, discrete = TRUE, option = "H") +
      xlab("A = log10(Intensity)") + ylab("M = log2(Ratio)") +
      coord_fixed(0.05) + theme_bw() + ylim(-ylim, ylim) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"), axis.text = element_text(size = 3L),
            strip.text.x = element_text(angle = 0, hjust = 0, vjust = 0.5, size = 7),
            strip.text.y = element_text(angle = 0, hjust = 0, vjust = 0.5, size = 7))
    #
    # Finish me: here would be a good point to generate individual MA plots
    MAfl <- paste0(dir, "/", MAttl)
    tmpFl <- tempfile(fileext = ".RDS")
    readr::write_rds(tmpDat, tmpFl)
    clusterExport(parClust,
                  list("plot", "tmpFl", "MAttl", "MAfl", "smplKl", "annot2"),
                  envir = environment())
    invisible(clusterCall(parClust, \() {
      assign("tmpDat", readr::read_rds(tmpFl), envir = .GlobalEnv)
      return()
    }))
    unlink(tmpFl)
    MAplotFls[smpls] <- parSapply(parClust, smpls, \(smpl) { #smpl <- smpls[1L] 
      dat <- tmpDat[which(tmpDat$`Parent sample` == smpl),]
      ann <- annot2[which(annot2$`Parent sample` == smpl),]
      plot$data <- dat
      ttl <- paste0(MAttl, " - ", smpl)
      plot <- plot + ggtitle(ttl) +
        geom_text(data = ann, aes(x = Amax, y = Y, label = Tag), hjust = 1, cex = 2L) +
        new_scale_color() +
        geom_smooth(aes(x = `A (mean log10 Intensity)`, y = `M (mean log2 FC)`, color = "GAM", linetype = "GAM"),
                    method = "gam", formula = y ~ s(x, bs = "cs")) +
        geom_smooth(aes(x = `A (mean log10 Intensity)`, y = `M (mean log2 FC)`, color = "LOESS", linetype = "LOESS"),
                    method = "loess", formula = y ~ x, se = FALSE) +
        scale_color_manual(name = "Method", values = c(LOESS = "green", GAM = "magenta")) +
        scale_linetype_manual(name = "Method", values = c(LOESS = "dotted", GAM = "dashed"))
      # Estimate plot dimensions
      gt <- ggplotGrob(plot)
      width  <- sum(gt$widths)
      height <- sum(gt$heights)
      a <- convertWidth(width, "in", valueOnly = TRUE)
      b <- convertHeight(height, "in", valueOnly = TRUE)
      #poplot(plot, 12L, 22L)
      fl <- paste0(MAfl, " - ", smpl)
      ggsave(paste0(fl, ".jpeg"), plot, width = 10L, height = 10L*b/a, units = "in")
      ggsave(paste0(fl, ".pdf"), plot, width = 10L, height = 10L*b/a, units = "in")
      return(fl)
    })
    # You could then visualize those individually in the parameters app as opposed to seeing a very zoomed out facet plot
    plot$data <- tmpDat
    plot <- if (length(wrpKl) == 1L) {
      plot + facet_wrap(as.formula(paste0(" ~ ", wrpKl2)), drop = TRUE)
    } else {
      # X and Y each have length 1!
      plot + facet_grid(as.formula(paste0("`", Y, "` ~ `", X, "`")), drop = TRUE)
    }
    plot <- plot + ggtitle(MAttl) +
      geom_text(data = annot2, aes(x = Amax, y = Y, label = Tag), hjust = 1, cex = 2L)
    # Estimate plot dimensions
    library(grid)
    gt <- ggplotGrob(plot)
    width <- sum(gt$widths)
    height <- sum(gt$heights)
    a <- convertWidth(width, "in", valueOnly = TRUE)
    b <- convertHeight(height, "in", valueOnly = TRUE)
    #poplot(plot, 12L, 22L)
    suppressMessages({
      ggsave(paste0(MAfl, ".jpeg"), plot, width = 10L, height = 10L*b/a, units = "in")
      ggsave(paste0(MAfl, ".pdf"), plot, width = 10L, height = 10L*b/a, units = "in")
    })
    ReportCalls <- AddPlot2Report(Space = FALSE, Title = MAttl)
  } else {
    msg <- paste0("Not enough valid data", c("", paste0(" for group ", grp))[(length(grps) > 1L) + 1L], " to draw an MA plot!")
    ReportCalls <- AddMsg2Report(Offset = TRUE, Space = FALSE, Warning = TRUE)
  }
}
ReportCalls$Calls[[LRepCalls]] <- append(ReportCalls$Calls[[LRepCalls]],
                                         "body_add_par(Report, \"\", style = \"Normal\")")
