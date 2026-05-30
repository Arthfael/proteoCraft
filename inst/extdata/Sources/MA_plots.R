###########################################################################################################
# MA plots                                                                                                #
###########################################################################################################
# This will create MA plots to check whether any variance vs intensity correction is required on the data (loess or vsn)
#
# What are the options here?
# We can have:
# - label-free, one experiment => skip
# - label-free, multiple MQ experiments (+/- fractions)
# - isobaric (+/- MQ experiments/fractions)
# Update:
ev2fr %<o% match(ev$"Raw file path", Frac.map$"Raw file") # Update again
#
if ((length(MQ.Exp) > 1L)||(LabelType == "Isobaric")) { # Should be always TRUE
  source(parSrc, local = FALSE)
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
      quntKol2 <- gsub(topattern(ev.ref["Original"]), "", quntKol)
      data2a <- as.data.table(data2[, quntKol])
      colnames(data2a) <- quntKol
      data2a$Group <- tmp
      data2a <- data2a[, lapply(.SD, sum, na.rm = TRUE), keyby = Group]
      colnames(data2a) <- c("Group", quntKol)
      data2a <- as.data.frame(data2a)
      data2a[, quntKol2] <- log10(data2a[, quntKol])
      data2a <- data2a[, which(!colnames(data2a) %in% quntKol)]
      data2a[, lsKl] <- data2[match(data2a$Group, tmp), lsKl]
      data2a$Group <- NULL
      data2 <- data2a
      smplKl <- "Label"
    }
    if (LabelType == "LFQ") {
      quntKol2 <- if (X == "Parent sample") { get("MQ.Exp") } else { get(X) }
      data2a <- data.table(Intensity = data2[[ev.col["Original"]]], Group = tmp)
      data2a <- data2a[, list(`log10(Intensity)` = sum(Intensity, na.rm = TRUE)),
                       keyby = Group]
      data2a$`log10(Intensity)` <- log10(data2a$`log10(Intensity)`)
      data2a <- as.data.frame(data2a)
      data2a[, lsKl] <- data2[match(data2a$Group, tmp), lsKl]
      data2a$Group <- NULL
      data2 <- spread(data2a, X, "log10(Intensity)")
      smplKl <- "Parent sample"
    }
    tmp <- data2[, quntKol2]
    avkol <- rowMeans(data2[, quntKol2], na.rm = TRUE)
    data2[, quntKol2] <- sweep(data2[, quntKol2], 1L, avkol, "-")/log10(2L)
    kol_ <- c("Modified sequence", Y)
    if ("PTM-enrich." %in% colnames(data2)) { kol_ <- c(kol_, "PTM-enrich.") }
    #data2 <- as.data.table(data2)
    data2 <- set_colnames(reshape2::melt(data2, id.vars = kol_),
                          c(kol_, X, "log2(Ratio)"))
    data2 <- as.data.frame(data2)
    data2$"Mean log10(Intensity)" <- avkol
    w1 <- is.finite(data2$"Mean log10(Intensity)")
    w2 <- is.finite(data2$"log2(Ratio)")
    w <- which(w1 & w2)
    if (length(w)) {
      data2 <- data2[w,]
      tst <- vapply(c(Y, X), \(x) { length(unique(data2[[x]])) }, 1L)
      wrpKl <- c(Y, X)
      wrpKl2 <- paste0("`", wrpKl, "`")
      if (1 %in% tst) {
        w <- which(tst > 1)
        wrpKl <- wrpKl[w]
        wrpKl2 <- wrpKl2[w]
        #if ((length(wrpKl) == 1)&&(grepl(" ", wrpKl))) { wrpKl <- paste0("`", wrpKl, "`") }
      } else {
        if ((tst[1L] >= tst[2L]*3) || (tst[1L] <= tst[2L]/3)) {
          wrpKl2 <- paste0("`", paste(wrpKl, collapse = " | "), "`")
          wrpKl <- paste(wrpKl, collapse = " | ")
          data2[[wrpKl]] <- as.factor(do.call(paste, c(data2[, c(X, Y)], sep = " | ")))
        }
      }
      lst <- lapply(wrpKl, \(k) { data2[[k]] })
      annot <- aggregate(data2$"log2(Ratio)", lst, \(x) {
        x <- x[which(is.finite(x))]
        c(paste0("Median: ", signif(median(x), 3L)), paste0("IQR: ", signif(IQR(x), 3)))
      })
      annot[, c("Median", "IQR")] <- do.call(as.data.frame, annot)
      annot$x <- NULL
      colnames(annot)[1L:length(wrpKl)] <- wrpKl
      filtFun <- \(x) { x[which(is.finite(x))] }
      annot$Amax <- max(filtFun(data2$"Mean log10(Intensity)"))*1.1
      annot$Amin <- min(filtFun(data2$"Mean log10(Intensity)"))*1.1
      annot$Mmax <- max(filtFun(data2$"log2(Ratio)"))*1.1
      annot$Mmin <- min(filtFun(data2$"log2(Ratio)"))*1.1
      annot2 <- annot[, c(wrpKl, "Amax", "Mmin", "Mmax")] 
      annot2 <- rbind(annot2, annot2)
      annot2$Tag <- c(annot$Median, annot$IQR)
      data2[[Y]] <- factor(data2[[Y]], levels = get(Y))
      ttl <- "MA plot - PSMs"
      if (length(grps) > 1L) { ttl <- paste0(ttl, " - ", grp) }
      ylim <- max(c(abs(c(annot$Mmax, annot$Mmin, (annot$Amax-annot$Amin)/4))))
      annot2$Y <- ylim*0.9
      w <- grep("^IQR: ", annot2$Tag)
      annot2$Y[w] <- -ylim*0.9
      if ("PTM-enrich." %in% colnames(data2)) {
        data2 <- data2[order(data2$"PTM-enrich."),]
        data2$"PTM-enrich." <- factor(data2$"PTM-enrich.", levels = c("", Modifs$`Full name`))
        data2$Alpha <- c(0.1, 1)[(data2$"PTM-enrich." != "")+1]
        myColors <- setNames(c("darkblue", "red"), c("Not-enriched/Flow-through", "Enriched"))
        plot <- ggplot(data2) +
          geom_scattermore(aes(x = `Mean log10(Intensity)`, y = `log2(Ratio)`, colour = `PTM-enrich.`#, alpha = Alpha
          ), shape = 16L, size = 1L#, alpha = 0.1
          )
        if (length(unique(data$`PTM-enrich.`)) == 1L) { plot <- plot + guides(colour = "none") }
      } else {
        plot <- ggplot(data2) + geom_scattermore(aes(x = `Mean log10(Intensity)`, y = `log2(Ratio)`,
                                                     colour = .data[[smplKl]]), size = 1L#, alpha = 0.1
        )
      }
      plot <- plot + scale_color_viridis(begin = 0.25, discrete = TRUE, option = "D") +
        geom_hline(yintercept = 0, colour = "grey") +
        xlab("A = mean log10(Intensity)") + ylab("M = log2(Ratio)") +
        geom_smooth(aes(x = `Mean log10(Intensity)`, y = `log2(Ratio)`), color = "purple", formula = "y ~ s(x, bs = 'cs')", # Note that this fails sometimes, may be a package version issue
                    linewidth = 0.1, linetype = "dashed", method = 'gam') +
        geom_text(data = annot2, aes(x = Amax, y = Y, label = Tag), hjust = 1, cex = 2L) +
        coord_fixed(log10(2)/log2(10)) + theme_bw() + ggtitle(ttl) + ylim(-ylim, ylim) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black"), axis.text = element_text(size = 3L),
              strip.text.x = element_text(angle = 0, hjust = 0, vjust = 0.5, size = 7),
              strip.text.y = element_text(angle = 0, hjust = 0, vjust = 0.5, size = 7))
      plot <- if (length(wrpKl) == 1L) {
        plot + facet_wrap(as.formula(paste0(" ~ ", wrpKl2)), drop = TRUE)
      } else {
        # X and Y each have length 1!
        plot + facet_grid(as.formula(paste0("`", Y, "` ~ `", X, "`")), drop = TRUE)
      }
      #poplot(plot, 12L, 22L)
      suppressMessages({
        ggsave(paste0(dir, "/", ttl, ".jpeg"), plot, width = 10L, height = 10L, units = "in")
        ggsave(paste0(dir, "/", ttl, ".pdf"), plot, width = 10L, height = 10L, units = "in")
      })
      ReportCalls <- AddPlot2Report(Space = FALSE)
    } else {
      msg <- paste0("Not enough valid data", c("", paste0(" for group ", grp))[(length(grps) > 1L) + 1L], " to draw an MA plot!")
      ReportCalls <- AddMsg2Report(Offset = TRUE, Space = FALSE, Warning = TRUE)
    }
  }
  ReportCalls$Calls[[LRepCalls]] <- append(ReportCalls$Calls[[LRepCalls]],
                                           "body_add_par(Report, \"\", style = \"Normal\")")
} else {
  stop("Uh, I think you have the wrong analysis pipeline here...\nwhere are my sample groups and replicates!?!?!")
}
