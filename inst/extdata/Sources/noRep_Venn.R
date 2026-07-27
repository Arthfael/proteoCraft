#### Venn diagrams (no-replicates workflow version)
setwd(wd)
packs <- c("ggplot2", "ggpolypath", "venn")
for (pck in packs) {
  if (!require(pck, character.only = TRUE)) {
    pak::pak(pck, upgrade = FALSE, ask = FALSE)
  }
}
for (pck in packs) {
  library(pck, character.only = TRUE)
}
#
if (!exists("plotly_Venn")) { plotly_Venn <- list() }
plotly_Venn %<o% plotly_Venn
#
HdrStlVenn <- createStyle(textDecoration = "bold", halign = "left", valign = "bottom", wrapText = TRUE,
                          numFmt = "TEXT", fontSize = 12L, textRotation = 60)
wb <- createWorkbook()
wbKount <- 0L
VennMx <- 7L
if (Venn_Obs) {
  dir <- paste0(wd, "/Venn diagrams")
  if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
  ref <- PG.int.cols["Original"]
  Grps <- unique(SamplesMap$`Ratios group`)
  for (grp in Grps) {
    sm <- SamplesMap[which(SamplesMap$`Ratios group` == grp),]
    Xp <- sm$Experiment
    test <- setNames(lapply(Xp, \(exp) {
      x <- PG[[paste0(ref, exp)]]
      which(is.finite(x))
    }), Xp)
    w <- which(lengths(test) > 0L)
    VennExp <- Xp[w]
    OK <- length(w) > 1L
    if (OK) {
      if (length(w) > VennMx) {
        msg <- paste0("Too many samples, select at least 2 and up to ", VennMx,
                      " to include in the Venn diagram (comma-separated):")
        if (length(Grps) > 1L) { msg <- paste0("Ratios group ", grp, ": ", msg) }
        tst <- !sm$Reference[match(VennExp, sm$Experiment)]
        tmp <- setNames(vapply(VennExp, \(x) { paste(c(x, rep(" ", 200L-nchar(x))), collapse = "") }, ""), VennExp)
        VennExp <- names(tmp)[match(dlg_list(tmp, tmp[which(tst)[1L:min(c(sum(tst, VennMx)))]], TRUE, msg)$res, tmp)]
        if (length(VennExp) < 2L) {
          msg <- "Skipping per-sample observations Venn diagrams"
          if (length(Grps) > 1L) { msg <- paste0(msg, " for ratios group ", grp) }
          msg <- paste0(msg, ": you should have selected at least 2 samples!")
          warning(msg)
          OK <- FALSE
        }
        if (length(VennExp) > VennMx) {
          msg <- "Skipping per-sample observations Venn diagrams"
          if (length(Grps) > 1L) { msg <- paste0(msg, " for ratios group ", grp) }
          msg <- paste0(msg, ": you should have selected ", VennMx,
                        " samples at most!")
          warning(msg)
          OK <- FALSE
        }
      }
    }
    if (OK) {
      ttl <- "Observations_LFQ_Venn_diagram_-_global"
      SheetNm <- "Samples composition"
      if (length(Grps) > 1L) {
        ttl <- paste0(ttl, " - group ", grp)
        SheetNm <- paste0(SheetNm, "_grp", grp)
      }
      test <- test[VennExp]
      AnalysisParam$"Venn diagram (LFC) - samples" <- list(VennExp)
      plot <- venn(test, ilabels = "counts", ellipse = TRUE, zcolor = "style", ggplot = TRUE)
      subTtl <- "Global, LFQ"
      plot <- plot + ggtitle("Venn diagram", subtitle = subTtl) +
        theme(plot.title = element_text(size = 15L), plot.subtitle = element_text(size = 10L))
      plot@layers$geom_path <- NULL
      plotLy <- ggplotly(plot, tooltip = "label")
      plotLy <- plotly::config(plotLy,
                               modeBarButtonsToRemove = c("select2d", "lasso2d"))
      plotly_Venn[[subTtl]] <- plotLy
      poplot(plot)
      suppressMessages({
        ggsave(paste0(dir, "/", ttl, ".jpg"), plot, dpi = 150L)
        ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 150L)
      })
      #system(paste0("open \"", dir, "/", ttl, ".jpg", "\""))
      wbKount <- wbKount+1L
      if (SheetNm %in% names(wb)) { removeWorksheet(wb, SheetNm) } 
      addWorksheet(wb, SheetNm)
      writeData(wb, SheetNm, PG[, c("id", "Leading protein IDs", "Genes")], 1L, 1L)
      l <- length(test)
      tmp <- sapply(names(test), \(smpl) {
        res <- rep("", nrow(PG))
        res[test[[smpl]]] <- "+"
        return(res)
      })
      writeData(wb, SheetNm, tmp, 4L, 1L)
      setRowHeights(wb, SheetNm, 1L, 120L)
      addStyle(wb, SheetNm, HdrStlVenn, 1L, 1L:(l+3L))
    } else { warning("Skipping LFC Venn diagrams: not enough valid samples!") }
  }
}
if (Venn_Ratios) {
  VennTypes <- ""
  if (RatiosThresh_2sided) { VennTypes <- c(VennTypes, "up", "down") }
  Grps <- unique(SamplesMap$`Ratios group`)
  for (grp in Grps) {
    sm <- SamplesMap[which(SamplesMap$`Ratios group` == grp),]
    VennExp <- sm$Experiment[which(!sm$Reference)]
    fc_filt <- FC_filt[VennExp]
    tst <- lengths(fc_filt)
    w <- which(tst == 0L)
    if (0L %in% tst) {
      l <- length(w)
      tmp <- names(fc_filt)[which(tst == 0L)]
      if (l > 1L) { tmp <- paste0(paste(tmp[1L:(l-1L)], collapse = ", "), " and ", tmp[l]) }
      warning(paste0("No filtered proteins for samples ", tmp, ", they will be skipped!"))
    }
    fc_filt <- fc_filt[which(tst > 0L)]
    VennExp <- names(fc_filt)
    OK <- length(VennExp) > 1L
    if (OK) {
      if (length(fc_filt) > VennMx) {
        msg <- paste0("Too many samples, select at least 2 and up to ", VennMx,
                      " to include in the Venn diagram (comma-separated):")
        if (length(Grps) > 1L) { msg <- paste0("Ratios group ", grp, ": ", msg) }
        VennExp <- dlg_list(VennExp, VennExp[1L:VennMx], TRUE, msg)$res
        if (length(VennExp) < 2L) {
          msg <- "Skipping ratios Venn diagrams"
          if (length(Grps) > 1L) { msg <- paste0(msg, " for ratios group ", grp) }
          msg <- paste0(msg, ": you should have selected at least 2 samples!")
          warning(msg)
          OK <- FALSE
        }
        if (length(VennExp) > VennMx) {
          msg <- "Skipping ratios Venn diagrams"
          if (length(Grps) > 1L) { msg <- paste0(msg, " for ratios group ", grp) }
          msg <- paste0(msg, ": you should have selected ", VennMx,
                        " samples at most!")
          warning(msg)
          OK <- FALSE
        }
      }
      if (OK) {
        fc_filt <- fc_filt[VennExp]
        AnalysisParam$"Venn diagram (ratios) - samples" <- list(VennExp)
        if (RatiosThresh_2sided) {
          updowntst <- setNames(lapply(VennExp, \(x) {
            #x <- VennExp[1L]
            rs <- sign(PG[fc_filt[[x]], paste0(PG.rat.cols["Original"], x)])
            w <- which(is.na(rs))
            rs[w] <- c(1L, -1L)[is.na(PG[fc_filt[[x]][w], paste0(PG.int.cols["Original"], x)])+1L]
            return(rs)
          }), VennExp)
        }
        setwd(wd); suppressWarnings(dir.create("Venn diagrams"))
        for (vt in VennTypes) { #vt <- ""
          ttl <- sub("_\\(\\)$", "", paste0("Ratios_Venn_diagram_-_global", "_(", vt, ")"))
          SheetNm <- paste0(c("Up/down", "Up", "Down")[match(vt, VennTypes)], "-reg. PGs")
          if (length(Grps) > 1L) {
            msg <- paste0(msg, " for ratios group ", grp)
            ttl <- paste0(ttl, " - group ", grp)
            SheetNm <- paste0(SheetNm, "_grp", grp)
          }
          subTtl <- rat.col
          flt <- fc_filt
          if (vt %in% c("up", "down")) {
            subTtl <- paste0(subTtl, ", ", vt)
            flt <- setNames(lapply(VennExp, \(x) {
              flt[[x]][which(updowntst[[x]] == c(1L, -1L)[match(vt, c("up", "down"))])]
            }), VennExp)
          }
          w <- which(lengths(flt) > 0L)
          if (length(w) > 1L) {
            flt <- flt[w]
            plot <- venn(flt, ilabels = "counts", ellipse = TRUE, zcolor = "style", ggplot = TRUE)
            plot <- plot + ggtitle("Venn diagram", subtitle = subTtl) +
              theme(plot.title = element_text(size = 15L), plot.subtitle = element_text(size = 10L))
            plot@layers$geom_path <- NULL
            plotLy <- ggplotly(plot, tooltip = "label")
            plotLy <- plotly::config(plotLy,
                                     modeBarButtonsToRemove = c("select2d", "lasso2d"))
            plotly_Venn[[subTtl]] <- plotLy
            poplot(plot)
            suppressMessages({
              ggsave(paste0(dir, "/", ttl, ".jpg"), plot, dpi = 150L)
              ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 150L)
            })
            #system(paste0("open \"", dir, "/", ttl, ".jpg", "\""))
            wbKount <- wbKount+1L
            if (SheetNm %in% names(wb)) { removeWorksheet(wb, SheetNm) } 
            addWorksheet(wb, SheetNm)
            writeData(wb, SheetNm, PG[, c("id", "Leading protein IDs", "Genes")], 1L, 1L)
            l <- length(flt)
            tmp <- sapply(names(flt), \(smpl) {
              res <- rep("", nrow(PG))
              res[flt[[smpl]]] <- "+"
              return(res)
            })
            writeData(wb, SheetNm, tmp, 4L, 1L)
            setRowHeights(wb, SheetNm, 1L, 120L)
            addStyle(wb, SheetNm, HdrStlVenn, 1L, 1L:(l+3L))
          } else { message(sub(" \\(\\)$", "", paste0("No overlaps to plot (", vt, ")"))) }
        }
      }
    } else {
      msg <- "Skipping ratios Venn diagrams"
      if (length(Grps) > 1L) { msg <- paste0(msg, " for ratios group ", grp) }
      msg <- paste0(msg, ": not enough valid samples!")
      warning(msg)
    }
  }
}
if (wbKount) { saveWorkbook(wb, paste0(wd, "/Venn diagrams/Venn diagrams.xlsx"), overwrite = TRUE) }
setwd(wd)
saveFun(plotly_Venn, paste0(dir, "/Venn_plotly.RDS"))
