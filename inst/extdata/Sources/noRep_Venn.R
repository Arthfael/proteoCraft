#### Venn diagrams (no-replicates workflow version)
setwd(wd)
library(ggplot2)
library(ggpolypath)
library(venn)
library(plotly)
library(htmlwidgets)
#
if (!exists("plotly_Venn")) { plotly_Venn <- list() }
#
HdrStlVenn <- openxlsx2::create_cell_style(num_fmt_id = "General",
                                           horizontal = "left",
                                           vertical = "bottom",
                                           text_rotation = 60,
                                           wrap_text = TRUE)
wb <- openxlsx2::wb_workbook()
wb$styles_mgr$add(HdrStlVenn, "Header_style")
wbKount <- 0L
myData <- PG
idKol <- "Leading protein IDs"
VennMx <- 7L
VennDir <- paste0(wd, "/Venn diagrams")
if (!dir.exists(VennDir)) { dir.create(VennDir, recursive = TRUE) }
Venn_fl %<o% paste0(VennDir, "/Venn_plotly.RDS")
if (Venn_Obs) {
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
      plot <- venn::venn(test, ilabels = "counts", ellipse = TRUE, zcolor = "style", ggplot = TRUE)
      subTtl <- "Global, LFQ"
      plot <- plot + ggtitle("Venn diagram", subtitle = subTtl) +
        theme(plot.title = element_text(size = 15L), plot.subtitle = element_text(size = 10L))
      plot@layers$geom_path <- NULL
      plotLy <- ggplotly(plot, tooltip = "label")
      plotLy$x$layout$xaxis$autorange <- TRUE
      plotLy$x$layout$yaxis$autorange <- TRUE
      plotLy <- htmlwidgets::onRender(plotLy, jsFun$global_auto_range)
      plotLy <- plotly::config(plotLy,
                               modeBarButtonsToRemove = c("select2d", "lasso2d"))
      plotLy <- plotly::plotly_build(plotLy)
      plotly_Venn[[subTtl]] <- plotLy
      poplot(plot)
      suppressMessages({
        ggsave(paste0(VennDir, "/", ttl, ".svg"), plot, dpi = 150L)
      })
      #system(paste0("open \"", VennDir, "/", ttl, ".jpg", "\""))
      wbKount <- wbKount+1L
      if (SheetNm %in% names(wb)) { removeWorksheet(wb, SheetNm) }
      if (SheetNm %in% wb_get_sheet_names(wb)) { wb <- openxlsx2::wb_remove_worksheet(wb, SheetNm) }
      wb <- openxlsx2::wb_add_worksheet(wb, SheetNm)
      l <- length(test)
      tmp <- sapply(VennExp, \(grp) {
        res <- rep("", nrow(myData))
        res[test[[grp]]] <- "+"
        return(res)
      })
      tmp <- as.data.frame(tmp)
      tmp[[idKol]] <- myData[[idKol]]
      tmp <- tmp[, c(idKol, VennExp)]
      wb <- openxlsx2::wb_add_data(wb, SheetNm, tmp, openxlsx2::wb_dims(2L, 2L))
      wb <- openxlsx2::wb_set_row_heights(wb, SheetNm, 2L, 120L)
      wb <- openxlsx2::wb_set_col_widths(wb, SheetNm, 1L, 3L)
      wb <- openxlsx2::wb_set_col_widths(wb, SheetNm, 1L, 10L)
      hdrDms <- openxlsx2::wb_dims(2L, 1L:(l+1L) + 1L)
      wb <- openxlsx2::wb_add_font(wb, SheetNm, hdrDms, size = 12L, bold = TRUE)
      wb <- openxlsx2::wb_set_cell_style(wb, SheetNm, hdrDms, wb$styles_mgr$get_xf_id("Header_style"))
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
      tmp <- .commAnd(names(fc_filt)[which(tst == 0L)])
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
        setwd(wd)
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
            plotLy$x$layout$xaxis$autorange <- TRUE
            plotLy$x$layout$yaxis$autorange <- TRUE
            plotLy <- htmlwidgets::onRender(plotLy, jsFun$global_auto_range)
            plotLy <- plotly::config(plotLy,
                                     modeBarButtonsToRemove = c("select2d", "lasso2d"))
            plotLy <- plotly::plotly_build(plotLy)
            plotly_Venn[[subTtl]] <- plotLy
            poplot(plot)
            suppressMessages({
              ggsave(paste0(VennDir, "/", ttl, ".svg"), plot, dpi = 150L)
            })
            #system(paste0("open \"", VennDir, "/", ttl, ".jpg", "\""))
            wbKount <- wbKount+1L
            if (SheetNm %in% wb_get_sheet_names(wb)) { wb <- openxlsx2::wb_remove_worksheet(wb, SheetNm) }
            wb <- openxlsx2::wb_add_worksheet(wb, SheetNm)
            l <- length(flt)
            tmp <- sapply(names(flt), \(smpl) {
              res <- rep("", nrow(myData))
              res[flt[[smpl]]] <- "+"
              return(res)
            })
            tmp <- as.data.frame(tmp)
            tmp[[idKol]] <- myData[[idKol]]
            tmp <- tmp[, c(idKol, names(flt))]
            wb <- openxlsx2::wb_add_data(wb, SheetNm, tmp, openxlsx2::wb_dims(2L, 2L))
            wb <- openxlsx2::wb_set_row_heights(wb, SheetNm, 2L, 120L)
            wb <- openxlsx2::wb_set_col_widths(wb, SheetNm, 1L, 3L)
            wb <- openxlsx2::wb_set_col_widths(wb, SheetNm, 1L, 10L)
            hdrDms <- openxlsx2::wb_dims(2L, 1L:(l+1L) + 1L)
            wb <- openxlsx2::wb_add_font(wb, SheetNm, hdrDms, size = 12L, bold = TRUE)
            wb <- openxlsx2::wb_set_cell_style(wb, SheetNm, hdrDms, wb$styles_mgr$get_xf_id("Header_style"))
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
if (wbKount) {
  vennXLfl <- paste0(wd, "/Venn diagrams/Venn diagrams.xlsx")
  openxlsx2::wb_save(wb, vennXLfl)
}
setwd(wd)
for (nm in names(plotly_Venn)) { #nm <- names(plotly_Venn)[1L]
  p <- plotly::plotly_build(plotly_Venn[[nm]])
  p$x$layout$xaxis$autorange <- TRUE
  p$x$layout$yaxis$autorange <- TRUE
  p <- htmlwidgets::onRender(p, jsFun$global_auto_range)
  plotly_Venn[[nm]] <- plotly::config(p,
                                      modeBarButtonsToRemove = c("select2d", "lasso2d"))
}
saveFun(plotly_Venn, Venn_fl)
