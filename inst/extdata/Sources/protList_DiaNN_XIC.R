# Peptide XIC extraction - core script
dianXDirs <- inDirs[which(SearchSoft == "DIANN")]
if (length(dianXDirs)) {
  xicDirs <- setNames(paste0(dianXDirs, "/report_xic"),
                      dianXDirs)
  xicDirs <- xicDirs[which(dir.exists(xicDirs))]
  dianXDirs <- names(xicDirs)
}
if (length(dianXDirs)) {
  xicFiles <- setNames(lapply(dianXDirs, \(dr) {
    list.files(xicDirs[dr], "\\.xic\\.parquet$", full.names = TRUE)
  }), dianXDirs)
  xicFiles <- xicFiles[which(lengths(xicFiles) > 0L)]
  dianXDirs <- names(xicFiles)
  xicDirs <- xicDirs[dianXDirs]
}
if (length(dianXDirs)) {
  library(arrow)
  for (inDir in dianXDirs) { #inDir <- dianXDirs[1L] #inDir <- dianXDirs[2L]
    xicDir <- xicDirs[inDir]
    XIC_fls <- xicFiles[[inDir]]
    #ms1Mob_fls <- list.files(xicDir, "\\.ms1_mobilogram\\.parquet$", full.names = TRUE)
    #ms2Mob_fls <- list.files(xicDir, "\\.ms2_mobilogram\\.parquet$", full.names = TRUE)
    source(parSrc, local = FALSE)
    if (dataType == "protList") {
      g <- grsep2(protlspep, ev$Proteins)
      xicOutDir <- paste0(wd, "/XIC")
    }
    if (dataType == "modPeptides") {
      g <- which(ev$`Modified sequence` %in% regPep)
      xicOutDir <- paste0(modDirs[1L], "/XIC")
    }
    if (!dir.exists(xicOutDir)) { dir.create(xicOutDir, recursive = TRUE) }
    myProts_pepList <- listMelt(strsplit(ev$Proteins[g], ";"),
                                ev$"Mod. seq. (DiaNN format)"[g])
    myProts_pepList <- aggregate(myProts_pepList$L1, list(myProts_pepList$value), unique)
    myProts_pepList <- setNames(myProts_pepList$x, myProts_pepList$Group.1)
    myProts <- names(myProts_pepList)
    u <- unique(unlist(myProts_pepList))
    tmp <- Frac.map$`Raw files name`
    clusterExport(parClust, list("tmp", "u"), envir = environment())
    XICs <- parLapply(parClust, XIC_fls, \(x) { #x <- XIC_fls[1L]
      res <- arrow::read_parquet(x)
      res$"Mod. seq." <- proteoCraft::gsub_Rep("[0-9]+$", "", res$pr)
      res <- res[which(res$"Mod. seq." %in% u),]
      nm <- gsub(".*/|\\.xic\\.parquet$", "", x)
      res$File <- nm
      res$"Seq_Run" <- do.call(paste, c(res[, c("pr", "File")], sep = ">>>"))
      res$File <- factor(res$File, levels = tmp)
      res <- res[which((!is.na(res$File)) & (res$feature != "index")),]
      return(res)
    })
    #View(XICs[[1L]])
    XICs <- plyr::rbind.fill(XICs)
    if (nrow(XICs)) {
      #View(XICs[1L:100L,])
      #
      m <- match(XICs$"Mod. seq.", ev$"Mod. seq. (DiaNN format)")
      myKol <- c("Proteins", "Sequence", "Modified sequence", "PEP", "Quantity Quality")
      XICs[, myKol] <- ev[m, myKol]
      ev$tmp <- ">>>"
      Boundaries <- do.call(paste0, c(ev[, c("Mod. seq. (DiaNN format)", "Charge", "tmp", "Raw file")]))
      ev$tmp <- NULL
      w <- which(Boundaries %in% XICs$Seq_Run)
      Boundaries <- data.frame(Seq_Run = Boundaries[w],
                               RT = ev$`Retention time`[w],
                               `RT (start)` = ev$`Retention time (start)`[w],
                               `RT (end)` = ev$`Retention time (end)`[w],
                               check.names = FALSE)
      Boundaries$File <- gsub(".*>>>", "", Boundaries$Seq_Run)
      tmpFl1 <- tempfile(fileext = ".rds")
      tmpFl2 <- tempfile(fileext = ".rds")
      clusterExport(parClust, list("tmpFl1", "tmpFl2", "myProts_pepList", "xicOutDir"), envir = environment())
      readr::write_rds(XICs, tmpFl1)
      readr::write_rds(Boundaries, tmpFl2)
      invisible(clusterCall(parClust, \() {
        assign("XICs", readr::read_rds(tmpFl1), envir = .GlobalEnv)
        assign("Boundaries", readr::read_rds(tmpFl2), envir = .GlobalEnv)
        return()
      }))
      unlink(tmpFl1)
      dirlist <- union(dirlist, paste0(xicOutDir, "/", myProts))
      invisible(parLapply(parClust, myProts, \(pr) { #pr <- myProts[1L]
        xicDir2 <- paste0(xicOutDir, "/", pr)
        if (!dir.exists(xicDir2)) { dir.create(xicDir2, recursive = TRUE) }
        pp <- myProts_pepList[[pr]]
        g <- which(XICs$`Mod. seq.` %in% pp)
        if (length(g)) {
          pkBnds <- Boundaries[which(Boundaries$Seq_Run %in% XICs$Seq_Run[g]),]
          lapply(pp, \(sq) { #sq <- pp[1L] #sq <- pp[11L]
            w <- which(XICs$`Mod. seq.` == sq)
            if (!length(w)) { return() } # I have no idea why this sometimes occurs, but it does!
            XIC <- XICs[w,]
            yMax <- aggregate(XIC$value, list(XIC$File), max)
            xMin <- min(XIC$rt)
            bnds <- pkBnds[which(pkBnds$Seq_Run %in% XIC$Seq_Run),]
            bnds$yMax <- yMax$x[match(bnds$File, yMax$Group.1)]
            wMS1 <- which(XIC$feature == "ms1")
            wMS2 <- which(XIC$feature != "ms1")
            aNNOt <- aggregate(XIC[, c("PEP", "Quantity Quality")], list(XIC$File), \(x) {
              signif(mean(x, na.rm = TRUE), 3L)
            })
            colnames(aNNOt)[1L] <- "File"
            aNNOt$PEP <- paste0("PEP = ", aNNOt$PEP)
            aNNOt$"Quantity Quality" <- paste0("Quantity Quality = ", aNNOt$"Quantity Quality")
            aNNOt$Text <- do.call(paste, c(aNNOt[, c("PEP", "Quantity Quality")], sep = "\n"))
            aNNOt$y <- yMax$x[match(aNNOt$File, yMax$Group.1)]*0.9
            plot <- ggplot2::ggplot() + ggplot2::scale_y_continuous(expand = c(0L, 10L))
            if (nrow(bnds)) {
              plot <- plot +
                ggplot2::geom_rect(data = bnds, ggplot2::aes(xmin = `RT (start)`, ymin = 0, xmax = `RT (end)`, ymax = yMax),
                                   fill = "lightblue", alpha = 0.2) +
                ggplot2::geom_vline(data = bnds, ggplot2::aes(xintercept = RT),
                                    color = "darkblue", linewidth = 0.5) +
                ggplot2::geom_vline(data = bnds, ggplot2::aes(xintercept = `RT (start)`),
                                    color = "darkblue", linewidth = 0.5, linetype = "dashed") +
                ggplot2::geom_vline(data = bnds, ggplot2::aes(xintercept = `RT (end)`),
                                    color = "darkblue", linewidth = 0.5, linetype = "dashed")
            }
            if (length(wMS1)) {
              plot <- plot +
                ggplot2::geom_line(data = XIC[wMS1,], ggplot2::aes(x = rt, y = value), color = "red",
                                   linewidth = 0.5)
            }
            if (length(wMS2)) {
              plot <- plot +
                ggplot2::geom_line(data = XIC[wMS2,], ggplot2::aes(x = rt, y = value, color = feature),
                                   linewidth = 0.3)
            }
            sq1 <- gsub("^_|_$", "", sq)
            sq2 <- gsub(":", "_", sq1)
            plot <- plot +
              ggplot2::geom_text(data = aNNOt, ggplot2::aes(label = Text, y = y), x = xMin, hjust = 0, vjust = 1, size = 3) +
              ggplot2::scale_colour_viridis_d() +
              ggplot2::facet_wrap(~File, scales = "free_y", ) + ggplot2::theme_bw() +
              ggplot2::scale_y_continuous(expand = c(0L, 0L)) +
              ggplot2::ggtitle(sq1, subtitle = pr) +
              ggplot2::xlab("Retention time") + ggplot2::ylab("Intensity")
            #proteoCraft::poplot(plot, 12L, 22L)
            #
            suppressMessages({
              ggplot2::ggsave(paste0(xicDir2, "/", sq2, ".jpeg"), plot, dpi = 450L, height = 10L, width = 15L)
              ggplot2::ggsave(paste0(xicDir2, "/", sq2, ".pdf"), plot, height = 10L, width = 10L)
            })
          })
        }
        return()
      }))
    }
  }
}
