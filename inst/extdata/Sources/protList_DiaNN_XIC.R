# DiaNN: draw XIC for proteins of interest 
if (length(protlspep)) { # XICs
  for (indir in inDirs[which(SearchSoft == "DIANN")]) { #indir <- inDirs[which(SearchSoft == "DIANN")][1L]
    xicDir <- paste0(indir, "/report_xic")
    if (dir.exists(xicDir)) {
      XIC_fls <- list.files(xicDir, "\\.xic\\.parquet$", full.names = TRUE)
      #ms1Mob_fls <- list.files(xicDir, "\\.ms1_mobilogram\\.parquet$", full.names = TRUE)
      #ms2Mob_fls <- list.files(xicDir, "\\.ms2_mobilogram\\.parquet$", full.names = TRUE)
      if (length(XIC_fls)) {
        library(arrow)
        source(parSrc, local = FALSE)
        g <- grsep2(protlspep, ev$Proteins)
        u <- ev$"Mod. seq. (DiaNN format)"[g]
        tmp <- Frac.map$`Raw files name`
        clusterExport(parClust, list("tmp", "g", "u"), envir = environment())
        XICs <- parLapply(parClust, XIC_fls, \(x) {
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
        clusterExport(parClust, list("tmpFl1", "tmpFl2"), envir = environment())
        readr::write_rds(XICs, tmpFl1)
        invisible(clusterCall(parClust, \(x) {
          assign("XICs", readr::read_rds(tmpFl1), envir = .GlobalEnv)
          return()
        }))
        unlink(tmpFl1)
        for (pr in protlspep) { #pr <- protlspep[1L]
          xicDir2 <- paste0(wd, "/XIC/", pr)
          if (!dir.exists(xicDir2)) { dir.create(xicDir2, recursive = TRUE) }
          dirlist <- unique(c(dirlist, xicDir2))
          g <- grsep2(pr, XICs$Proteins)
          if (length(g)) {
            pkBnds <- Boundaries[which(Boundaries$Seq_Run %in% XICs$Seq_Run[g]),]
            u <- unique(XICs$"Modified sequence"[g])
            clusterExport(parClust, list("g", "xicDir2", "pr"), envir = environment())
            readr::write_rds(pkBnds, tmpFl2)
            invisible(clusterCall(parClust, \(x) {
              assign("pkBnds", readr::read_rds(tmpFl2), envir = .GlobalEnv)
              return()
            }))
            unlink(tmpFl2)
            invisible(parLapply(parClust, u, \(sq) { #sq <- u[1L] #sq <- u[2L]
              XIC <- XICs[g,]
              sq2 <- gsub("^_|_$", "", sq)
              ppXIC <- XIC[which(XIC$"Modified sequence" == sq),]
              yMax <- aggregate(ppXIC$value, list(ppXIC$File), max)
              xMin <- min(ppXIC$rt)
              bnds <- pkBnds[which(pkBnds$Seq_Run %in% ppXIC$Seq_Run),]
              bnds$yMax <- yMax$x[match(bnds$File, yMax$Group.1)]
              wMS1 <- which(ppXIC$feature == "ms1")
              wMS2 <- which(ppXIC$feature != "ms1")
              aNNOt <- aggregate(ppXIC[, c("PEP", "Quantity Quality")], list(ppXIC$File), \(x) { signif(mean(x, na.rm = TRUE), 3L) })
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
                  ggplot2::geom_line(data = ppXIC[wMS1,], ggplot2::aes(x = rt, y = value), color = "red",
                                     linewidth = 0.5)
              }
              if (length(wMS2)) {
                plot <- plot +
                  ggplot2::geom_line(data = ppXIC[wMS2,], ggplot2::aes(x = rt, y = value, color = feature),
                                     linewidth = 0.3)
              }
              plot <- plot +
                ggplot2::geom_text(data = aNNOt, ggplot2::aes(label = Text, y = y), x = xMin, hjust = 0, vjust = 1, size = 3) +
                ggplot2::scale_colour_viridis_d() +
                ggplot2::facet_wrap(~File, scales = "free_y", ) + ggplot2::theme_bw() +
                ggplot2::scale_y_continuous(expand = c(0L, 0L)) +
                ggplot2::ggtitle(paste0(pr, " ", sq2)) +
                ggplot2::xlab("Retention time") + ggplot2::ylab("Intensity")
              #proteoCraft::poplot(plot, 12L, 22L)
              #
              suppressMessages({
                ggplot2::ggsave(paste0(xicDir2, "/", sq2, ".jpeg"), plot, dpi = 450L, height = 10L, width = 15L)
                ggplot2::ggsave(paste0(xicDir2, "/", sq2, ".pdf"), plot, height = 10L, width = 10L)
              })
              return()
            }))
          }
        }
      }
    }
  }
}
