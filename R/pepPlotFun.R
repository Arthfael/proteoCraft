#' pepPlotFun
#' 
#' @description 
#' Useful function for checking before/after peptides normalisation.
#' 
#' @param df1 Quantitative data prior to normalisation.
#' @param df2 Normalized quantitative data.
#' @param ttl Title.
#' @param dstDir Destination directory.
#' @param save Logical, should we save the plot (default = TRUE).
#' @param xpMap Experiment map, default = Exp.map
#' @param VPAL Volcano.plots.Aggregate.Level object, default = Volcano.plots.Aggregate.Level
#' 
#' @returns
#' This function does not return anything.
#' 
#' @export

pepPlotFun <- function(df1,
                         df2,
                         ttl,
                         dstDir,
                         save = TRUE,
                         xpMap = Exp.map,
                         VPAL = Volcano.plots.Aggregate.Level) {
  tst1 <- df1
  tst2 <- df2
  colnames(tst1) <- sub(".+ - ", "", colnames(tst1))
  colnames(tst2) <- sub(".+ - ", "", colnames(tst2))
  tst1 <- dfMelt(tst1)
  tst2 <- dfMelt(tst2)
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
  plot <- ggplot2::ggplot(tst1) +
    ggplot2::geom_density(stat = "density", alpha = 0.1,
                          ggplot2::aes(x = value, colour = variable, fill = variable)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::facet_grid(Norm~Group) + ggplot2::ggtitle(ttl) + ggplot2::theme_bw()
  if (grepl("ratio", ttl, ignore.case = TRUE)) { ntrcpt <- 0 } else {
    ntrcpt <- median(tst1$value[which((tst1$Group != "References")&
                                        (tst1$Norm == "Original"))])
  }
  plot <- plot + ggplot2::geom_vline(xintercept = ntrcpt, linetype = "dashed")
  print(plot) # This type of QC plot does not need to pop up, the side panel is fine
  if (save) {
    suppressMessages({
      ggplot2::ggsave(paste0(dstDir, "/", ttl, ".jpeg"), plot, dpi = 150)
      ggplot2::ggsave(paste0(dstDir, "/", ttl, ".pdf"), plot, dpi = 150)
    })
  }
  return()
}
