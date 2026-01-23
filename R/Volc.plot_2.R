#' Volc.plot_2
#'
#' @description 
#' Another function for volcano plots, meant to be simpler than the main, awfully complex volcano.plot() function.
#' This has not been used for a long time so may prove unreliable.
#' Candidate for deletion.
#' 
#' @param data The dataframe that contains the values.
#' @param ids The name of the column of ID values to display.
#' @param x The name of the column with values to display on the X axis.
#' @param y The name of the column with values to display on the Y axis.
#' @param a The optional name of the column which will be used to calculate the Alpha values of dots. NULL by default.
#' @param N Optional. The number of replicates that were used to calculate P values (Y axis).
#' @param title The title of the plot.
#' @param v.lines Vertical threshold (positive value: this is assumed to be symmetrical).
#' @param h.lines Horizontal thresholds.
#' @param v.lines_min Minimal value for vertical threshold. NULL by default.
#' @param labels Should labels be displayed (for dots beyond thresholds)? Default = FALSE.
#' @param displ.limit Maximum number of labels to display. Default = 50. Set to NULL to turn off.
#' @param save FALSE by default. Set to format you want pictures saved as.
#' @param print TRUE by default. But can be set to FALSE because sometimes you just want to save a plot for later without printing it.
#' @param col_return FALSE by default. Should a column of significance be added to the output?
#' 
#' @examples
#' test <- Volc.plot_2(data = PG, ids = "Short ID", x = a1, y = a2,
#'                     N = length(grep(a1, colnames(PG), fixed = TRUE)), title = "Volcano plot",
#'                     v.lines = 2^v, h.lines = h)
#'

Volc.plot_2 <- function(data, ids, x, y, a = NULL, N, title, v.lines, h.lines, v.lines_min = NULL, labels = FALSE,
                        displ.limit = 50, save = FALSE, print = TRUE, col_return = FALSE) {
  stopifnot(v.lines > 0)
  if (is.null(a)) {
    tmp <- data[, c(ids, x, y)]; colnames(tmp) <- c("Names", "X", "Y"); Alpha <- FALSE
  } else {
    tmp <- data[, c(ids, x, y, a)]; colnames(tmp) <- c("Names", "X", "Y", "Alpha"); Alpha <- TRUE
  }
  tmp$Colour <- ""
  if (!is.null(v.lines)) {
    a <- v.lines
    if (!is.null(v.lines_min)) { a <- max(a, v.lines_min) }
    tmp$Colour[which(tmp$X > log2(a))] <- "Increased"
    tmp$Colour[which(tmp$X < -log2(a))] <- "Decreased" 
  }
  h.lines <- sort(unlist(h.lines), decreasing = FALSE)
  tmp$Colour[which(tmp$Y < -log10(max(h.lines)))] <- ""
  tmp$Names[which(tmp$Colour == "")] <- ""
  myColors <- c("red","blue","black")
  names(myColors) <- c("Increased","Decreased","")
  colScale <- ggplot2::scale_colour_manual(name = "colour", values = myColors)
  filter <- which(proteoCraft::is.all.good(tmp$Y, mode = "logical"))
  tmp <- tmp[filter,]
  xm <- max(abs(tmp$X))
  ym <- max(tmp$Y)
  test <- which(tmp$Names != "")
  if ((!is.null(displ.limit))&&(displ.limit > 0)) {
    if (length(test) > displ.limit) {
      a <- abs(tmp$X)
      a[which(tmp$Names == "")] <- 0
      b <- length(a) + 1 - rank(a, ties.method = "average")
      tmp$Names[which(b > displ.limit)] <- ""
    }
  }
  if (!Alpha) {
    plot1 <- ggplot2::ggplot(tmp) +
      ggplot2::geom_point(ggplot2::aes(x = X, y = Y, colour = Colour)) +
      ggplot2::theme_bw() +
      ggplot2::coord_fixed(ratio = 2*xm/ym, xlim = c(-xm, xm)*1.05, ylim = c(0, ym)*1.05,
                           expand = FALSE) +
      colScale +
      ggplot2::xlab("log2(Ratio)") + ggplot2::ylab("-log10(P value)") +
      ggplot2::ggtitle(title)
  } else {
    plot1 <- ggplot2::ggplot(tmp) +
      ggplot2::geom_point(ggplot2::aes(x = X, y = Y, colour = Colour, alpha = Alpha)) +
      ggplot2::theme_bw() +
      ggplot2::coord_fixed(ratio = 2*xm/ym, xlim = c(-xm, xm)*1.05, ylim = c(0, ym)*1.05,
                           expand = FALSE) +
      colScale +
      ggplot2::xlab("log2(Ratio)") + ggplot2::ylab("-log10(P value)") +
      ggplot2::ggtitle(title)
  }
  if (!missing(N)) {
    plot1 <- plot1 +
      ggplot2::geom_text(x = -xm, y = ym*1.025, label = paste0("N = ", N), hjust = 0)
  }
  if (!is.null(v.lines)) {
    plot1 <- plot1 +
      ggplot2::geom_vline(xintercept = -log2(v.lines), color = "blue") +
      ggplot2::geom_vline(xintercept = log2(v.lines), color = "red") +
      ggplot2::geom_text(x = -log2(v.lines)-0.02*xm, y = 2*ym/3,
                         label = paste0("log2(Ratio) = ", round(-log2(v.lines),2)),
                         cex = 3, hjust = 0, colour = "blue", angle = 90) +
      ggplot2::geom_text(x = log2(v.lines)+0.015*xm, y = 2*ym/3,
                         label = paste0("log2(Ratio) = ", round(log2(v.lines), 2)),
                         cex = 3, hjust = 0, colour = "red", angle = 90)
    if ((!is.null(v.lines_min)&&(v.lines < v.lines_min))) {
      plot1 <- plot1 +
        ggplot2::geom_vline(xintercept = -log2(v.lines_min), color = "darkblue") +
        ggplot2::geom_vline(xintercept = log2(v.lines_min), color = "darkred")
    }
  }
  h.lines <- proteoCraft::is.all.good(h.lines)
  if (length(h.lines) > 0) {
    for (i in 1:length(h.lines)) {
      plot1 <- plot1 +
        ggplot2::geom_hline(yintercept = -log10(h.lines[i]),
                            colour = "orange") +
        ggplot2::geom_text(x = -xm*0.95, y = -log10(h.lines[i])+0.01*ym,
                           label = names(h.lines[i]), cex = 3, hjust = 0,
                           colour = "orange")
    }
  }
  plot2 <- plot1 +
    ggrepel::geom_text_repel(ggplot2::aes(x = X, y = Y, label = Names, colour = Colour),
                             cex = 2.5)
  if (print) {
    if (labels) { proteoCraft::poplot(plot2) } else { proteoCraft::poplot(plot1) }
  }
  if ((length(save) > 1)||(save != FALSE)) {
    t <- gsub("/|:|\\*|\\?|<|>|\\|", "-", title)
    save <- unique(gsub("^jpg$", "jpeg", gsub("^\\.", "", tolower(save))))
    for (s in save) {
      if (s %in% c("jpeg", "tiff", "png", "bmp")) {
        ggplot2::ggsave(paste0(t, ".", s), plot1,
                        dpi = 300, width = 10, height = 10, units = "in")
        ggplot2::ggsave(paste0(t, "_labelled.", s), plot2,
                        dpi = 300, width = 10, height = 10, units = "in")
      } else {
        ggplot2::ggsave(paste0(t, ".", s), plot1)
        ggplot2::ggsave(paste0(t, "_labelled.", s), plot2)
      }
    }
  }
  plot1 <- ggplotify::as.ggplot(ggplotify::as.grob(plot1))
  plot2 <- ggplotify::as.ggplot(ggplotify::as.grob(plot2))
  res <- list("Volcano plot" = plot1,
              "Labelled volcano plot" = plot2)
  if (col_return) {
    col <- rep("", nrow(data))
    filter2 <- which(tmp$Colour != "")
    h <- -log10(h.lines)
    names(h) <- gsub("^Threshold-", "_", names(h))
    tmp2 <- apply(tmp[filter2, c("Y", "Colour")], 1, function(x) {
      x <- unlist(x)
      t <- which(as.numeric(h) <= as.numeric(x[1]))
      if (length(t) > 0) { r <- paste0(x[2], names(h)[min(t)]) } else { r <- "" }
      return(r)
    })
    col[filter[filter2]] <- tmp2
    res[["Significance"]] <- col
  }
  return(res)
}
