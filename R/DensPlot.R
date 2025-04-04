#' DensPlot
#'
#' @description
#' Just a small wrapper function for density plots, because I am tired of always typing the same stuff.
#' 
#' @param data The dataframe that contains the values.
#' @param X The name of the column of values whose density distribution is to be displayed.
#' @param Series The name of the column that maps values to series.
#' @param Col The series' colours. If there are more than one, must be a vector of colour names, named as series names, of the same length as the number of series. If unprovided default ggplot2 colours will be used instead.
#' @param Stat One of "count" or "density" (density).
#' @param log Should we log or not? If TRUE, the base is taken to be 10 by default; else you can instead specify the log.
#' @param print TRUE by default. But can be set to FALSE because sometimes you just want to save a plot for later without printing it.
#' @param FUN The optional function whose curve should be overlaid on the plot.
#' @param FUN.color The color of the curve of the optional function above. Default = "red".
#' @param title The title of the plot. Default = "Density plot"
#' @param save FALSE by default. Set to format you want pictures saved as.
#' @param return Should we return the plot? Default = FALSE.
#' 
#' @examples
#' Kolours <- colorRampPalette(c("blue", "red"))(length(unique(ev$Raw.file)))
#' names(Kolours) <- unique(ev$Raw.file)
#' DensPlot(ev, "Intensity", "Raw.file", Col = Kolours)
#' 
#' @export

DensPlot <- function(data,
                     X,
                     Series = NULL,
                     Col = "black",
                     Stat = "density",
                     log = FALSE,
                     print = TRUE,
                     FUN = NULL,
                     FUN.color = "red",
                     title = "Density plot",
                     save = FALSE,
                     return = FALSE) {
  if (is.logical(log)&&log) {
    log <- 10
  }
  if (!is.logical(log)) {
    X1 <- paste0("log", log, X)
    data[[X1]] <- log(data[[X]], log)
  } else {
    X1 <- X
  }
  data$X <- data[[X1]]
  if (missing("Series")) {
    plot <- ggplot2::ggplot(data) + ggplot2::geom_density(stat = Stat, ggplot2::aes(x = X), colour = Col) +
      ggplot2::xlab(X1)
  } else {
    data$Y <- data[[Series]]
    plot <- ggplot2::ggplot(data) + ggplot2::geom_density(stat = Stat, ggplot2::aes(x = X, colour = Y)) +
      ggplot2::xlab(X1)
    if (!missing(Col)) {
      if (length(Col) != length(unique(Col))) {
        stop("The number of series and of series colours differ!")
      } else {
        plot <- plot + ggplot2::scale_colour_manual(name = "colour", values = Col)
      }
    }
  }
  plot <- plot + ggplot2::ggtitle(title)
  if (!missing(FUN)) { plot <- plot + ggplot2::stat_function(fun = FUN, colour = FUN.color) }
  if (print) { proteoCraft::poplot(plot) }
  if ((length(save) > 1)||(save != FALSE)) {
    t <- gsub("/|:|\\*|\\?|<|>|\\|", "-", title)
    save <- unique(gsub("^jpg$", "jpeg", gsub("^\\.", "", tolower(save))))
    for (s in save) {
      if (s %in% c("jpeg", "tiff", "png", "bmp")) {
        ggplot2::ggsave(filename = paste0(t, ".", s),
                        plot = plot, dpi = 300, width = 10, height = 10, units = "in")
      } else {
        ggplot2::ggsave(filename = paste0(t, ".", s), plot = plot)
      }
    }
  }
  if (return) {
    plot <- ggplotify::as.ggplot(ggplotify::as.grob(plot))
    return(plot)
  }
}
