#' multiplot
#' 
#' @description 
#' I did not actually create this function but nicked it from the R cookbook, so credit goes to them:
#' http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)
#' Currently unused, thus candidate for deletion.
#' 
#' @param plotlist The plots to plot. Duh!
#' @param cols Number of columns. The number of rows will be inferred from this.
#' @param layout A layout matrix. Usually assumed to be null so created by the function.
#' @param save FALSE by default. Set to format you want pictures saved as.
#' @param filename Required only if save = TRUE. Extensions will be ignored as they are provided by the save argument.
#' 
#' @export

multiplot <- function(..., plotlist = NULL, cols = 1, layout = NULL, save = FALSE, filename = NULL) {
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots <- length(plots)
  rows <- ceiling(numPlots/cols)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)), ncol = cols, nrow = rows)
  }
  if (numPlots == 1) {
    if (.Platform$OS.type == "windows") { windows(width = width, height = height) }
    if (.Platform$OS.type == "unix") { x11(width = width, height = height) }
    print(plots[[1]])
    if ((length(save) > 1)||(save != FALSE)) {
      filename <- gsub("/|:|\\*|\\?|<|>|\\|", "-", filename)
      save <- gsub("^\\.", "", save)
      save <- unique(gsub("jpg", "jpeg", gsub("^\\.", "", tolower(save))))
      for (s in save) {
        f <- paste(filename, ".", s, sep = "")
        dev <- get(s)
        if (s %in% c("tiff", "png", "bmp")) {
          dev(filename = f, units = "in", width = cols*10, height = rows*10, res = 300)
        } else {
          if (s == "jpeg") {
            dev(filename = f, units = "in", width = cols*10, height = rows*10, res = 300, compression = "none")
          } else { dev(filename = f) }
        }
        grid::grid.newpage()
        print(plots[[1]])
        dev.off()
      }
    }
  } else {
    # Set up the page
    if (.Platform$OS.type == "windows") { windows(width = width, height = height) }
    if (.Platform$OS.type == "unix") { x11(width = width, height = height) }
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow(layout),
                                                                 ncol(layout))))
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = grid::viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
    }
    if ((length(save) > 1)||(save != FALSE)) {
      filename <- gsub(":|\\*|\\?|<|>|\\||/", "-", filename)
      save <- gsub("^\\.", "", save)
      save <- unique(gsub("jpg", "jpeg", gsub("^\\.", "", tolower(save))))
      for (s in save) {
        f <- paste(filename, ".", s, sep = "")
        dev <- get(s)
        if (s %in% c("tiff", "png", "bmp")) {
          dev(filename = f, units = "in", width = cols*10, height = rows*10, res = 300)
        } else {
          if (s == "jpeg") {
            dev(filename = f, units = "in", width = cols*10, height = rows*10, res = 300, quality = 100)
          } else { dev(filename = f) }
        }
        grid::grid.newpage()
        grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow(layout),
                                                                     ncol(layout))))
        # Make each plot, in the correct location
        for (i in 1:numPlots) {
          # Get the i,j matrix positions of the regions that contain this subplot
          matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
          print(plots[[i]],
                vp = grid::viewport(layout.pos.row = matchidx$row,
                                    layout.pos.col = matchidx$col))
        }
        dev.off()
      }
    }
  }
}


