# Save plotly plots
cat(" -> Saving plotly plots...\n")
if (!dir.exists(dr)) {
  dr2 <- paste0(wd, "/", dr)
  if (dir.exists(dr2)) {
    # At this stage if it doesn't exist then it was a sub-directory
    dr <- dr2;rm(dr2)
  }
} 
# Courtesy of chatGPT:
#  - Pre-extract compact JSON specs
l <- length(myPlotLys)
if (l) {
  tmpFls <- paste0(dr, "/tmp", 1:l, ".RDS")
  invisible(lapply(1:l, function(i) { #i <- 1
    x <- myPlotLys[[i]]
    if ("Render" %in% names(x)) { # Solution for when a plotly is buggy: pre-render it before exporting
      pl <- x$Render
    } else {
      pl <- plotly::plotly_build(x$Plot)
    }
    # keep as plain JSON string (compact, cheap to ship)
    pl <- list(data = pl$x$data,
               layout = pl$x$layout)
    #format(object.size(b), "MB")
    readr::write_rds(pl, tmpFls[i])
    return()
  }))
  plot_Ttls <- vapply(myPlotLys, function(x) { x$Ttl }, "")
  plot_Paths <- paste0(gsub("[^A-Za-z0-9\\. ]", "_", plot_Ttls), ".html")
  #
  invisible(parallel::clusterCall(parClust, function() {
    library(plotly)
    library(htmlwidgets)
    library(jsonlite)
  }))
  parallel::clusterExport(parClust, list("dr", "plot_Paths", "tmpFls"), envir = environment())
  save_widget_from_def <- function(i) {
    curDir <- getwd() 
    setwd(dr)
    def <- readr::read_rds(tmpFls[i])
    w <- plotly::plot_ly()
    w$x$data   <- def$data
    w$x$layout <- def$layout
    pth <- file.path(dr, plot_Paths[i])
    htmlwidgets::saveWidget(w, pth, selfcontained = TRUE)
    setwd(curDir)
    return(pth)
  }
  tst <- parallel::parLapply(parClust, 1:l, save_widget_from_def)
  unlink(tmpFls)
  setwd(wd)
  # - Done!
}
