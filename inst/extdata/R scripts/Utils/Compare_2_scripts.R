# Compare 2 versions of the same script to identify divergences and harmonize them 
options(stringsAsFactors = FALSE)
require(proteoCraft)
require(rstudioapi)
require(reshape2)
require(magrittr)
require(openxlsx)

wrkDr <- dirname(rstudioapi::getActiveDocumentContext()$path)
Sys.sleep(1L) # Required because sometimes R will attempt to set the word directory too quickly and fail.
setwd(wrkDr)
checkPlot <- FALSE

homePath <- paste0(normalizePath(Sys.getenv("HOME"), winslash = "/"), "/R/proteoCraft")
dlft <- read.xlsx(paste0(homePath, "/Default_locations.xlsx"))

defltdir <- if ((exists("dir1"))&&(length(dir1))&&(nchar(dir1))&&(dir.exists(dir1))) { dir1 } else {
  dlft$Path[match("Search folder", dlft$Folder)]
}
f1 <- rstudioapi::selectFile("Select first script", path = gsub("/+", "/", paste0(defltdir, "/*.R")))
dir1 <- gsub("/$", "", dirname(f1))
setwd(dir1)
f2 <- rstudioapi::selectFile("Select second script", path = paste0(dir1, "/*.R"))
dir2 <- gsub("/$", "", dirname(f2))
#setwd(dir2)
setwd(wrkDr)
cat(paste0("Last modified dates:\n - File 1: ", file.info(f1)$mtime, "\n - File 2: ", file.info(f2)$mtime, "\n"))

pat <- "^ *$| *$"   # Remove all useless lines
# If you want to remove comments:
pat <- paste0(pat, "| *#.*")

isOpen1 <- isOpen2 <- FALSE
checkPlotXprs <- expression({
  require(ggplot2)
  g1 <- which(fl1 != ""); g2 <- which(fl2 != "")
  l1 <- length(g1); l2 <- length(g2)
  wa <- which(c(l1, l2) == max(c(l1, l2)))[1L]
  wb <- rev(which(c(l1, l2) == min(c(l1, l2))))[1L]
  fla <- get(paste0("fl", 1L:2L)[wa])
  flb <- get(paste0("fl", 1L:2L)[wb])
  fla <- fla[g1]
  flb <- flb[g2]
  test <- lapply(fla, function(x) { which(flb == x) })
  test <- listMelt(test, NULL, c("Y", "X"))
  Xlab <- basename(get(paste0("f", 1L:2L)[wa]))
  Ylab <- basename(get(paste0("f", 1L:2L)[wb]))
  if (Xlab == Ylab) {
    Xlab <- dirname(get(paste0("f", 1L:2L)[wa]))
    Ylab <- dirname(get(paste0("f", 1L:2L)[wb]))
  }
  plot <- ggplot(test) + geom_tile(aes(x = X, y = Y), fill = "white", height = 1L, width = 1L) +
    theme(panel.background = element_rect(fill = "black"),
          plot.background = element_rect(fill = "transparent", colour = NA),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          plot.margin = unit(c(0L, 0L, 0L, 0L), "null"),
          panel.spacing = unit(c(0L, 0L, 0L, 0L), "null"),
          axis.line = element_blank(),
          legend.position = "none",
          axis.ticks.length = unit(0L, "null"),
          legend.margin = margin(0L, 0L, 0L, 0L, "cm")) +
    xlab(Xlab) + ylab(Ylab) +
    coord_fixed() +
    scale_x_continuous(expand = c(0L, 0L), limits = range(test$X)) +
    scale_y_continuous(expand = c(0L, 0L), limits = range(test$Y))
  poplot(plot, 12L, 20L)
})
checkScriptXprs <- expression({
  fl1 <- readLines(f1)
  fl2 <- readLines(f2)
  fl1 <- gsub(pat, "", fl1)
  fl2 <- gsub(pat, "", fl2)
  w1 <- which(fl1 != "")
  w2 <- which(fl2 != "")
  w <- suppressWarnings(which(fl1[w1] != fl2[w2]))
  if (length(w) >= i) {
    # Let's open both documents in this session
    #if (!isOpen1) {
      rstudioapi::documentOpen(f1, w1[w[i]])
      #isOpen1 <- TRUE
    #}
    #if (!isOpen2) {
      rstudioapi::documentOpen(f2, w2[w[i]])
      #isOpen2 <- TRUE
    #}
    cat(paste0("Discrepancy at line ", paste(unique(c(w1[w[i]], w2[w[i]])), collapse = "/"),
               ":\n\n - file 1: ", fl1[w1][w[i]], "\n\n - file 2: ", fl2[w2][w[i]], "\n"))
  } else {
    print("Both scripts are identical!")
  }
})

## 
#system(paste0("open \"", f1, "\""))
#system(paste0("open \"", f2, "\""))
i <- 1L
#
eval(checkScriptXprs)
#
if (checkPlot) { eval(checkPlotXprs) }
