options(stringsAsFactors = FALSE)
require(proteoCraft)
require(rstudioapi)
require(reshape2)
require(magrittr)
require(ggplot2)
require(openxlsx)

wrkDr <- dirname(rstudioapi::getActiveDocumentContext()$path)
Sys.sleep(1) # Required because sometimes R will attempt to set the word directory too quickly and fail.
setwd(wrkDr)


homePath <- paste0(normalizePath(Sys.getenv("HOME"), winslash = "/"), "/R/proteoCraft")
dlft <- read.xlsx(paste0(homePath, "/Default_locations.xlsx"))

if ((exists("dir1"))&&(length(dir1))&&(nchar(dir1))&&(dir.exists(dir1))) { defltdir <- dir1 } else { defltdir <- dlft$Path[match("Search folder", dlft$Folder)] }
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


## 
#system(paste0("open \"", f1, "\""))
#system(paste0("open \"", f2, "\""))
i <- 1
fl1 <- readLines(f1); fl2 <- readLines(f2); fl1 <- gsub(pat, "", fl1); fl2 <- gsub(pat, "", fl2); w <- suppressWarnings(which(fl1 != fl2)); if (length(w) > 0) { cat(paste0("Line ", w[i], ":\n\n - file 1: ", fl1[w[i]], "\n\n - file 2: ", fl2[w[i]], "\n")) } else { system(paste0("open \"", f1, "\"")); system(paste0("open \"", f2, "\"")); print("Both scripts are identical!") }
#

g1 <- which(fl1 != ""); g2 <- which(fl2 != "")
l1 <- length(g1); l2 <- length(g2)
wa <- which(c(l1, l2) == max(c(l1, l2)))[1]
wb <- rev(which(c(l1, l2) == min(c(l1, l2))))[1]
fla <- get(paste0("fl", 1:2)[wa])
flb <- get(paste0("fl", 1:2)[wb])
fla <- fla[g1]
flb <- flb[g2]
test <- lapply(fla, function(x) { which(flb == x) })
test <- listMelt(test, NULL, c("Y", "X"))
Xlab <- basename(get(paste0("f", 1:2)[wa]))
Ylab <- basename(get(paste0("f", 1:2)[wb]))
if (Xlab == Ylab) {
  Xlab <- dirname(get(paste0("f", 1:2)[wa]))
  Ylab <- dirname(get(paste0("f", 1:2)[wb]))
}
plot <- ggplot(test) + geom_tile(aes(x = X, y = Y), fill = "white", height = 1, width = 1) +
  theme(panel.background = element_rect(fill = "black"),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        panel.spacing = unit(c(0, 0, 0, 0), "null"),
        axis.line = element_blank(),
        legend.position = "none",
        axis.ticks.length = unit(0, "null"),
        legend.margin = margin(0, 0, 0, 0, "cm")) +
  xlab(Xlab) + ylab(Ylab) +
  coord_fixed() +
  scale_x_continuous(expand = c(0, 0), limits = range(test$X)) +
  scale_y_continuous(expand = c(0, 0), limits = range(test$Y))
poplot(plot, 12, 20)
