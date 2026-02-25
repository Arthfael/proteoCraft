# PDB files
dir <- "...My_Dataset"
fls <- list.files(dir, "\\.pdb$", full.names = TRUE)

for (fl in fls) { #fl <- fls[1]
  txt <- readLines(fl)
  model <- txt[grep("^ATOM ", txt)[1]:length(txt)]
  no <- grep("^ATOM ", model, invert = TRUE)
  if (length(no)) {
    model <- model[1:(no[1]-1)]
  }
  tmp <- strsplit(gsub("^ATOM +", "", model), " +")
  tst <- sapply(tmp, length)
  tst <- aggregate(1:length(tst), list(tst), list)
  colnames(tst) <- c("Entries", "Rows")
  tst$nRows <- sapply(tst$Rows, length)
  if (nrow(tst) > 1) {
    n <- sum(tst$nRows[which(tst$Entries != 11)])
    i <- (n > 1)+1
    warning(paste0("Invalid PDB file, ", n, " model row", c("", "s")[i], " do", c("es", "")[i], " not contain the expected number (11) of entries!\n(See \"tst\" object for row ind", c("ex", "ices")[i], "...)\nTrying alternate parsing method..."))
    dat <- data.frame(V1 = gsub("^ +| +$", "", substr(model, 5, 11)),
                      V2 = gsub("^ +| +$", "", substr(model, 12, 17)),
                      V3 = gsub("^ +| +$", "", substr(model, 18, 20)),
                      V4 = gsub("^ +| +$", "", substr(model, 21, 22)),
                      V5 = gsub("^ +| +$", "", substr(model, 23, 26)),
                      V6 = gsub("^ +| +$", "", substr(model, 27, 38)),
                      V7 = gsub("^ +| +$", "", substr(model, 39, 46)),
                      V8 = gsub("^ +| +$", "", substr(model, 47, 54)),
                      V9 = gsub("^ +| +$", "", substr(model, 55, 60)),
                      V10 = gsub("^ +| +$", "", substr(model, 61, 66)),
                      V11 = gsub("^ +| +$", "", substr(model, 67, nchar(model))))
  } else {
    dat <- as.data.frame(t(sapply(tst, unlist)))
  }
  colnames(dat) <- c("Atom serial number", "Branch indicator", "Residue", "Chain identifier", "Residue sequence number",
                     "X", "Y", "Z", "Occupancy", "Temperature factor (B-factor)", "Element symbol")
  dat$X <- as.numeric(dat$X)
  dat$Y <- as.numeric(dat$Y)
  dat$Z <- as.numeric(dat$Z)
  dat$`Residue sequence number` <- as.integer(dat$`Residue sequence number`)
  dat$`Temperature factor (B-factor)` <- as.factor(dat$`Temperature factor (B-factor)`)
  dat2 <- dat[which(dat$`Branch indicator` == "CA"),] # CA = carbon alpha I guess?
  L <- nrow(dat2)
  rg <- 1:(L-1)
  dat3 <- dat2[rg+1,]
  dat3$`Residue sequence number` <- dat3$`Residue sequence number`-0.5
  dat3$X <- (dat2$X[rg]+dat2$X[rg+1])/2
  dat3$Y <- (dat2$Y[rg]+dat2$Y[rg+1])/2
  dat3$Z <- (dat2$Z[rg]+dat2$Z[rg+1])/2
  dat3$Residue <- ""
  dat3 <- rbind(dat2, dat3)
  dat3 <- dat3[order(dat3$`Residue sequence number`),]
  require(plotly)
  require(htmlwidgets)
  my3dplotly <- plot_ly(dat3, x = ~X, y = ~Y, z = ~Z, type = "scatter3d", mode = "lines",
                        text = ~Residue, opacity = 1, hoverinfo = "text",
                        line = list(width = 6, color = ~`Temperature factor (B-factor)`, reverscale = FALSE))
  wdBckp <- getwd()
  setwd(dir) # Somehow the self-contained argument of saveWidget only works if the current work directory is the same as the location of the destination html 
  saveWidget(my3dplotly,
             gsub("\\.pdb$", ".html", fl),
             selfcontained = TRUE)
  setwd(wdBckp)
}
