#' cov3D
#'
#' @description
#' This function can be used to create a 3d coverage map of identified pepides overlaid over a protein model in PDB format.
#' 
#' @param pdb A PDB file, either already read into R or a path to an existing such file.\cr The first model in the file will be detected, parsed, and the collection of peptides of interest will be mapped onto the model.  
#' @param peptides The collection of peptide sequences to map on the protein sequence(s).
#' @param colscale Name of a color scale acceptable by plotly. The default is "viridis" if asRatios = FALSE and "plasma" if asRatios = TRUE
#' @param path Where to save our httml file.
#' @param intensities Values to map to the "colscale" argument so that peptides can be printed with different colors, e.g. mapped to abundance.
#' @param display Logical: should we open the plot? TRUE by default.
#' @param asRatios Is the data fed to the intensities argument actually ratios data (default = FALSE)? If TRUE, this changes the default color scale from "viridis" to "plasma"
#' 
#' @examples
#' cov3D(pdb, peptides)
#' 
#' @export

cov3D <- function(pdb,
                  peptides,
                  colscale,
                  path,
                  intensities = NULL,
                  display = TRUE,
                  asRatios = FALSE) {
  TESTING <- FALSE
  #aRmel::DefArg(aRmel::cov3D) ;TESTING = TRUE
  #pdb = pdbFl; peptides = tmpDat$Group.1; path = paste0(Par_dir, "/FLAG_KCC2_coverage (", pdbNm, ").html"); intensities = tmpDat$x
  #pdb = pdbFl; peptides = tmpDat$Group.1; path = paste0(Par_dir, "/FLAG_KCC2_coverage - ", nm, " (", pdbNm, ").html"); intensities = tmpDat$x
  #
  if (TESTING) {
    # Note:
    # This is not a perfect alternative to missing but will work in most cases, unless x matches a function imported by a package 
    misFun <- function(x) { return(!exists(deparse(substitute(x)))) }
  } else { misFun <- missing }
  #
  wdBckp <- getwd()
  #
  if ((misFun(asRatios))||(!is.logical(asRatios))||(is.na(asRatios))) {
    asRatios <- FALSE
  }
  if ((misFun(colscale))||(!is.character(colscale))) {
    colscale <- c("viridis", "plasma")[asRatios+1]
  }
  #
  stopifnot(!misFun(pdb),
            !misFun(peptides))
  #
  #pdb <- readLines(pdb)
  if ((length(pdb) == 1)&&(file.exists(pdb))) { pdb <- readLines(pdb) }
  model <- pdb[grep("^ATOM ", pdb)[1]:length(pdb)]
  #nm <- gsub("^TITLE +| +$", "", grep("^TITLE +", pdb, value = TRUE)[1])
  no <- grep("^ATOM ", model, invert = TRUE)
  if (length(no)) { model <- model[1:(no[1]-1)] }
  tmp <- strsplit(gsub("^ATOM +", "", model), " +")
  tst <- sapply(tmp, length)
  tst <- aggregate(1:length(tst), list(tst), list)
  colnames(tst) <- c("Entries", "Rows")
  tst$nRows <- sapply(tst$Rows, length)
  if (nrow(tst) > 1) {
    n <- sum(tst$nRows[which(tst$Entries != 11)])
    i <- (n > 1)+1
    warning(paste0("Parsing method 1 failed: ", n, " model row", c("", "s")[i], " do", c("es", "")[i], " not contain the expected number (11) of entries!\nTrying parsing method 2..."))
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
  dat2 <- dat[which(dat$`Branch indicator` == "CA"),] # CA = carbon alpha it seems?
  #
  # Compute coverage
  #peptides <- aggregate(ev$Intensity, list(ev$"Modified sequence"), sum, na.rm = TRUE)
  #intensities <- peptides$x
  #peptides <- peptides$Group.1
  aa321 <- data.frame(aa1 = as.character(aRmel::AA_table$AA),
                      aa3 = c("ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU", "MET", "ASN", "PYL", "PRO", "GLN", "ARG", "SER", "THR", "SEC", "VAL", "TRP", "TYR"))
  dat2$AA <- aa321$aa1[match(dat2$Residue, aa321$aa3)]
  dat2$I2L <- gsub("I", "L", dat2$AA)
  if ((misFun(intensities))||(is.null(intensities))) { intensities <- rep(1, length(peptides)) }
  pepTbl <- data.frame(seq = peptides,
                       int = intensities,
                       I2Lpep = gsub("^_|_$", "", gsub("I", "L", peptides)))
  pepTbl$I2Lpep <- setNames(aRmel::annot_to_tabl(pepTbl$I2Lpep), peptides)
  pepTbl$matches <- setNames(lapply(pepTbl$I2Lpep, function(x) {
    l <- nrow(x)
    w <- sapply(1:l, function(y) {
      which(dat2$I2L == x$Sequence[y]) + 1 - y
    })
    w <- unlist(w)
    w <- aggregate(w, list(w), length)
    return(w$Group.1[which(w$x == l)])
  }), peptides)
  dat2$Intensity <- 0
  pepTbl$L <- sapply(pepTbl$I2Lpep, nrow)
  for (i in 1:nrow(pepTbl)) {#i <- 2
    m <- pepTbl$matches[[i]]
    if (length(m)) {
      rg <- m:(m + pepTbl$L[i] - 1)
      dat2$Intensity[rg] <- dat2$Intensity[rg] +  pepTbl$int[i]
    }
  }
  L <- nrow(dat2)
  # Label N/C-termini
  dat2$AA[1] <- paste0("N-term. ", dat2$AA[1])
  dat2$AA[L] <- paste0("C-term. ", dat2$AA[L])
  dat2$Type <- "AA"
  #
  # Now create the skeleton (twice more rows -1, because we change color mid-peptide bond)
  rg <- 1:(L-1)
  dat3 <- dat2[rg+1,]
  dat3$`Residue sequence number` <- dat3$`Residue sequence number`-0.5
  dat3$X <- (dat2$X[rg]+dat2$X[rg+1])/2
  dat3$Y <- (dat2$Y[rg]+dat2$Y[rg+1])/2
  dat3$Z <- (dat2$Z[rg]+dat2$Z[rg+1])/2
  dat3$Residue <- ""
  dat3 <- rbind(dat2, dat3)
  dat3 <- dat3[order(dat3$`Residue sequence number`),]
  #
  modsTst <- sapply(pepTbl$I2Lpep, function(x) {
    x <- x$Annotations
    return(which(x != ""))
  })
  w <- which((sapply(modsTst, length) > 0)&(sapply(pepTbl$matches, length) > 0))
  if (length(w)) {
    modsTst <- setNames(lapply(w, function(x) {
      y <- list(Tbl = pepTbl$I2Lpep[[x]],
                Match = pepTbl$matches[[x]])
      y$Tbl$Pos <- 1:nrow(y$Tbl)
      y$Tbl <- y$Tbl[modsTst[[x]], , drop = FALSE]
      return(y)
    }), peptides[w])
    modsTst <- setNames(lapply(modsTst, function(x) { #x <- modsTst[[1]]
      y <- lapply(x$Match, function(y) {
        z <- x$Tbl
        z$Pos <- z$Pos + y - 1
        #return(z[, c("Annotations", "Pos"), drop = FALSE])
        return(z)
      })
      y <- do.call(rbind, y)
      return(y)
    }), names(modsTst))
    modsTst <- do.call(rbind, modsTst)
    modsTst <- aggregate(modsTst$Annotations, list(modsTst$Pos, modsTst$Sequence), unique)
    modsTst <- modsTst[order(modsTst$x),]
    modsTst$Pos <- do.call(paste, c(modsTst[, c("Group.2", "Group.1")], sep = ""))
    modsTst <- aggregate(modsTst$x, list(modsTst$Pos), function(x) { paste(unique(x), collapse = " / ") })
    modsTst$Pos <- as.integer(gsub("^[A-Z]", "", modsTst$Group.1))
    modsTst[, c("Residue", "AA", "X", "Y", "Z")] <- dat2[modsTst$Pos, c("Residue", "AA", "X", "Y", "Z")]
    modsTst$PTM <- do.call(paste, c(modsTst[, c("AA", "x")], sep = " "))
    modsTst$Type <- "mod. AA"
    wN <- which(!1:L %in% modsTst$Pos)
    dat2 <- dat2[wN,] # We will only plot normal amino acids here if they do not figure already in the modified amino acid object (which is plotted differently)
  }
  #
  symb <- setNames(c("circle", "diamond"), c("AA", "mod. AA"))
  #
  lst2 <- list(size = 6,
               sizemode = "area",
               colorbar = list(xanchor = "left",
                               title = "Intensity"), 
               showscale = FALSE)
  lst3 <- list(width = 6,
               showscale = FALSE)
  #
  w0 <- which((is.na(dat2$Intensity))|(dat2$Intensity == 0))
  dat2$Intensity[w0] <- NA
  w0 <- which((is.na(dat3$Intensity))|(dat3$Intensity == 0))
  dat3$Intensity[w0] <- NA
  my3dplotly <- plotly::plot_ly(symbols = symb)
  my3dplotly <- plotly::add_trace(my3dplotly, data = dat3, x = ~X, y = ~Y, z = ~Z, type = "scatter3d", mode = "lines",
                                  color = ~Intensity, colors = colscale, line = lst3 <- list(width = 6, showscale = FALSE),
                                  opacity = 1, hoverinfo = "none", inherit = FALSE)
  #my3dplotly
  my3dplotly <- plotly::add_trace(my3dplotly, data = dat2, x = ~X, y = ~Y, z = ~Z, type = "scatter3d", mode = "markers",
                                  color = ~Intensity, colors = colscale,
                                  text = ~AA, opacity = 1, hoverinfo = "text", marker = lst2,
                                  symbol = ~Type, inherit = FALSE)
  #my3dplotly
  if (length(w)) {
    lst4 <- list(width = 8,
                 reverscale = FALSE,
                 size = 8,
                 sizemode = "area",
                 showscale = FALSE,
                 colorbar = list(title = "PTM",
                                 orientation = "h",
                                 xanchor = "center",
                                 yanchor = "bottom",
                                 x = 0))
    my3dplotly <- plotly::add_trace(my3dplotly, data = modsTst, x = ~X, y = ~Y, z = ~Z, type = "scatter3d", mode = "markers",
                                    color = ~PTM, text = ~PTM, opacity = 1, hoverinfo = "text", marker = lst4,
                                    symbol = ~Type, inherit = FALSE)
    #my3dplotly
  }
  #
  #path <- paste0(wdBckp, "test.html")
  if (!misFun(path)) {
    dir <- dirname(path)
    if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
    setwd(dir) # Somehow the self-contained argument of saveWidget only works if the current work directory is the same as the location of the destination html 
    htmlwidgets::saveWidget(my3dplotly, path, selfcontained = TRUE)
    setwd(wdBckp)
    if (display) { system(paste0("open \"", path, "\"")) }
  } else {
    if (display) { print(my3dplotly) }
  }
}
