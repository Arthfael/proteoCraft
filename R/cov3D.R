#' cov3D
#'
#' @description
#' This function can be used to create a 3d coverage map of identified pepides overlaid over a protein model in PDB format.
#' 
#' @param pdb A PDB file, either already read into R or a path to an existing such file.\cr The first model in the file will be detected, parsed, and the collection of peptides of interest will be mapped onto the model.  
#' @param peptides The collection of peptide sequences to map on the protein sequence(s).
#' @param colscale Name of a color scale acceptable by plotly. The default is "viridis" if asRatios = FALSE and "plasma" if asRatios = TRUE
#' @param path Where to save our html file.
#' @param intensities Values to map to the "colscale" argument so that peptides can be printed with different colors, e.g. mapped to abundance.
#' @param display Logical: should we open the plot? TRUE by default.
#' @param asRatios Is the data fed to the intensities argument actually ratios data (default = FALSE)? If TRUE, this changes the default color scale from "viridis" to "plasma"
#' @param I_eq_L Should we consider I and L identical? Currently, by default, TRUE for both DIA and DDA: see https://github.com/vdemichev/DiaNN/discussions/1631
#' @param ttl Title of the plot.
#'
#' @returns
#' This function does not return anything.
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
                  asRatios = FALSE,
                  I_eq_L = TRUE,
                  ttl = NULL) {
  TESTING <- FALSE
  #DefArg(cov3D) ;TESTING = TRUE
  #pdb = pdbFl; peptides = tmpDat$Group.1; path = paste0(Par_dir, "/FLAG_KCC2_coverage (", pdbNm, ").html"); intensities = tmpDat$x
  #pdb = pdbFl; peptides = tmpDat$Group.1; path = paste0(Par_dir, "/FLAG_KCC2_coverage - ", nm, " (", pdbNm, ").html"); intensities = tmpDat$x
  #pdb = pdbFl; peptides = tmpDat$Group.1; path = paste0(Par_dir, "/FLAG_KCC2_coverage - ", nm, " (", pdbNm, ").html"); intensities = tmpDat$x
  #pdb = fl; peptides = seq1; path = pth; ttl = nm; intensities = int1[[x]]; display = FALSE
  #pdb = fl; peptides = seq1; path = pth; ttl = nm; intensities = intVect[[x]][grs]; display = FALSE
  misFun <- if (TESTING) {
    # Note:
    # This is not a perfect alternative to missing but will work in most cases, unless x matches a function imported by a package 
    \(x) { return(!exists(deparse(substitute(x)))) }
  } else { missing }
  #
  wdBckp <- getwd()
  #
  if ((misFun(I_eq_L))||(!is.logical(I_eq_L))||(is.na(I_eq_L))) {
    I_eq_L <- if ((exists("isDIA"))&&(is.logical(isDIA))&&(!is.na(isDIA))) {
      !isDIA
    } else { TRUE }
  }
  if ((misFun(asRatios))||(!is.logical(asRatios))||(is.na(asRatios))) {
    asRatios <- FALSE
  }
  if ((misFun(colscale))||(!is.character(colscale))) {
    colscale <- c("viridis", "plasma")[asRatios+1L]
  }
  #
  stopifnot(!misFun(pdb),
            !misFun(peptides))
  #
  #pdb <- readLines(pdb)
  if ((length(pdb) == 1L)&&(file.exists(pdb))) { pdb <- readLines(pdb) }
  model <- pdb[grep("^ATOM ", pdb)[1L]:length(pdb)]
  #nm <- gsub("^TITLE +| +$", "", grep("^TITLE +", pdb, value = TRUE)[1L])
  no <- grep("^ATOM ", model, invert = TRUE)
  if (length(no)) { model <- model[1L:(no[1L]-1L)] }
  # tmp <- strsplit(gsub("^ATOM +", "", model), " +")
  # tst <- lengths(tmp)
  # tst <- aggregate(1L:length(tst), list(tst), list)
  # colnames(tst) <- c("Entries", "Rows")
  # tst$nRows <- lengths(tst$Rows)
  #if (tst$nRows > 0L) {
    # n <- sum(tst$nRows[which(tst$Entries != 11L)])
    # i <- (n > 1L)+1L
    # warning(paste0("Parsing method 1 failed: ", n, " model row", c("", "s")[i], " do", c("es", "")[i], " not contain the expected number (11) of entries!\nTrying parsing method 2..."))
    dat <- data.frame(V1 = gsub("^ +| +$", "", substr(model, 5L, 11L)),
                      V2 = gsub("^ +| +$", "", substr(model, 12L, 17L)),
                      V3 = gsub("^ +| +$", "", substr(model, 18L, 20L)),
                      V4 = gsub("^ +| +$", "", substr(model, 21L, 22L)),
                      V5 = gsub("^ +| +$", "", substr(model, 23L, 26L)),
                      V6 = gsub("^ +| +$", "", substr(model, 27L, 38L)),
                      V7 = gsub("^ +| +$", "", substr(model, 39L, 46L)),
                      V8 = gsub("^ +| +$", "", substr(model, 47L, 54L)),
                      V9 = gsub("^ +| +$", "", substr(model, 55L, 60L)),
                      V10 = gsub("^ +| +$", "", substr(model, 61L, 66L)),
                      V11 = gsub("^ +| +$", "", substr(model, 67L, nchar(model))))
  # } else {
  #    dat <- as.data.frame(t(sapply(tst, unlist)))
  # }
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
  aa321 <- data.frame(aa1 = as.character(AA_table$AA),
                      aa3 = c("ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU", "MET", "ASN", "PYL", "PRO", "GLN", "ARG", "SER", "THR", "SEC", "VAL", "TRP", "TYR"))
  dat2$AA <- aa321$aa1[match(dat2$Residue, aa321$aa3)]
  dat2$I2L <- dat2$AA
  if (I_eq_L) {
    dat2$I2L <- gsub("I", "L", dat2$AA)
  }
  if ((misFun(intensities))||(is.null(intensities))||(length(intensities) != length(peptides))) {
    intensities <- rep(1, length(peptides))
  }
  pepTbl <- data.frame(seq = peptides,
                       int = intensities,
                       I2Lpep = gsub("^_|_$", "", peptides))
  if (I_eq_L) {
    pepTbl$I2Lpep <- gsub("I", "L", pepTbl$I2Lpep)
  }
  pepTbl$I2Lpep <- setNames(annot_to_tabl(pepTbl$I2Lpep), peptides)
  pepTbl$I2Lpep <- setNames(lapply(pepTbl$I2Lpep, \(x) { x[which(x$Sequence != "_"),] }), peptides)
  pepTbl$matches <- setNames(lapply(pepTbl$I2Lpep, \(x) { #x <- pepTbl$I2Lpep[[1L]]
    rg <- 1L:nrow(x)
    l <- length(rg)
    w <- sapply(rg, \(y) {
      which(dat2$I2L == x$Sequence[y]) + 1L - y
    })
    w <- unlist(w)
    w <- aggregate(w, list(w), length)
    return(w$Group.1[which(w$x == l)])
  }), peptides)
  # To check matches:
  # pepTbl$tst <- vapply(1:nrow(pepTbl), \(i) { #i <- 1L
  #   l <- nrow(pepTbl$I2Lpep[[i]])
  #   unique(unlist(lapply(pepTbl$matches[[i]], \(m) { paste(dat2$I2L[m:(m+l-1)], collapse = "") })))
  # }, "")
  #View(pepTbl[, c("seq", "tst")])
  dat2$Intensity <- 0
  pepTbl$L <- vapply(pepTbl$I2Lpep, nrow, 1L)
  for (i in 1L:nrow(pepTbl)) { #i <- 2L
    m <- pepTbl$matches[[i]]
    if (length(m)) {
      rg <- m:(m + pepTbl$L[i] - 1L)
      dat2$Intensity[rg] <- dat2$Intensity[rg] +  pepTbl$int[i]
    }
  }
  L <- nrow(dat2)
  # Label N/C-termini
  dat2$AA[1L] <- paste0("N-term. ", dat2$AA[1L])
  dat2$AA[L] <- paste0("C-term. ", dat2$AA[L])
  dat2$Type <- "AA"
  #
  # Now create the skeleton (twice more rows -1, because we change color mid-peptide bond)
  rg <- 1L:(L-1L)
  rg_1 <- rg+1L
  dat3 <- dat2[rg_1,]
  dat3$`Residue sequence number` <- dat3$`Residue sequence number`-0.5
  dat3$X <- (dat2$X[rg]+dat2$X[rg_1])/2
  dat3$Y <- (dat2$Y[rg]+dat2$Y[rg_1])/2
  dat3$Z <- (dat2$Z[rg]+dat2$Z[rg_1])/2
  dat3$Residue <- ""
  dat3 <- rbind(dat2, dat3)
  dat3 <- dat3[order(dat3$`Residue sequence number`),]
  #
  modsTst <- sapply(pepTbl$I2Lpep, \(x) {
    x <- x$Annotations
    return(which(x != ""))
  })
  w <- which((lengths(modsTst) > 0L)&(lengths(pepTbl$matches) > 0L))
  if (length(w)) {
    modsTst <- setNames(lapply(w, \(x) {
      y <- list(Tbl = pepTbl$I2Lpep[[x]],
                Match = pepTbl$matches[[x]])
      y$Tbl$Pos <- 1L:nrow(y$Tbl)
      y$Tbl <- y$Tbl[modsTst[[x]], , drop = FALSE]
      return(y)
    }), peptides[w])
    modsTst <- setNames(lapply(modsTst, \(x) { #x <- modsTst[[1L]]
      y <- lapply(x$Match, \(y) {
        z <- x$Tbl
        z$Pos <- z$Pos + y - 1L
        #return(z[, c("Annotations", "Pos"), drop = FALSE])
        return(z)
      })
      y <- do.call(rbind, y)
      return(y)
    }), names(modsTst))
    modsTst <- do.call(rbind, modsTst)
    modsTst <- aggregate(modsTst$Annotations, list(modsTst$Pos, modsTst$Sequence), unique)
    modsTst2 <- listMelt(modsTst$x, 1L:nrow(modsTst), c("x", "Row"))
    modsTst2[, c("Group.1", "Group.2")] <- modsTst[modsTst2$Row, c("Group.1", "Group.2")]
    modsTst <- modsTst2
    modsTst <- modsTst[order(modsTst$x),]
    modsTst$Pos <- do.call(paste, c(modsTst[, c("Group.2", "Group.1")], sep = ""))
    modsTst <- aggregate(modsTst$x, list(modsTst$Pos), \(x) { paste(unique(x), collapse = " / ") })
    modsTst$Pos <- as.integer(gsub("^[A-Z]", "", modsTst$Group.1))
    modsTst[, c("Residue", "AA", "X", "Y", "Z")] <- dat2[modsTst$Pos, c("Residue", "AA", "X", "Y", "Z")]
    modsTst$PTM <- do.call(paste, c(modsTst[, c("AA", "x")], sep = " "))
    modsTst$Type <- "mod. AA"
    wN <- which(!1L:L %in% modsTst$Pos)
    dat2 <- dat2[wN,] # We will only plot normal amino acids here if they do not figure already in the modified amino acid object (which is plotted differently)
  }
  #
  symb <- setNames(c("circle", "diamond"), c("AA", "mod. AA"))
  #
  lst2 <- list(size = 6L,
               sizemode = "area",
               colorbar = list(xanchor = "left",
                               title = "Intensity"), 
               showscale = FALSE)
  lst3 <- list(width = 6L,
               showscale = FALSE)
  #
  w0 <- which((is.na(dat2$Intensity))|(dat2$Intensity == 0))
  dat2$Intensity[w0] <- NA
  dat2$AA <- do.call(paste, c(dat2[, c("AA", "Residue sequence number")], sep = " "))
  w0 <- which((is.na(dat3$Intensity))|(dat3$Intensity == 0))
  dat3$Intensity[w0] <- NA
  my3dplotly <- plotly::plot_ly(symbols = symb)
  my3dplotly <- plotly::add_trace(my3dplotly, data = dat3, x = ~X, y = ~Y, z = ~Z, type = "scatter3d", mode = "lines",
                                  color = ~Intensity, colors = colscale, line = lst3 <- list(width = 6L, showscale = FALSE),
                                  opacity = 1L, hoverinfo = "none", inherit = FALSE)
  #my3dplotly
  my3dplotly <- plotly::add_trace(my3dplotly, data = dat2, x = ~X, y = ~Y, z = ~Z, type = "scatter3d", mode = "markers",
                                  color = ~Intensity, colors = colscale,
                                  text = ~AA, opacity = 1L, hoverinfo = "text", marker = lst2,
                                  symbol = ~Type, inherit = FALSE, showlegend = FALSE)
  #my3dplotly
  if (length(w)) {
    lst4 <- list(width = 8L,
                 reverscale = FALSE,
                 size = 8L,
                 sizemode = "area",
                 showscale = FALSE,
                 colorbar = list(title = "PTM",
                                 orientation = "h",
                                 xanchor = "center",
                                 yanchor = "bottom",
                                 x = 0))
    my3dplotly <- plotly::add_trace(my3dplotly, data = modsTst, x = ~X, y = ~Y, z = ~Z,
                                    type = "scatter3d", mode = "markers",
                                    color = ~PTM, text = ~PTM, opacity = 1L, hoverinfo = "text",
                                    marker = lst4, symbol = ~Type, inherit = FALSE, showlegend = FALSE)
    #my3dplotly
  }
  if ((!misFun(ttl))&&(is.character(ttl))&&(nchar(ttl))) {
    my3dplotly <- plotly::layout(my3dplotly,
                                 title = list(text = ttl, x = 0.5, y = 0.95, xanchor = "center", yanchor = "bottom"))
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
  return()
}
