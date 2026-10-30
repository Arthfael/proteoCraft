
################################################################
#                                                              #
#                       Thermo raw files                       #
#      Get MS2 scans around peak apex for targeted method      #
#                                                              #
################################################################

# Packages
requireNamespace("proteoCraft")

# Library paths
RPath <- as.data.frame(library()$results)
RPath <- normalizePath(RPath$LibPath[match("proteoCraft", RPath$Package)], winslash = "/")
libPath <- paste0(RPath, "/proteoCraft")

# Default locations 
dflt <- "B:/group/mspecgrp/_Archive/Acquired_data"
if (!dir.exists(dflt)) { dflt <- "B:/group/lsfgrp/Mass_Spec/Acquired_data_v2" }
if (!dir.exists(dflt)) { dflt <- "B:/archive/lsfgrp/MS/Acquired_data" }
if (!dir.exists(dflt)) { dflt <- "D:/Data" }
if (!dir.exists(dflt)) { dflt <- "D:/groups_temp" }

# Create cluster
parSrc <- paste0(libPath, "/extdata/Sources/make_check_Cluster.R")
#rstudioapi::documentOpen(parSrc)
source(parSrc)

# Some useful functions
cleanRawNm <- \(rawFileNm) { sub("\\.raw$", "", rawFileNm, ignore.case = TRUE) }
environment(cleanRawNm) <- .GlobalEnv
parallel::clusterExport(parClust, "cleanRawNm", envir = environment())

# Install/load rawrr
rawrrSrc <- paste0(libPath, "/extdata/Sources/install_rawrr.R")
#rstudioapi::documentOpen(rawrrSrc)
source(rawrrSrc)
getInt <- !inherits(rawrr_tst, "try-error")
#
if (!getInt) { stop("Check your rawrr installation!") }
# Function to get an XIC
xicFun <- \(fl, mass, tol, filter) { #fl <- fls[1L]
  #mass <- 134.0472
  #filter <- "FTMS - p ESI Full ms2 328.0452"
  x <- try(rawrr::readChromatogram(fl, type = "xic", mass, tol, filter), silent = TRUE)
  if (!inherits(x, "try-error")) {
    res <- list(Outcome = TRUE,
                Output = data.frame("Raw file" = fl,
                                    "Raw file name" = cleanRawNm(basename(fl)),
                                    "M/Z" = mass,
                                    "Tolerance" = tol,
                                    "Retention time" = as.numeric(x[[1L]]$times),
                                    "Intensity" = as.numeric(x[[1L]]$intensities),
                                    check.names = FALSE))
  } else { res <- list(Outcome = FALSE,
                       Error = x) }
  return(res)
}
environment(xicFun) <- .GlobalEnv

# Select Thermo Raw files
# Raw files
if (exists("fls")) {
  fls <- grep("\\.raw$", fls, ignore.case = TRUE, value = TRUE)
  fls <- fls[file.exists(fls)]
} else { fls <- c() }
l <- length(fls)
if (l) {
  msg <- paste0("There are ", l, " input Raw files already present in environment from a previous run:")
  opt <- c("Remove them                                                                                                             ",
           "Reprocess them all                                                                                                      ",
           "Change the selection                                                                                                    ")
  myChoice <- svDialogs::dlg_list(opt, opt[1L], title = msg)$res
  if (myChoice == opt[1L]) { fls <- c() }
}
l <- length(fls)
if (!l) {
  wd <- rstudioapi::selectDirectory(path = dflt)
  #wd <- unique(dirname(fls))
  #filt <- matrix(data = c("Thermo raw file", "*.raw;*.RAW"), ncol = 2, dimnames = list("raw file"))
  #fls <- "B:/archive/lsfgrp/MS/Acquired data/frimlgrp/LCMS_JiFrATeplova1_m"
  allRawFls <- list.files(wd, "*\\.raw$", FALSE, TRUE, TRUE, TRUE)
  allRawFls <- setNames(allRawFls, sub(topattern(paste0(wd, "/")), "", allRawFls))
  age <- setNames(vapply(allRawFls, \(x) { file.info(x)$mtime }, 1), allRawFls)
  allRawFls <- allRawFls[order(age, decreasing = TRUE)]
  blnks <- grep("blank", allRawFls, invert = TRUE)
  fls <- svDialogs::dlg_list(names(allRawFls), names(allRawFls)[blnks], TRUE, "Select raw files to analyse")$res
  fls <- allRawFls[match(fls, names(allRawFls))]
  names(fls) <- NULL
  #fls <- normalizePath(choose.files(paste0(dflt, "/*.raw"), filters = filt), winslash = "/")
  #fls <- fls[order(file.info(fls)$mtime, decreasing = FALSE)]
  if (getInt) {
    # Here we try to obtain more realistic file modified times than from Windows, using the header
    tst <- parallel::parSapply(parClust, fls, \(fl) { #fl <- fls[1L] #fl <- B:/group/lsfgrp/Mass_Spec/Acquired_data_v2/frimlgrp/LCMS_JiFrHAi3_m/saturation_20250410102029.raw
      rs <- try({ rawrr::readFileHeader(fl)$`Creation date` }, silent = TRUE)
      if (inherits(rs, "try-error")) {
        #rs <- NA
        rs <- format(file.info(fl)$mtime, "%d/%m/%Y %H:%M:%S", tz = Sys.timezone())
      }
      return(rs)
    })
    tst <- as.POSIXct(tst, tz = Sys.timezone(), format = "%d/%m/%Y %H:%M:%S")
    fls <- fls[order(tst, decreasing = FALSE)]
  }
}

# Extract full chromatogram for quantifier ion of interest
tol <- 20 # In ppm!
#tol <- as.numeric(svDialogs::dlg_input("Enter mass tolerance (ppm)", tol)$res)
MS2s <- parallel::parLapply(parClust, fls, \(fl) { #fl <- fls[1L]
  x <- rawrr::readIndex(fl)
  x <- x[x$MSOrder == "Ms2",]
  return(unique(sub("@.*", "", unique(x$scanType))))
})
ms2s <- unique(unlist(MS2s))
if (length(ms2s) == 1L) {
  cat(paste0("Single MS2 filter detected: \"", ms2s, "\"\n"))
} else {
  tmp <- vapply(ms2s, \(x) { paste(c(x, rep(" ", 250L-nchar(x))), collapse = "") }, "")
  tmp <- svDialogs::dlg_list(tmp, tmp[1L], TRUE, "Select filter(s) for which you want to get an XIC")$res
  ms2s <- sub(" *$", "", tmp)
}
dflts <- c("134.0472", "150.0421", NA)
masses <- setNames(lapply(ms2s, \(ms2) { #ms2 <- ms2s[1L]
  msg <- paste0(ms2, " XIC: enter M/Z to extract (if multiple M/Zs, separate them with \" \")")
  m <- match(ms2, c("FTMS - p ESI Full ms2 328.0452",
                    "FTMS - p ESI Full ms2 344.0402"))
  if (is.na(m)) { m <- 3L }
  dflt <- dflts[m]
  mass <- as.numeric(unlist(strsplit(svDialogs::dlg_input(msg, dflt)$res, " ")))
  return(mass[!is.na(mass)])
}), ms2s)
refRT <- lapply(masses, \(x) { NA })
w <- which(vapply(refRT, is.na, TRUE))
while (length(w)) {
  for (ms2 in ms2s[w]) {
    dflt <- refRT[[ms2]]
    if (is.na(dflt)) { dflt <- 6 }
    msg <- paste0(ms2, ": enter expected retention time (in minutes)")
    refRT[[ms2]] <- as.numeric(svDialogs::dlg_input(msg, dflt)$res)
  }
  w <- which(vapply(refRT, is.na, TRUE))
}

# XICs
fltms2 <- listMelt(masses, names(masses), c("Mass", "Filter"))
xicNms <- do.call(paste, c(fltms2[, c("Mass", "Filter")], sep = " for "))
parallel::clusterExport(parClust, list("tol", "ms2s", "fltms2", "xicFun"), envir = environment())
xic <- setNames(lapply(1L:nrow(fltms2), \(i) {
  setNames(parallel::parLapply(parClust, fls, \(fl) { xicFun(fl, fltms2$Mass[[i]], tol, fltms2$Filter[[i]]) }), fls)
}), xicNms)
w <- setNames(lapply(xicNms, \(x) { #x <- names(xic)[1L]
  which(vapply(xic[[x]], \(y) { y$Outcome }, TRUE))
}), xicNms)
xic <- setNames(lapply(xicNms, \(x) { #x <- names(xic)[1L]
  setNames(lapply(fls[w[[x]]], \(fl) { xic[[x]][[fl]]$Output }), fls[w[[x]]])
}), xicNms)
w <- which(setNames(lapply(xicNms, length), xicNms) > 0)
xic <- xic[w]
nms <- names(xic)
nms2 <- setNames(vapply(strsplit(nms, " for FTMS - p ESI Full ms2 "), \(x) {
  paste(rev(unlist(x)), collapse = "_")
}, ""), nms)
peakRTs <- list()
xicDir <- paste0(wd, "/XICs")
if (!dir.exists(xicDir)) { dir.create(xicDir, recursive = TRUE) }
for (nm in names(xic)) { #nm <- names(xic)[1L]
  ms2 <- sub(".* for ", "", nm)
  if (length(xic[[nm]])) {
    peakRTs[[nm]] <- list()
    for (fl in names(xic[[nm]])) { #fl <- names(xic[[nm]])[1L]
      flNm <- gsub(".*/|\\.raw$", "", fl, ignore.case = TRUE)
      dat <- xic[[nm]][[fl]]
      dat <- dat[, c("Retention time", "Intensity")]
      tst <- abs(dat$`Retention time` - refRT[[ms2]])
      w <- which(tst == min(tst))[1L]
      rg <- (-50L:50L)+w
      rgMx <- max(dat$Intensity[rg])
      w2 <- rg[dat$Intensity[rg] == rgMx]
      peakRTs[[nm]][[fl]] <- dat$`Retention time`[w2]
      # Save
      write.csv(dat, paste0(xicDir, "/XIC_", flNm, "_-_", nms2[nm], ".csv"), row.names = FALSE)
    }
  }
}
#
# Get averaged MS2 for each file
ms2Dir <- paste0(wd, "/MS2s")
if (!dir.exists(ms2Dir)) { dir.create(ms2Dir, recursive = TRUE) }
for (nm in names(xic)) { #nm <- names(xic)[1L]
  ms2 <- sub(".* for ", "", nm)
  cat(" ->", ms2, "\n")
  for (fl in names(xic[[nm]])) { #fl <- names(xic[[nm]])[1L]
    cat("    -", fl, "\n")
    flNm <- gsub(".*/|\\.raw$", "", fl, ignore.case = TRUE)
    # Load file index
    ind <- rawrr::readIndex(fl)
    ind$scanType <- sub("@.*", "", ind$scanType)
    ind <- ind[ind$scanType == ms2,]
    ind$Diff <- abs(ind$StartTime - peakRTs[[nm]][[fl]])
    w <- which(ind$Diff == min(ind$Diff))[1L]
    # Get +/- 5 scans
    w <- w+(-5L:5L)
    parallel::clusterExport(parClust, "fl", envir = environment())
    MS2scns <- parallel::parLapply(parClust, ind$scan[w], \(x) {
      rawrr::readSpectrum(fl, x)[[1L]]
    })
    MS2scns <- lapply(MS2scns, \(x) { data.frame(mz = x$mZ,
                                                 Intensity = x$intensity,
                                                 Scan = paste0("scan ", x$scan))
    })
    MS2scns <- do.call(rbind, MS2scns)
    # Save summed scan
    sumMS2scns <- aggregate(MS2scns$Intensity, list(round(MS2scns$mz, 3L)), sum)
    colnames(sumMS2scns) <- c("m/z", "Intensity")
    write.csv(sumMS2scns, paste0(ms2Dir, "/MS2_", flNm, "_-_", nms2[nm], ".csv"), row.names = FALSE)
    # Plot
    MS2scns <- rbind(MS2scns,
                     data.frame(mz = sumMS2scns$`m/z`,
                                Intensity = sumMS2scns$Intensity,
                                Scan = "Summed"))
    MS2scns$Scan <- factor(MS2scns$Scan, levels = unique(MS2scns$Scan))
    plot <- ggplot2::ggplot(MS2scns) +
      ggplot2::geom_line(ggplot2::aes(x = mz, y = Intensity, group = Scan, color = Scan)) +
      ggplot2::facet_wrap(~Scan, scale = "free") + ggplot2::theme_bw() +
      ggplot2::ggtitle(nm, subtitle = paste0("File = ", gsub(".*/|\\.raw", "", fl, ignore.case = TRUE)))
    # Facet scale should at least be:
    #  - "free_x" or "free" if averaging
    #  - "free" if summing
    # but summing gives the best S/N ratio.
    poplot(plot, 12L, 22L)
    suppressMessages({
      ggplot2::ggsave(paste0(ms2Dir, "/MS2_", flNm, "_-_", nms2[nm], ".jpeg"), plot,
                      dpi = 300L, units = "in", width = 10L)
      ggplot2::ggsave(paste0(ms2Dir, "/MS2_", flNm, "_-_", nms2[nm], ".pdf"), plot,
                      dpi = 300L, units = "in", width = 10L)
    })
  }
}
