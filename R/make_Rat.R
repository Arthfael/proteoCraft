#' make_Rat
#'
#' @description
#' Internal function to output consistent ratios and intensity values.
#' 
#' @param myData Input data.
#' @param Priority If set to "int", ratios calculated will reflect the real ratios of intensities. If set to "rat", intensities will be recalculated to reflect ratios.
#' @param int.log Input intensities log scale, default = 10
#' @param rat.log Input ratios log scale, default = 2
#' @param refGroups Aggregate defining ratio groups.
#' @param experiment.map Experiment map.
#' @param experiment.map_col Name of the sample column in the experiment map.
#' @param int.root Root of intensity column names.
#' @param rat.root Root of ratios column names.
#' @param verbose Logical, default = TRUE
#' 
#' @returns
#' A list of two objects:
#'  - Ratios_to_Refs: map of ratios to references,
#'  - Data: data supplemented with new ratios and intensities columns
#' 
#' @export

make_Rat <- function(myData = pep,
                     Priority = "int",
                     int.log = 10,
                     rat.log = 2,
                     refGroups,
                     experiment.map,
                     experiment.map_col = "Ref.Sample.Aggregate",
                     int.root,
                     rat.root,
                     verbose = TRUE) {
  #proteoCraft::DefArg(make_Rat)
  #proteoCraft::DefArg(make_Rat, silent = FALSE)
  if (is.logical(int.log)) {
    if (!is.na(int.log)&&int.log) {
      logTrans <- TRUE
      warning("Assuming default log10 intensities.")
      int.log <- 10
    } else {
      int.log <- FALSE
      logTrans <- FALSE
    }
  }
  if (is.numeric(int.log)) {
    stopifnot(int.log > 0)
    logTrans <- TRUE
  }
  stopifnot(is.numeric(rat.log),
            rat.log > 0)
  if (logTrans) {
    log_Int2Rat <- log(rat.log, int.log)
  }
  rat.2.ref <- lapply(refGroups$values, function(i) { #i <- refGroups$values[1]
    j <- setNames(unlist(strsplit(i, "___")), refGroups$names)
    temp <- lapply(refGroups$names, function(x) {
      if (j[[x]] == "NA") {
        return(which((is.na(experiment.map[[x]]))|(experiment.map[[x]] == j[[x]])))
      }
      return(which(experiment.map[[x]] == j[[x]]))
    })
    temp2 <- sort(unique(unlist(temp)))
    test <- vapply(temp2, function(x) { sum(vapply(temp, function(y) { x %in% unlist(y) }, TRUE)) }, 1)
    temp2 <- temp2[which(test == length(temp))]
    temp3 <- experiment.map[temp2,]
    temp3 <- temp3[which(temp3$Reference),]
    smpls <- temp3[[experiment.map_col]]
    if (!length(smpls)) { return() }
    kol <- paste0(int.root, smpls)
    w <- which(kol %in% colnames(myData))
    if (!length(w)) { return() }
    # Values in log
    val <- as.numeric(apply(myData[, kol[w], drop = FALSE], 1, function(x) {
      if (logTrans) {
        proteoCraft::log_ratio_av(x) # To avoid double calculations, we keep its original base here (from INTENSITIES!)
      } else {
        proteoCraft::log_ratio_av(log(x, rat.log)) # Log with base from RATIOS!
      }
    }))
    nm <- paste0(int.root, i, ".REF")
    map <- data.frame(Name = nm)
    map$Source <- list(smpls)
    return(list(Values = val,
                Name = nm,
                Map = map))
  })
  rat.2.ref <- rat.2.ref[which(vapply(rat.2.ref, function(x) { "Values" %in% names(x) }, TRUE))]
  nms <- vapply(rat.2.ref, function(x) { x$Name }, "")
  #
  # Reference values
  vals <- lapply(rat.2.ref, function(x) { x$Values })
  vals <- as.data.frame(do.call(cbind, vals))
  colnames(vals) <- nms
  myData[, nms] <- vals # vals is always log scaled
  if (logTrans) {
    # In this case, the data in vals is log but base from INTENSITIES! => we need to adjust base for ratios calculations!
    vals <- vals/log_Int2Rat
  } else {
    # Here base is RATIOS and we need to delog
    myData[, nms] <- rat.log^myData[, nms]
  }
  # Now vals has always log base from RATIOS
  #
  # Mappings of ratios to references
  rat.2.ref <- lapply(rat.2.ref, function(x) { x$Map })
  rat.2.ref <- plyr::rbind.fill(rat.2.ref)
  #
  # Re-calculate individual ratios/intensities to relevant reference,
  # If the reference is an average, we will also calculate individual ref ratios to it;
  # this will be useful further down the line.
  tmpRes <- lapply(refGroups$values, function(i) { #i <- refGroups$values[1]
    j <- set_names(unlist(strsplit(i, "___")), refGroups$names)
    temp <- lapply(refGroups$names, function(x) {
      if (j[[x]] == "NA") {
        return(which((is.na(experiment.map[[x]]))|(experiment.map[[x]] == j[[x]])))
      }
      return(which(experiment.map[[x]] == j[[x]]))
    })
    temp2 <- sort(unique(unlist(temp)))
    test <- vapply(temp2, function(x) { sum(vapply(temp, function(y) { x %in% unlist(y) }, TRUE)) }, 1)
    temp2 <- temp2[which(test == length(temp))]
    temp3 <- experiment.map[temp2,]
    # Get reference
    k <- j[which(names(j) %in% refGroups$names)]
    b <- paste0(int.root, paste(k, collapse = "___"), ".REF")
    rfVct <- vals[[b]] # Already has the correct RATIOS log base
    #
    A <- a <- temp3[[experiment.map_col]]
    if (length(which(temp3$Reference)) == 1) {
      # If there is only one ref in the group, remove it as there is no point calculating a ratio to itself,
      # or re-calculating its expression to itself.
      a <- temp3[which(!temp3$Reference), experiment.map_col]
    }
    a <- a[which(paste0(int.root, a) %in% colnames(myData))]
    if (!length(a)) {
      warning(paste0(paste0("There are no ", c("ratios", "expression values")[match(Priority, c("int", "rat"))],
                            " to calculate for level ", i)))
      return()
    }
    intKl <- paste0(int.root, a)
    ratKl <- paste0(rat.root, a)
    if (Priority == "int") { # Here we calculate log ratios so they precisely reflect log expression
      nuDat <- myData[, intKl, drop = FALSE]
      # Adjust log base
      if (logTrans) {
        nuDat <- nuDat/log_Int2Rat
      } else {
        nuDat <- log(nuDat, rat.log)
      }
      #
      nuDat <- as.data.frame(sweep(nuDat, 1, rfVct, "-"))
      colnames(nuDat) <- ratKl
      #myData[, ratKl] <- nuDat
    }
    if (Priority == "rat") { # Here we re-calculate log expression to reflect log ratios
      nuDat <- vals[, ratKl, drop = FALSE] # Base RATIOS
      nuDat <- sweep(nuDat, 1, rfVct, "+")
      # Adjust log base
      if (logTrans) {
        nuDat <- nuDat*log_Int2Rat
      } else {
        nuDat <- rat.log^(nuDat)
      }
      #
      nuDat <- as.data.frame(nuDat)
      colnames(nuDat) <- intKl
      #myData[, intKl] <- nuDat
    }
    msgs <- c()
    for (a1 in A) { #a1 <- A[1]
      tmp <- proteoCraft::is.all.good(myData[[paste0(int.root, a1)]])
      if (!int.log) { tmp <- tmp[which(tmp > 0)] }
      emd <- signif(median(tmp), 3)
      esd <- signif(sd(tmp), 3)
      if (logTrans) {
        msg <- paste0("-> Sample ", proteoCraft::cleanNms(a1), ":\n",
                      paste0("   - log", int.log, "(expr.): median = ", emd, ", SD = ", esd, "\n"))
      } else {
        msg <- paste0("-> Sample ", proteoCraft::cleanNms(a1), ":\n",
                      paste0("   - expr.: median = ", emd, ", SD = ", esd, "\n"))
      }
      if (paste0(rat.root, a1) %in% colnames(myData)) {
        tmp <- proteoCraft::is.all.good(myData[[paste0(rat.root, a1)]])
        rmd <- signif(median(tmp), 3)
        rsd <- signif(sd(tmp), 3)
        msg <- paste0(msg, paste0("   - log", rat.log, "(ratio): median = ", rmd, ", SD = ", rsd, "\n"))
      } # else { msg <- paste0(msg, "     (no ratios calculated)\n") }
      msgs <- c(msgs, msg)
    }
    myRes <- list(Data = nuDat,
                  Ref = b,
                  Samples = a,
                  Message = msgs)
    return(myRes)
  })
  #
  nuData <- as.data.frame(do.call(cbind, lapply(tmpRes, function(x) { x$Data })))
  #
  tmp <- setNames(lapply(tmpRes, function(x) { x$Samples }),
                  vapply(tmpRes, function(x) { x$Ref }, ""))
  rat.2.ref$Used_by <- lapply(rat.2.ref$Name, function(x) { tmp[[x]] })
  #
  if (verbose) {
    tmp <- paste(unlist(lapply(tmpRes, function(x) { x$Message })), collapse = "")
    cat(tmp, "\n")
  }
  return(list(Data = nuData,
              Ratios_to_Refs = rat.2.ref))
}
