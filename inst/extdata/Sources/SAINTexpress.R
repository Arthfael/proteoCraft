################
# SAINTexpress #
################
if (saintExprs) { saintExprs <- "Target" %in% colnames(Exp.map) }
if (saintExprs) {
  SaintRoot <- "C:/SAINTexpress"
  SaintDir <- paste0(SaintRoot, "/SAINTexpress_v3.6.3__2018-03-09")
  if (!dir.exists(SaintDir)) { dir.create(SaintDir, recursive = TRUE) }
  SaintEx <- paste0(SaintDir, "/Precompiled_binaries/Windows64/SAINTexpress-int.exe")
  if (!file.exists(SaintEx)) {
    url <- "https://download.sourceforge.net/saint-apms/SAINTexpress_v3.6.3__2018-03-09.tar.gz"
    packs <- c("curl")
    for (pack in packs) {
      if (!suppressMessages(require(pack, character.only = TRUE))) { install.packages(pack, update = FALSE) }
      require(pack, character.only = TRUE)
    }
    cran_req <- unique(c(cran_req, packs))
    destFl <- paste0(SaintRoot, "/SAINTexpress_v3.6.3__2018-03-09.tar.gz")
    kount <- 0L
    while ((!kount)||((kount < 5L)&&(inherits(tst, "try-error")))) {
      tst <- try(download.file(url, destFl), silent = TRUE)
      kount <- kount+1L
    }
    if (inherits(tst, "try-error")) { saintExprs <- FALSE } else {
      gunzip(destFl)
      utils::untar(gsub("\\.gz$", "", destFl), exdir = SaintRoot)
      unlink(destFl)
      unlink(gsub("\\.gz$", "", destFl))
      saintExprs <- file.exists(SaintEx)
    }
  }
}
if (saintExprs) {
  #
  msg <- "Running SAINTexpress analysis..."
  ReportCalls <- AddMsg2Report(Space = FALSE)
  #
  subDr <- "Reg. analysis/SAINTexpress"
  saintDir <- paste0(gsub("/+$", "", wd), "/", subDr)
  if (!dir.exists(saintDir)) { dir.create(saintDir, recursive = TRUE) }
  dirlist <- unique(c(dirlist, saintDir))
  fcRt <- "log2(FC) - "
  fdrRt <- "-log10(BFDR) - "
  kol_A <- c("AvgP", "MaxP", "TopoAvgP", "TopoMaxP", "SaintScore")
  kol_B <- c("OddsScore", "BFDR", "boosted_by")
  #
  # Copy the SAINTexpress manual
  #url <- "https://raw.githubusercontent.com/bornea/APOSTL/master/wk_images/SAINTexpress-manual.pdf"
  #destfile <- paste0(saintDir, "/SAINT express manual.pdf")
  #tst <- try(download.file(url, destfile, "curl"), silent = TRUE)
  #if (inherits(tst, "try-error")) { try(download.file(url, destfile, "wget"), silent = TRUE) }
  # Description:
  # http://saint-apms.sourceforge.net/Main.html
  fl <- system.file("extdata", "SAINTexpress-manual.pdf", package = "proteoCraft")
  try(file.copy(fl, saintDir), silent = TRUE)
  #file.exists(paste0(saintDir, "/", basename(fl)))
  #
  Indiq <- c("T", "C")
  #
  # Translate contrasts into Bait structure
  whSingle <- which(!myContrasts$isDouble)
  m <- match(Exp.map$Ref.Sample.Aggregate, rownames(expMap)) # just to be safe...
  EM <- expMap[m,]
  myContrasts$A_full <- Exp.map[match(myContrasts$A, EM[[VPAL$limmaCol]]), VPAL$column]
  myContrasts$B_full <- Exp.map[match(myContrasts$B, EM[[VPAL$limmaCol]]), VPAL$column]
  allContr <- aggregate(myContrasts$A_full[whSingle], list(myContrasts$B_full[whSingle]), list)
  colnames(allContr) <- c("Ref", "Treat")
  allContr$Contr <- apply(allContr[, c("Ref", "Treat")], 1L, \(x) { #x <- allContr[1L, c("Ref", "Treat")]
    A <- unlist(x[2L])
    B <- unlist(x[1L])
    vapply(A, \(y) { myContrasts$Contrast[which((myContrasts$A_full == y)&(myContrasts$B_full == B))] }, "")
  })
  #
  # Prey and GO tables are not Target specific
  mtch <- listMelt(strsplit(PG$"Leading protein IDs", ";"), PG$id, c("Protein", "PG id"))
  mtch <- set_colnames(aggregate(mtch$`PG id`, list(mtch$Protein), unique), c("Protein", "PG ids"))
  Prey <- data.frame(Protein = mtch$Protein)
  Prey[, c("Sequence", "Gene")] <- db[match(gsub("^CON__", "", Prey$Protein), gsub("^CON__", "", db$`Protein ID`)),
                                      c("Sequence", "Gene")]
  Prey$Length <- nchar(Prey$Sequence)
  w <- which(is.na(Prey$Length))
  if (length(w)) {
    Prey$Length[w] <- median(Prey$Length, na.rm = TRUE)
  }
  Prey$Sequence <- NULL
  Prey <- Prey[, c("Protein", "Length", "Gene")]
  w <- which((is.na(Prey$Gene))|(Prey$Gene == ""))
  l <- length(w)
  if (l) { Prey$Gene[w] <- paste0("Dummy_Gene", seq_len(l)) }
  #
  # Prepare output table
  allSAINTs %<o% Prey
  # - match to PGs
  tmp1 <- listMelt(strsplit(PG$"Leading protein IDs", ";"), PG$id, ColNames = c("Prot", "PG"))
  allSAINTs$PG_id <- as.integer(tmp1$PG[match(allSAINTs$Protein, tmp1$Prot)])
  #w <- which(is.na(allSAINTs$PG_id)) # should be length 0 because we are working from "Leading protein IDs" when creating Prey/allSAINTs
  allSAINTs$"Rel. av. log10 abundance" <- PG$"Rel. av. log10 abundance"[match(allSAINTs$PG_id, PG$id)] # Not exactly, but good enough for now (this is using PG-level data for a Protein-level table)
  # Ultimately, we should also introduce a proteins-(not PG!)-level quantitative table...
  tmp <- listMelt(strsplit(pep$Proteins, ";"), seq_len(nrow(pep)), c("Protein", "Row"))
  tmp <- aggregate(tmp$Row, list(tmp$Protein), length)
  allSAINTs$"Rel. log10(Peptides count)" <- tmp$x[match(allSAINTs$Protein, tmp$Group.1)]
  #
  if (Annotate) {
    w <- which(db$`Protein ID` %in% allSAINTs$Protein)
    GO <- listMelt(strsplit(db$`GO-ID`[w], ";"), db$`Protein ID`[w], c("L1", "value"))
  }
  #
  clusterExport(parClust,
                list("Exp", "VPAL", "RG", "Annotate", "wd", "saintDir", "fcRt", "SaintEx", "kol_A", "kol_B",
                     "cleanNms", "topattern", "is.all.good"),
                envir = environment())
  #
  SAINT_list <- lapply(1L:nrow(allContr), \(ii) { #ii <- 1L
    w0 <- which(Exp.map[[VPAL$column]] %in% allContr$Ref[ii])
    w1 <- which(Exp.map[[VPAL$column]] %in% allContr$Treat[[ii]])
    em <- Exp.map[c(w0, w1), c("Ref.Sample.Aggregate", VPAL$column, RG$column, "Target")]
    em$Reference <- c(rep(TRUE, length(w0)),
                      rep(FALSE, length(w1)))
    em$Contrast <- NA_character_
    w0 <- which(em$Reference)
    w1 <- which(!em$Reference)
    em$Contrast[w1] <- allContr$Contr[[ii]][match(em[w1, VPAL$column], allContr$Treat[[ii]])]
    em$IP_name <- gsub("\\.", "", cleanNms(em$Ref.Sample.Aggregate, rep = ""))
    #
    Bait <- data.frame(IP_name = em$IP_name,
                       Bait = em$Target,
                       Indicator = Indiq[em$Reference+1L])
    #Bait <- Bait[which((!is.na(Bait$Bait))&(Bait$Bait %in% Prey$Protein)),]
    Bait$Bait[which((is.na(Bait$Bait))|(!Bait$Bait %in% Prey$Protein))] <- "CONTROL"
    # if (("CONTROL" %in% Bait$Bait)&&(!"CONTROL" %in% Prey$Protein)) {
    #   # Add a dummy CONTROL protein if any co-IPs do not have a bait (typically IP- isotype controls)
    #   Prey <- rbind(Prey,
    #                 data.frame(Protein = "CONTROL",
    #                            Length = 10000, # Cannot be NA
    #                            Gene = "dummyGene"))
    # }
    #
    kol0 <- paste0(Prot.Expr.Root, em$Ref.Sample.Aggregate[w0])
    kol1 <- paste0(Prot.Expr.Root, em$Ref.Sample.Aggregate[w1])
    kol <- c(kol0, kol1)
    kol0k <- paste0("Evidences count - ", em$Ref.Sample.Aggregate[w0])
    kol1k <- paste0("Evidences count - ", em$Ref.Sample.Aggregate[w1])
    kolk <- c(kol0k, kol1k)
    klnms <- gsub("\\.", "", cleanNms(em$Ref.Sample.Aggregate[c(w0, w1)], rep = ""))
    tmp <- PG[, c("id", kol, kolk)]
    clusterExport(parClust,
                  list("tmp", "kol", "kolk", "em"),
                  envir = environment())
    Interact <- allSAINTs[, "Protein", drop = FALSE]
    # We are removing:
    #  - values not based on PSMs (e.g. when limpa used, or imputed values)
    #  - values in the same row with nearly identical value (down to 6 decimals)
    Interact[, klnms] <- as.data.frame(t(parSapply(parClust, mtch$`PG ids`, \(x) {
      #x <- mtch$`PG ids`[[1L]]
      #x <- mtch$`PG ids`[[grep(prot.list[1L], mtch$Protein)]]
      m <- match(unlist(x), tmp$id)
      if (!length(m)) { stop(x) }
      x1 <- x2 <- round(as.matrix(tmp[m, kol, drop = FALSE]), 6L)
      # - remove values based on no PSMs
      y <- as.matrix(tmp[m, kolk, drop = FALSE])
      w <- which(y == 0L, arr.ind = TRUE)
      x2[w] <- NA_real_
      # - also remove (in case we used limpa) values from the same row and with the same value as removed values
      x2 <- lapply(1L:nrow(x1), \(i) {
        rmv <- x1[i, which(is.na(x2[i,]) & !is.na(x1[i,]))]
        keep <- x2[i,]
        keep[which(keep %in% rmv)] <- NA_real_
        return(keep)
      })
      x2 <- do.call(rbind, x2)
      #
      x2 <- as.numeric(apply(x2, 2L, \(x) { mean(is.all.good(x), na.rm = TRUE) }))
      return(x2)
    })))
    Interact <- dfMelt(Interact, c("Protein", "IP_name", "Intensity"), "Protein")
    Interact <- Interact[which(!is.na(Interact$Intensity)),]
    Interact$IP_name <- as.character(Interact$IP_name)
    Interact$Intensity <- 10L^Interact$Intensity
    Interact$Bait <- Bait$Bait[match(Interact$IP_name, Bait$IP_name)]
    Interact <- Interact[, c("IP_name", "Bait", "Protein", "Intensity")]
    Grps <- unique(em$Contrast[which(!em$Reference)])
    #
    tmpFl1 <- tempfile(fileext = ".rds")
    tmpFl2 <- tempfile(fileext = ".rds")
    readr::write_rds(Interact, tmpFl1)
    readr::write_rds(Prey, tmpFl2)
    exports <- list("Bait", "Grps", "tmpFl1", "tmpFl2")
    if (Annotate) {
      tmpFl3 <- tempfile(fileext = ".rds")
      readr::write_rds(GO, tmpFl3)
      exports <- append(exports, "tmpFl3")
    } 
    clusterExport(parClust,
                  exports,
                  envir = environment())
    invisible(clusterCall(parClust, \() {
      assign("Interact", readr::read_rds(tmpFl1), envir = .GlobalEnv)
      assign("Prey", readr::read_rds(tmpFl2), envir = .GlobalEnv)
      if (Annotate) {
        assign("GO", readr::read_rds(tmpFl3), envir = .GlobalEnv)
      }
      return()
    }))
    unlink(tmpFl1)
    unlink(tmpFl2)
    if (Annotate) {
      unlink(tmpFl3)
    }
    saint_Lst <- setNames(parLapply(parClust, Grps, \(grp) { #grp <- Grps[1L]
      grpMtch <- match(grp, Grps)
      dr <- paste0(saintDir, "/", grp)
      if (!dir.exists(dr)) { dir.create(dr, recursive = TRUE) }
      setwd(dr)
      #
      em_ <- em[which(em$Reference | em$Contrast == grp),]
      Bait2 <- Bait[match(em_$IP_name, Bait$IP_name),]
      Interact2 <- Interact[which(Interact$IP_name %in% Bait2$IP_name),]
      Prey2 <- Prey[which(Prey$Protein %in% c(Bait2$Bait, Interact2$Protein)),]
      Interact2 <- Interact2[which(Interact2$Protein %in% c(Bait2$Bait, Prey$Protein)),]
      baitFl <- paste0(dr, "/", grpMtch, "_tempBait.txt")
      preyFl <- paste0(dr, "/", grpMtch, "_tempPrey.txt")
      interFl <- paste0(dr, "/", grpMtch, "_tempInteract.txt")
      data.table::fwrite(Bait2, baitFl, quote = FALSE, col.names = FALSE, row.names = FALSE, eol = "\n", sep = "\t", na = "NA")
      data.table::fwrite(Prey2, preyFl, quote = FALSE, col.names = FALSE, row.names = FALSE, eol = "\n", sep = "\t", na = "NA")
      data.table::fwrite(Interact2, interFl, quote = FALSE, col.names = FALSE, row.names = FALSE, eol = "\n", sep = "\t", na = "NA")
      goFl <- paste0(dr,  "/", grpMtch, "_tempGO.txt")
      if (Annotate) {
        GO2 <- GO[which(GO$value %in% Prey$Protein),]
        GO2 <- aggregate(GO2$value, list(GO2$L1), \(x) { paste(unique(x), collapse = " ") })
        data.table::fwrite(GO2, goFl, quote = FALSE, col.names = FALSE, row.names = FALSE, eol = "\n", sep = "\t", na = "NA")
      }
      fls <- c(SaintEx, c(interFl, preyFl, baitFl, goFl)[seq_len(3L+Annotate)])
      w <- which(!file.exists(fls))
      stopifnot(length(w) == 0L)
      #fls[w]
      fls2 <- paste0("\"", fls, "\"")
      # cmd <- paste0(c(fls2[1L],
      #                 paste0("-L", sum(Exp.map$Reference)),
      #                 fls2[2L:length(fls2)]), collapse = " ")
      lstFl <- paste0(dr, "/list.txt") # Important: it saves the list in the local directory!!!
      if (file.exists(lstFl)) { unlink(lstFl) }
      cmd <- paste0(c(fls2[1L],
                      paste0("-L", sum(em_$Reference)),
                      fls2[2L:length(fls2)]), collapse = " ")
      #cat(cmd)
      #writeClipboard(cmd)
      #cat("   ", grp2, "\n")
      system(cmd)
      #cat("    -> done\n")
      rs <- list(Outcome = file.exists(lstFl))
      if (rs$Outcome) {
        # Read and process results
        tmpDF <- data.table::fread(lstFl, integer64 = "numeric", check.names = FALSE, data.table = FALSE)
        tmpDF <- tmpDF[order(tmpDF$SaintScore, decreasing = TRUE),]
        #View(tmpDF)
        #View(tmpDF[grep(prot.list[1L], tmpDF$Prey),])
        for (k in kol_A) {
          w <- which(tmpDF[[k]] == ".")
          if (length(w)) { tmpDF[w, k] <- NA }
          tmpDF[[k]] <- as.numeric(tmpDF[[k]])
        }
        tmpDF[, paste0(kol_A, " - ", grp)] <- tmpDF[, kol_A]
        for (k in kol_B) {
          tmpDF[[paste0(k, " - ", grp)]] <- tmpDF[[k]]
        }
        tmpDF[[paste0(fcRt, grp)]] <- log2(tmpDF$FoldChange)
        tmpDF <- tmpDF[, which(!colnames(tmpDF) %in% c(kol_A, kol_B))]
        tmpDF$"Av. log10 abundance" <- log10(tmpDF$AvgIntensity)
        rs$Table <- tmpDF
      } else { warning(paste0("SAINTexpress analysis failed for group ", grp2)) }
      setwd(wd) # Important to release the subfolder
      return(rs)
    }), Grps)
  })
  #
  m <- match(allSAINTs$Protein, db$`Protein ID`)
  kol <- c("Potential contaminant", "Common Name", "Gene")
  allSAINTs[, kol] <- db[m, kol]
  #
  SUCCESS <- FALSE
  for (ii in seq_along(SAINT_list)) {
    wOK <- which(vapply(SAINT_list[[ii]], \(x) { x$Outcome }, TRUE))
    if (length(wOK)) {
      Grps <- names(SAINT_list[[ii]])[wOK]
      for (grp in Grps) {
        SUCCESS <- TRUE
        kol <- paste0(c(fcRt, paste0(c(kol_A, kol_B), " - ")), grp)
        allSAINTs[, kol] <- NA
        w <- which(allSAINTs$Protein %in% SAINT_list[[ii]][[grp]]$Table$Prey)
        m <- match(allSAINTs$Protein[w], SAINT_list[[ii]][[grp]]$Table$Prey)
        allSAINTs[w, kol] <- SAINT_list[[ii]][[grp]]$Table[m, kol]
      }
    }
  }
  #View(allSAINTs)
  if (SUCCESS) {
    ArbThr <- data.frame(yintercept = -log10(BH.FDR),
                         slope = 0,
                         xintercept = NA_real_,
                         colour = colorRampPalette(c("orange", "red"))(length(BH.FDR)),
                         label = paste0(BH.FDR*100, "% FDR"))
    k1 <- grep("^BFDR - ", colnames(allSAINTs), value = TRUE)
    k2 <- gsub("^BFDR - ", fdrRt, k1)
    tmp1 <- as.matrix(allSAINTs[, k1])
    tmp2 <- -log10(tmp1)
    w <- which(tmp1 == 0, arr.ind = TRUE)
    nr <- nrow(w)
    if (nr) {
      mx <- max(c(is.all.good(as.numeric(tmp2))+0.2, 2))
      ArbThr <- rbind(ArbThr,
                      data.frame(yintercept = mx-0.1,
                                 slope = 0,
                                 xintercept = NA,
                                 colour = "#F000FF",
                                 label = "^ Infinite value (BFDR = 0)"))
      msg <- paste0("      Replacing ", nr, " infinite -log10(BFDR) value", c("", "s")[(nr > 1L)+1L],
                    " with a high but finite value of ", mx, " for the purpose of plotting them!")
      ReportCalls <- AddMsg2Report(Space = FALSE)
      tmp2[w] <- mx
    }
    allSAINTs[, k2] <- as.data.frame(tmp2)
    #
    Parma <- Param
    Parma$Plot.labels <- "Common Name"
    Parma$Plot.threshold.metrics <- Parma$Plot.threshold.values <- Parma$Plot.threshold.tests <- Parma$Plot.threshold.colours <- ""
    # Volcano plot
    setwd(wd)
    volcPlot_args2 <- volcPlot_args
    volcPlot_args2$Prot <- allSAINTs
    volcPlot_args2$Proteins.col <- "Protein"
    volcPlot_args2$X.root <- fcRt
    volcPlot_args2$Y.root <- fdrRt
    volcPlot_args2$parameters <- Parma
    volcPlot_args2$title <- "SAINTexpress volcano plot "
    volcPlot_args2$subfolder <- subDr
    volcPlot_args2$plotly_labels <- setNames(c("Common Name", "Protein", "Gene"),
                                             c("Protein name", "Protein ID", "Gene"))
    volcPlot_args2$arbitrary.lines <- ArbThr
    volcPlot_args2$labels <- "thresholds"
    volcPlot_args2$cl <- parClust
    # For testing:
    #DefArg(Volcano.plot);TESTING <- TRUE
    #invisible(lapply(names(volcPlot_args2), \(x) { assign(x, volcPlot_args2[[x]], envir = .GlobalEnv); return() }))
    tempVPip <- do.call(Volcano.plot, volcPlot_args2)
    #
    # Save plotly plots
    dr <- saintDir
    myPlotLys <- tempVPip$`Plotly plots`
    Src <- paste0(libPath, "/extdata/Sources/save_Plotlys.R")
    #rstudioapi::documentOpen(Src)
    source(Src, local = FALSE)
    #
    VP_list <- tempVPip
    insrt <- ""
    Src <- paste0(libPath, "/extdata/Sources/thresholds_Excel.R")
    #rstudioapi::documentOpen(Src)
    source(Src, local = FALSE)
    #
    # Folder cleanup
    dirs <- list.dirs(saintDir)
    dirs <- setdiff(dirs, saintDir)
    unlink(dirs, force = TRUE)
    dirs <- list.dirs(saintDir)
    dirs <- setdiff(dirs, saintDir)
    if (length(dirs)) {
      for (dr in dirs) {
        shell(paste0("RMDIR /S /Q \"", dr, "\""), mustWork = FALSE)
      }
    }
    #
    # Regulated columns
    # Here we do not use the ones from Volcano.plot because the logic is different than usual
    nms <- allContr$Contr
    regkol <- paste0("Regulated - ", nms)
    Reg_filters$"SAINTexpress" <- list()
    if ("con" %in% filter_types) {
      Reg_filters$"SAINTexpress"$"By condition" <- list()
    }
    allSAINTs[, regkol] <- ""
    for (nm in nms) { #nm <- nms[1L]
      twoSided <- !myContrasts$`Up-only`[match(nm, myContrasts$Contrast)]
      thresh <- tempVPip$Thresholds$Absolute[[nm]]
      up <- thresh$Value[match("up", thresh$Levels)]
      down <- thresh$Value[match("down", thresh$Levels)]
      k1 <- paste0("BFDR - ", nm)
      k2 <- paste0("log2(FC) - ", nm)
      k3 <- paste0("Regulated - ", nm)
      allSAINTs[which(allSAINTs[[k1]] <= max(BH.FDR)), k3] <- "too small FC" # base level
      for (f in rev(BH.FDR)) {
        allSAINTs[which((allSAINTs[[k1]] <= f)&(allSAINTs[[k2]] >= up)), k3] <- paste0("up, FDR = ", f*100,"%")
        if (twoSided) {
          allSAINTs[which((allSAINTs[[k1]] <= f)&(allSAINTs[[k2]] <= down)), k3] <- paste0("down, FDR = ", f*100,"%")
        }
      }
      # Create SAINTexpress-based filters
      if ("con" %in% filter_types) {
        up <- grep("^up|^Specific", unique(unlist(allSAINTs[, regkol])), value = TRUE)
        down <- grep("^down|^Anti-specific", unique(unlist(allSAINTs[, regkol])), value = TRUE)
        #
        flt <- list(Columns = k3,
                    Filter_up = sort(which(allSAINTs[[k3]] %in% up)),
                    Filter_down = sort(which(allSAINTs[[k3]] %in% down)),
                    Filter = sort(which(allSAINTs[[k3]] %in% c(up, down))))
        # We want to translate it into rows from PG
        u <- unique(unlist(flt[paste0("Filter", c("_up", "_down", ""))]))
        m <- match(allSAINTs$PG_id[u], PG$id)
        flt$PG_Filter_up <- m[match(flt$Filter_up, u)]
        flt$PG_Filter_down <- m[match(flt$Filter_down, u)]
        flt$PG_Filter <- m[match(flt$Filter, u)]
        Reg_filters$"SAINTexpress"$"By condition"[[nm]] <- flt
      }
    }
    #
    # Z-scored clustering heatmaps of regulated proteins
    clustersTest <- try({
      clustMode <- "SAINTexpress"
      dataType <- "PG"
      clstSrc <- paste0(libPath, "/extdata/Sources/cluster_Heatmap_Main.R")
      #rstudioapi::documentOpen(clstSrc)
      source(clstSrc, local = FALSE)
    }, silent = TRUE) # Allowed to fail, but with a warning!
    if (inherits(clustersTest, "try-error")) {
      warning("Could not draw heatmap for SAINTexpress results!")
    }
    #
    # Write final table
    data.table::fwrite(allSAINTs, paste0(saintDir, "/SAINTexpress interactions list.tsv"),
                       row.names = FALSE, sep = "\t", na = "NA")
    #
    # Mat-meth text
    l <- length(DatAnalysisTxt)
    DatAnalysisTxt[l] <- paste0(DatAnalysisTxt[l],
                                " Data was also tested with the SAINTexpress algorithm - which directly outputs FDR values - using GO annotations as boosting information.")
    #
  }
  # "Master, don't forget the Athenians!"
  dlg_message("TO DO: add SAINTq for PSM-level analysis of DIA data!", "ok")
}
setwd(wd)
