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
    kount <- 0
    while ((!kount)||((kount < 5)&&("try-error" %in% class(tst)))) {
      tst <- try(download.file(url, destFl), silent = TRUE)
      kount <- kount+1
    }
    if ("try-error" %in% class(tst)) { saintExprs <- FALSE } else {
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
  saintDir <- gsub("/+", "/", paste0(wd, "/", subDr))
  if (!dir.exists(saintDir)) { dir.create(saintDir, recursive = TRUE) }
  dirlist <- unique(c(dirlist, saintDir))
  fcRt <- "log2(FC) - "
  fdrRt <- "-log10(BFDR) - "
  #
  # Copy the SAINTexpress manual
  #url <- "https://raw.githubusercontent.com/bornea/APOSTL/master/wk_images/SAINTexpress-manual.pdf"
  #destfile <- paste0(saintDir, "/SAINT express manual.pdf")
  #tst <- try(download.file(url, destfile, "curl"), silent = TRUE)
  #if ("try-error" %in% class(tst)) { try(download.file(url, destfile, "wget"), silent = TRUE) }
  # Description:
  # http://saint-apms.sourceforge.net/Main.html
  fl <- system.file("extdata", "SAINTexpress-manual.pdf", package = "proteoCraft")
  try(file.copy(fl, saintDir), silent = TRUE)
  #file.exists(paste0(saintDir, "/", basename(fl)))
  #
  Indiq <- c("T", "C")
  if (Mirror.Ratios) { Indiq <- rev(Indiq) } # Deprecate me please!!!
  # Prey and GO tables are not Target specific
  mtch <- listMelt(strsplit(PG$"Protein IDs", ";"), PG$id, c("Protein", "PG id"))
  mtch <- set_colnames(aggregate(mtch$`PG id`, list(mtch$Protein), unique), c("Protein", "PG ids"))
  Interact <- Prey <- data.frame(Protein = mtch$Protein)
  Prey[, c("Sequence", "Gene")] <- db[match(gsub("^CON__", "", Prey$Protein), gsub("^CON__", "", db$`Protein ID`)), c("Sequence", "Gene")]
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
  if (Annotate) {
    GO <- listMelt(strsplit(PG$`GO-ID`, ";"), PG$id)
    GO$Proteins <- PG$`Protein IDs`[match(GO$L1, PG$id)]
    GO <- listMelt(strsplit(GO$Proteins, ";"), GO$value)
  }
  Bait <- data.frame(IP_name = gsub("\\.", "", cleanNms(Exp.map$Ref.Sample.Aggregate, rep = "")),
                     Bait = Exp.map$Target,
                     Indicator = Indiq[Exp.map$Reference+1])
  #Bait <- Bait[which((!is.na(Bait$Bait))&(Bait$Bait %in% Prey$Protein)),]
  Bait$Bait[which((is.na(Bait$Bait))|(!Bait$Bait %in% Prey$Protein))] <- "CONTROL"
  # if (("CONTROL" %in% Bait$Bait)&&(!"CONTROL" %in% Prey$Protein)) {
  #   # Add a dummy CONTROL protein if any co-IPs do not have a bait (typically IP- isotype controls)
  #   Prey <- rbind(Prey,
  #                 data.frame(Protein = "CONTROL",
  #                            Length = 10000, # Cannot be NA
  #                            Gene = "dummyGene"))
  # }
  kol <- paste0(Prot.Expr.Root, Exp.map$Ref.Sample.Aggregate)
  klnms <- gsub("\\.", "", cleanNms(Exp.map$Ref.Sample.Aggregate, rep = ""))
  tmp <- PG[, c("id", kol)]
  clusterExport(parClust, list("tmp", "kol"), envir = environment())
  Interact[, klnms] <- as.data.frame(t(parSapply(parClust, mtch$`PG ids`, function(x) {
    m <- match(unlist(x), tmp$id)
    if (!length(m)) { stop(x) }
    x <- tmp[m, kol, drop = FALSE]
    x <- as.numeric(apply(x, 2, function(x) { mean(proteoCraft::is.all.good(x)) }))
    return(x)
  })))
  Interact <- reshape::melt(Interact, id.vars = "Protein")
  colnames(Interact) <- c("Protein", "IP_name", "Intensity")
  Interact$IP_name <- as.character(Interact$IP_name)
  Interact$Intensity <- 10^Interact$Intensity
  Interact <- Interact[which(is.all.good(Interact$Intensity, 2)),]
  Interact <- Interact[which(Interact$Intensity > 0),]
  #Interact <- Interact[which(Interact$IP_name %in% Bait$IP_name),]
  Interact$Bait <- Bait$Bait[match(Interact$IP_name, Bait$IP_name)]
  Interact <- Interact[, c("IP_name", "Bait", "Protein", "Intensity")]
  #data.table::fwrite(Bait, baitFl, quote = FALSE, col.names = FALSE, row.names = FALSE, eol = "\n", sep = "\t", na = "NA")
  #data.table::fwrite(Prey, preyFl, quote = FALSE, col.names = FALSE, row.names = FALSE, eol = "\n", sep = "\t", na = "NA")
  #data.table::fwrite(Interact, interFl, quote = FALSE, col.names = FALSE, row.names = FALSE, eol = "\n", sep = "\t", na = "NA")
  Grps <- unique(Exp.map[which(!Exp.map$Reference), VPAL$column])
  Parma <- Param
  Parma$Plot.labels <- "Common Name"
  Parma$Plot.threshold.metrics <- Parma$Plot.threshold.values <- Parma$Plot.threshold.tests <- Parma$Plot.threshold.colours <- ""
  readr::write_rds(Interact, paste0(saintDir, "/Interact.RDS"))
  readr::write_rds(Prey, paste0(saintDir, "/Prey.RDS"))
  if (Annotate) { readr::write_rds(GO, paste0(saintDir, "/GO.RDS")) } 
  kol1 <- c("AvgP", "MaxP", "TopoAvgP", "TopoMaxP", "SaintScore")
  kol2 <- c("OddsScore", "BFDR", "boosted_by")
  clusterExport(parClust,
                list("Exp.map", "Exp", "VPAL", "RG", "Bait", "Annotate", "SaintEx", "saintDir", "wd", "fcRt", "kol1", "kol2"),
                envir = environment())
  invisible(clusterCall(parClust, function() {
    Interact <<- readr::read_rds(paste0(saintDir, "/Interact.RDS"))
    Prey <<- readr::read_rds(paste0(saintDir, "/Prey.RDS"))
    if (Annotate) {
      GO2 <<- readr::read_rds(paste0(saintDir, "/GO.RDS"))
    } 
    return()
  }))
  unlink(paste0(saintDir, "/Interact.RDS"))
  unlink(paste0(saintDir, "/Prey.RDS"))
  if (Annotate) {
    unlink(paste0(saintDir, "/GO.RDS"))
  }
  saintst <- setNames(parLapply(parClust, Grps, function(grp) { #grp <- Grps[1]
    grp2 <- proteoCraft::cleanNms(grp)
    dr <- paste0(saintDir, "/", grp2)
    if (!dir.exists(dr)) { dir.create(dr, recursive = TRUE) }
    setwd(dr)
    ratgrp <- unique(Exp.map[match(grp, Exp.map[[VPAL$column]]), RG$column])
    m <- Exp.map[which(Exp.map[[RG$column]] == ratgrp),]
    m <- m[which((m[[VPAL$column]] == grp)|(m$Reference)),]
    Bait2 <- Bait[which(Bait$IP_name %in% gsub("\\.", "",
                                               proteoCraft::cleanNms(m$Ref.Sample.Aggregate, rep = ""))),]
    Interact2 <- Interact[which(Interact$IP_name %in% Bait2$IP_name),]
    Prey2 <- Prey[which(Prey$Protein %in% c(Bait2$Bait, Interact2$Protein)),]
    baitFl <- paste0(dr, "/tempBait.txt")
    preyFl <- paste0(dr, "/tempPrey.txt")
    interFl <- paste0(dr, "/tempInteract.txt")
    lstFl <- paste0(dr, "/list.txt") # Important: it saves the list in the local directory!!!
    data.table::fwrite(Bait2, baitFl, quote = FALSE, col.names = FALSE, row.names = FALSE, eol = "\n", sep = "\t", na = "NA")
    data.table::fwrite(Prey2, preyFl, quote = FALSE, col.names = FALSE, row.names = FALSE, eol = "\n", sep = "\t", na = "NA")
    data.table::fwrite(Interact2, interFl, quote = FALSE, col.names = FALSE, row.names = FALSE, eol = "\n", sep = "\t", na = "NA")
    goFl <- paste0(dr, "/tempGO.txt")
    if (Annotate) {
      GO2 <- GO2[which(GO2$value %in% Prey$Protein),]
      GO2 <- aggregate(GO2$value, list(GO2$L1), function(x) { paste(unique(x), collapse = " ") })
      data.table::fwrite(GO2, goFl, quote = FALSE, col.names = FALSE, row.names = FALSE, eol = "\n", sep = "\t", na = "NA")
    }
    fls <- c(SaintEx, c(interFl, preyFl, baitFl, goFl)[seq_len(3+Annotate)])
    w <- which(!file.exists(fls))
    stopifnot(length(w) == 0)
    #fls[w]
    fls2 <- paste0("\"", fls, "\"")
    # cmd <- paste0(c(fls2[1],
    #                 paste0("-L", sum(Exp.map$Reference)),
    #                 fls2[2:length(fls2)]), collapse = " ")
    if (file.exists(lstFl)) { unlink(lstFl) }
    cmd <- paste0(c(fls2[1],
                    paste0("-L", sum(m$Reference)),
                    fls2[2:length(fls2)]), collapse = " ")
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
      for (k in kol1) {
        w <- which(tmpDF[[k]] == ".")
        if (length(w)) { tmpDF[w, k] <- NA }
        tmpDF[[k]] <- as.numeric(tmpDF[[k]])
      }
      tmpDF[, paste0(kol1, " - ", grp)] <- tmpDF[, kol1]
      for (k in kol2) {
        tmpDF[[paste0(k, " - ", grp)]] <- tmpDF[[k]]
      }
      tmpDF[[paste0(fcRt, grp)]] <- log2(tmpDF$FoldChange)
      tmpDF <- tmpDF[, which(!colnames(tmpDF) %in% c(kol1, kol2))]
      tmpDF$"Av. log10 abundance" <- log10(tmpDF$AvgIntensity)
      rs$Table <- tmpDF
    } else { warning(paste0("SAINTexpress analysis failed for group ", grp2)) }
    setwd(wd) # Important to release the subfolder
    return(rs)
  }), Grps)
  setwd(wd)
  saintst <- saintst[which(vapply(saintst, function(x) { x$Outcome }, TRUE))]
  l <- length(saintst)
  if (l) {
    nms <- names(saintst)
    saintst <- setNames(lapply(saintst, function(x) { x$Table }), nms)
    msg <- "   -> Reading results\n"
    ReportCalls <- AddMsg2Report(Space = FALSE)
    allSAINTs %<o% data.frame(Protein = unique(unlist(lapply(saintst, function(x) { x$Prey }))))
    for (i in seq_len(l)) {
      kol <- paste0(c("log2(FC)", kol1, kol2), " - ", nms[i])
      tmp <- saintst[[nms[i]]][, c("Prey", kol)]
      allSAINTs[, kol] <- NA
      w <- which(allSAINTs$Protein %in% tmp$Prey)
      allSAINTs[w, kol] <- tmp[match(allSAINTs$Protein[w], tmp$Prey), kol]
    }
    m <- match(allSAINTs$Protein, db$`Protein ID`)
    kol <- c("Potential contaminant", "Common Name", "Gene")
    allSAINTs[, kol] <- db[m, kol]
    tmp <- listMelt(strsplit(PG$`Protein IDs`, ";"), PG$id)
    tmp <- tmp[which(tmp$value %in% allSAINTs$Protein),]
    #
    tmp$"Rel. av. log10 abundance" <- PG$"Rel. av. log10 abundance"[match(as.numeric(tmp$L1), PG$id)] # Not exactly, but good enough for now (this is using PG-level data for a Protein-level table)
    allSAINTs$"Rel. av. log10 abundance" <- tmp$"Rel. av. log10 abundance"[match(allSAINTs$Protein, tmp$value)]
    tmp <- listMelt(strsplit(pep$Proteins, ";"), seq_len(nrow(pep)), c("Protein", "Row"))
    #
    tmp <- aggregate(tmp$Row, list(tmp$Protein), length)
    allSAINTs$"Rel. log10(Peptides count)" <- tmp$x[match(allSAINTs$Protein, tmp$Group.1)]
    # Since here we will be plotting FDR, not P-values,
    # we can use our parameter BH FDR thresholds as they are
    ArbThr <- data.frame(yintercept = -log10(BH.FDR),
                         slope = 0,
                         xintercept = NA,
                         colour = colorRampPalette(c("orange", "red"))(length(BH.FDR)),
                         label = paste0(BH.FDR*100, "% FDR"))
    #
    k1 <- grep("^BFDR - ", colnames(allSAINTs), value = TRUE)
    k2 <- gsub("^BFDR - ", fdrRt, k1)
    tmp1 <- as.matrix(allSAINTs[, k1])
    tmp2 <- -log10(tmp1)
    w <- which(tmp1 == 0, arr.ind = TRUE)
    nr <- nrow(w)
    if (nr) {
      mx <- max(is.all.good(as.numeric(tmp2)))
      ArbThr <- rbind(ArbThr,
                      data.frame(yintercept = mx+0.1,
                                 slope = 0,
                                 xintercept = NA,
                                 colour = "#F000FF",
                                 label = "^ Infinite value (BFDR = 0)"))
      rp <- mx+0.2
      msg <- paste0("      Replacing ", nr, " infinite -log10(BFDR) values with a high but finite value of ", rp, " for the purpose of plotting them!")
      ReportCalls <- AddMsg2Report(Space = FALSE)
      tmp2[w] <- rp
    }
    allSAINTs[, k2] <- tmp2
    #
    # Match to PGs
    tmp1 <- listMelt(setNames(strsplit(PG$`Leading protein IDs`, ";"), PG$id), ColNames = c("Prot", "PG"))
    tmp2 <- listMelt(setNames(strsplit(PG$`Protein IDs`, ";"), PG$id), ColNames = c("Prot", "PG"))
    allSAINTs$PG_id <- tmp1$PG[match(allSAINTs$Protein, tmp1$Prot)]
    w <- which(is.na(allSAINTs$PG_id))
    if (length(w)) {
      allSAINTs$PG_id[w] <- tmp2$PG[match(allSAINTs$Protein[w], tmp2$Prot)]
    }
    allSAINTs$PG_id <- as.integer(allSAINTs$PG_id)
    #
    data.table::fwrite(allSAINTs, paste0(saintDir, "/SAINTexpress interactions list.tsv"),
                       row.names = FALSE, sep = "\t", na = "NA")
    #
    # Volcano plot
    labKol <- setNames(c("Common Name", "Protein", "Gene"),
                       c("Protein name", "Protein ID", "Gene"))
    tempVPip <- Volcano.plot(Prot = allSAINTs, Proteins.col = "Protein",
                             mode = "custom",
                             experiments.map = Exp.map,
                             X.root = fcRt,
                             Y.root = fdrRt,
                             aggregate.map = Aggregate.map,
                             aggregate.name = VPAL$aggregate,
                             aggregate.list = Aggregate.list, parameters = Parma,
                             save = c("jpeg", "pdf"), labels = "thresholds",
                             Ref.Ratio.values = Ref.Ratios,
                             Ref.Ratio.method = paste0("obs", RefRat_Mode),
                             ratios.FDR = as.numeric(Param$Ratios.Contamination.Rates),
                             arbitrary.lines = ArbThr,
                             proteins = prot.list, proteins_split = protsplit,
                             return = TRUE, return.plot = TRUE,
                             title = "SAINTexpress volcano plot ",
                             subfolder = subDr,
                             subfolderpertype = FALSE, Symmetrical = TwoSided,
                             Alpha = "Rel. log10(Peptides count)",
                             Size = "Rel. av. log10 abundance", Size.max = 2,
                             plotly = create_plotly, plotly_local = create_plotly_local,
                             plotly_labels = labKol,
                             cl = parClust)
    #
    # Save plotly plots
    dr <- paste0(wd, "/", subDr)
    myPlotLys <- tempVPip$`Plotly plots`
    Src <- paste0(libPath, "/extdata/R scripts/Sources/save_Plotlys.R")
    #rstudioapi::documentOpen(Src)
    source(Src, local = FALSE)
    #
    VP_list <- tempVPip
    insrt <- ""
    Src <- paste0(libPath, "/extdata/R scripts/Sources/thresholds_Excel.R")
    #rstudioapi::documentOpen(Src)
    source(Src, local = FALSE)
    #
    # Folder cleanup
    unlink(paste0(saintDir, "/", c("Interact", "Prey"), ".RDS"))
    if (Annotate) { unlink(paste0(saintDir, "/GO.RDS")) }
    unlink(paste0(saintDir, "/", cleanNms(nms)))
    #
    # Regulated columns
    # Here we do not use the ones from Volcano.plot because the logic is different than usual
    regkol <- paste0("Regulated - ", nms)
    Reg_filters$"SAINTexpress" <- list()
    if ("con" %in% filter_types) {
      Reg_filters$"SAINTexpress"$"By condition" <- list()
    }
    allSAINTs[, regkol] <- ""
    for (nm in nms) { #nm <- nms[1]
      thresh <- tempVPip$Thresholds$Absolute[[nm]]
      up <- thresh$Value[match("up", thresh$Levels)]
      if (TwoSided) { down <- thresh$Value[match("down", thresh$Levels)] }
      k1 <- paste0("BFDR - ", nm)
      k2 <- paste0("log2(FC) - ", nm)
      k3 <- paste0("Regulated - ", nm)
      allSAINTs[which(allSAINTs[[k1]] <= max(BH.FDR)), k3] <- "too small FC" # base level
      for (f in rev(BH.FDR)) {
        allSAINTs[which((allSAINTs[[k1]] <= f)&(allSAINTs[[k2]] >= up)), k3] <- paste0("up, FDR = ", f*100,"%")
        if (TwoSided) {
          allSAINTs[which((allSAINTs[[k1]] <= f)&(allSAINTs[[k2]] <= down)), k3] <- paste0("down, FDR = ", f*100,"%")
        }
      }
      # Create SAINTexpress-based filters
      if ("con" %in% filter_types) {
        up <- grep("^up|^Specific", unique(unlist(allSAINTs[, regkol])), value = TRUE)
        down <- grep("^down|^Anti-specific", unique(unlist(allSAINTs[, regkol])), value = TRUE) # The "anti-specific" part will only become relevant if in future I disconnect symmetry and/or adding specific tags from pull-down experiments
        Reg_filters$"SAINTexpress"$"By condition"[[nm]] <- list(Columns = k3,
                                                                Filter_up = sort(which(allSAINTs[[k3]] %in% up)),
                                                                Filter_down = sort(which(allSAINTs[[k3]] %in% down)),
                                                                Filter = sort(which(allSAINTs[[k3]] %in% c(up, down))))
      }
    }
    #
    # Z-scored clustering heatmaps of regulated proteins
    clustersTest <- try({
      clustMode <- "SAINTexpress"
      clstSrc <- paste0(libPath, "/extdata/R scripts/Sources/cluster_Heatmap_Main.R")
      #rstudioapi::documentOpen(clstSrc)
      source(clstSrc, local = FALSE)
    }, silent = TRUE) # Allowed to fail, but with a warning!
    if ("try-error" %in% class(clustersTest)) {
      warning("Could not draw heatmap for SAINTexpress results!")
    }
    #
    # Mat-meth text
    DatAnalysisTxt <- paste0(DatAnalysisTxt, " Data was also tested with the SAINTexpress algorithm - which directly outputs FDR values - using GO annotations as boosting information.")
    #
  }
  # "Master, don't forget the Athenians!"
  dlg_message("TO DO: add SAINTq for PSM-level analysis of DIA data!", "ok")
}
setwd(wd)
