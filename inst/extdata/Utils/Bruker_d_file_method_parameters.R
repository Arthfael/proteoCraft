FILES <- c()
Moar <- TRUE
if ((!exists("dflt"))||(!dir.exists(dflt))) { dflt <- "...Search_Folder" }
if ((!exists("dflt"))||(!dir.exists(dflt))) { dflt <- "D:/" }
while (Moar) {
  msg <- c("Select input Bruker .d folder.", "Select another Bruker .d folder, or escape to stop adding more.")[(length(FILES) > 0L)+1L]
  fl <- rstudioapi::selectDirectory(msg, path = dflt)
  Moar <- length(fl) > 0L
  if (Moar) {
    dflt <- sub("/[^/]+$", "", fl)
    FILES <- unique(c(FILES, fl))
  }
}
path <- paste0(normalizePath(Sys.getenv("HOME"), winslash = "/"), "/R")
if (length(FILES)) {
  require(XML)
  require(RSQLite)
  # LC method
  tstLC <- try({
    LC_method <- paste0(FILES, sub("\\.d$", ".m/hystar.method", sub(".+_", "/", FILES)))
    LC_method.exists <- file.exists(LC_method)
    w <- which(LC_method.exists)
    LC_meth <- lapply(LC_method, \(fl) { #fl <- LC_method[1L]
      x <- xmlToList(fl)
      lc <- x$LCMethodData$ModuleMethods$ModuleMethodData$text
      lc <- xmlToList(lc)
      #lc$ModuleMethodData$Method$AdvancedSettings
      return(lc)
    })
    LC <- data.frame(Files = FILES)
    LC$NANO <- gsub(" +", " ",
                    gsub("Bruker *", "",
                         sub("\t.*", "",
                             vapply(LC_meth, \(x) { x$HyStarMethodData$ModuleName }, "")), ignore.case = TRUE))
    w <- which(LC$NANO == "nanoElute")
    LC$Column <- sapply(LC_meth, \(x) { x$ModuleMethodData$Method$SeparatorName })
    LC$Trap <- sapply(LC_meth, \(x) { x$ModuleMethodData$Method$TrapName })
    LC$usesTrap <- as.logical(toupper(sapply(LC_meth, \(x) { x$ModuleMethodData$Method$UsesTrapColumn })))
    LC$OvenT <- sapply(LC_meth, \(x) {
      c(NA, x$ModuleMethodData$Method$OvenTemperature)[as.logical(toupper(x$ModuleMethodData$Method$IsSetTemperature))+1L]
    })
    LC$Gradient <- lapply(LC_meth, \(x) {
      nm <- names(x$ModuleMethodData$Method$Gradient)
      grd <- sapply(1L:length(nm), \(y) { x$ModuleMethodData$Method$Gradient[[y]] })
      grd <- t(as.data.frame(grd))
      rownames(grd) <- NULL
      return(grd)
    })
    LC$Flow <- sapply(LC_meth, \(x) { x$HyStarMethodData$MainFlow })
    LC$RunTime <- sapply(LC_meth, \(x) { x$HyStarMethodData$NoStandardMethodData })
    kolFl <- paste0(path, "/proteoCraft/LC_columns.xlsx")
    if (!file.exists(kolFl)) { proteoCraft::Configure() }
    allKolumns <- openxlsx::read.xlsx(kolFl, check.names = FALSE)
    if (!"Description" %in% colnames(allKolumns)) { allKolumns$Description <- "" }
    kolDescr <- \(Colonnes) {
      apply(Colonnes[, c("Name", "Material", "Length (cm)", "ID (µm)", "Particles size (µm)", "Vendor", "P/N"), drop = FALSE], 1L, \(x) {
        #x <- Colonnes[1L, c("Name", "Material", "Length (cm)", "ID (µm)", "Particles size (µm)", "Vendor", "P/N")]
        x <- gsub("^ +| +$", "", as.character(unlist(x)))
        dimz <- c(x[3L], x[4L])
        w <- which(!is.na(dimz))
        dimzTst <- length(w) > 0L
        dimz <- if (dimzTst) { paste(paste0(sub("\\.0+$", "", dimz), c(" cm", " µm ID"))[w], collapse = " * ") } else { NA }
        pn <- c(x[6L], x[7L])
        w <- which(!is.na(pn))
        pnTst <- length(w) > 0L
        if (pnTst) { pn <- paste(paste0(c("", "P/N "), pn)[w], collapse = " ") } else { pn <- NA }
        pnTst <- !is.na(x[7L])
        if (pnTst) { x[7L] <- paste0("P/N ", x[7L]) }
        prtcls <- c(x[5L], x[2L])
        w <- which(!is.na(prtcls))
        prtclsTst <- length(w) > 0L
        if (prtclsTst) { prtcls <- paste(paste0(prtcls, c(" µm", "-coated particles"))[w], collapse = " ") } else { prtcls <- NA }
        res <- paste0(x[1L], " column (", paste(c(prtcls, dimz, pn)[c(prtclsTst, dimzTst, pnTst)], collapse = ", "), ")")
        res <- gsub("column column", "column", res, ignore.case = TRUE)
        return(res)
      })
    }
    #w <- which((nchar(allKolumns$Description) == 0L) | is.na(allKolumns$Description))
    w <- 1L:nrow(allKolumns)
    allKolumns$Description[w] <- kolDescr(allKolumns[w,])
    Kolumns <- unique(c(allKolumns$Description[allKolumns$Function == "Analytical"], "None (direct infusion)", "Add new..."))
    preKolumns <- unique(c(allKolumns$Description[allKolumns$Function == "Trap"],  "None (direct injection)", "Add new..."))
    LC$Par <- lapply(LC_meth, \(x) {
      nm <- names(x$ModuleMethodData$Metho$AdvancedSettings$Parameters)
      pr <- sapply(1L:length(nm), \(y) { x$ModuleMethodData$Metho$AdvancedSettings$Parameters[[y]] })
      colnames(pr) <- pr[1L,]
      pr <- pr[2L:nrow(pr),]
      return(as.data.frame(pr))
    })
    kol <- c("Files", "NANO", "Column", "Trap", "usesTrap", "OvenT", "Gradient", "Flow", "RunTime")
    #kol %in% colnames(uBrMeth)
    LC_txt  <- setNames(apply(LC[, kol], 1L, \(x) {
      #x <- uBrMeth[1L, c("Samples", "NANO", "Column", "Trap", "usesTrap", "OvenT", "Gradient", "Flow", "RunTime")]
      NANO <- x[[2L]]
      Moult <- length(x[[1L]]) > 1L
      NANOMAKER <- "unknown manufacturer"
      if (NANO %in% c("RSLC Nano", "Vanquish NEO")) { NANOMAKER <- "ThermoFisher Scientific" }
      if (NANO %in% c(paste0("nanoElute", c("", paste0(" ", 1L:2L))),
                      paste0("Bruker nanoElute", c("", paste0(" ", 1L:2L))))) { NANOMAKER <- "Bruker Daltonics" }
      tstNano <- tolower(substring(NANO, 1L, 1L)) %in% c("a", "e", "i", "o", "u")
      COLUMN <- x[[3L]]
      tstCol <- tolower(substring(COLUMN, 1L, 1L)) %in% c("a", "e", "i", "o", "u")
      USESTRAP <- x[[5L]]
      PRECOLUMN <- x[[4L]]
      ovenTEMP <- x[[6L]]
      tstPreCol <- tolower(substring(PRECOLUMN, 1L, 1L)) %in% c("a", "e", "i", "o", "u")
      FLOWVALUE <- x[[8L]]
      Grad <- as.data.frame(x[[7L]])
      Grad$Time <- as.numeric(Grad$Time)/60
      Grad$Mix <- round(as.numeric(Grad$Mix)*100, 2L)
      GRADLENGTH <- max(Grad$Time)
      Grad.Baseline <- min(Grad$Mix)
      Grad.Plateau <- max(Grad$Mix)
      wGrdPlt <- which(Grad$Mix == Grad.Plateau)
      Grad.PlateauL <- round(Grad$Time[rev(wGrdPlt)[1L]] - Grad$Time[wGrdPlt[1L]], 0L)
      GradStrt <- Grad$Time[Grad$Mix == Grad.Baseline][1L]
      tmp1 <- Grad$Time[wGrdPlt[1L]]
      tmp0 <- Grad$Time[which(Grad$Time == tmp1)-1L]
      GradEnd <- c(tmp0, tmp1)[(tmp1 - tmp0 > 1L)+1L]
      GRADLENGTH <- GradEnd-GradStrt
      Grad <- Grad[which(Grad$Time == GradStrt)[1L]:which(Grad$Time == GradEnd),]
      Grad <- paste(apply(Grad[, c("Time", "Mix")], 1L, \(x) { paste0(x[[1L]], " min, ", x[[2L]], "%") }), collapse = "; ")
      Grad <- paste0(Grad, ", followed immediately by a ", Grad.PlateauL, " min plateau at ", Grad.Plateau, "%")
      LCtxt <- paste0(c("The s", "S")[Moult+1L], "ample", c(" was", "s were")[Moult+1L], " analysed by LC-MS/MS on a",
                      c("", "n")[tstNano+1L], " ", NANO, " nano-HPLC (",
                      NANOMAKER, ").")
      if (COLUMN != "None (direct infusion)") {
        if (USESTRAP) {
          LCtxt <- paste0(LCtxt, ", concentrated over a", c("", "n")[tstPreCol+1L], " ", PRECOLUMN, ", then ")
        } else { LCtxt <- paste0(LCtxt, ". ", c("The s", "S")[Moult2+1L], "ample", c(" was", "s were")[Moult2+1L], " ") }
        LCtxt <- paste0(LCtxt, "bound to a", c("", "n")[tstCol+1L], " ", COLUMN, " heated at ", ovenTEMP, "°C and eluted over the following ", GRADLENGTH,
                        " min gradient: solvent A, CHECK!; solvent B, CHECK!; ",
                        "constant ", FLOWVALUE, " nL/min flow; ",
                        "B percentage: ", Grad, ".")
      } else { LCtxt <- paste0(LCtxt, " by direct infusion.") }
      return(LCtxt)
    }), FILES)
  }, silent = TRUE)
  if (!inherits(tstLC, "try-error")) { print(LC_txt) } else {
    warning(tstLC)
  }
  #
  # Get MS method
  MS_method <- paste0(FILES, "/analysis.tdf")
  MS_method.exists <- file.exists(MS_method)
  w <- which(MS_method.exists)
  MS_meth <- setNames(lapply(MS_method[w], \(fl) { #fl <- MS_method[w[1L]]
    #x <- readr::read_lines(fl)
    #x[1L:100L]
    SQltTDF <- dbConnect(drv = RSQLite::SQLite(), dbname = fl)
    #Obj <- setNames(lapply(nms, \(nm) { suppressWarnings(dbGetQuery(SQltTDF, statement = paste0("SELECT * FROM '", nm, "'"))) }), nms)
    #for (nm in nms) { if (nrow(Obj[[nm]])) { Obj[[nm]] <- cbind(rep(nm, nrow(Obj[[nm]])), Obj[[nm]]); colnames(Obj[[nm]])[1L] <- "Table_name" } }
    #for (nm in nms) { if (nrow(Obj[[nm]])) { View(Obj[[nm]]) } }
    GlobMetDat <- suppressWarnings(dbGetQuery(SQltTDF, statement = paste0("SELECT * FROM '", "GlobalMetadata", "'")))
    PropsDat <- suppressWarnings(dbGetQuery(SQltTDF, statement = paste0("SELECT * FROM '", "Properties", "'")))
    PropsDefDat <- suppressWarnings(dbGetQuery(SQltTDF, statement = paste0("SELECT * FROM '", "PropertyDefinitions", "'")))
    FrmDat <- try(suppressWarnings(dbGetQuery(SQltTDF, statement = paste0("SELECT * FROM '", "PasefFrameMsMsInfo", "'"))), silent = TRUE)
    if (inherits(FrmDat, "try-error")) { FrmDat <- NA }
    DIADat1 <- try(suppressWarnings(dbGetQuery(SQltTDF, statement = paste0("SELECT * FROM '", "DiaFrameMsMsWindows", "'"))), silent = TRUE)
    if (inherits(DIADat1, "try-error")) { DIADat1 <- NA } else {
      DIADat1$CollisionEnergy <- signif(DIADat1$CollisionEnergy, 2L) # Do we really need more precision?!?!
    }
    #DIADat2 <- suppressWarnings(dbGetQuery(SQltTDF, statement = paste0("SELECT * FROM '", "DiaFrameMsMsInfo", "'")))
    #DIADat3 <- suppressWarnings(dbGetQuery(SQltTDF, statement = paste0("SELECT * FROM '", "DiaFrameMsMsWindowGroups", "'")))
    PRMDat <- try(suppressWarnings(dbGetQuery(SQltTDF, statement = paste0("SELECT * FROM '", "PrmFrameMsMsInfo", "'"))), silent = TRUE)
    if (inherits(PRMDat, "try-error")) { PRMDat <- NA }
    dbDisconnect(SQltTDF)
    return(list(Method = list(Global = GlobMetDat,
                              Properties = PropsDat,
                              Definitions = PropsDefDat,
                              Frames = FrmDat,
                              DIA = DIADat1,
                              PRM = PRMDat),
                Name = GlobMetDat$Value[GlobMetDat$Key == "MethodName"]))
  }), FILES[w])
  MS_method_parsed <- setNames(lapply(MS_meth, \(x) { #x <- MS_meth[[1L]]
    fl <- paste0(path, "/proteoCraft/", x$Name)
    GlobMetDat <- x$Method$Global
    PropsDat <- x$Method$Properties
    PropsDefDat <- x$Method$Definitions
    FrmDat <- x$Method$Frames
    DIADat1 <- x$Method$DIA
    PRMDat <- x$Method$PRM
    INSTR <- GlobMetDat$Value[GlobMetDat$Key == "InstrumentName"]
    INSTRMAKER <- GlobMetDat$Value[GlobMetDat$Key == "InstrumentVendor"]
    if (INSTRMAKER == "Bruker") { INSTRMAKER <- "Bruker Daltonics"}
    #View(GlobMetDat)
    #View(PropsDat)
    PropsDat$Property <- PropsDefDat$PermanentName[match(PropsDat$Property, PropsDefDat$Id)]
    #View(PropsDefDat)
    #View(FrmDat)
    MssRng <- sub("\\.$", "", sub("0+$", "", GlobMetDat$Value[match(c("MzAcqRangeLower", "MzAcqRangeUpper"), GlobMetDat$Key)]))
    MobRng <- sub("\\.$", "", sub("0+$", "", GlobMetDat$Value[match(c("OneOverK0AcqRangeLower", "OneOverK0AcqRangeUpper"), GlobMetDat$Key)]))
    if (!is.logical(FrmDat)) {
      isoWdths <- data.frame(Width = c(min(FrmDat$IsolationWidth), max(FrmDat$IsolationWidth)))
      isoWdths$MZ <- round(c(max(FrmDat$IsolationMz[FrmDat$IsolationWidth == isoWdths$Width[1L]]),
                             min(FrmDat$IsolationMz[FrmDat$IsolationWidth == isoWdths$Width[2L]])))
      CEs <- data.frame(CE = round(c(min(FrmDat$CollisionEnergy), max(FrmDat$CollisionEnergy))),
                        IM = as.numeric(MobRng))
      #require(ggplot2)
      #plot1 <- ggplot(FrmDat) + geom_point(aes(x = IsolationMz, y = IsolationWidth, colour = Frame)) + theme_bw()
      #poplot(plot1) # -> Isolation width is strictly a function of m/z and based on the table in the timsControl GUI...
      #plot2 <- ggplot(FrmDat) + geom_point(aes(x = IsolationMz, y = CollisionEnergy, colour = Frame)) + theme_bw()
      #poplot(plot2) # ... but Collision Energy isn't.
      #plot2a <- ggplot(FrmDat) + geom_point(aes(x = IsolationMz, y = CollisionEnergy, colour = ScanNumBegin)) + theme_bw()
      #poplot(plot2a) # ... but Collision Energy isn't.
    } else { isoWdths <- CEs <- NA }
    #
    props <- c("FocusPreTOF_Lens1_TransferTime_Set", "FocusPreTOF_Lens1_PrePulseStorageTime_Set",
               "TOF_DetectorTof_HighSensitivity_Enabled",
               "Mode_IonPolarity", "Mode_ScanMode",
               "IMS_Cycle_RampTime", "IMS_Cycle_AccumulationTime",
               "MSMS_Pasef_NumRampsPerCycle", "MSMS_Pasef_TotalCycleTime",
               "MSMS_Pasef_AllowedCharge_Minimum", "MSMS_Pasef_AllowedCharge_Maximum",
               #"MSMSAuto_PreferChargeRangeLow", "MSMSAuto_PreferChargeRangeHigh", "MSMS_DefaultChargeState",
               "MSMS_Pasef_IntensityThreshold", "MSMS_Pasef_Scheduler_TargetIntensity",
               "MSMS_Pasef_ExclusionReleaseTime", "MSMS_Pasef_ExclusionReconsiderActive", "MSMS_Pasef_ExclusionReconsiderRatio",
               "MSMS_Pasef_ExclusionWindowMassWidth", "MSMS_Pasef_ExclusionWindow1PerK0Width")
    PropsDF <- PropsDat[PropsDat$Property %in% props,]
    PropsDF$Frame <- as.integer(PropsDF$Frame)
    PropsDF <- aggregate(PropsDF$Frame, list(PropsDF$Property, PropsDF$Value), \(x) { paste0(min(x), "-", max(x)) })
    colnames(PropsDF) <- c("Property", "Value", "Frames")
    PropsDF <- PropsDF[match(props, PropsDF$Property),]
    PropsDF[, c("Unit", "Name", "Group", "Explanation")] <- PropsDefDat[match(PropsDF$Property, PropsDefDat$PermanentName), c("DisplayDimension", "DisplayName", "DisplayGroupName", "DisplayValueText")]
    w <- which(PropsDF$Explanation != "")
    PropsDF$Explanation <- strsplit(PropsDF$Explanation, ";")
    PropsDF$Explanation[w] <- lapply(PropsDF$Explanation[w], \(x) { #x <- PropsDF$Explanation[w[1L]]
      x <- unlist(x)
      x <- t(as.data.frame(strsplit(x, ":")))
      return(x)
    })
    PropsDF$Value <- as.character(PropsDF$Value)
    PropsDF$Value[w] <- apply(PropsDF[w, c("Value", "Explanation")], 1L, \(x) { x[[2L]][match(x[[1L]], x[[2L]][, 1L]), 2L] })
    g <- lapply(c("Minimum", "Maximum"), \(x) { sub(paste0(" ", x, "$"), "", grep(paste0(" ", x, "$"), PropsDF$Name, value = TRUE)) })
    u <- unique(unlist(g))
    u <- u[vapply(u, \(x) { (x %in% g[[1L]]) & (x %in% g[[2L]]) }, TRUE)]
    if (length(u)) {
      tmp <- PropsDF[PropsDF$Name == paste0(u, " Minimum"),]
      tmp$Name <- sub(" Minimum$", "", tmp$Name)
      tmp$Value <- apply(tmp[, c("Value", "Name", "Frames")], 1L, \(x) {
        paste0(x[[1L]], "-", PropsDF$Value[(PropsDF$Name == paste0(x[[2L]], " Maximum")) & (PropsDF$Frames == x[[3L]])])
      })
      tmp$Property <- gsub(".Minimum", "", tmp$Property)
      for (v in u) {
        wY <- which(PropsDF$Name %in% paste0(v, " ", c("Minimum", "Maximum")))
        if (min(wY) > 1L) { tmp <- rbind(PropsDF[1L:(min(wY)-1L),], tmp) }
        if (max(wY) < nrow(PropsDF)) { tmp <- rbind(tmp, PropsDF[(max(wY)+1L):nrow(PropsDF),]) }
        PropsDF <- tmp
      }
    }
    PropsDF2 <- aggregate(PropsDF[, c("Value", "Frames")], list(PropsDF$Name), list)
    colnames(PropsDF2)[1L] <- "Name"
    PropsDF2[, c("Group", "Unit", "Property")] <- PropsDF[match(PropsDF2$Name, PropsDF$Name), c("Group", "Unit", "Property")]
    PropsDF2 <- PropsDF2[match(PropsDF$Property, PropsDF2$Property),]
    AuF <- unique(PropsDF$Frames)
    PropsDF2$Frame_group <- lapply(PropsDF2$Frames, \(x) { match(x, AuF) })
    PropsDF2$Text <- ""
    PropsDF2$Text <- apply(PropsDF2[, c("Name", "Value", "Unit", "Frame_group")], 1L, \(x) {
      l <- length(x[[2L]])
      res <- NA
      if (l > 1L) {
        x[[2L]] <- .commAnd(x[[2L]])
        x[[4L]] <- paste0("(frame groups ",  .commAnd(x[[4L]]), ", resp.)")
      } else {
        if ((x[[1L]] == "PASEF Current Intens/Previous Intens")&&(x[[2L]] == "Off")) { res <- "" }
      }
      if (is.na(res)) {
        res <- paste0(tolower(x[[1L]]), " = ", x[[2L]])
        if (nchar(x[[3L]])) { res <- paste0(res, " ", x[[3L]]) }
        if (l > 1L) { res <- paste0(res, " ", x[[4L]]) }
      }
      return(res)
    })
    PropsDF2$Text <- sub("^(ims cycle|pasef|set) ", "", PropsDF2$Text)
    PropsDF2$Text <- sub("^number of ", "", PropsDF2$Text)
    PropsDF2$Text <- sub("prepulsestorage ", "pre-pulse storage ", PropsDF2$Text)
    PropsDF2$Text <- sub("current intens/previous intens ", "current/previous intensity ratio ", PropsDF2$Text)
    PropsDF2$Text <- sub("= Off$", "= off", sub("= On$", "= on", PropsDF2$Text))
    PropsDF2$Text <- sub("= Off$", "= off", sub("= On$", "= on", PropsDF2$Text))
    Types <- "TIMS"
    if (grepl("PASEF", PropsDF2$Text[PropsDF2$Name == "Scan Mode"], ignore.case = TRUE)) { Types <- c(Types, "PASEF") }
    Txt <- paste(sapply(Types, \(x) {
      paste0(x, " parameters: ", paste0(PropsDF2$Text[PropsDF2$Group == x], collapse = ", "))
    }), collapse = "; ")
    Txt <- paste0(paste(PropsDF2$Text[!PropsDF2$Group %in% Types], collapse = ", "), "; ", Txt)
    if (length(AuF) > 1L) {
      AuF <- apply(data.frame(N = 1L:length(AuF), Fr = AuF), 1L, paste, collapse = " = ")
      Txt <- paste0("Frame ranges: ", AuF, "; ", Txt)
    }
    DIAtst <- grepl("dia", PropsDF2$Text[PropsDF2$Name == "Scan Mode"], ignore.case = TRUE)
    PRMtst <- grepl("prm", PropsDF2$Text[PropsDF2$Name == "Scan Mode"], ignore.case = TRUE)
    if ((DIAtst == 0L) && is.logical(FrmDat) && is.na(FrmDat)) {
      Txt <- paste0("Isolation: ", isoWdths$Width[1L], " Th for target precursors up to ",  isoWdths$MZ[1L],
                    " Th, then linear increase up to a maximum of ",
                    isoWdths$Width[2L], " Th at ",  isoWdths$MZ[2L], " Th; ",
                    Txt)
    }
    Txt <- paste0("MS method: M/Z range = ", paste(MssRng, collapse = "-"), " Th, ion mobility range = ", paste(MobRng, collapse = "-"), " 1/K0; ", Txt, ".")
    Txt <- gsub("cm2", "cm²", Txt)
    if (DIAtst) {
      DIAfl <- sub("\\.m$", " - DIA_windows.csv", fl)
      message(paste0("Writing DIA windows to ", DIAfl))
      write.csv(DIADat1, DIAfl, row.names = FALSE)
      #plot <- ggplot(DIADat1, aes(xmin = IsolationMz-IsolationWidth/2, xmax = IsolationMz+IsolationWidth/2,
      #                            ymin = ScanNumBegin, ymax = ScanNumEnd, colour = WindowGroup)) + theme_bw() +
      #  geom_rect(alpha = 0.1) + scale_y_reverse() + xlab("M/Z") + ylab("Scan number")
      #poplot(plot)
      tmp1 <- paste0("scan ranges: ", paste(apply(DIADat1[, c("ScanNumBegin", "ScanNumEnd")], 1L, paste, collapse = "-"), collapse = "/"))
      tmp2 <- paste0("isolation M/Z: ", paste(round(DIADat1$IsolationMz, 3L), collapse = "/")) # That is already wayyyyyy more precise than we need!!!
      tmp3 <- paste0("collision energy: ", paste(DIADat1$CollisionEnergy, collapse = "/"))
      if (length(unique(DIADat1$IsolationWidth)) == 1L) {
        DIATxt <- paste0("DIA windows scheme: fixed width = ", unique(DIADat1$IsolationWidth), ", ", tmp1, ", ", tmp2, ", ", tmp3, ".")
      } else {
        tmp4 <- paste0("isolation width: ", paste(DIADat1$IsolationWidth, collapse = "/"))
        DIATxt <- paste0("DIA windows scheme: ", tmp1, ", ", tmp2, ", ", tmp4, ", ", tmp3, ".")
      }
      Txt <- paste0(Txt, " ", DIATxt)
    } else { DIADat1 <- DIAfl <- NA }
    if (PRMtst) {
      PRMfl <- sub("\\.m$", " - PRM_inclusion_list.csv", fl)
      message(paste0("Writing PRM inclusion list to ", PRMfl))
      write.csv(PRMDat, PRMfl, row.names = FALSE)
      Txt <- paste0(Txt, " PRM INCLUSION LIST EXPORTED AT \"", PRMfl, "\"")
    } else { PRMDat <- PRMfl <- NA }
    res <- list("Instrument" = INSTR, "Vendor" = INSTRMAKER, "M/Z range" = MssRng, "IM range" = MobRng, "Isolation widths" = isoWdths,
                "Properties" = PropsDF2, "Text" = Txt,
                "DIA?" = DIAtst, "DIA windows" = DIADat1,
                "PRM?" = PRMtst, "PRM inclusion list" = PRMDat, "PRM list path" = PRMfl)
    return(res)
  }), FILES[w])
  print(sapply(names(MS_method_parsed), \(fl) { MS_method_parsed[[fl]]$`M/Z range` }))
  print(sapply(names(MS_method_parsed), \(fl) { MS_method_parsed[[fl]]$`IM range` }))
  print(sapply(names(MS_method_parsed), \(fl) { MS_method_parsed[[fl]]$Text }))
  writeClipboard(sapply(names(MS_method_parsed), \(fl) { MS_method_parsed[[fl]]$Text }))
  #a <- MS_meth[[1L]]$Method$Properties
  #b <- MS_meth[[1L]]$Method$Definitions
  #a$PermanentName <- b$PermanentName[match(a$Property, b$Id)]
  #grep("Polygon", b$PermanentName, value = TRUE, ignore.case = TRUE)
  #grep("Range", b$PermanentName, value = TRUE, ignore.case = TRUE)
  #grep("Mass", b$PermanentName, value = TRUE, ignore.case = TRUE)
  #grep("MSMS", b$PermanentName, value = TRUE, ignore.case = TRUE)
  #a$Value[a$PermanentName == "MSMSAuto_ListPreferMassStart"][1L]
  #a$Value[a$PermanentName == "Quadrupole_IsolationMass_Set"][1L]
  #a$Value[a$PermanentName == "MSMS_Pasef_Masses"][1L]
  #unique(a$Value[a$PermanentName == "MSMS_ListCollisionEnergy"])
  #unique(a$Value[a$PermanentName == "MSMSAuto_ListExcludeMassStart"])
  #sort(unique(a$Value[a$PermanentName == "MSMS_Pasef_Masses"]))
  #sort(unique(a$Value[a$PermanentName == "MSMSManual_ListIsolationMass"]))
  #grep("Quadrupole", b$PermanentName, value = TRUE)
  #w <- grep("^IMS_PolygonFilter", a$PermanentName)
  #tst <- aggregate(a$Value[w], list(a$PermanentName[w]), unique)
  #grep("Mass", b$PermanentName, value = TRUE, ignore.case = TRUE)
  #w <- which(a$PermanentName == "IMS_PolygonFilter_Type")
  #tst <- aggregate(a$Value[w], list(a$PermanentName[w]), unique)
  #tst$x
  #length(unique(a$Value))
}
