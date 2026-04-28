#' MatMet_LCMS
#'
#' @description
#' A function to write the LC-MS part of Materials & Methods, pulling information from the methods file.
#' 
#' Note to self:
#' !!! When editing function, always save with UTF-8 encoding!!!
#' 
#' @param RawFiles The raw files from which to attempt to extract the method.
#' @param LocalRoot What is the root of the folder structure in the local folder. If the raw files are not local but were moved to an archive, this root will be replaced by that of the archive, see ArchiveRoot. Default = "D:/groups_temp/"
#' @param ArchiveRoot If the raw files are not local but were moved to an archive, where would that be? Default = "Q:/MS/Acquired data/"
#' @param ScanHdsMnLoc ScanHeadsman path. Default = "C:/ScanHeadsman-1.2.20200730"; so far for us 1.2 works, 1.3 does not.
#' @param Columns Path to a table of LC columns to pick from.
#' @param N.clust A limit on the number of vCPUs to use. If left as NULL (default), uses the number of available clusters - 1, to a minimum of 1.
#' @param N.reserved Default = 1. Number of reserved vCPUs the function is not to use. Note that for obvious reasons it will always use at least one.
#' @param cl Already have a cluster handy? Why spend time making a new one, which on top of that may invalidate the old one. Just pass it along!
#' 
#' @details
#' Only tested on Bruker timsTOF HT with nanoElute 2 and Thermo Q-Exactive HT with Ultimate 3000 or Vanquish, so will probably fail on different machines from the same vendors, and obviously also on machines from other vendors.
#' 
#' @returns
#' A list with two entries:
#' - "Text": the text for each LCMS method
#' - "Instruments": all LC and MS instrument models.
#' 
#' @examples
#' MatMet_Text_list <- MatMet_LCMS()
#' 
#' @export

MatMet_LCMS <- function(ScanHdsMnLoc = "C:/ScanHeadsman-1.2.20200730", # Should eventually be shifted to a value in the default locations table
                        LocalRoot = "...Search_Folder/",
                        ArchiveRoot = "...Archive/Acquired data/",
                        RawFiles = rawFiles,
                        Columns = "default",
                        N.clust,
                        N.reserved = 1L,
                        cl) {
  TESTING <- FALSE
  #DefArg(MatMet_LCMS);TESTING <- TRUE
  #RawFiles <- rawFiles <-...
  #
  misFun <- if (TESTING) {
    # Note:
    # This is not a perfect alternative to missing but will work in most cases, unless x matches a function imported by a package 
    function(x) { return(!exists(deparse(substitute(x)))) }
  } else { missing }
  #
  # Create cluster
  # Create cluster
  stopCl <- FALSE
  if ((is.null(cl))||(!inherits(cl, "cluster"))) {
    dc <- parallel::detectCores()
    if (misFun(N.reserved)) { N.reserved <- 1L }
    nMax <- max(c(dc - N.reserved, 1L))
    if (misFun(N.clust)) { N.clust <- nMax } else {
      if (N.clust > nMax) {
        warning("More cores specified than allowed, I will ignore the specified number! You should always leave at least one free for other processes, see the \"N.reserved\" argument.")
        N.clust <- nMax
      }
    }
    cat("     Making fresh cluster...\n")
    cl <- parallel::makeCluster(N.clust, type = "SOCK")
    stopCl <- TRUE
  }
  N.clust <- length(cl)
  #
  Moult <- length(RawFiles) > 1L
  #
  # LC columns
  colChar <- c("Name",
               "Class",
               "Vendor",
               "Length (cm)",
               "ID (µm)",
               "Particles size (µm)",
               "Pore size (Å)",
               "Material",
               "Type",
               "Description",
               "P/N",
               "Function")
  colTst <- (Columns != "default")&(file.exists(Columns))
  libPath <- as.data.frame(library()$results)
  libPath <- normalizePath(libPath$LibPath[match("proteoCraft", libPath$Package)], winslash = "/")
  proteoPath <- paste0(libPath, "/proteoCraft")
  myPath <- pkgPath <- paste0(normalizePath(Sys.getenv("HOME"), winslash = "/"), "/R/proteoCraft")
  if (!colTst) {
    Columns <- paste0(pkgPath, "/LC_columns.xlsx")
    if (!file.exists(Columns)) {
      dfltLocFl <- paste0(pkgPath, "/Default_locations.xlsx")
      dfltLoc <- openxlsx2::read_xlsx(dfltLocFl)
      tmpPath <- dfltLoc$Path[match("Temporary folder", dfltLoc$Folder)]
      myPath <- tmpPath
      Columns <- paste0(tmpPath, "/LC_columns.xlsx")
    }
    if (!file.exists(Columns)) {
      Columns <- paste0(proteoPath, "/extdata/LC_columns.xlsx")
    }
  }
  Columns2 <- paste0(myPath, "/LC_columns.xlsx")
  #cat(Columns, "\n")
  #cat(Columns2, "\n")
  colTst <- file.exists(Columns)
  if (colTst) {
    allKolumns <- openxlsx2::read_xlsx(Columns)
  } else {
    # Should never happen now that we distribute "LC_columns.xlsx" with the package
    Kolumns <- c("EasySpray C18 column (2 µm particle size, 75 µm ID x 50 cm length, ThermoFisher Scientific P/N ES903)",
                 "200 cm C18 µPAC column (micro-Pillar Array Column, PharmaFluidics P/N 5525031518210B)",
                 "50 cm C18 µPAC GEN2 column (2nd generation micro-Pillar Array Column, PharmaFluidics P/N COL-nano050G2B)")
    preKolumns <- c("Acclaim PepMap C18 pre-column (5 µm particle size, 0.3 mm ID x 5 mm length, ThermoFisher Scientific P/N 160454)",
                    "C18 µPAC trapping-column (PharmaFluidics P/N 55250200018001)")
    allKolumns <- data.frame(Name = c(Kolumns, preKolumns),
                             Class = c(Kolumns, preKolumns),
                             Material = "C18",
                             "Vendor" = c("ThermoFisher Scientific", "PharmaFluidics", "PharmaFluidics", "ThermoFisher Scientific", "PharmaFluidics"),
                             "Length (cm)" = c(50, 200, 50, 0.5, NA),
                             "ID (µm)" = c(75, NA, NA, 300, NA),
                             "Particles size (µm)" = c(2, NA, NA, 5, NA),
                             "Pore size (Å)" = NA,
                             "P/N" = c("ES903", "5525031518210B", "COL-nano050G2B", "160454", "55250200018001"),
                             Type = c("Analytical", "Analytical", "Analytical", "Trap",  "Trap"),
                             Description = "",
                             check.names = FALSE)
  }
  if (!"Description" %in% colnames(allKolumns)) { allKolumns$Description <- "" }
  kolDescr <- function(Colonnes) {
    apply(Colonnes[, c("Name", "Material", "Length (cm)", "ID (µm)", "Particles size (µm)", "Vendor", "P/N"),
                   drop = FALSE], 1L, \(x) {
      #x <- Colonnes[1, c("Name", "Material", "Length (cm)", "ID (µm)", "Particles size (µm)", "Vendor", "P/N")]
      x <- gsub("^ +| +$", "", as.character(unlist(x)))
      dimz <- c(x[3L], x[4L])
      w <- which(!is.na(dimz))
      dimzTst <- length(w) > 0
      dimz <- if (dimzTst) { paste(paste0(gsub("\\.0+$", "", dimz), c(" cm", " µm ID"))[w], collapse = " * ") } else { NA }
      pn <- c(x[6L], x[7L])
      w <- which(!is.na(pn))
      pnTst <- length(w) > 0
      pn <- if (pnTst) { paste(paste0(c("", "P/N "), pn)[w], collapse = " ") } else { NA }
      pnTst <- !is.na(x[7L])
      if (pnTst) { x[7L] <- paste0("P/N ", x[7L]) }
      prtcls <- c(x[5L], x[2L])
      w <- which(!is.na(prtcls))
      prtclsTst <- length(w) > 0L
      prtcls <- if (prtclsTst) { paste(paste0(prtcls, c(" µm", "-coated particles"))[w], collapse = " ") } else { NA }
      res <- paste0(x[1L], " column (", paste(c(prtcls, dimz, pn)[which(c(prtclsTst, dimzTst, pnTst))], collapse = ", "), ")")
      res <- gsub("column column", "column", res, ignore.case = TRUE)
      return(res)
    })
  }
  #w <- which((nchar(allKolumns$Description) == 0)|(is.na(allKolumns$Description)))
  w <- 1L:nrow(allKolumns)
  allKolumns$Description[w] <- kolDescr(allKolumns[w,])
  Kolumns <- unique(c(allKolumns$Description[which(allKolumns$Function == "Analytical")], "None (direct infusion)", "Add new..."))
  preKolumns <- unique(c(allKolumns$Description[which(allKolumns$Function == "Trap")],  "None (direct injection)", "Add new..."))
  AddKolKount  <- 0L
  #View(allKolumns)
  #
  MethodsTbl <- data.frame(Raw.file = gsub("\\\\", "/", RawFiles))
  MethodsTbl$Folder <- dirname(RawFiles)
  w <- which(MethodsTbl$Folder == ".")
  if (length(w)) { MethodsTbl$Folder[w] <- getwd() }
  MethodsTbl$Name <- basename(MethodsTbl$Raw.file)
  MethodsTbl$Ext <- vapply(basename(MethodsTbl$Raw.file), \(x) {
    x <- if (!grepl("\\.", x)) { NA } else { rev(unlist(strsplit(x, "\\.")))[1L] }
    return(x)
  }, "")
  w <- which(is.na(MethodsTbl$Ext))
  if (length(w)) {
    l <- (length(w)>1)+1
    warning(paste0(c("One", "Some")[l]," MS file", c("'s", "s'")[l], " extension", c("", "s")[l], " could not be detected, skipping:", paste("\n -> ", MethodsTbl$Raw.file[w], collapse = ""), "\n"))
    MethodsTbl <- MethodsTbl[which(!is.na(MethodsTbl$Ext)),]
  }
  # For now only Thermo raw files and Bruker d files (folders) are supported
  w <- which(!tolower(MethodsTbl$Ext) %in% c("d", "raw"))
  if (length(w)) { 
    l <- (length(w)>1L)+1L
    warning(paste0("Only Thermo .raw or Bruker .d MS files are supported, skipping file", c("", "s")[l], ":",
                   paste("\n -> ", MethodsTbl$Raw.file[w], collapse = ""), "\n"))
    MethodsTbl <- MethodsTbl[which(MethodsTbl$Ext %in% c("d", "raw")),]
  }
  #
  Ext <- unique(MethodsTbl$Ext)
  if (length(Ext) > 1L) { warning("Were you aware that this dataset is a mixture of Thermo .raw files and Bruker .d folders?") }
  # Vector of material and methods text for each LC/MS method
  all_LCMS_txts <- c()
  # List of all LC and MS instruments used in the dataset, without connection to methods
  all_LCMS_instr <- list(LC = c(),
                         MS = c())
  # Bruker .d folders
  if ("d" %in% tolower(Ext)) {
    library(XML)
    library(RSQLite)
    ADflt <- "MS-grade H₂O + 0.1% formic acid"
    BDflt <- "100% acetonitrile + 0.1% formic acid"
    BrMeth <- MethodsTbl[which(tolower(MethodsTbl$Ext) == "d"),]
    BrMeth$FlInfo <- paste0(BrMeth$Raw.file, "/SampleInfo.xml")
    w <- which(file.exists(BrMeth$FlInfo))
    BrMeth[, c("LC_method_fl", "MS_method_fl")] <- NA
    BrMeth[w, c("LC_method_fl", "MS_method_fl")] <- as.data.frame(t(as.data.frame(sapply(BrMeth$FlInfo[w], \(fl) {
      #fl <- BrMeth$FlInfo[w][1L]
      x <- xmlToList(fl)
      res <- if (".attrs" %in% names(x$Sample)) { c(x$Sample$.attrs["Method"], x$Sample$.attrs["MS_Method"]) } else {
        c(x$Sample["Method"], x$Sample["MS_Method"])
      }
      return(gsub("\\?.+", "", gsub("\\\\", "/", res)))
    }))))
    # I dislike .d files because they are folders...
    # ... but...
    # ... I like .d folders because they contain a human-readable method for LC!
    # ... and a parsable (though not directly human-readable) method for MS, with all (but one!) parameters of interest.
    #
    # Another thing: those files embed MS methods but these are not the original theoretical method, those are realised methods.
    # Thus they will be slightly different for files run with the same method. This is why it makes sense to check method unicity using the path to the original method,
    # not the method as embedded in the file.
    #
    #
    BrMeth$LCMS_method_fls <- NA
    BrMeth$LCMS_method_fls[w] <- apply(BrMeth[w, c("LC_method_fl", "MS_method_fl")], 1L, \(x) {
      list(LC = x[[1L]],
           MS = x[[2L]])
    })
    tmp <- unique(BrMeth$LCMS_method_fls)
    BrRoots <- "Bruker"
    if (length(tmp) > 1L) { BrRoots <- paste0(BrRoots, " meth#", 1L:length(tmp)) }
    uBrMeth <- data.frame(MethodID = BrRoots,
                          Column = NA, Trap = NA, SolvA = NA, SolvB = NA)
    uBrMeth$Method <- tmp
    BrMeth$MethodID <- sapply(BrMeth$LCMS_method_fls, \(x) {
      uBrMeth$MethodID[which(vapply(uBrMeth$Method, \(y) { identical(x, y) }, TRUE))]
    })
    if (length(tmp) > 1L) {
      tmp2 <- paste(sapply(uBrMeth$MethodID, \(x) { paste(c(paste0(" - ", x, " files:"),
                                                            paste0("   > ", BrMeth$Raw.file[which(BrMeth$MethodID == x)])), collapse = "\n") }),
                   collapse = "\n\n")
      msg <- paste0(length(tmp), " different Bruker LC-MS/MS acquisition methods detected:\n\n", tmp2, "\n")
      cat(msg)
    }
    # Get LC methods
    BrMeth$Run_ID <- as.integer(gsub(".*_|\\.d$", "", BrMeth$Raw.file))
    BrMeth$LC_method <- paste0(BrMeth$Raw.file, "/", as.character(BrMeth$Run_ID), ".m/hystar.method")
    BrMeth$LC_method.exists <- file.exists(BrMeth$LC_method)
    w <- which(BrMeth$LC_method.exists)
    uBrMeth$LC_method <- BrMeth$LC_method[w][match(uBrMeth$MethodID, BrMeth$MethodID[w])]
    uBrMeth$LC_meth <- lapply(uBrMeth$LC_method, \(fl) { #fl <- uBrMeth$LC_method[1L]
      if (is.na(fl)) { return() }
      x <- XML::xmlToList(fl)
      lc <- x$LCMethodData$ModuleMethods$ModuleMethodData$text
      lc <- XML::xmlToList(lc)
      return(lc)
      #lc$ModuleMethodData$Method$AdvancedSettings
    })
    # Get MS method
    BrMeth$MS_method <- paste0(BrMeth$Raw.file, "/analysis.tdf")
    BrMeth$MS_method.exists <- file.exists(BrMeth$MS_method)
    w <- which(BrMeth$MS_method.exists)
    uBrMeth$MS_method <- BrMeth$MS_method[w][match(uBrMeth$MethodID, BrMeth$MethodID[w])]
    uBrMeth$MS_meth <- lapply(uBrMeth$MS_method, \(fl) { #fl <- BrMeth$MS_method[w[1L]]
      #x <- readLines(fl)
      #x[1L:100L]
      SQltTDF <- dbConnect(drv = RSQLite::SQLite(), dbname = fl)
      #nms <- dbListTables(SQltTDF)
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
                  Name = GlobMetDat$Value[which(GlobMetDat$Key == "MethodName")]))
    })
    # Put them together
    uBrMeth$Method <- apply(uBrMeth[, c("LC_meth", "MS_meth")], 1L, \(x) {
      list(LC = x[[1L]],
           MS = x[[2L]]) })
    # Parse them and create text
    uBrMeth$Samples <- lapply(uBrMeth$MethodID, \(x) { BrMeth$Raw.file[which(BrMeth$MethodID == x)] })
    uBrMeth$NANO <- gsub(" +", " ", gsub("Bruker *", "", gsub("\t.*", "",
                                                              vapply(uBrMeth$Method, \(x) { x$LC$HyStarMethodData$ModuleName }, "")
                                                              ), ignore.case = TRUE))
    w <- which(uBrMeth$NANO == "nanoElute")
    uBrMeth$Root <- apply(uBrMeth[, c("MethodID", "NANO")], 1L, paste, collapse = " - ")
    opt <- paste0(1L:2L, paste(rep(" ", 200L), collapse = ""))
    for (i in w) {
      uBrMeth$NANO[i] <- paste0("nanoElute",
                                c("", " 2"))[match(gsub(" ", "",
                                                        svDialogs::dlg_list(opt, opt[2L],
                                                                            title = paste0(uBrMeth$Root[i], " - confirm LC version:"))$res),
                                                   as.character(1L:2L))]
    }
    uBrMeth$Root <- apply(uBrMeth[, c("MethodID", "NANO")], 1L, paste, collapse = " - ")
    uBrMeth$Column <- vapply(uBrMeth$Method, \(x) { x$LC$ModuleMethodData$Method$SeparatorName }, "")
    uBrMeth$usesTrap <- as.logical(toupper(vapply(uBrMeth$Method, \(x) { x$LC$ModuleMethodData$Method$UsesTrapColumn }, "")))
    w <- which(uBrMeth$usesTrap)
    uBrMeth$Trap[w] <- vapply(uBrMeth$Method[w], \(x) { x$LC$ModuleMethodData$Method$TrapName }, "")
    uBrMeth$OvenT <- sapply(uBrMeth$Method, \(x) {
      c(NA, x$LC$ModuleMethodData$Method$OvenTemperature)[as.logical(toupper(x$LC$ModuleMethodData$Method$IsSetTemperature))+1L]
    })
    uBrMeth$Gradient <- lapply(uBrMeth$Method, \(x) {
      nm <- names(x$LC$ModuleMethodData$Method$Gradient)
      grd <- sapply(1L:length(nm), \(y) { x$LC$ModuleMethodData$Method$Gradient[[y]] })
      grd <- t(as.data.frame(grd))
      rownames(grd) <- NULL
      return(grd)
    })
    uBrMeth$Flow <- sapply(uBrMeth$Method, \(x) { x$LC$HyStarMethodData$MainFlow })
    uBrMeth$RunTime <- sapply(uBrMeth$Method, \(x) { x$LC$HyStarMethodData$NoStandardMethodData })
    for (i in 1L:nrow(uBrMeth)) {
      kol <- paste0("\" ", uBrMeth$Column[i], " \"")
      kol <- svDialogs::dlg_list(unique(c(kol, Kolumns)), kol,
                                 title = paste0(uBrMeth$Root[i], " - select column used:"))$res
      if (kol == "Add new...") {
        colCharDf <- data.frame(Description = "Analytical")
        for (kk in colChar[which(colChar != "Function")]) {
          colCharDf[[kk]] <- svDialogs::dlg_input(paste0("Enter value for field \"", kk, "\""), "")$res
        }
        kol <- colCharDf$Name
        allKolumns <- plyr::rbind.fill(allKolumns, colCharDf)
        rownames(allKolumns) <- NULL
        tst <- apply(allKolumns, 1L, paste, collapse = "|")
        tst <- aggregate(1L:length(tst), list(tst), min)
        allKolumns <- allKolumns[sort(tst$x),]
        Kolumns <- unique(c(allKolumns$Description[which(allKolumns$Function == "Analytical")],
                            "None (direct infusion)",
                            "Add new..."))
        AddKolKount <- AddKolKount + 1L
        kol <- tmp$Description
      }
      uBrMeth$Column[i] <- kol
      if (uBrMeth$usesTrap[i]) {
        kol <- paste0("\" ", uBrMeth$Trap[i], " \"")
        kol <- svDialogs::dlg_list(unique(c(kol, preKolumns)), kol,
                                   title = paste0(uBrMeth$Root[i],  " - select trap used:"))$res
        if (kol == "Add new...") {
          colCharDf <- data.frame(Description = "Trap")
          for (kk in colChar[which(colChar != "Function")]) {
            colCharDf[[kk]] <- svDialogs::dlg_input(paste0("Enter value for field \"", kk, "\""), "")$res
          }
          kol <- colCharDf$Name
          allKolumns <- plyr::rbind.fill(allKolumns, colCharDf)
          rownames(allKolumns) <- NULL
          tst <- apply(allKolumns, 1L, paste, collapse = "|")
          tst <- aggregate(1L:length(tst), list(tst), min)
          allKolumns <- allKolumns[sort(tst$x),]
          preKolumns <- unique(c(allKolumns$Description[which(allKolumns$Function == "Trap")],
                                 "None (direct injection)",
                                 "Add new..."))
          AddKolKount <- AddKolKount + 1L
          kol <- tmp$Description
        }
        uBrMeth$Trap[i] <- kol
        uBrMeth$SolvA[i] <- svDialogs::dlg_input(paste0(uBrMeth$Root[i],
                                                        " confirm solvent A composition:"),
                                                 ADflt)$res
        uBrMeth$SolvB[i] <- svDialogs::dlg_input(paste0(uBrMeth$Root[i],
                                                        " confirm solvent B composition:"),
                                                 BDflt)$res
      }
    }
    uBrMeth$Par <- lapply(uBrMeth$Method, \(x) {
      nm <- names(x$LC$ModuleMethodData$Metho$AdvancedSettings$Parameters)
      pr <- sapply(1L:length(nm), \(y) { x$LC$ModuleMethodData$Metho$AdvancedSettings$Parameters[[y]] })
      colnames(pr) <- pr[1L,]
      pr <- pr[2L:nrow(pr),]
      return(as.data.frame(pr))
    })
    uBrMeth$MS_method_parsed <- lapply(uBrMeth$MS_meth, \(x) { #x <- uBrMeth$MS_meth[[1L]]
      fl <- paste0(getwd(), "/", x$Name)
      GlobMetDat <- x$Method$Global
      PropsDat <- x$Method$Properties
      PropsDefDat <- x$Method$Definitions
      FrmDat <- x$Method$Frames
      DIADat1 <- x$Method$DIA
      PRMDat <- x$Method$PRM
      INSTR <- GlobMetDat$Value[which(GlobMetDat$Key == "InstrumentName")]
      INSTRMAKER <- GlobMetDat$Value[which(GlobMetDat$Key == "InstrumentVendor")]
      if (INSTRMAKER == "Bruker") { INSTRMAKER <- "Bruker Daltonics"}
      #View(GlobMetDat)
      #View(PropsDat)
      PropsDat$Property <- PropsDefDat$PermanentName[match(PropsDat$Property, PropsDefDat$Id)]
      #View(PropsDefDat)
      #View(FrmDat)
      MssRng <- gsub("\\.$", "", gsub("0+$", "", GlobMetDat$Value[match(c("MzAcqRangeLower", "MzAcqRangeUpper"), GlobMetDat$Key)]))
      MobRng <- gsub("\\.$", "", gsub("0+$", "", GlobMetDat$Value[match(c("OneOverK0AcqRangeLower", "OneOverK0AcqRangeUpper"), GlobMetDat$Key)]))
      if (!is.logical(FrmDat)) {
        isoWdths <- data.frame(Width = c(min(FrmDat$IsolationWidth), max(FrmDat$IsolationWidth)))
        isoWdths$MZ <- round(c(max(FrmDat$IsolationMz[which(FrmDat$IsolationWidth == isoWdths$Width[1L])]),
                               min(FrmDat$IsolationMz[which(FrmDat$IsolationWidth == isoWdths$Width[2L])])))
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
      PropsDF <- PropsDat[which(PropsDat$Property %in% props),]
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
      g <- lapply(c("Minimum", "Maximum"), \(x) { gsub(paste0(" ", x, "$"), "", grep(paste0(" ", x, "$"), PropsDF$Name, value = TRUE)) })
      u <- unique(unlist(g))
      u <- u[which(sapply(u, \(x) { x %in% g[[1L]] & x %in% g[[2L]] }))]
      if (length(u)) {
        tmp <- PropsDF[which(PropsDF$Name == paste0(u, " Minimum")),]
        tmp$Name <- gsub(" Minimum$", "", tmp$Name)
        tmp$Value <- apply(tmp[, c("Value", "Name", "Frames")], 1L, \(x) {
          paste0(x[[1L]], "-", PropsDF$Value[which((PropsDF$Name == paste0(x[[2L]], " Maximum"))&(PropsDF$Frames == x[[3L]]))])
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
          x[[2L]] <- paste0(paste(x[[2L]][1L:(l-1L)], collapse = ", "), " and ", x[[2L]][l])
          x[[4L]] <- paste0("(frame groups ",  paste(x[[4L]][1L:(l-1L)], collapse = ", "), " and ",
                            x[[4L]][l], ", resp.)")
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
      PropsDF2$Text <- gsub("^ims cycle |^pasef |^set ", "", PropsDF2$Text)
      PropsDF2$Text <- gsub("^number of ", "", PropsDF2$Text)
      PropsDF2$Text <- gsub("prepulsestorage ", "pre-pulse storage ", PropsDF2$Text)
      PropsDF2$Text <- gsub("current intens/previous intens ", "current/previous intensity ratio ", PropsDF2$Text)
      PropsDF2$Text <- gsub("= Off$", "= off", gsub("= On$", "= on", PropsDF2$Text))
      PropsDF2$Text <- gsub("= Off$", "= off", gsub("= On$", "= on", PropsDF2$Text))
      Types <- "TIMS"
      if (grepl("PASEF", PropsDF2$Text[which(PropsDF2$Name == "Scan Mode")], ignore.case = TRUE)) { Types <- c(Types, "PASEF") }
      Txt <- paste(sapply(Types, \(x) {
        paste0(x, " parameters: ", paste0(PropsDF2$Text[which(PropsDF2$Group == x)], collapse = ", "))
      }), collapse = "; ")
      Txt <- paste0(paste(PropsDF2$Text[which(!PropsDF2$Group %in% Types)], collapse = ", "), "; ", Txt)
      if (length(AuF) > 1L) {
        AuF <- apply(data.frame(N = 1L:length(AuF), Fr = AuF), 1L, paste, collapse = " = ")
        Txt <- paste0("Frame ranges: ", AuF, "; ", Txt)
      }
      DIAtst <- grepl("dia", PropsDF2$Text[which(PropsDF2$Name == "Scan Mode")], ignore.case = TRUE)
      PRMtst <- grepl("prm", PropsDF2$Text[which(PropsDF2$Name == "Scan Mode")], ignore.case = TRUE)
      if ((!DIAtst)&&(is.logical(FrmDat))&&(is.na(FrmDat))) {
        Txt <- paste0("Isolation: ", isoWdths$Width[1L], " Th for target precursors up to ",
                      isoWdths$MZ[1L], " Th, then linear increase up to a maximum of ",
                      isoWdths$Width[2L], " Th at ",  isoWdths$MZ[2L], " Th; ",
                      Txt)
      }
      Txt <- paste0("MS method: M/Z range = ", paste(MssRng, collapse = "-"), " Th, ion mobility range = ", paste(MobRng, collapse = "-"), " 1/K0; ", Txt, ".")
      Txt <- gsub("cm2", "cm²", Txt)
      if (DIAtst) {
        DIAfl <- gsub("\\.m$", " - DIA_windows.csv", fl)
        cat("Writing DIA windows to", DIAfl, "\n")
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
        PRMfl <- gsub("\\.m$", " - PRM_inclusion_list.csv", fl)
        cat("Writing PRM inclusion list to", PRMfl, "\n")
        write.csv(PRMDat, PRMfl, row.names = FALSE)
        Txt <- paste0(Txt, " PRM INCLUSION LIST EXPORTED AT \"", PRMfl, "\"")
      } else { PRMDat <- PRMfl <- NA }
      res <- list("Instrument" = INSTR, "Vendor" = INSTRMAKER, "M/Z range" = MssRng, "IM range" = MobRng, "Isolation widths" = isoWdths,
                  "Properties" = PropsDF2, "Text" = Txt,
                  "DIA?" = DIAtst, "DIA windows" = DIADat1,
                  "PRM?" = PRMtst, "PRM inclusion list" = PRMDat, "PRM list path" = PRMfl)
      return(list(Method = res,
                  MS = INSTR,
                  MS_vendor = INSTRMAKER))
    })
    uBrMeth$MS_txt <- vapply(uBrMeth$MS_method_parsed, \(x) { x$Method$Text }, "")
    uBrMeth$MS <- vapply(uBrMeth$MS_method_parsed, \(x) { x$MS }, "")
    kol <- c("Samples", "NANO", "Column", "Trap", "usesTrap", "OvenT", "Gradient", "Flow", "RunTime", "SolvA", "SolvB", "MS_method_parsed")
    #kol %in% colnames(uBrMeth)
    uBrMeth$LC_method_parsed  <- apply(uBrMeth[, kol], 1L, \(x) {
      #x <- uBrMeth[1, c("Samples", "NANO", "Column", "Trap", "usesTrap", "OvenT", "Gradient", "Flow", "RunTime", "SolvA", "SolvB", "MS_method_parsed")]
      #x <- uBrMeth[1, kol]
      NANO <- x[[2L]]
      Moult <- length(x[[1L]]) > 1L
      INSTR <- x[[12L]]$Instrument
      if (is.null(INSTR)) { INSTR <- x[[12L]]$Method$Instrument }
      if (is.null(INSTR)) { INSTR <- x[[12L]]$MS }
      INSTRMAKER <- x[[12L]]$Vendor
      if (is.null(INSTRMAKER)) { INSTRMAKER <- x[[12L]]$Method$Vendor }
      NANOMAKER <- "unknown manufacturer"
      if (NANO %in% c("RSLC Nano", "Vanquish NEO")) { NANOMAKER <- "ThermoFisher Scientific" }
      if (NANO %in% c(paste0("nanoElute", c("", paste0(" ", 1L:2L))),
                      paste0("Bruker nanoElute", c("", paste0(" ", 1L:2L))))) { NANOMAKER <- "Bruker Daltonics" }
      tstNano <- tolower(substring(NANO, 1L, 1L)) %in% c("a", "e", "i", "o", "u")
      tstInstr <- tolower(substring(INSTR, 1L, 1L)) %in% c("a", "e", "i", "o", "u")
      COLUMN <- x[[3L]]
      tstCol <- tolower(substring(COLUMN, 1L, 1L)) %in% c("a", "e", "i", "o", "u")
      USESTRAP <- x[[5L]]
      PRECOLUMN <- x[[4L]]
      ovenTEMP <- x[[6L]]
      tstPreCol <- tolower(substring(PRECOLUMN, 1L, 1L)) %in% c("a", "e", "i", "o", "u")
      SolvA <- x[[10L]]
      SolvB <- x[[11L]]
      FLOWVALUE <- x[[8L]]
      Grad <- as.data.frame(x[[7L]])
      Grad$Time <- as.numeric(Grad$Time)/60
      Grad$Mix <- round(as.numeric(Grad$Mix)*100, 2L)
      GRADLENGTH <- max(Grad$Time)
      Grad.Baseline <- min(Grad$Mix)
      Grad.Plateau <- max(Grad$Mix)
      Grad.PlateauL <- round(Grad$Time[rev(which(Grad$Mix == Grad.Plateau))[1L]] - Grad$Time[which(Grad$Mix == Grad.Plateau)[1L]], 0)
      GradStrt <- Grad$Time[which(Grad$Mix == Grad.Baseline)[1L]]
      tmp1 <- Grad$Time[which(Grad$Mix == Grad.Plateau)[1L]]
      tmp0 <- Grad$Time[which(Grad$Time == tmp1)-1L]
      GradEnd <- c(tmp0, tmp1)[(tmp1 - tmp0 > 1L)+1L]
      GRADLENGTH <- GradEnd-GradStrt
      Grad <- Grad[which(Grad$Time == GradStrt)[1L]:which(Grad$Time == GradEnd),]
      Grad <- paste(apply(Grad[, c("Time", "Mix")], 1L, \(x) { paste0(x[[1L]], " min, ", x[[2L]], "%") }), collapse = "; ")
      Grad <- paste0(Grad, ", followed immediately by a ", Grad.PlateauL, " min plateau at ", Grad.Plateau, "%")
      LC_txt <- paste0(c("The s", "S")[Moult+1L], "ample", c(" was", "s were")[Moult+1L], " analyzed by LC-MS/MS on a",
                     c("", "n")[tstNano+1L], " ", NANO, " nano-HPLC (",
                     NANOMAKER, ") coupled with a", c("", "n")[tstInstr+1L], " ", INSTR, " (", INSTRMAKER, ")")
      if (COLUMN != "None (direct infusion)") {
        LC_txt <- if (USESTRAP) {
          paste0(LC_txt, ", concentrated over a", c("", "n")[tstPreCol+1L], " ", PRECOLUMN, ", then ")
        } else { paste0(LC_txt, ". ", c("The s", "S")[Moult2+1L], "ample", c(" was", "s were")[Moult2+1L], " ") }
        LC_txt <- paste0(LC_txt, "bound to a", c("", "n")[tstCol+1L], " ", COLUMN, " heated at ", ovenTEMP, "°C and eluted over the following ", GRADLENGTH,
                         " min gradient: solvent A, ", SolvA, "; ", "solvent B, ", SolvB, "; ",
                         "constant ", FLOWVALUE, " nL/min flow; ",
                         "B percentage: ", Grad, ".")
      } else { LC_txt <- paste0(LC_txt, " by direct infusion.") }
      return(list(Method = LC_txt,
                  LC = NANO,
                  LC_vendor = NANOMAKER))
    })
    uBrMeth$LC_txt <- vapply(uBrMeth$LC_method_parsed, \(x) { x$Method }, "")
    uBrMeth$LC <- vapply(uBrMeth$LC_method_parsed, \(x) { x$LC }, "")
    tmp <- uBrMeth[, c("LC_txt", "MS_txt")]
    #tmp <- tmp[c(1L, 2L, 2L),]
    #tmp$LC_txt[2L] <- tmp$LC_txt[1L]
    #tmp$MS_txt[3L] <- tmp$MS_txt[1L]
    tst <- vapply(c("LC_txt", "MS_txt"), \(x) { length(unique(tmp[[x]])) }, 1L)
    for (i in c("LC", "MS")) {
      k <- paste0(i, "_txt")
      if (tst[[k]] > 1L) {
        m <- match(tmp[[k]], unique(tmp[[k]]))
        tmp[[k]] <- paste0(i, " method #", m, ": ", tmp[[k]])
      }
    }
    tst <- tst < nrow(tmp)
    if (sum(tst)) {
      if (tst["LC_txt"]) {
        tmp <- aggregate(tmp$MS_txt, list(tmp$LC_txt), unique)
        colnames(tmp) <- c("LC_txt", "MS_txt")
        if (inherits(tmp$MS_txt, "list")) {
          l <- lengths(tmp$MS_txt)
          w <- which(l > 1L)
          tmp$MS_txt[w] <- vapply(tmp$MS_txt[w], paste, "", collapse = "\n")
          tmp$MS_txt <- unlist(tmp$MS_txt)
        }
      }
      if (tst["MS_txt"]) {
        tmp <- aggregate(tmp$LC_txt, list(tmp$MS_txt), unique)[, 2L:1L]
        colnames(tmp) <- c("LC_txt", "MS_txt")
        if (inherits(tmp$LC_txt, "list")) {
          l <- lengths(tmp$LC_txt)
          w <- which(l > 1L)
          tmp$LC_txt[w] <- vapply(tmp$LC_txt[w], paste, "", collapse = "\n")
          tmp$LC_txt <- unlist(tmp$LC_txt)
        }
      }
    }
    nr <- nrow(tmp)
    if (nr > 1L) {
      tmp$LC_txt <- paste0("Method #", 1L:nr, ":\n", tmp$LC_txt)
    }
    brLCMS_txt <- do.call(paste, c(tmp[, c("LC_txt", "MS_txt")], sep = "\n"))
    brLCMS_txt <- paste(brLCMS_txt, collapse = "\n")
    #cat(brLCMS_txt)
    # # Vector of material and methods text for each LC/MS method
    all_LCMS_txts <- c(all_LCMS_txts, brLCMS_txt)
    #
    all_LCMS_instr$LC <- unique(c(all_LCMS_instr$LC, uBrMeth$LC))
    all_LCMS_instr$MS <- unique(c(all_LCMS_instr$MS, uBrMeth$MS))
  }
  if ("raw" %in% tolower(Ext)) {
    rawMeth <- MethodsTbl[which(tolower(MethodsTbl$Ext) == "raw"),]
    rawMeth$Method <- gsub("\\.[^\\.]+$", ".methods.txt", rawMeth$Raw.file)
    rawMeth$Exists <- file.exists(gsub("\\.[^\\.]+$", ".raw", rawMeth$Raw.file)) # Because the file could be another format
    rawMeth$Method.exists <- file.exists(rawMeth$Method)
    w <- which(!rawMeth$Exists)
    if (length(w)) {
      RootPat <- topattern(LocalRoot)
      wRY <- w[grep(RootPat, rawMeth$Folder[w])]
      wRN <- w[grep(RootPat, rawMeth$Folder[w], invert = TRUE)]
      if (length(wRY)) {
        rawMeth$Folder[wRY] <- gsub(RootPat, ArchiveRoot, rawMeth$Folder[wRY])
        rawMeth$Method[wRY] <- gsub(RootPat, ArchiveRoot, rawMeth$Method[wRY])
        rawMeth$Raw.file[wRY] <- gsub(RootPat, ArchiveRoot, rawMeth$Raw.file[wRY])
      }
      if (length(wRN)) {
        rawMeth$Folder[wRN] <- ArchiveRoot
        rawMeth$Method[wRN] <- paste0(ArchiveRoot, "/", basename(rawMeth$Method[wRN]))
        rawMeth$Raw.file[wRY] <- paste0(ArchiveRoot, "/", basename(rawMeth$Raw.file[wRY]))
      }
      rawMeth$Exists[w] <- file.exists(rawMeth$Raw.file[w])
      rawMeth$Method.exists[w] <- file.exists(rawMeth$Method[w])
    }
    w1 <- which(rawMeth$Method.exists)
    w2 <- which((!rawMeth$Method.exists)&(rawMeth$Exists))
    if (length(c(w1, w2))) {
      if (length(w2)) {
        ScanHdsMnTst <- FALSE
        if (!misFun(ScanHdsMnLoc)) {
          ScanHdsMnTst <- file.exists(paste0(ScanHdsMnLoc, "/ScanHeadsman.exe"))
        }
        if (!ScanHdsMnTst) {
          if (!misFun(ScanHdsMnLoc)) {
            warning("Invalid ScanHdsMnLoc argument, no \"ScanHeadsman.exe\" found there!")
          }
          tst <- grep("ScanHeadsman-1\\.2", list.dirs("C:/", recursive = FALSE), value = TRUE)
          ScanHdsMnTst <- length(tst) > 0
          if (ScanHdsMnTst) { ScanHdsMnLoc <- normalizePath(tst[1L], winslash = "/") } else {
            warning("No valid version of ScanHeadsman-1.2 could be located in the system!")
          }
        }
        if (ScanHdsMnTst) {
          parallel::clusterExport(cl, list("ScanHdsMnLoc"), envir = environment())
          f0 <- function(rwfl) { #rwfl <- rawMeth$Raw.file[w2][1L]
            cmd <- paste0("\"", ScanHdsMnLoc, "/ScanHeadsman.exe\" \"", rwfl, "\" -n -m=1 -t=",
                          1)
            #cat(paste0(cmd, "\n"))
            system(cmd)
          }
          environment(f0) <- .GlobalEnv
          tst2 <- parallel::parSapply(cl, rawMeth$Raw.file[w2], f0)
          w2 <- w2[which(!tst2)]
        }
      }
      w <- unique(order(c(w1, w2)))
      f0 <- function(methfl) { #methfl <- rawMeth$Method[w][1L]
        con <- file(methfl, encoding = "UTF-16")
        xmeth <- readLines(con, skipNul = TRUE, encoding = "UTF-8")
        return(list(xmeth))
      }
      environment(f0) <- .GlobalEnv
      methods <- setNames(parallel::parSapply(cl, rawMeth$Method[w], f0), rawMeth$Raw.file[w])
      tmp <- sapply(methods, paste, collapse = "___")
      tst <- aggregate(1L:length(methods), list(tmp), list)
      tst$Raws <- lapply(tst$Group.1, \(x) { names(tmp)[which(tmp == x)] })
      ThRoots <- "Thermo"
      if (nrow(tst) > 1L) {
        ThRoots <- paste0(ThRoots, " meth#", 1L:nrow(tst))
        tmp <- paste(sapply(1L:nrow(tst), \(x) { paste(c(paste0(" - ", ThRoots[x], " files:"),
                                                                    paste0("   > ", tst$Raws[[x]])), collapse = "\n") }),
                     collapse = "\n\n")
        msg <- paste0(nrow(tst), " different Thermo LC-MS/MS acquisition methods detected:\n\n", tmp, "\n")
        cat(msg)
      }
      meths <- lapply(tst$x, \(method) { methods[[method[[1L]]]] })
      thMeth <- lapply(seq(meths), \(i_meth) {
        LC_txt <- MS_txt <- "TEMPLATE"
        MSok <- LCok <- FALSE
        TXT <- c()
        #View(data.frame(Method = meth[which(meth != "")]))
        meth <- meths[[i_meth]]
        if ((length(meth))&&(sum(nchar(meth)))) {
          g <- grep("Method 1 is ", meth)
          lcMeth <- meth[1L:(g-1L)]
          msMeth <- meth[g:length(meth)]
          INSTR <- gsub("^ +Method of ", "", msMeth[3L])
          NANO <- gsub(" +", " ", gsub("Thermo", "", gsub("_", " ", gsub("^Instrument: | on .*", "", grep("^Instrument: ", lcMeth, value = TRUE, ignore.case = TRUE))), ignore.case = TRUE))
          ADflt <- "MS-grade H₂O + 0.1% formic acid"
          BDflt <- "80% acetonitrile in H₂O + 0.08% formic acid"
          if (NANO %in% c("RSLC Nano", "Vanquish NEO")) {
            NANOMAKER <- "ThermoFisher Scientific"
            if (NANO == "RSLC Nano") {
              NANO <- "Ultimate 3000 RSLC_Nano"
              FloNomPat <- "PumpModule\\.NC_Pump\\.Flow\\.Nominal"
              GradBPat <- "PumpModule\\.NC_Pump\\.%B\\.Value"
            }
            if (NANO == "Vanquish NEO") {
              FloNomPat <- "Neo\\.PumpModule\\.Pump\\.Flow\\.Nominal"
              GradBPat <- "Neo\\.PumpModule\\.Pump\\.%B\\.Value"
            }
            g <- grep(" *\\[min\\]", lcMeth)
            RnLns <- gsub(" *\\[min\\].*", "", lcMeth[g])
            w <- which(!is.na(suppressWarnings(as.numeric(RnLns))))
            RnLns <- as.numeric(RnLns[w]) ; g <- g[w]
            RnStrtLn <- min(RnLns)
            RnEndLn <- max(RnLns)
            g <- c(g, length(lcMeth))
            tmp <- lapply(1L:length(RnLns), \(x) {
              lcMeth[g[x]:(g[x+1L]-1L)]
            })
            Grad <- data.frame(Time = RnLns)
            Grad$Flow.Nominal <- sapply(tmp, \(x) { #x <- tmp[[1L]]
              suppressWarnings(as.numeric(gsub(paste0(" *", FloNomPat, ": * | *\\[.+$"), "",
                                               grep(paste0(" *", FloNomPat, ": *"), x, value = TRUE))))*1000L
            })
            Grad$B.Value <- sapply(tmp, \(x) { #x <- tmp[[1L]]
              suppressWarnings(as.numeric(gsub(paste0(" *", GradBPat, ": * | *\\[.+$"), "",
                                               grep(paste0(" *", GradBPat, ": *"), x, value = TRUE))))
            })
            Grad <- Grad[which(sapply(Grad$Flow.Nominal, length) == 1L),]
            Grad$Flow.Nominal <- as.numeric(Grad$Flow.Nominal)
            Grad$B.Value <- as.numeric(Grad$B.Value)
            Grad <- aggregate(Grad[, c("Flow.Nominal", "B.Value")], list(Grad$Time), unique)
            colnames(Grad)[1L] <- "Time"
            FLOWVALUE <- unique(Grad$Flow.Nominal)
            kol <- svDialogs::dlg_list(Kolumns, title = paste0(ThRoots[i_meth],
                                                               " - select column used:"))$res
            if (kol == "Add new...") {
              cat(colChar, "\n")
              kol <- svDialogs::dlg_input("Fill in new column's characteristics (\"|\"-separated):",
                                          colChar)$res
              kol <- unlist(strsplit(kol, " *\\| *"))
              l <- length(kol)
              if (l < 7L) { kol <- c(kol, rep(NA, 7L-l)) }
              kol <- kol[1L:7L]
              tmp <- as.data.frame(t(as.data.frame(kol)))
              colnames(tmp) <- c("Name", "Material", "Length (cm)", "ID (µm)", "Particles size (µm)", "Vendor", "P/N")
              tmp$Function <- "Analytical"
              tmp$Description <- kolDescr(tmp)
              k <- colnames(allKolumns)
              k <- k[which(!k %in% colnames(tmp))]
              tmp[, k] <- NA
              tmp[1L, which(tmp[1L,] == "NA")] <- NA
              allKolumns <- rbind(allKolumns, tmp)
              rownames(allKolumns) <- NULL
              tst <- apply(allKolumns, 1L, paste, collapse = "|")
              tst <- aggregate(1L:length(tst), list(tst), min)
              allKolumns <- allKolumns[sort(tst$x),]
              Kolumns <- unique(c(allKolumns$Description[which(allKolumns$Function == "Analytical")],
                                  "None (direct infusion)", "Add new..."))
              AddKolKount <- AddKolKount + 1L
              kol <- tmp$Description
            }
            COLUMN <- kol
            tstCol <- tolower(substring(COLUMN, 1L, 1L)) %in% c("a", "e", "i", "o", "u")
            if (COLUMN != "None (direct infusion)") {
              kol <- svDialogs::dlg_list(preKolumns, preKolumns[1L],
                                         title = paste0(ThRoots[i_meth],
                                                        " - select trap used:"))$res
              if (kol == "Add new...") {
                cat(colChar, "\n")
                kol <- svDialogs::dlg_input("Fill in new trap's characteristics (\"|\"-separated):",
                                            colChar)$res
                tmp <- as.data.frame(t(as.data.frame(unlist(strsplit(kol, " *\\| *")))))
                colnames(tmp) <- c("Name", "Material", "Length (cm)", "ID (µm)", "Particles size (µm)", "Vendor", "P/N")
                tmp$Function <- "Trap"
                k <- colnames(allKolumns)
                k <- k[which(!k %in% colnames(tmp))]
                tmp[,k] <- NA
                tmp$Description <- kolDescr(tmp)
                allKolumns <- rbind(allKolumns, tmp)
                rownames(allKolumns) <- NULL
                tst <- apply(allKolumns, 1L, paste, collapse = "|")
                tst <- aggregate(1L:length(tst), list(tst), min)
                allKolumns <- allKolumns[sort(tst$x),]
                preKolumns <- unique(c(allKolumns$Description[which(allKolumns$Function == "Trap")],
                                       "None (direct injection)", "Add new..."))
                AddKolKount <- AddKolKount + 1L
                kol <- tmp$Description
              }
              PRECOLUMN <- kol
              tstPreCol <- tolower(substring(PRECOLUMN, 1L, 1L)) %in% c("a", "e", "i", "o", "u")
              SolvA <- svDialogs::dlg_input(paste0(ThRoots[i_meth],
                                                   " - confirm nanoLC solvent A composition:"),
                                            ADflt)$res
              SolvB <- svDialogs::dlg_input(paste0(ThRoots[i_meth],
                                                   " - confirm nanoLC solvent B composition:"),
                                            BDflt)$res
              stopifnot(length(FLOWVALUE) == 1L) # Unexpected situation, currently un-handled
              Grad.Baseline <- min(Grad$B.Value)
              Grad.Plateau <- max(Grad$B.Value)
              Grad.PlateauL <- round(Grad$Time[rev(which(Grad$B.Value == Grad.Plateau))[1L]] - Grad$Time[which(Grad$B.Value == Grad.Plateau)[1L]], 0)
              GradStrt <- Grad$Time[which(Grad$B.Value == Grad.Baseline)[1L]]
              tmp1 <- Grad$Time[which(Grad$B.Value == Grad.Plateau)[1L]]
              tmp0 <- Grad$Time[which(Grad$Time == tmp1)-1L]
              GradEnd <- c(tmp0, tmp1)[(tmp1 - tmp0 > 1)+1L]
              GRADLENGTH <- GradEnd-GradStrt
              Grad <- Grad[which(Grad$Time == GradStrt)[1L]:which(Grad$Time == GradEnd),]
              Grad <- paste(apply(Grad[, c("Time", "B.Value")], 1L, \(x) { paste0(x[[1L]], " min, ", x[[2L]], "%") }), collapse = "; ")
              Grad <- paste0(Grad, ", followed immediately by a ", Grad.PlateauL, " min plateau at ", Grad.Plateau, "%")
            }
            tstNano <- tolower(substring(NANO, 1L, 1L)) %in% c("a", "e", "i", "o", "u")
            tstInstr <- tolower(substring(INSTR, 1L, 1L)) %in% c("a", "e", "i", "o", "u")
            Moult2 <- length(RawFiles) > 1L
            LCok <- TRUE
          } else {
            warning("This function can currently only process Thermo Dionex Ultimate 3000 RSLC_Nano or Waters nanoACQUITY methods!")
          }
          if (grepl("Q[-, ]Exactive", INSTR, ignore.case = TRUE)) {
            INSTRMAKER <- "ThermoFisher Scientific"
            POL <- gsub("^ *Polarity *| *$", "", grep("^Polarity", msMeth, value = TRUE))
            POL <- c("+", "-")[match(POL, c("positive", "negative"))]
            FWHM <- gsub("^ *Chrom\\. peak width \\(FWHM\\) *| *s *$", "", grep("Chrom\\. peak width \\(FWHM\\)", msMeth, value = TRUE))
            msMeth2 <- msMeth[grep(" +Experiments?$", msMeth)[1L]:length(msMeth)]
            MS1 <- data.frame(Start = grep("FULL MS*", msMeth2))
            DIAtst <- c("dd-MS² / dd-SIM", "DIA") %in% msMeth2
            if (!sum(DIAtst)) {
              warning("Message from function:\nSorry, I can currently only deal with DDA or DIA methods.")
            }
            if (sum(DIAtst) == 2L) {
              warning("Message from function:\nThis LC-MS method seems to combine both DDA and DIA scans, I cannot make sense of it in my current code.\nPlease rewrite me!")
            }
            if (sum(DIAtst) == 1L) {
              MS2Type <- c("dd-MS² / dd-SIM", "DIA")[which(DIAtst)]
              MS2 <- data.frame(Start = which(msMeth2 == MS2Type)[1L],
                                End = which(msMeth2 == "                                     Setup")-2L)
              MS1$End <- MS2$Start-2L
              MS1$Data <- apply(MS1[, c("Start", "End")], 1L, \(x) { list(msMeth2[x[[1L]]:x[[2L]]]) })
              MS2$Data <- apply(MS2[, c("Start", "End")], 1L, \(x) { list(msMeth2[x[[1L]]:x[[2L]]]) })
              MSs <- setNames(lapply(c("MS1", "MS2"), \(x) { get(x) }), c("MS1", "MS2"))
              SPECTTYPE <- sapply(c("MS1", "MS2"), \(x) { #x <- "MS1"
                ln <- grep("Spectrum data type", unlist(MSs[[x]]$Data))
                res <- if (length(ln) == 1L) { gsub("^ *Spectrum data type *| *$", "", unlist(MSs[[x]]$Data)[ln]) } else { NA }
                return(res)
              })
              RES <- sapply(c("MS1", "MS2"), \(x) { #x <- "MS1"
                ln <- grep("Resolution", unlist(MSs[[x]]$Data))
                stopifnot(length(ln) == 1L)
                return(gsub("^ *Resolution *| *$", "", unlist(MSs[[x]]$Data)[ln]))
              })
              AGC <- sapply(c("MS1", "MS2"), \(x) { #x <- "MS1"
                ln <- grep("AGC target *[1-9][0-9]*e[1-9][0-9]*", unlist(MSs[[x]]$Data))
                stopifnot(length(ln) == 1L)
                return(gsub("^ *AGC target *| *$", "", unlist(MSs[[x]]$Data)[ln]))
              })
              USCNS <- sapply(c("MS1", "MS2"), \(x) { #x <- "MS1"
                ln <- grep("Microscans", unlist(MSs[[x]]$Data))
                stopifnot(length(ln) == 1L)
                return(as.integer(gsub("^ *Microscans *| *$", "", unlist(MSs[[x]]$Data)[ln])))
              })
              MAXIT <- sapply(c("MS1", "MS2"), \(x) { #x <- "MS1"
                ln <- grep("^ *Maximum IT*", unlist(MSs[[x]]$Data))
                stopifnot(length(ln) == 1L)
                return(gsub("^ *Maximum IT *| *$", "", unlist(MSs[[x]]$Data)[ln]))
              })
              TOPN <- gsub("^ *Loop count *| *$", "", grep("Loop count", unlist(MSs$MS2$Data), value = TRUE))
              ISOWIN <- gsub("^ *Isolation window *| *$", "", grep("Isolation window", unlist(MSs$MS2$Data), value = TRUE))
              ISOWINOFFSET <- gsub("^ *Isolation offset *| *$", "", grep("Isolation offset", unlist(MSs$MS2$Data), value = TRUE))
              ISOWINOFFSET <- if (ISOWINOFFSET == "0.0 m/z") { "(no offset)"} else { paste0("with ", ISOWINOFFSET, " offset") }
              F1MSS <- gsub("^ *Fixed first mass *| *$", "", grep("Fixed first mass", unlist(MSs$MS2$Data), value = TRUE))
              F1MSStst <- suppressWarnings(as.numeric(F1MSS))
              F1MSS <- if ((!is.na(F1MSStst))&&(F1MSStst)) { paste0(F1MSS, " m/z fixed first mass, ") } else { "" }
              NCE <- gsub("^ *\\(N\\)CE */ *stepped *\\(N\\)CE +nce: *| *$", "", grep("\\(N\\)CE */ *stepped \\(N\\)CE +nce:", unlist(MSs$MS2$Data), value = TRUE))
              SCRNG <- sapply(c("MS1", "MS2"), \(x) { #x <- "MS1"
                ln <- grep("^ *Scan range*", unlist(MSs[[x]]$Data))
                #stopifnot(length(ln) == 1L)
                return(gsub("^ *Scan range *| *$", "", unlist(MSs[[x]]$Data)[ln]))
              })
              if (DIAtst[[1L]]) {
                # DDA method
                ExclZ <- gsub("^ *Charge exclusion *| *$", "", grep("Charge exclusion", unlist(MSs$MS2$Data), value = TRUE))
                ExclZ <- unlist(strsplit(ExclZ, ", +"))
                tst1 <- "unassigned" %in% ExclZ
                EXCLZ <- "excluding"
                if (tst1) {
                  #EXCLZ <- paste0(EXCLZ, " unassigned charges, charges ")
                  ExclZ <- ExclZ[which(ExclZ != "unassigned")]
                }
                ExclZ[which(nchar(ExclZ) == 1L)] <- paste0(ExclZ[which(nchar(ExclZ) == 1L)], POL)
                tst2 <- (length(grep("^>", ExclZ)) == 1L)&&(grep("^>", ExclZ) == length(ExclZ))
                if (tst2) { ExclZ <- ExclZ[1L:(length(ExclZ)-1L)] }
                if (length(ExclZ) > 1L) { ExclZ <- c(paste(ExclZ[1L:(length(ExclZ)-1L)], collapse = ", "), ExclZ[length(ExclZ)]) }
                ExclZ <- paste(ExclZ, collapse = c(" and ", ", ")[(tst1|tst2)+1L])
                EXCLZ <- paste0("excluding charges ", ExclZ)
                if (tst1||tst2) {
                  EXCLZ <- if (tst1+tst2) { paste0(EXCLZ, " and higher or unassigned,") } else {
                    paste0(EXCLZ, " and ", c("higher", "unassigned")[which(c(tst1, tst2))], ",")
                  }
                }
                DYNEXCL <- as.numeric(gsub("^ *Dynamic exclusion *|  *s *$", "", grep("Dynamic exclusion", unlist(MSs$MS2$Data), value = TRUE)))
                MS_txt <- paste0(MS_txt, "up to ", TOPN, " MS2s per cycle. ",
                                 "DDA MS2 parameters: ")
                if (!is.na(SPECTTYPE["MS2"])) { MS_txt <- paste0(MS_txt, SPECTTYPE["MS2"], " mode, ") }
                MS_txt <- paste0(MS_txt, USCNS["MS2"], " microscan", c("", "s")[(USCNS["MS2"] > 1L)+1L], ", ",
                                 RES["MS2"], " resolution, ",
                                 "AGC target ", AGC["MS2"], ", ",
                                 MAXIT["MS2"], " maximum IT, ",
                                 ISOWIN, " isolation window ", ISOWINOFFSET, ", ", 
                                 F1MSS,
                                 SCRNG["MS1"], ", ",
                                 "NCE ", NCE, ", ",
                                 EXCLZ, " ",
                                 DYNEXCL, "s dynamic exclusion.")
              }
              if (DIAtst[[2L]]) {
                # DIA method
                InclStrt <- grep(" +Mass +Formula +Species +CS +Polarity +Start +End +\\(N\\)CE +MSX +ID +Comment", msMeth2)
                InclHdr <- msMeth2[InclStrt]
                InclUnit <- msMeth2[InclStrt+1L]
                InclMn <- InclStrt+2L
                InclRg <- grep(" *[0-9]+\\.*[0-9]* +((Posi)|(Nega))tive +", msMeth2[InclMn:length(msMeth2)]) + InclMn - 1L
                InclMx <- max(InclRg[which(sapply(InclRg, \(x) {
                  sum(!InclMn:x %in% InclRg) == 0L
                }))])
                InclWind <- msMeth2[InclMn:InclMx]
                brks <- unlist(strsplit(InclHdr, ""))
                l <- length(brks)
                brks <- which(vapply(2L:l, \(x) { (brks[x] == " ")&(brks[x-1] != " ") }, TRUE))
                brks <- c(0L, brks)
                l <- length(brks)
                InclWind <- as.data.frame(sapply(2L:l, \(x) {
                  gsub("^ +", "", substr(InclWind, brks[x-1]+1L, brks[x]))
                }))
                colnames(InclWind) <- sapply(2L:l, \(x) {
                  df <- data.frame(x = gsub("^ +", "", substr(InclHdr, brks[x-1L]+1L, brks[x])),
                                   y = gsub("^ +", "", substr(InclUnit, brks[x-1L]+1L, brks[x])))
                  x <- gsub(" +$", "", do.call(paste, c(df, sep = " ")))
                })
                InclWind$`Mass [m/z]` <- as.numeric(InclWind$`Mass [m/z]`)
                InclWind$`Start [min]` <- as.numeric(InclWind$`Start [min]`)
                InclWind$`End [min]` <- as.numeric(InclWind$`End [min]`)
                InclWind <- InclWind[order(InclWind$`Start [min]`, InclWind$`Mass [m/z]`, decreasing = FALSE),]
                NWind <- nrow(InclWind)
                #if (length(InclWind) == TOPN) {
                InclWindL <- InclWind$`Mass [m/z]`-as.numeric(gsub(" m/z$", "", ISOWIN))/2
                InclWindR <- InclWind$`Mass [m/z]`+as.numeric(gsub(" m/z$", "", ISOWIN))/2
                OLtst <- unique(round(InclWindR[1L:(NWind-1L)] - InclWindL[2L:NWind], 3L))
                stopifnot(length(OLtst) == 1L)
                if (OLtst == 0L) { OLtxt <- "no overlap" }
                if (OLtst > 0L) { OLtxt <- paste0(OLtst, " m/z overlap") }
                if (OLtst < 0L) { OLtxt <- paste0(OLtst, " m/z non-covered gap between adjacent windows") }
                #}
              }
              MSok <- TRUE
            }
          } else {
            warning("This function currently can only process Q-Exactive type methods!")
          }
          if (LCok) {
            LC_txt <- paste0(c("The s", "S")[Moult2+1L], "ample", c(" was", "s were")[Moult2+1L], " analyzed by LC-MS/MS on a",
                             c("", "n")[tstNano+1L], " ", NANO, " nano-HPLC (",
                             NANOMAKER, ") coupled with a", c("", "n")[tstInstr+1L], " ", INSTR, " (", INSTRMAKER, ")")
            if (COLUMN != "None (direct infusion)") {
              LC_txt <- if (PRECOLUMN != "None (direct injection)") {
                paste0(LC_txt, ", concentrated over a", c("", "n")[tstPreCol+1L], " ", PRECOLUMN, ", then ")
              } else { paste0(LC_txt, ". ", c("The s", "S")[Moult2+1L], "ample", c(" was", "s were")[Moult2+1L], " ") }
              LC_txt <- paste0(LC_txt, "bound to a", c("", "n")[tstCol+1L], " ", COLUMN, " and eluted over the following ", GRADLENGTH,
                               " min gradient: solvent A, ", SolvA, "; ", "solvent B, ", SolvB, "; ",
                               "constant ", FLOWVALUE, " nL/min flow; ",
                               "B percentage: ", Grad, ".")
            } else { LC_txt <- paste0(LC_txt, " by direct infusion.") }
          }
          if (MSok) {
            MS_txt <- paste0("Mass spectra were acquired in positive mode with a Data ", c("D", "Ind")[which(DIAtst)], "ependent Acquisition (D",
                             c("D", "I")[which(DIAtst)], "A) method: chromatographic peak width = ",
                             "FWHM ", FWHM, " s, ",
                             "MS1 parameters: ")
            if (!is.na(SPECTTYPE["MS1"])) { MS_txt <- paste0(MS_txt, SPECTTYPE["MS1"], " mode, ") }
            MS_txt <- paste0(MS_txt,
                             USCNS["MS1"], " microscan", c("", "s")[(USCNS["MS1"] > 1L)+1L], ", ",
                             RES["MS1"], " resolution, ",
                             "AGC target ", AGC["MS1"], ", ",
                             MAXIT["MS1"], " maximum IT, ",
                             SCRNG["MS1"], "; ")
            if (DIAtst[[1L]]) {
              MS_txt <- paste0(MS_txt, "up to ", TOPN, " MS2s per cycle. ",
                               "DDA MS2 parameters: ")
              if (!is.na(SPECTTYPE["MS2"])) { MS_txt <- paste0(MS_txt, SPECTTYPE["MS2"], " mode, ") }
              MS_txt <- paste0(MS_txt, USCNS["MS2"], " microscan", c("", "s")[(USCNS["MS2"] > 1L)+1L], ", ",
                               RES["MS2"], " resolution, ",
                               "AGC target ", AGC["MS2"], ", ",
                               MAXIT["MS2"], " maximum IT, ",
                               ISOWIN, " isolation window ", ISOWINOFFSET, ", ", 
                               F1MSS,
                               SCRNG["MS1"], ", ",
                               "NCE ", NCE, ", ",
                               EXCLZ, " ",
                               DYNEXCL, "s dynamic exclusion.")
            }
            if (DIAtst[[2L]]) {
              MS_txt <- paste0(MS_txt, "DIA scans: ",
                               c("", paste0(TOPN, " MS2 scans per cycle, "))[(TOPN != NWind)+1L],
                               NWind, " windows of ", ISOWIN, " width per cycle covering the range from ", min(InclWindL), " to ",
                               max(InclWindR), " m/z (", OLtxt, "), ")
              if (!is.na(SPECTTYPE["MS2"])) { MS_txt <- paste0(MS_txt, "spectra acquired in ", SPECTTYPE["MS2"], " mode, ") }
              MS_txt <- paste0(MS_txt,
                               USCNS["MS2"], " microscan", c("", "s")[(USCNS["MS2"] > 1L)+1L], ", at ",
                               RES["MS2"], " resolution; ",
                               "AGC target ", AGC["MS2"], ", ",
                               MAXIT["MS2"], " maximum IT, ",
                               F1MSS,
                               "NCE ", NCE, ".")
            }
          }
          
          TXT <- paste0(LC_txt, " ", MS_txt)
          #writeClipboard(TXT)
        }
        return(list(Method = data.frame(LC_txt = LC_txt,
                                        MS_txt = MS_txt),
                    LC = NANO,
                    LC_vendor = NANOMAKER,
                    MS = INSTR,
                    MS_vendor = INSTRMAKER))
      })
      thMeth_LC <- vapply(thMeth, \(x) { x$LC }, "")
      thMeth_MS <- vapply(thMeth, \(x) { x$MS }, "")
      all_LCMS_instr$LC <- unique(c(all_LCMS_instr$LC, thMeth_LC))
      all_LCMS_instr$MS <- unique(c(all_LCMS_instr$MS, thMeth_MS))
      thLCMS_txt <- lapply(thMeth, \(x) { x$Method })
      thLCMS_txt <- do.call(rbind, thLCMS_txt)
      tmp <- thLCMS_txt[, c("LC_txt", "MS_txt")]
      #tmp <- tmp[c(1L, 2L, 2L),]
      #tmp$LC_txt[2L] <- tmp$LC_txt[1L]
      #tmp$MS_txt[3L] <- tmp$MS_txt[1L]
      tst <- vapply(c("LC_txt", "MS_txt"), \(x) { length(unique(tmp[[x]])) }, 1L)
      for (i in c("LC", "MS")) {
        k <- paste0(i, "_txt")
        if (tst[[k]] > 1L) {
          m <- match(tmp[[k]], unique(tmp[[k]]))
          tmp[[k]] <- paste0(i, " method #", m, ": ", tmp[[k]])
        }
      }
      tst <- tst < nrow(tmp)
      if (sum(tst)) {
        if (tst["LC_txt"]) {
          tmp <- aggregate(tmp$MS_txt, list(tmp$LC_txt), unique)
          colnames(tmp) <- c("LC_txt", "MS_txt")
          if (inherits(tmp$MS_txt, "list")) {
            l <- lengths(tmp$MS_txt)
            w <- which(l > 1L)
            tmp$MS_txt[w] <- vapply(tmp$MS_txt[w], paste, "", collapse = "\n")
            tmp$MS_txt <- unlist(tmp$MS_txt)
          }
        }
        if (tst["MS_txt"]) {
          tmp <- aggregate(tmp$LC_txt, list(tmp$MS_txt), unique)[, 2L:1L]
          colnames(tmp) <- c("LC_txt", "MS_txt")
          if (inherits(tmp$LC_txt, "list")) {
            l <- lengths(tmp$LC_txt)
            w <- which(l > 1L)
            tmp$LC_txt[w] <- vapply(tmp$LC_txt[w], paste, "", collapse = "\n")
            tmp$LC_txt <- unlist(tmp$LC_txt)
          }
        }
      }
      nr <- nrow(tmp)
      if (nr > 1L) {
        tmp$LC_txt <- paste0("Method #", 1L:nr, ":\n", tmp$LC_txt)
      }
      thLCMS_txt <- do.call(paste, c(tmp[, c("LC_txt", "MS_txt")], sep = "\n"))
      thLCMS_txt <- paste(thLCMS_txt, collapse = "\n")
    } else {
      warning("All files and their methods are missing, at least one file or one method is required!")
      thLCMS_txt <- "TEMPLATE"
    }
    all_LCMS_txts <- c(all_LCMS_txts, thLCMS_txt)
  }
  if (AddKolKount) {
    require(openxlsx)
    Head <- createStyle(textDecoration = c("bold", "underline"))
    wb <- createWorkbook()
    addWorksheet(wb, "Sheet1")
    writeData(wb, "Sheet1", allKolumns, 1L, 1L)
    addStyle(wb, "Sheet1", Head, 1L, 1L:ncol(allKolumns))
    tst <- vapply(1L:ncol(allKolumns), \(x) {
      x1 <- ceiling(nchar(colnames(allKolumns)[x])*1.2)
      x2 <- ceiling(nchar(allKolumns[[x]])*1.2)
      return(max(c(x1, min(c(60L, x2), na.rm = TRUE))))
    }, 1L)
    setColWidths(wb, "Sheet1", 1L:ncol(allKolumns), tst)
    try(saveWorkbook(wb, Columns2, TRUE), silent = TRUE) # We do not want the function to fail if this fails, it isn't worth the trouble
    #openwd(dirname(Columns2))
  }
  #
  if (stopCl) { parallel::stopCluster(cl) }
  return(list(Text = all_LCMS_txts,
              Instruments = all_LCMS_instr))
}
