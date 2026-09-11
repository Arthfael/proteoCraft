#' load_Bckp
#'
#' @description
#' This function:
#'  - cleans the existing environment,
#'  - reloads a data analysis backup,
#'  - sets the proper working directory,
#'  - loads all necessary packages and
#'  - (if necessary) re-creates the original parallel cluster.
#' 
#' @param bckp The backup file to reload. If missing, opens a popup selection window.
#' @param startDir Here you can specify the directory in which the selection window will start.
#' @param clean If set to TRUE (default), will first remove all objects existing in the current environment.
#' @param loadPack Should we load packages? Default = TRUE
#'
#' @export

load_Bckp <- function(backup,
                      startDir,
                      clean = TRUE,
                      loadPack = TRUE) {
  # Cleanup workspace here
  TESTING <- FALSE
  #DefArg(load_Bckp);TESTING <- TRUE
  if (clean) { suppressWarnings(rm(list = setdiff(ls(), c("TESTING", "clean", "loadPack")), envir = .GlobalEnv)) }
  misFun <- if (TESTING) {
    # Note:
    # This is not a perfect alternative to missing but will work in most cases, unless x matches a function imported by a package
    \(x) { return(!exists(deparse(substitute(x)))) }
  } else { missing }
  #
  RPath <- as.data.frame(library()$results)
  RPath <- normalizePath(RPath$LibPath[match("proteoCraft", RPath$Package)], winslash = "/")
  libPath <- paste0(RPath, "/proteoCraft")
  homePath <- paste0(normalizePath(Sys.getenv("HOME"), winslash = "/"), "/R/proteoCraft")
  #
  wd_test <- \() {
    if (!exists("wd", .GlobalEnv)) { return(FALSE) }
    wd <- get("wd", .GlobalEnv)
    return((length(wd) == 1L) && (!is.na(wd)) && is.character(wd) && dir.exists(wd))
  }
  if (misFun(backup)) {
    if (misFun(startDir)) {
      if (wd_test()) {
        defltdir <- wd
      } else {
        dlft <- openxlsx::read.xlsx(paste0(homePath, "/Default_locations.xlsx"))
        defltdir <- dlft$Path[match("Temporary folder", dlft$Folder)]
      }
    } else { defltdir <- startDir }
    #
    bckp <- rstudioapi::selectFile("Select backup file to load...",
                                   path = defltdir,
                                   filter = "Backup file (*.RDS | *.RData)") # RData is legacy: I used the wrong extension, we have been saving RDS for a while now!
  } else {
    bckp <- backup
  }
  if (grepl("^~", bckp)) { bckp <- gsub("^~", Sys.getenv("R_USER"), bckp) }
  if ((!nchar(bckp)) || (length(bckp) != 1L) || (!is.character(bckp))) {
    warning("\"backup\" must be a single length 1 character path to a valid proteoCraft backup file!")
    return()
  }
  if ((!file.exists(bckp))) {
    warning("The specifid \"backup\" file does not exist!")
    return()
  }
  wdExisted <- exists("wd", .GlobalEnv) && is.character(wd) && (length(wd) == 1L) && (!is.na(wd)) && dir.exists(wd)
  bckpDeerayktoray <- dirname(bckp)
  #
  # Now, I have recently switched to a different, faster way of saving backups using parallelization.
  # The bit below is meant to allow some form of backwards compatibility with older backups.
  tst <- "Didnae work, matey!"
  inst <- as.data.frame(installed.packages())
  cat("Reloading backup file\n   ", bckp, "\n...\nplease hold...\n")
  tst <- try(loadFun(bckp, tryClassic = TRUE), silent = TRUE)
  if (inherits(tst, "try-error") || (is.character(tst) && (length(tst) == 1L) && (tst == "Didnae work, matey!"))) {
    stop("Backup re-loading failed!")
  }
  assign("backupFile", bckp, envir = .GlobalEnv)
  if (!exists(".obj")) { .obj %<o% ".obj"  }
  assign(".obj", union(.obj, "backupFile"), envir = .GlobalEnv) # Exception: this one we keep always at the far end!
  #
  if ((!exists("wd")) || (!wd_test())) {
    if (!wdExisted) {
      warning("Invalid work directory loaded from backup file, using its parent directory instead!")
    } else {
      warning("Setting work directory to the backup file's parent directory.")
    }
    wd <- bckpDeerayktoray
  }
  tst <- try(setwd(wd), silent = TRUE)
  #
  # Important
  # #########
  # Sometimes, when the package (proteoCraft) has been updated, and a backup is then loaded using this function,
  # it appears like old versions of the package are being loaded too.
  # The chunk below is an attempt to fix the issue, but may not work.
  #
  sessInf <- sessionInfo()
  tmp <- c(names(sessInf$loadedOnly),
           names(sessInf$otherPkgs))
  if (length(tmp)) {
    #lapply(tmp, require, character.only = TRUE) # This was on in previous versions but may be leftover code from testing... no idea why it was on here
    # It may have been the cause of issues I had when running the function multiple times in one session, where packages could be loaded in a fresh session but not a 
    # fresh one...
    tmp <- paste0('package:', tmp)
    try(invisible(lapply(tmp, detach, character.only = TRUE, unload = TRUE)), silent = TRUE)
  }
  allCores <- parallel::detectCores()
  maxCores <- max(c(round(allCores*0.95)-1L, 1L)) # New slightly more conservative default
  if (loadPack) {
    if (exists("cran_req")) {
      for (pack in cran_req) {
        loadTst <- try({ library(pack, character.only = TRUE) }, silent = TRUE)
        if (inherits(loadTst, "try-error")) { warning(loadTst) }
      }
      if ("data.table" %in% cran_req) { setDTthreads(threads = maxCores) }
    }
    if (exists("bioc_req")) {
      for (pack in bioc_req) {
        loadTst <- try({ library(pack, character.only = TRUE) }, silent = TRUE)
        if (inherits(loadTst, "try-error")) { warning(loadTst) }
      }
    }
    library("proteoCraft", character.only = TRUE)
  }
  # Re-create parallel cluster
  usePar <- FALSE
  if (!sum(!c("parClust", "N.clust") %in% .obj)) {
    #
    # Check N.clust parameter validity
    N.clust <- try(as.integer(N.clust), silent = TRUE)
    if ((!exists("N.clust")) || (!is.numeric(N.clust)) || is.na(N.clust) || (N.clust < 1L) || (N.clust > maxCores)) {
      N.clust <- maxCores
    }
    assign("N.clust", N.clust, envir = .GlobalEnv) # Export
    if (exists("cran_req") && ("data.table" %in% cran_req)) { data.table::setDTthreads(threads = N.clust) }
    #
    # Do we already have a compatible cluster running?
    if (exists("parClust") && exists("N.clust") && inherits(parClust, "cluster")) {
      currNodes <- as.integer(sub(" .*", "", sub("socket cluster with ", "", capture.output(parClust))))
      if (currNodes != N.clust) { parallel::stopCluster(parClust) } # Delete if there is one but is is incompatible
    }
    # If not, create it:
    a <- 1
    tst <- try(parallel::clusterExport(parClust, "a", envir = environment()), silent = TRUE)
    if (inherits(tst, "try-error")) {
      if (exists("parClust")) { try(parallel::stopCluster(parClust), silent = TRUE) }
      parClust <- parallel::makeCluster(N.clust, type = "SOCK")
      assign("parClust", parClust, envir = .GlobalEnv) # Export cluster
    }
    usePar <- TRUE
  }
  #
  cat(paste0("Backup \"", bckp, "\" loaded, work directory set and packages loaded.\n"))
  if (exists(".obj") && (sum(!.obj %in% c("backupFile", ".obj", "%<o%", "%<c%")))) {
    if (exists("ScriptPath")) {
      cat("Analysis script used ---> ", ScriptPath, "\n")
      if (file.exists(ScriptPath)) {
        try({
          # This code rests on the assumption that we follow the following pattern when calling sourced sub-scripts:
          # .*Src ((<-)|(%<[co]%)) # ... define path to source
          # source(.*Src) # source call
          #
          # Process script and sources and identify which remanent objects they create
          allSrcs <- data.frame(Path = c(ScriptPath,
                                         list.files(paste0(libPath, "/extdata/Sources"), ".R$", full.names = TRUE, recursive = TRUE),
                                         list.files(paste0(libPath, "/extdata/Pepper"), ".R$", full.names = TRUE, recursive = TRUE)))
          allSrcs$Name <- basename(allSrcs$Path)
          allSrcs$Code <- lapply(allSrcs$Path, \(x) { sub(" *#.*", "", readLines(x)) })
          scrptCode <- sub("^ *", "", sub(" *#.*", "", allSrcs$Code[[1L]]))
          lScrpt <- length(scrptCode)
          w_bckpCalls <- c(1L, grep("source\\( *bckpSrc", scrptCode), lScrpt)
          # Identify other scripts sourced by each script
          allSrcs$g <- lapply(allSrcs$Code, \(x) { grep("source\\(", x) }) # ... first where the source call occurs...
          allSrcs$gSrc <- lapply(allSrcs$Code, \(x) { grep("Src *((<-)|(%<[co]%)) *", x) }) #... then where...
          allSrcs$allSrcs <- lapply(1L:nrow(allSrcs), \(x) { allSrcs$Code[[x]][allSrcs$gSrc[[x]]] }) # ... and how the path to which source is called is defined
          # Case where we assign a source to a non-generic object name:
          # specSrcDefs = source paths which are not assigned to the recycled "Src" object but assigned specific, often remanent object names
          specSrcDefs <- sort(unique(sub(" ((<-)|(*%<[co]%)) *", " <- ", sub("^ *", "", unlist(allSrcs$allSrcs)))))
          specSrcDefs <- grep("^Src", specSrcDefs, value = TRUE, invert = TRUE)
          specSrcDefs <- data.frame(obj = sub(".* +", "", sub(" *((<-)|(%<[co]%)) *.*", "", specSrcDefs)),
                                    src = sub(".*/", "", sub("\\.R\".*", "", specSrcDefs)))
          specSrcDefs <- specSrcDefs[which(specSrcDefs$src == make.names(specSrcDefs$src)),]
          specSrcDefs$src <- paste0(specSrcDefs$src, ".R")
          specSrcDefs <- specSrcDefs[which((specSrcDefs$src != basename(ScriptPath)) & (specSrcDefs$src %in% allSrcs$Name)),]
          #
          allSrcs$subSrcs <- lapply(1L:nrow(allSrcs), \(x) { #x <- 1L #x <- x+1L #x <- 53L #x <- 75L #x <- 110L
            g <- allSrcs$g[[x]]
            lg <- length(g)
            if (!lg) { return() }
            gSrc <- allSrcs$gSrc[[x]]
            locSrcs <- length(gSrc)
            if (locSrcs) {
              srcDefs <- allSrcs$allSrcs[[x]]
              srcDefs <- data.frame(gSrc = gSrc,
                                    obj = sub(".* +", "", sub(" *((<-)|(%<[co]%)) *.*", "", srcDefs)),
                                    src = sub(".*/", "", sub("\\.R\".*", "", srcDefs)))
              srcDefs <- srcDefs[which(srcDefs$src == make.names(srcDefs$src)),]
              locSrcs <- nrow(srcDefs)
            }
            if (locSrcs) {
              srcDefs$src <- paste0(srcDefs$src, ".R")
              srcDefs <- srcDefs[which((srcDefs$src != basename(ScriptPath)) & (srcDefs$src %in% allSrcs$Name)),]
              locSrcs <- nrow(srcDefs)
            }
            Srcs <- sub("(, *local *= *((FALSE)|(TRUE)|(F)|(T)))? *\\).*", "", sub(".*source\\( *", "", allSrcs$Code[[x]][g]))
            srcNms <- rep(NA, lg)
            for (i in 1L:lg) { #i <- 1L
              Src <- Srcs[i]
              j <- c()
              if (locSrcs) {
                j <- which((srcDefs$gSrc < g[i]) & (srcDefs$obj == Src))
              }
              if (length(j)) { srcNms[i] <- srcDefs$src[max(j)] } else {
                m <- match(Src, specSrcDefs$obj)
                if (!is.na(m)) { srcNms[i] <- specSrcDefs$src[m] } # else { print(Src) }
              }
            }
            w <- which(!is.na(srcNms))
            return(data.frame(i = g[w],
                              source = srcNms[w]))
          })
          # Replace code in parent scripts
          allSrcs$Code_lst <- lapply(allSrcs$Code, as.list)
          lSubSrcs <- vapply(allSrcs$subSrcs, \(x) { return( if (!is.data.frame(x)) { 0L } else { nrow(x) } ) }, 1L)
          wh1 <- which(lSubSrcs == 0L)
          wh2 <- which(lSubSrcs > 0L)
          while (length(wh2)) {
            tmp <- lapply(wh2, \(i) { #i <- wh2[1L] #i <- wh2[4L]
              Srcs <- allSrcs$subSrcs[[i]]
              Code <- allSrcs$Code_lst[[i]]
              w1 <- which(Srcs$source %in% allSrcs$Name[wh1])
              w2 <- which(!Srcs$source %in% allSrcs$Name[wh1])
              if (length(w1)) {
                Code[Srcs$i[w1]] <- allSrcs$Code_lst[match(Srcs$source[w1], allSrcs$Name)]
                Srcs <- Srcs[w2,]
              }
              return(list(sources = Srcs,
                          code = Code))
            })
            allSrcs$Code_lst[wh2] <- lapply(tmp, \(x) { x$code })
            allSrcs$subSrcs[wh2] <- lapply(tmp, \(x) { x$sources })
            lSubSrcs <- vapply(allSrcs$subSrcs, \(x) { return( if (!is.data.frame(x)) { 0L } else { nrow(x) } ) }, 1L)
            wh1 <- which(lSubSrcs == 0L)
            wh2 <- which(lSubSrcs > 0L)
          }
          code_lst <- allSrcs$Code_lst[[1L]]
          code_lst <- listMelt(code_lst, 1L:length(code_lst), c("code", "row"))
          #
          # Analyse objects created
          # Identify remanent objects created by each source
          w_obj1 <- grep("%<[co]%", code_lst$code)
          w_obj2 <- grep(".obj <- ((union)|(c)|(unique\\(c))\\(", code_lst$code)
          obj <- list()
          if (length(w_obj1)) {
            obj1 <- sub(".* +", "", sub(" *%<[co]%.*", "", code_lst$code[w_obj1]))
            w <- which(obj1 == make.names(obj1)) # Check that we only parsed valid variable names
            obj1 <- data.frame(obj = obj1[w],
                               row = code_lst$row[w_obj1[w]])
            obj$v1 <- obj1
          }
          if (length(w_obj2)) {
            obj2 <- sub("\\).*", "", sub(".*, *", "", code_lst$code[w_obj2])) # Sometimes the left element is an expression: we cannot handle those cases 
            w <- which(obj2 == make.names(obj2)) # Check that we only parsed valid variable names
            obj2 <- data.frame(obj = obj2[w],
                               row = code_lst$row[w_obj2[w]])
            obj$v2 <- obj2
          }
          obj <- do.call(rbind, obj)
          obj <- obj[order(obj$row),]
          obj$prediction <- vapply(obj$row, \(x) { min(w_bckpCalls[which(w_bckpCalls > x)]) }, 1L)
          cat(paste0("\n   FYI, the last remanent object created before this backup was made is \"", .obj[1L], "\""))
          w <- which(.obj %in% obj$obj)
          if (length(w)) {
            pred <- lapply(.obj[w], \(x) { unique(obj$prediction[which(obj$obj == x)]) })
            pred1 <- pred[[1L]]         # Prediction from 1st object in .obj (last added)
            pred2 <- max(unlist(pred))  # Prediction from all objects
            tst <- (pred2 %in% pred1)
            if (pred2 < lScrpt) {
              rg <- (pred2+1L):lScrpt
              rg <- rg[which(scrptCode[rg] != "")]
              if (length(rg)) { pred2 <- rg[1L] }
            }
            if (pred2 > rev(w_bckpCalls)[2L]) {
              msg <- paste0("   ", c("However, b", "B")[tst + 1L],
                            "ackup analysis suggests that this backup had reached the end of the analysis, so there should be nothing more to run...\n   But maybe you want to re-run some parts without starting from scratch?\n   (opening script...)\n")
              cat(msg)
              suppressWarnings(rstudioapi::documentOpen(ScriptPath))
            } else {
              cat(paste0(c("   However, it seems that at some point the script had also been run with this data beyond that point, thus we suggest starting execution from",
                           "   Backup analysis suggests starting execution from")[tst + 1L], pred2, ")\n"))
              suppressWarnings(rstudioapi::documentOpen(ScriptPath, pred2))
            }
          }
        }, silent = TRUE)
      } else {
        cat(" ... but it appears the file doesn't exist anymore...\n")
        cat(paste0("   FYI, the last object listed in .obj is \"", rev(.obj)[1L], "\".\n"))
      }
    } else {
      cat(paste0("   FYI, the last object listed in .obj is \"", rev(.obj)[1L], "\".\n"))
    }
  }
  if (exists("mySeed")) { set.seed(mySeed) }
  cat("\nYou're good to go!\n")
  return()
}
