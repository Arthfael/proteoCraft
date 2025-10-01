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
  if (clean) { rm(list = ls(), envir = .GlobalEnv) }
  #
  TESTING <- FALSE
  #proteoCraft::DefArg(proteoCraft::load_Bckp); TESTING <- TRUE
  if (TESTING) {
    # Note:
    # This is not a perfect alternative to missing but will work in most cases, unless x matches a function imported by a package
    misFun <- function(x) { return(!exists(deparse(substitute(x)))) }
  } else { misFun <- missing }
  #
  homePath <- paste0(normalizePath(Sys.getenv("HOME"), winslash = "/"), "/R/proteoCraft")
  #
  if (misFun(backup)) {
    if (misFun(startDir)) {
      if ((exists("wd"))&&("character" %in% class(wd))&&(dir.exists(wd))) {
        defltdir <- wd
      } else {
        dlft <- openxlsx::read.xlsx(paste0(homePath, "/Default_locations.xlsx"))
        defltdir <- dlft$Path[match("Temporary folder", dlft$Folder)]
      }
    } else { defltdir <- startDir }
    #bckp <- normalizePath(choose.files(paste0(defltdir, "/*.RData"), multi = FALSE), winslash = "/")
    bckp <- rstudioapi::selectFile("Select backup file to load...",
                                   path = paste0(defltdir, "/*.RData"),
                                   filter = "RData file (*.RData)")
  } else {
    bckp <- backup
  }
  if ((!nchar(bckp))||(length(bckp) != 1)||(!"character" %in% class(bckp))) {
    warning("\"backup\" must be a single length = 1 character path to a valid proteoCraft backup file!")
    return()
  }
  if ((!file.exists(bckp))) {
    warning("The specifid \"backup\" file does not exist!")
    return()
  }
  wdExisted <- exists("wd", .GlobalEnv)
  bckpDeerayktoray <- dirname(bckp)
  #
  #bckp <- "~/R/proteoCraft/AN_GNRGFL1_5637917142/Backup.RData"
  # Now, I have recently switched to a different, faster way of saving backups using parallelization.
  # The bit below is meant to allow some form of backwards compatibility with older backups.
  tst <- "Didnae work, matey!"
  inst <- as.data.frame(installed.packages())
  tst <- try(proteoCraft::loadFun(bckp), silent = TRUE)
  if (("try-error" %in% class(tst))||(("character" %in% class(tst))&&(length(tst) == 1)&&(tst == "Didnae work, matey!"))) {
    stop("Backup re-loading failed!")
  }
  #
  if ((exists("wd"))&&(!dir.exists(wd))) {
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
  maxCores <- max(c(round(allCores*0.95)-1, 1)) # New slightly more conservative default
  if (loadPack) {
    if (exists("cran_req")) {
      for (pack in cran_req) { library(pack, character.only = TRUE) }
      if ("data.table" %in% cran_req) { setDTthreads(threads = maxCores) }
    }
    if (exists("bioc_req")) {
      for (pack in bioc_req) { library(pack, character.only = TRUE) }
    }
    library("proteoCraft", character.only = TRUE)
  }
  # Re-create parallel cluster
  usePar <- FALSE
  if ((exists(".obj"))&&(sum(c("parClust", "N.clust") %in% .obj) == 2)) {
    #
    # Check N.clust parameter validity
    N.clust <- try(as.integer(N.clust), silent = TRUE)
    if ((!exists("N.clust"))||(!is.numeric(N.clust))||(is.na(N.clust))||(N.clust < 1)||(N.clust > maxCores)) {
      N.clust <- maxCores
    }
    N.clust <<- N.clust # Export
    if ((exists("cran_req"))&&("data.table" %in% cran_req)) { data.table::setDTthreads(threads = N.clust) }
    #
    # Do we already have a compatible cluster running?
    if ((exists("parClust"))&&(exists("N.clust"))&&("cluster" %in% class(parClust))) {
      currNodes <- gsub(" .*", "", gsub("socket cluster with ", "", capture.output(parClust)))
      currNodes <- as.integer(currNodes)
      if (currNodes != N.clust) { parallel::stopCluster(parClust) }
    }
    # If not, create it:
    a <- 1
    tst <- try(parallel::clusterExport(parClust, "a", envir = environment()), silent = TRUE)
    if ("try-error" %in% class(tst)) {
      if (exists("parClust")) { try(parallel::stopCluster(parClust), silent = TRUE) }
      parClust <- parallel::makeCluster(N.clust, type = "SOCK")
    }
    usePar <- TRUE
  }
  #
  cat(paste0("Backup \"", bckp, "\" loaded, work directory set and packages loaded.\n"))
  if (exists(".obj")) {
    if (exists("ScriptPath")) {
      cat("Analysis script used ---> ", ScriptPath, "\n")
      if (file.exists(ScriptPath)) {
        scrpt <- readLines(ScriptPath)
        scrpt <- data.frame(call = scrpt)
        scrpt$row <- 1:nrow(scrpt)
        scrpt$listCall <- as.list(scrpt$call)
        g1 <- grep("^ *Src <- ", scrpt$call)
        g2 <- grep("^ *source\\(Src\\)", scrpt$call)
        l <- length(g1)
        while (l) {
          if ((l == length(g2))
              &&(sum(g1 > g2) == 0)
              &&(sum(g2[1:(l-1)] > g1[2:l]) == 0)) {
            srcFls <- sapply(gsub("^ *Src <- ", "", scrpt$call[g1]), function(x) {
              x <- eval(parse(text = x))
              if (is.null(x)) { x <- "" }
              return(x)
            })
            wY <- which((nchar(srcFls) > 0)&(file.exists(srcFls)))
            wN <- which((nchar(srcFls) == 0)|(!file.exists(srcFls)))
            srcs <- suppressWarnings(lapply(srcFls[wY], readLines))
            scrpt$listCall[g1] <- ""
            scrpt$listCall[g2[wY]] <- srcs
            scrpt$listCall[g2[wN]] <- ""
            scrpt <- listMelt(scrpt$listCall, scrpt$row, c("call", "row"))
            scrpt$listCall <- as.list(scrpt$call)
            g1 <- grep("^ *Src <- ", scrpt$call)
            g2 <- grep("^ *source\\(Src\\)", scrpt$call)
            l <- length(g1)
          } else { l <- 0 }
        }
        g0 <- grep("^ *saveImgFun\\(((\"Backup\\.RData\")|(BckUpFl))\\)", scrpt$call)
        if (length(g0)) {
          ghash <- grep("^ *#", scrpt$call, invert = TRUE)
          f0 <- function(x) { #x <- .obj[4]
            rs <- NA
            y <- c(grep(paste0("^ *", x, " *%<(o|c)%"), scrpt$call),
                   grep(paste0("^ *\\.obj *<- *unique\\(c\\(", x, "\\)\\)"), scrpt$call))
            if (length(y)) {
              y <- scrpt$row[min(y)]
              rs <- scrpt$row[g0][which(scrpt$row[g0] > y)][1]
            }
            return(rs)
          }
          #environment(f0) <- .GlobalEnv
          ok <- FALSE
          if (usePar) {
            parallel::clusterExport(parClust, list("g0", "scrpt"), envir = environment())
            tst <- try(setNames(parSapply(parClust, .obj, f0), .obj), silent = TRUE)
            #
            # NB:
            #   Recent change, now object names are added to .obj left, not right:
            #   > .obj <- unique(c(objNm, .obj)))
            #   thus we now have two ways to detect how far the script progressed:
            #    - Using all names in .obj, we can figure out how far we have gone at any run of the script from these inputs assuming the backup wasn't discarded.
            #    - Using just the first name in .obj for which we have a non NA value, we can figure out what was the last run chunk.
            #
            ok <- !("try-error" %in% class(tst))
          } else {
            tst <- try(setNames(sapply(.obj, f0), .obj), silent = TRUE)
          }
          if (!ok) { tst <- setNames(vapply(.obj, f0, as.integer(1)), .obj) }
          tst <- tst[which((vapply(names(tst), exists, TRUE))&(!is.na(tst)))]
          if (length(tst)) {
            m1 <- tst[1]
            o1 <- names(tst)[1]
            tst <- tst[order(tst, decreasing = TRUE)]
            m2 <- tst[1]
            o2 <- names(tst)[1]
            rs1 <- scrpt$row[ghash][which(scrpt$row[ghash] > m1)][1]
            rs2 <- scrpt$row[ghash][which(scrpt$row[ghash] > m2)][1]
            # rs1 should always be smaller or equal to rs2
            msg <- paste0("\n   FYI, the last remanent object which was added chronologically to .obj was \"", o1, "\"")
            if (o1 != o2) {
              msg <- paste0(msg,
                            ",\n   but it seems that at some point the script had also been run with this data beyond this point\n   (up to the creation of remanent object ",
                            o2, ")\n")
            } else {
              msg <- paste0(msg, ",\n   which seems to also be how far the script has been run with this data at any point.\n")
            }
            cat(msg)
            #system(paste0("open \"", ScriptPath, "\""))
            if (rs1 >= max(scrpt$row[g0])) {
              cat("\n   Backup analysis suggests that this backup had reached the end of the analysis, so there should be nothing more to run...\nBut maybe you want to re-run some parts without starting from scratch?\n")
              cat("   (opening script...)\n")
              rstudioapi::documentOpen(ScriptPath)
            } else {
              if (rs1 < rs2) {
                cat(paste0("\n   -> We thus suggest starting execution from either row ", rs1, " or ", rs2, "...\n"))
                cat(paste0("   (opening script at the line ", rs1, ")\n"))
              } else {
                cat(paste0("\n   -> We thus suggest starting execution from row ", rs1, "...\n"))
                cat("   (opening script at the corresponding line...)\n")
              }
              rstudioapi::documentOpen(ScriptPath, line = rs1)
            }
            #system(paste0("open \"", ScriptPath, "\""))
          }
        }
      } else {
        cat(" ... but it appears the file doesn't exist anymore...\n")
      }
    } else {
      cat(paste0("   FYI, the last object listed in .obj is \"", rev(.obj)[1], "\".\n"))
    }
  }
  if (exists("mySeed")) { set.seed(mySeed) }
  if (usePar) { parClust <<- parClust } # Export cluster
  cat("\nYou're good to go!\n")
}
