#' Configure
#'
#' @description
#' Configuration function to run once at least after first installation, and once in a while after updating the package.
#' 
#' @param updateOntologies Should we update ontologies? This is a very slow process, so false by default.
#' 
#' @details
#' This function:
#'  - creates default LC column definitions
#'  - creates a default sample solvents
#'  - creates and edits with user input default folder locations\
#'  - check for a Python installation and if one is found, pip installs sdrf-pipelines
#'  - moves some files from the package to its subfolder in HOME
#'  - downloads some ontologies used by some scripts (mainly the SDRF editor)
#' 
#' @returns
#' This function does not return anything.
#' 
#' @examples
#' proteoCraft::Configure()
#' 
#' @export

Configure <- function(updateOntologies = FALSE) { #updateOntologies = TRUE
  libPath <- as.data.frame(library()$results)
  libPath <- normalizePath(libPath$LibPath[match("proteoCraft", libPath$Package)], winslash = "/")
  proteoPath <- paste0(libPath, "/proteoCraft")
  homePath <- paste0(normalizePath(Sys.getenv("HOME"), winslash = "/"), "/R/proteoCraft")
  if (!dir.exists(homePath)) {
    dir.create(homePath, recursive = TRUE)
    cat(paste0("Created package HOME in \"", homePath, "\"\n"))
  } else {
    paste0("Existing package HOME found in \"", homePath, "\"\n")
  }
  KolDef1 <- paste0(proteoPath, "/extdata/LC_columns.xlsx")
  KolDef2 <- paste0(homePath, "/LC_columns.xlsx")
  if ((file.exists(KolDef1)) && (!file.exists(KolDef2))) {
    tst <- try(file.copy(KolDef1, homePath, overwrite = FALSE), silent = TRUE)
    print(tst)
    if (file.exists(KolDef2)) { cat(" - Created default LC columns .xlsx table in HOME...\n") }
  } else {
    cat(" - Default LC columns Excel table already found in HOME.\n")
  }
  locDirsFl1 <- paste0(proteoPath, "/extdata/Default_locations.xlsx")
  locDirsFl2 <- paste0(homePath, "/Default_locations.xlsx")
  writeDefltLoc <- function(locFile = locDirsFl2, dirTbl = locDirs) {
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, "Default_folders")
    openxlsx::writeData(wb, "Default_folders", dirTbl)
    openxlsx::addStyle(wb, "Default_folders",
                       openxlsx::createStyle(textDecoration = c("bold", "underline")),
                       1L,
                       1L:ncol(dirTbl),
                       stack = TRUE)
    openxlsx::addStyle(wb, "Default_folders",
                       openxlsx::createStyle(fontName = "Consolas", textDecoration = "italic"),
                       1L:nrow(dirTbl) + 1L,
                       match("Path", colnames(dirTbl)), stack = TRUE)
    openxlsx::addStyle(wb, "Default_folders", openxlsx::createStyle(textDecoration = "italic"), 
                       1L:nrow(dirTbl) + 1L, match("Help", colnames(dirTbl)), 
                       stack = TRUE)
    openxlsx::setColWidths(wb,
                           "Default_folders",
                           match(c("Folder", "Path", "Help"), colnames(dirTbl)),
                           c(25L, 60L, 150L))
    openxlsx::saveWorkbook(wb, locFile, overwrite = TRUE)
  }
  if ((file.exists(locDirsFl1))&&(!file.exists(locDirsFl2))) {
    tst <- try(file.copy(locDirsFl1, homePath, overwrite = FALSE), silent = TRUE)
    print(tst)
    if (file.exists(locDirsFl2)) {
      cat(" - Created default directories locations .xlsx table in HOME...\n")
    }
  } else {
    cat(" - Default directories locations Excel table already found in HOME")
    tmpDir <- openxlsx2::read_xlsx(locDirsFl1)
    locDirs <- openxlsx2::read_xlsx(locDirsFl2)
    w <- which(!tmpDir$Folder %in% locDirs$Folder)
    if (length(w)) {
      cat("\n   ... but we will append newly created folder definitions.\n")
      locDirs <- rbind(locDirs, tmpDir[w, ])
      writeDefltLoc(locDirsFl2, locDirs)
    } else { cat(".\n") }
  }
  locDirs <- openxlsx2::read_xlsx(locDirsFl2)
  locDirs <- locDirs[which(!is.na(locDirs$Folder)), ]
  w <- which((locDirs$Path == "")|(!file.exists(locDirs$Path)))
  if (length(w)) {
    for (i in w) {
      locDirs$Path[i] <- if (locDirs$Folder[i] == "Python") {
        rstudioapi::selectFile(paste0("Select ", locDirs$Folder[i], " executable (.exe file)"),
                               path = "C:/PROGRA~1/*.exe")
      } else {
        rstudioapi::selectDirectory(gsub("( folder)+", " folder", paste0("Select ", locDirs$Folder[i], " folder")))
      }
    }
    writeDefltLoc(locDirsFl2, locDirs)
  }
  g <- grep("^~/", locDirs$Path)
  if (length(g)) {
    locDirs$Path[g] <- base::path.expand(locDirs$Path[g])
    writeDefltLoc(locDirsFl2, locDirs)
  }
  #
  fl0 <- paste0(proteoPath, "/extdata/Sample_solvents.txt")
  fl1 <- paste0(homePath, "/Sample_solvents.txt")
  if (file.exists(fl1)) {
    opt0 <- readLines(fl0)
    opt1 <- readLines(fl1)
    opt1 <- unique(c(opt1, opt0))
    opt1 <- opt1[which(nchar(opt1) > 0L)]
    opt1 <- opt1[which(!is.na(opt1))]
    write(opt1, fl1)
  } else {
    try(file.copy(fl0, homePath, overwrite = FALSE), silent = TRUE)
  }
  scrpts2 <- c("Regulation analysis - detailed script", "Regulation analysis - detailed script_pepOnly", 
               "No replicates analysis - detailed script", "Histones_analysis_-_DiaNN_FragPipe_Skyline_or_alphaDIA_input", 
               "Reload_renv_from_lock_file")
  extDr1 <- paste0(proteoPath, "/extdata/R scripts")
  fls <- paste0(extDr1, "/", scrpts2, ".R")
  for (fl in fls) {
    tst <- try(file.copy(fl, homePath, overwrite = TRUE), 
               silent = TRUE)
    print(tst)
  }
  cat("Updated analysis scripts in HOME...")
  pyPath <- c()
  if ("Python" %in% locDirs$Folder) {
    try(.pyConfig(), silent = TRUE)
  }
  #
  locScriptsDir <- gsub("/[^/]+$", "/proteoCraft_localScripts", locDirs$Path[match("Temporary folder", locDirs$Folder)])
  if (!dir.exists(locScriptsDir)) { dir.create(locScriptsDir) }
  write(c("Save here any script (or a windows shortcut to any script) which you want to run automatically at the beginning of the main analysis scripts.",
          "For instance, as the developer, I keep in this folder a script to source my local working versions of the package's functions,",
          "so that whenever I am testing the workflows I use the latest versions, not those from the currently installed version of the package.",
          "You could do something similar if you decide to modify this package/fix bugs (but PLEASE report them too!)",
          "",
          "These scripts will be run in alphabetic order! If order is important, make sure that file names provide the correct one!",
          ""),
        paste0(locScriptsDir, "/README.txt"))
  #
  ontFls <- c("Tissues", "Modifications", "Cell_types", "Diseases", 
              "MS_quant_meth", "MS_label_meth", "MS_acq_meth", "MS_models", 
              "DevStages_Fly", "DevStages_Zebrafish", "DevStages_Mouse", 
              "DevStages_Human", "DevStages_Worms")
  ontFls <- setNames(paste0(homePath, "/", ontFls, ".csv"),
                     ontFls)
  if ((sum(!file.exists(ontFls)))||(updateOntologies)) {
    if (!require(rols, quietly = TRUE)) { pak::pak("rols") }
    require(rols)
    if (updateOntologies) { cat("Updating ontology files...\n") }
    get_all_terms <- function(ontology, size = 500L) {
      all_terms <- list()
      page <- 0L
      repeat {
        res <- httr::GET(sprintf("https://www.ebi.ac.uk/ols4/api/ontologies/%s/terms", ontology),
                         query = list(size = size, page = page))
        httr::stop_for_status(res)
        dat <- httr::content(res, as = "parsed", simplifyDataFrame = TRUE)
        terms <- dat$`_embedded`$terms
        if (!length(terms)) 
          break
        w <- which(vapply(colnames(terms), function(x) { is.data.frame(terms[[x]]) }, TRUE))
        while (length(w)) {
          b <- lapply(w, function(i) { terms[[i]] })
          b <- do.call(cbind, b)
          terms <- terms[, -w]
          terms[, colnames(b)] <- b
          w <- which(vapply(colnames(terms), function(x) { is.data.frame(terms[[x]]) }, TRUE))
        }
        page <- page + 1L
        all_terms[[page]] <- terms
      }
      return(do.call(plyr::rbind.fill, all_terms))
    }
    ontoDF <- data.frame(Term = c("bto", "mod", "CL", "DOID"),
                         Name = c("Tissues", "Modifications", "Cell_types", "Diseases"))
    for (i in 1L:nrow(ontoDF)) { #i <- 1L
      fl <- ontFls[ontoDF$Name[i]]
      #print(fl)
      if ((!file.exists(fl))||(updateOntologies)) {
        my_terms <- get_all_terms(ontoDF$Term[i], size = 1000L)
        my_terms <- my_terms[, c("label", "description", "synonyms", "iri", "short_form")]
        w <- which(vapply(colnames(my_terms), function(x) { is.list(my_terms[[x]]) }, TRUE))
        if (length(w)) { for (i in w) { my_terms[[i]] <- vapply(my_terms[[i]], function(x) { c(unlist(x), "")[1L] }, "") } }
        data.table::fwrite(my_terms, fl, row.names = FALSE)
        rm(fl)
      }
    }
    if ((!file.exists(ontFls["MS_quant_meth"]))||
        (!file.exists(ontFls["MS_label_meth"]))||
        (!file.exists(ontFls["MS_acq_meth"]))||
        (updateOntologies)) {
      pride <- rols::olsOntology("pride")
      prideTrms <- rols::olsTerms(pride)
      prideDigDeep <- function(x, myTrms) {
        trm_x <- myTrms[[x]]
        trms_x <- rols::children(trm_x)
        rsx <- data.frame(ID = trm_x@obo_id,
                          Name = trm_x@label, 
                          Description = paste(unlist(trm_x@description), collapse = "; "))
        if ((!is.null(trms_x))&&(length(trms_x))) {
          rsi <- lapply(1L:length(trms_x), function(i) {
            trm_i <- trms_x[[i]]
            trms_ij <- rols::children(trm_i)
            rs_i <- data.frame(ID = trm_i@obo_id, Name = trm_i@label, 
                               Description = paste(unlist(trm_i@description), collapse = "; "))
            if ((!is.null(trms_ij)) && (length(trms_ij))) {
              rs_j <- lapply(1L:length(trms_ij), function(j) {
                trm_j <- trms_ij[[j]]
                trms_ijk <- rols::children(trm_j)
                rs_j <- data.frame(ID = trm_j@obo_id, 
                                   Name = trm_j@label,
                                   Description = paste(unlist(trm_j@description), collapse = "; "))
                if ((!is.null(trms_ijk)) && (length(trms_ijk))) {
                  rs_k <- lapply(1L:length(trms_ijk), 
                                 function(k) {
                                   trm_k <- trms_ijk[[k]]
                                   data.frame(ID = trm_k@obo_id, 
                                              Name = trm_k@label,
                                              Description = paste(unlist(trm_k@description), collapse = "; "))
                                 })
                  rs_k <- do.call(rbind, rs_k)
                  rs_j <- rbind(rs_j, rs_k)
                }
                return(rs_j)
              })
              rs_j <- do.call(rbind, rs_j)
              rs_i <- rbind(rs_i, rs_j)
            }
            return(rs_i)
          })
          rsi <- do.call(rbind, rsi)
          rsx <- rbind(rsx, rsi)
        }
        return(rsx)
      }
      # Labelling methods
      labMeth <- prideTrms[["PRIDE:0000514"]]
      labMethChildren <- rols::children(labMeth)
      labMethChildren <- setNames(labMethChildren@x, vapply(labMethChildren@x, function(x) { x@label }, ""))
      labMethChildren <- lapply(1L:length(labMethChildren), prideDigDeep, 
                                myTrms = labMethChildren)
      labMethChildren <- do.call(rbind, labMethChildren)
      labMethChildren <- aggregate(1L:nrow(labMethChildren), list(labMethChildren$Name), function(x) {
        data.frame(ID = labMethChildren$ID[x[1L]],
                   Description = labMethChildren$Description[x[1L]])
      })
      colnames(labMethChildren)[1L] <- "Name"
      labMethChildren[, c("ID", "Description")] <- do.call(rbind, labMethChildren$x)
      labMethChildren$x <- NULL
      if ((!file.exists(ontFls["MS_label_meth"]))||(updateOntologies)) {
        write.csv(labMethChildren,
                  ontFls["MS_label_meth"], 
                  row.names = FALSE)
      }
      # Quantification methods
      quantMeth <- prideTrms[["PRIDE:0000309"]]
      quantMethChildren <- rols::children(quantMeth)
      quantMethChildren <- setNames(quantMethChildren@x, vapply(quantMethChildren@x, function(x) { x@label }, ""))
      quantMethChildren <- lapply(1L:length(quantMethChildren), prideDigDeep, 
                                  myTrms = quantMethChildren)
      quantMethChildren <- do.call(rbind, quantMethChildren)
      quantMethChildren <- aggregate(1L:nrow(quantMethChildren), list(quantMethChildren$Name), function(x) {
        data.frame(ID = quantMethChildren$ID[x[1L]],
                   Description = quantMethChildren$Description[x[1L]])
      })
      colnames(quantMethChildren)[1L] <- "Name"
      quantMethChildren[, c("ID", "Description")] <- do.call(rbind, quantMethChildren$x)
      quantMethChildren$x <- NULL
      #quantMethChildren <- quantMethChildren[which(!quantMethChildren$ID %in% labMethChildren$ID),]
      # (there is overlap with the labelling method term, not sure whether we should segregate them or not)
      if ((!file.exists(ontFls["MS_quant_meth"]))||(updateOntologies)) {
        write.csv(quantMethChildren,
                  ontFls["MS_quant_meth"], 
                  row.names = FALSE)
      }
      # Quantification methods
      acqMeth <- prideTrms[["PRIDE:0000659"]]
      acqMethChildren <- rols::children(acqMeth)
      acqMethChildren <- setNames(acqMethChildren@x, vapply(acqMethChildren@x, function(x) { x@label }, ""))
      acqMethChildren <- lapply(1L:length(acqMethChildren), prideDigDeep, 
                                myTrms = acqMethChildren)
      acqMethChildren <- do.call(rbind, acqMethChildren)
      acqMethChildren <- aggregate(1L:nrow(acqMethChildren), list(acqMethChildren$Name), function(x) {
        data.frame(ID = acqMethChildren$ID[x[1L]],
                   Description = acqMethChildren$Description[x[1L]])
      })
      colnames(acqMethChildren)[1L] <- "Name"
      acqMethChildren[, c("ID", "Description")] <- do.call(rbind, acqMethChildren$x)
      acqMethChildren$x <- NULL
      if ((!file.exists(ontFls["MS_acq_meth"]))||(updateOntologies)) {
        write.csv(acqMethChildren,
                  ontFls["MS_acq_meth"], 
                  row.names = FALSE)
      }
    }
    if ((!file.exists(ontFls["MS_models"]))||(updateOntologies)) {
      MS_onto <- rols::olsOntology("ms")
      MS_trms <- rols::olsTerms(MS_onto)
      instrMod <- MS_trms[["MS:1000031"]]
      vendors <- rols::children(instrMod)
      vendorsLab <- lapply(vendors@x, function(x) {
        x@label
      })
      vendorsInstr <- lapply(1L:length(vendorsLab), function(x) { #x <- 1L
        trm <- MS_trms[[names(vendorsLab)[[x]]]]
        models <- rols::children(trm)
        if (is.null(models)) { return(NULL) }
        rs <- vapply(models@x, function(y) {
          y@label
        }, "")
        rs2 <- setNames(lapply(1L:length(rs), function(y) {
          try({
            trm2 <- MS_trms[[names(rs)[y]]]
            mods <- rols::children(trm2)
            mods <- vapply(mods@x, function(z) {
              z@label
            }, "")
          }, silent = TRUE)
        }), names(rs))
        w <- which(vapply(rs2, inherits, TRUE, "try-error"))
        if (length(w)) {
          rs2[names(rs)[w]] <- rs[names(rs)[w]]
        }
        rs <- unlist(rs2)
        names(rs) <- gsub("^MS:[0-9]+\\.", "", names(rs))
        rs <- rs[grep("^MS:[0-9]+$", names(rs))]
        rs <- sort(rs)
        if (!length(rs)) { return(NULL) }
        rs <- data.frame(Vendor = gsub(" instrument model$", "", vendorsLab[[x]]),
                         Instrument = rs)
        return(rs)
      })
      vendorsInstr <- do.call(rbind, vendorsInstr)
      vendorsInstr$Term <- rownames(vendorsInstr)
      rownames(vendorsInstr) <- NULL
      write.csv(vendorsInstr, ontFls["MS_models"], row.names = FALSE)
    }
    devStgFun <- function(Ontology, parentTerm, fileName) {
      #Ontology = "mmusdv"; parentTerm = "MmusDv_0000000"
      onto <- rols::olsOntology(Ontology)
      trms <- rols::olsTerm(onto, parentTerm)
      dvStgs <- rols::children(trms)
      dvStgs <- setNames(dvStgs@x, vapply(dvStgs@x, 
                                          function(x) {
                                            x@label
                                          }, ""))
      dvStgsa <- lapply(1L:length(dvStgs), function(x) {
        trm <- dvStgs[[x]]
        trms <- rols::children(trm)
        if ((!is.null(trms)) && (length(trms))) {
          rs <- lapply(1L:length(trms), function(i) {
            trm <- trms[[i]]
            rsi <- data.frame(ID = trm@obo_id, Name = trm@label, 
                              Description = paste(unlist(trm@description), collapse = "; "))
            return(rsi)
          })
          rs <- do.call(rbind, rs)
        } else {
          rs <- data.frame(ID = trm@obo_id, Name = trm@label, 
                           Description = paste(unlist(trm@description), collapse = "; "))
        }
        return(rs)
      })
      dvStgsa <- do.call(rbind, dvStgsa)
      write.csv(dvStgsa, ontFls[[fileName]], row.names = FALSE)
      return()
    }
    ontoDF2 <- data.frame(Ontology = c("fbdv", "zfs", "mmusdv", "hsapdv", "wbls"),
                          Term = c("FBdv_00005259", "ZFS_0100000", "MmusDv_0000000", "HsapDv_0000000", "WBls_0000075"),
                          FileName = paste0("DevStages_", c("Fly", "Zebrafish", "Mouse", "Human", "Worms")))
    ontoDF2$File <- ontFls[ontoDF2$FileName]
    ontoDF2$Do <- (!file.exists(ontoDF2$File))|(updateOntologies)
    w <- which(ontoDF2$Do)
    if (length(w)) {
      for (i in w) {
        try(devStgFun(ontoDF2$Ontology[i], ontoDF2$Term[i], ontoDF2$FileName[i]), silent = TRUE)
      }
    }
  }
  cat("Done\n")
  return()
}
