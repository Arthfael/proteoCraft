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

Configure <- function(updateOntologies = FALSE) {
  libPath <- as.data.frame(library()$results)
  libPath <- normalizePath(libPath$LibPath[match("proteoCraft", libPath$Package)], winslash = "/")
  proteoPath <- paste0(libPath, "/proteoCraft")
  homePath <- paste0(normalizePath(Sys.getenv("HOME"), winslash = "/"), "/R/proteoCraft")
  #
  # Create home directory
  if (!dir.exists(homePath)) {
    dir.create(homePath, recursive = TRUE)
    cat(paste0("Created package HOME in \"", homePath, "\"\n"))
  } else { paste0("Existing package HOME found in \"", homePath, "\"\n") }
  #
  # Default LC column definitions
  KolDef1 <- paste0(proteoPath, "/extdata/LC_columns.xlsx")
  KolDef2 <- paste0(homePath, "/LC_columns.xlsx")
  if ((file.exists(KolDef1))&&(!file.exists(KolDef2))) {
    tst <- try(file.copy(KolDef1, homePath, overwrite = FALSE), silent = TRUE)
    print(tst)
    if (file.exists(KolDef2)) { cat(" - Created default LC columns .xlsx table in HOME...\n") }
  } else { cat(" - Default LC columns Excel table already found in HOME.\n") }
  #
  # Default locations table
  locDirsFl1 <- paste0(proteoPath, "/extdata/Default_locations.xlsx")
  locDirsFl2 <- paste0(homePath, "/Default_locations.xlsx")
  writeDefltLoc <- function(locFile = locDirsFl2, dirTbl = locDirs) {
    # Should be rewritten with openxlsx2!!!
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, "Default_folders")
    openxlsx::writeData(wb, "Default_folders", dirTbl)
    openxlsx::addStyle(wb, "Default_folders",
                       openxlsx::createStyle(textDecoration = c("bold", "underline")),
                       1, 1:ncol(dirTbl), stack = TRUE)
    openxlsx::addStyle(wb, "Default_folders",
                       openxlsx::createStyle(fontName = "Consolas", textDecoration = "italic"),
                       1:nrow(dirTbl) + 1, match("Path", colnames(dirTbl)), stack = TRUE)
    openxlsx::addStyle(wb, "Default_folders",
                       openxlsx::createStyle(textDecoration = "italic"), 1:nrow(dirTbl) + 1,
                       match("Help", colnames(dirTbl)), stack = TRUE)
    openxlsx::setColWidths(wb, "Default_folders", match(c("Folder", "Path", "Help"), colnames(dirTbl)), c(25, 60, 150))
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
      locDirs <- rbind(locDirs, tmpDir[w,])
      writeDefltLoc(locDirsFl2, locDirs)
    } else { cat(".\n") }
  }
  #
  locDirs <- openxlsx2::read_xlsx(locDirsFl2)
  locDirs <- locDirs[which(!is.na(locDirs$Folder)),] # Precaution...
  #openxlsx2::xl_open(locDirsFl2)
  w <- which((locDirs$Path == "")|((!dir.exists(locDirs$Path)))&(!file.exists(locDirs$Path)))
  if (length(w)) {
    for (i in w) {
      if (locDirs$Folder[i] == "Python") {
        locDirs$Path[i] <- rstudioapi::selectFile(paste0("Select ", locDirs$Folder[i], " executable (.exe file)"),
                                                  path = "C:/PROGRA~1/*.exe")
      } else {
        locDirs$Path[i] <- rstudioapi::selectDirectory(gsub("( folder)+", " folder",
                                                            paste0("Select ", locDirs$Folder[i], " folder")))
      }
    }
    writeDefltLoc(locDirsFl2, locDirs)
  }
  # Sample solvent definitions - used by MatMet_WetLab()
  fl0 <- paste0(proteoPath, "/extdata/Sample_solvents.txt")
  fl1 <- paste0(homePath, "/Sample_solvents.txt")
  if (file.exists(fl1)) {
    opt0 <- readLines(fl0)
    opt1 <- readLines(fl1)
    opt1 <- unique(c(opt1, opt0))
    opt1 <- opt1[which(nchar(opt1) > 0)]
    opt1 <- opt1[which(!is.na(opt1))]
    write(opt1, fl1)
  } else {
    try(file.copy(fl0, homePath, overwrite = FALSE), silent = TRUE)
  }
  #
  # Also copy main data analysis scripts to home
  scrpts2 <- c(#"Regulation analysis - master script",
               "Regulation analysis - detailed script",
               "Regulation analysis - detailed script_pepOnly",
               "No replicates analysis - detailed script",
               "Histones_analysis_-_DiaNN_FragPipe_Skyline_or_alphaDIA_input",
               "Reload_renv_from_lock_file")
  extDr1 <- paste0(proteoPath, "/extdata/R scripts")
  fls <- paste0(extDr1, "/", scrpts2, ".R")
  for (fl in fls) {
    tst <- try(file.copy(fl, homePath, overwrite = TRUE), silent = TRUE)
    print(tst)
  }
  cat("Updated analysis scripts in HOME...")
  #
  # Check that we have python installed
  pyPath <- c()
  if ("Python" %in% locDirs$Folder) {
   pyTest <- !as.logical(try({
      #pyPath0 <- locDirs$Path[match("Python", locDirs$Folder)]
      #cat(pyPath <- gsub("/", "\\\\", pyPath0))
      # Update pip
      cmd <- "pip install --upgrade pip"
      system(cmd, show.output.on.console = FALSE)
      # Install sdrf-pipelines
      cmd <- "pip install sdrf-pipelines"
      system(cmd, show.output.on.console = FALSE)
    }, silent = TRUE))
  }
  #
  # Download ontologies relevant for writing an SDRF file
  #   This sometimes fails with `Error: C stack usage SOMEABSURDLYLARGENUMBER is too close to the limit` when other packages are loaded,
  #   I suspect because of a conflict between internally used functions sharing a same name and not called safely with package::function()
  #   Solution: run this on a cluster without loading said packages!
  #
  # The next step will be to build a shiny app (one more!) into SDRF_4_PRIDE.R to allow selectized choice of valid term(s) for each ontology.
  #
  ontFls <- c("Tissues", #"Species",
              "Modifications", "Cell_types", "Diseases",
              #"PRIDE",
              "MS_Quant_meth", "MS_label_meth", "MS_acq_meth", "MS_models",
              "DevStages_Fly", "DevStages_Zebrafish", "DevStages_Mouse", "DevStages_Human", "DevStages_Worms")
  ontFls <- setNames(paste0(homePath, "/", ontFls, ".csv"), ontFls)
  if ((sum(!file.exists(ontFls)))||(updateOntologies)) {
    if (!require(rols)) { pak::pkg_install("rols") }
    if (updateOntologies) { cat("Updating ontology files...\n") }
    tmpCl <- parallel::makeCluster(1, "SOCK")
    get_all_terms <- function(ontology, size = 500) {
      all_terms <- list()
      page <- 0L
      repeat {
        res <- httr::GET(
          sprintf("https://www.ebi.ac.uk/ols4/api/ontologies/%s/terms", ontology),
          query = list(size = size, page = page)
        )
        httr::stop_for_status(res)
        dat <- httr::content(res, as = "parsed", simplifyDataFrame = TRUE)
        terms <- dat$`_embedded`$terms
        if (length(terms) == 0) break
        w <- which(vapply(colnames(terms), function(x) { "data.frame" %in% class(terms[[x]]) }, TRUE))
        while (length(w)) {
          b <- lapply(w, function(i) {
            a <- terms[[i]]
            return(a)
          })
          b <- do.call(cbind, b)
          terms <- terms[, -w]
          terms[, colnames(b)] <- b
          w <- which(vapply(colnames(terms), function(x) { "data.frame" %in% class(terms[[x]]) }, TRUE))
        }
        page <- page + 1L
        all_terms[[page]] <- terms
      }
      all_terms <- do.call(plyr::rbind.fill, all_terms)
    }
    parallel::clusterExport(tmpCl, list("homePath", "get_all_terms", "updateOntologies"), envir = environment())
    invisible(parallel::clusterCall(tmpCl, function() {
      library(rols)
      library(httr)
      library(jsonlite)
      return()
    }))
    # - Species: too large, takes forever and not worth it: we already usually have it (taxID is basically it)
    # if ((!file.exists(ontFls["Species"]))||(updateOntologies)) {
    #   parallel::clusterCall(tmpCl, function() {
    #     # Get NCBITaxon root terms from OLS4
    #     sp_terms <- get_all_terms("NCBITaxon", size = 1000)
    #     sp <- sp_terms[, c("label", "description", "synonyms", "iri", "short_form")]
    #     w <- which(vapply(colnames(sp), function(x) { "list" %in% class(sp[[x]]) }, TRUE))
    #     if (length(w)) {
    #       for (i in w) {
    #         #sp[[i]] <- vapply(sp[[i]], paste, "", collapse = ";")
    #         sp[[i]] <- vapply(sp[[i]], function(x) { x <- c(unlist(x), "")[1] }, "") # We don't need all that clutter here
    #       }
    #     }
    #     write.csv(sp, ontFls["Species"], row.names = FALSE)
    #   })
    # }
    # - Tissue
    if ((!file.exists(ontFls["Tissues"]))||(updateOntologies)) {
      parallel::clusterCall(tmpCl, function() {
        #ol <- rols::Ontology("bto") # Unfortunately this is not currently supported
        # Get BTO root terms from OLS4
        bto_terms <- get_all_terms("bto", size = 1000)
        #btoTerms <- unique(bto_terms$label)
        #btoTerms <- btoTerms[which((!is.na(btoTerms))&(nchar(btoTerms) > 0))]
        bto <- bto_terms[, c("label", "description", "synonyms", "iri", "short_form")]
        w <- which(vapply(colnames(bto), function(x) { "list" %in% class(bto[[x]]) }, TRUE))
        if (length(w)) {
          for (i in w) {
            #bto[[i]] <- vapply(bto[[i]], paste, "", collapse = ";")
            bto[[i]] <- vapply(bto[[i]], function(x) { x <- c(unlist(x), "")[1] }, "") # We don't need all that clutter here
          }
        }
        write.csv(bto, ontFls["Tissues"], row.names = FALSE)
        return()
      })
    }
    # - Modifications
    if ((!file.exists(ontFls["Modifications"]))||(updateOntologies)) {
      parallel::clusterCall(tmpCl, function() {
        # Get MOD root terms from OLS4
        mod_terms <- get_all_terms("mod", size = 1000)
        mod <- mod_terms[, c("label", "description", "synonyms", "iri", "short_form")]
        w <- which(vapply(colnames(mod), function(x) { "list" %in% class(mod[[x]]) }, TRUE))
        if (length(w)) {
          for (i in w) {
            #mod[[i]] <- vapply(mod[[i]], paste, "", collapse = ";")
            mod[[i]] <- vapply(mod[[i]], function(x) { x <- c(unlist(x), "")[1] }, "") # We don't need all that clutter here
          }
        }
        write.csv(mod, ontFls["Modifications"], row.names = FALSE)
        return()
      })
    }
    # - Cell type
    if ((!file.exists(ontFls["Cell_types"]))||(updateOntologies)) {
      parallel::clusterCall(tmpCl, function() {
        # Get Cell Type root terms from OLS4
        cl_terms <- get_all_terms("CL", size = 1000)
        cl <- cl_terms[, c("label", "description", "synonyms", "iri", "short_form")]
        w <- which(vapply(colnames(cl), function(x) { "list" %in% class(cl[[x]]) }, TRUE))
        if (length(w)) {
          for (i in w) {
            #cl[[i]] <- vapply(cl[[i]], paste, "", collapse = ";")
            cl[[i]] <- vapply(cl[[i]], function(x) { x <- c(unlist(x), "")[1] }, "") # We don't need all that clutter here
          }
        }
        write.csv(cl, ontFls["Cell_types"], row.names = FALSE)
        return()
      })
    }
    # - Disease
    if ((!file.exists(ontFls["Diseases"]))||(updateOntologies)) {
      parallel::clusterCall(tmpCl, function() {
        # Get Disease root terms from OLS4
        dis_terms <- get_all_terms("DOID", size = 1000)
        dis <- dis_terms[, c("label", "description", "synonyms", "iri", "short_form")]
        w <- which(vapply(colnames(dis), function(x) { "list" %in% class(dis[[x]]) }, TRUE))
        if (length(w)) {
          for (i in w) {
            #dis[[i]] <- vapply(dis[[i]], paste, "", collapse = ";")
            dis[[i]] <- vapply(dis[[i]], function(x) { x <- c(unlist(x), "")[1] }, "") # We don't need all that clutter here
          }
        }
        write.csv(dis, ontFls["Diseases"], row.names = FALSE)
        return()
      })
    }
    # - PRIDE // obsolete
    # if ((!file.exists(ontFls["PRIDE"]))||(!file.exists(ontFls["MS_Quant_meth"]))||(updateOntologies)) {
    #   parallel::clusterCall(tmpCl, function() {
    #     pride_terms <- get_all_terms("PRIDE", size = 1000)
    #     pride <- pride_terms[, c("label", "description", "synonyms", "iri", "short_form")]
    #     w <- which(vapply(colnames(pride), function(x) { "list" %in% class(pride[[x]]) }, TRUE))
    #     if (length(w)) {
    #       for (i in w) {
    #         #pride[[i]] <- vapply(pride[[i]], paste, "", collapse = ";")
    #         pride[[i]] <- vapply(pride[[i]], function(x) { x <- c(unlist(x), "")[1] }, "") # We don't need all that clutter here
    #       }
    #     }
    #     get_term_children <- function(children_url) {
    #       res <- httr::GET(children_url)
    #       httr::stop_for_status(res)
    #       dat <- jsonlite::fromJSON(httr::content(res, "text"), flatten = TRUE)
    #       dat$`_embedded`$terms
    #     }
    #     quantMeth <- list()
    #     k <- 0L
    #     parent <- pride_terms[which(pride_terms$label %in% c("Proteomics data acquisition method", "Quantification method")),]
    #     repeat {
    #       if (!"_links.children.href" %in% colnames(parent)) {
    #         # first loop, dealing with the modified output of get_all_terms()
    #         children_url <- parent[, grep("^href\\.[0-9]+$", colnames(parent), value = TRUE)]
    #         wChildr <- grep("children", children_url[1,])
    #         children_url <- unique(children_url[, wChildr])
    #       } else {
    #         # unmodified results
    #         children_url <- parent$"_links.children.href"
    #       }
    #       childr <- lapply(children_url, function(url) { get_term_children(url) })
    #       childr <- do.call(plyr::rbind.fill, childr)
    #       childr[, c("label", "obo_id")]
    #       childr <- childr[grep("[Gg]el-based", childr$label, invert = TRUE),]
    #       if (!nrow(childr)) break
    #       k <- k + 1L
    #       wN <- which(!childr$has_children)
    #       wY <- which(childr$has_children)
    #       quantMeth[[k]] <- childr[, c("label", "iri")]
    #       if (!length(wY)) break
    #       parent <- childr[wY,]
    #     }
    #     quantMeth <- do.call(plyr::rbind.fill, quantMeth)
    #     if ((!file.exists(ontFls["PRIDE"]))||(updateOntologies)) {
    #       write.csv(pride, ontFls["PRIDE"], row.names = FALSE)
    #     }
    #     if ((!file.exists(ontFls["MS_Quant_meth"]))||(updateOntologies)) {
    #       write.csv(quantMeth, ontFls["MS_Quant_meth"], row.names = FALSE)
    #     }
    #     return()
    #   })
    # }
    # - PRIDE: MS labelling and quantitative methods
    if ((!file.exists(ontFls["MS_Quant_meth"]))||(!file.exists(ontFls["MS_label_meth"]))||(!file.exists(ontFls["MS_acq_meth"]))||(updateOntologies)) {
      # For this rols works, yay!
      parallel::clusterCall(tmpCl, function() {
        ol <- rols::Ontology("pride")
        ol@config$fileLocation <- "http://purl.obolibrary.org/obo/pride/releases/2025-02-17/pride.owl" # HARD FIX
        prideDigDeep <- function(x, myTrms) { #x <- 5
          #myTrms = labTrms
          trm_x <- myTrms[[x]]
          trms_x <- rols::children(trm_x)
          rsx <- data.frame(ID = trm_x@obo_id,
                            Name = trm_x@label,
                            Description = paste(unlist(trm_x@description), collapse = "; "))
          if ((!is.null(trms_x))&&(length(trms_x))) {
            rsi <- lapply(1:length(trms_x), function(i) { #i <- 1
              trm_i <- trms_x[[i]]
              trms_ij <- rols::children(trm_i)
              rs_i <- data.frame(ID = trm_i@obo_id,
                                 Name = trm_i@label,
                                 Description = paste(unlist(trm_i@description), collapse = "; "))
              if ((!is.null(trms_ij))&&(length(trms_ij))) {
                rs_j <- lapply(1:length(trms_ij), function(j) { #j <- 1
                  trm_j <- trms_ij[[j]]
                  trms_ijk <- rols::children(trm_j)
                  rs_j <- data.frame(ID = trm_j@obo_id,
                                     Name = trm_j@label,
                                     Description = paste(unlist(trm_j@description), collapse = "; "))
                  if ((!is.null(trms_ijk))&&(length(trms_ijk))) {
                    rs_k <- lapply(1:length(trms_ijk), function(k) { #k <- 1
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
        topLabTrm <- rols::Term(ol, "PRIDE_0000514")
        labTrms <- rols::children(topLabTrm)  # get all levels
        labTrms <- setNames(labTrms@x, vapply(labTrms@x, function(x) { x@label }, ""))
        labTrmsa <- lapply(1:length(labTrms), prideDigDeep, myTrms = labTrms)
        labTrmsa <- do.call(rbind, labTrmsa)
        labTrmsa <- aggregate(1:nrow(labTrmsa), list(labTrmsa$Name), function(x) {
          data.frame(ID = labTrmsa$ID[x[1]],
                     Description = labTrmsa$Description[x[1]])
        })
        colnames(labTrmsa)[1] <- "Name"
        labTrmsa[, c("ID", "Description")] <- do.call(rbind, labTrmsa$x)
        labTrmsa$x <- NULL
        if ((!file.exists(ontFls["MS_label_meth"]))||(updateOntologies)) {
          write.csv(labTrmsa, ontFls["MS_label_meth"], row.names = FALSE)
        }
        # Quantitative methods
        topMethTrm <- rols::Term(ol, "PRIDE_0000309")
        methTrms <- rols::children(topMethTrm)  # get all levels
        methTrms <- setNames(methTrms@x, vapply(methTrms@x, function(x) { x@label }, ""))
        methTrmsa <- lapply(1:length(methTrms), prideDigDeep, myTrms = methTrms)
        methTrmsa <- do.call(rbind, methTrmsa)
        methTrmsa <- aggregate(1:nrow(methTrmsa), list(methTrmsa$Name), function(x) {
          data.frame(ID = methTrmsa$ID[x[1]],
                     Description = methTrmsa$Description[x[1]])
        })
        colnames(methTrmsa)[1] <- "Name"
        methTrmsa[, c("ID", "Description")] <- do.call(rbind, methTrmsa$x)
        methTrmsa$x <- NULL
        if ((!file.exists(ontFls["MS_quant_meth"]))||(updateOntologies)) {
          write.csv(methTrmsa, ontFls["MS_quant_meth"], row.names = FALSE)
        }
        # Acquisition methods
        topMethTrm <- rols::Term(ol, "PRIDE_0000659")
        methTrms <- rols::children(topMethTrm)  # get all levels
        methTrms <- setNames(methTrms@x, vapply(methTrms@x, function(x) { x@label }, ""))
        methTrmsb <- lapply(1:length(methTrms), prideDigDeep, myTrms = methTrms)
        methTrmsb <- do.call(rbind, methTrmsb)
        methTrmsb <- aggregate(1:nrow(methTrmsb), list(methTrmsb$Name), function(x) {
          data.frame(ID = methTrmsb$ID[x[1]],
                     Description = methTrmsb$Description[x[1]])
        })
        colnames(methTrmsb)[1] <- "Name"
        methTrmsb[, c("ID", "Description")] <- do.call(rbind, methTrmsb$x)
        methTrmsb$x <- NULL
        if ((!file.exists(ontFls["MS_acq_meth"]))||(updateOntologies)) {
          write.csv(methTrmsb, ontFls["MS_acq_meth"], row.names = FALSE)
        }
      })
    }
    #
    # - MS instrument ontology
    if ((!file.exists(ontFls["MS_models"]))||(updateOntologies)) {
      # For this rols works, yay!
      parallel::clusterCall(tmpCl, function() {
        ol <- rols::Ontology("ms")
        instrMod <- rols::Term(ol, "MS:1000031")
        vendors <- rols::children(instrMod)  # get all levels
        vendorsLab <- lapply(vendors@x, function(x) {
          x@label
        })
        vendorsInstr <- lapply(1:length(vendorsLab), function(x) { #x <- 2 #x <- 5 #x <- 7 #x <- 10
          trm <- rols::Term(ol, names(vendorsLab)[[x]])
          models <- rols::children(trm)
          if (!is.null(models)) {
            rs <- vapply(models@x, function(y) { y@label }, "")
            rs2 <- setNames(lapply(1:length(rs), function(y) { #y <- g[1]
              try({
                trm2 <- rols::Term(ol, names(rs)[y])
                mods <- rols::children(trm2)
                mods <- vapply(mods@x, function(z) { z@label }, "")
              }, silent = TRUE)
            }), names(rs))
            w <- which(vapply(rs2, function(x) { "try-error" %in% class(x) }, TRUE))
            if (length(w)) { rs2[names(rs)[w]] <- rs[names(rs)[w]] }
            rs <- unlist(rs2)
            names(rs) <- gsub("^MS:[0-9]+\\.", "", names(rs))
            rs <- rs[grep("^MS:[0-9]+$", names(rs))]
            rs <- sort(rs)
            if (length(rs)) {
              rs <- data.frame(Vendor = gsub(" instrument model$", "", vendorsLab[[x]]),
                               Instrument = rs)
              return(rs)
            } else { return(NULL) }
          } else { return(NULL) }
        })
        vendorsInstr <- do.call(rbind, vendorsInstr)
        vendorsInstr$Term <- rownames(vendorsInstr)
        rownames(vendorsInstr) <- NULL
        write.csv(vendorsInstr, ontFls["MS_models"], row.names = FALSE)
        return()
      })
    }
    #
    # Development stage ontologies
    # For now only Zebrafish (Danio rerio) and FlyBase (Drosophila melanogaster) seem to be supported by the SDRF format... weird...
    # I'll still get human, mouse and WormBase (C. Elegans + others) too.
    devStgFun <- function(Ontology, parentTerm, fileName) {
      ol <- rols::Ontology(Ontology)
      dvStgs <- rols::Term(ol, parentTerm)
      dvStgs2 <- rols::children(dvStgs)  # get all levels
      dvStgs2 <- setNames(dvStgs2@x, vapply(dvStgs2@x, function(x) { x@label }, ""))
      dvStgs2a <- lapply(1:length(dvStgs2), function(x) { #x <- 1 #x <- 98
        trm <- dvStgs2[[x]]
        trms <- rols::children(trm)
        if ((!is.null(trms))&&(length(trms))) {
          rs <- lapply(1:length(trms), function(i) { #i <- 1
            trm <- trms[[i]]
            rsi <- data.frame(ID = trm@obo_id,
                              Name = trm@label,
                              Description = paste(unlist(trm@description), collapse = "; "))
            return(rsi)
          })
          rs <- do.call(rbind, rs)
        } else {
          rs <- data.frame(ID = trm@obo_id,
                           Name = trm@label,
                           Description = paste(unlist(trm@description), collapse = "; "))
        }
        return(rs)
      })
      dvStgs2a <- do.call(rbind, dvStgs2a)
      write.csv(dvStgs2a, ontFls[[fileName]], row.names = FALSE)
      return()
    }
    parallel::clusterExport(tmpCl, "devStgFun", envir = environment())
    # Fly
    if ((!file.exists(ontFls[["DevStages_Fly"]]))||(updateOntologies)) {
      parallel::clusterCall(tmpCl, function() { devStgFun("fbdv", "FBdv_00005259", "DevStages_Fly") })
    }
    # Zebrafish
    if ((!file.exists(ontFls[["DevStages_Zebrafish"]]))||(updateOntologies)) {
      parallel::clusterCall(tmpCl, function() { devStgFun("zfs", "ZFS_0100000", "DevStages_Zebrafish") })
    }
    # Mouse
    if ((!file.exists(ontFls[["DevStages_Mouse"]]))||(updateOntologies)) {
      parallel::clusterCall(tmpCl, function() { devStgFun("mmusdv", "MmusDv_0000000", "DevStages_Mouse") })
    }
    # Human
    if ((!file.exists(ontFls[["DevStages_Human"]]))||(updateOntologies)) {
      parallel::clusterCall(tmpCl, function() { devStgFun("hsapdv", "HsapDv_0000000", "DevStages_Human") })
    }
    # C. elegans
    if ((!file.exists(ontFls[["DevStages_Worms"]]))||(updateOntologies)) {
      parallel::clusterCall(tmpCl, function() { devStgFun("wbls", "WBls_0000075", "DevStages_Worms") })
    }
    parallel::stopCluster(tmpCl)
  }
  #
  #
  cat("Done\n")
  return()
}
