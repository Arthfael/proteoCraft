#' Configure
#'
#' @description
#' Configuration function which will make sure some files from the package are moved to its subfolder in HOME upon first installation.
#' 
#' @param updateOntologies Should we update ontologies? This is a very slow process, so false by default.
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
  locDef1 <- paste0(proteoPath, "/extdata/Default_locations.xlsx")
  locDef2 <- paste0(homePath, "/Default_locations.xlsx")
  if ((file.exists(locDef1))&&(!file.exists(locDef2))) {
    tst <- try(file.copy(locDef1, homePath, overwrite = FALSE), silent = TRUE)
    print(tst)
    if (file.exists(locDef2)) {
      cat(" - Created default directories locations .xlsx table in HOME...\n")
    }
  } else {
    cat(" - Default directories locations Excel table already found in HOME")
    tmp1 <- openxlsx::read.xlsx(locDef1)
    tmp2 <- openxlsx::read.xlsx(locDef2)
    w <- which(!tmp1$Folder %in% tmp2$Folder)
    if (length(w)) {
      cat("\n   ... but we will append newly created folder definitions.\n")
      tmp2 <- rbind(tmp2, tmp1[w,])
      wb <- openxlsx::createWorkbook()
      openxlsx::addWorksheet(wb, "Default_folders")
      openxlsx::writeData(wb, "Default_folders", tmp2)
      openxlsx::addStyle(wb, "Default_folders",
                         openxlsx::createStyle(textDecoration = c("bold", "underline")),
                         1, 1:ncol(tmp2), stack = TRUE)
      openxlsx::addStyle(wb, "Default_folders",
                         openxlsx::createStyle(fontName = "Consolas", textDecoration = "italic"),
                         1:nrow(tmp2) + 1, match("Path", colnames(tmp2)), stack = TRUE)
      openxlsx::addStyle(wb, "Default_folders",
                         openxlsx::createStyle(textDecoration = "italic"), 1:nrow(tmp2) + 1,
                         match("Help",colnames(tmp2)), stack = TRUE)
      openxlsx::setColWidths(wb, "Default_folders", match(c("Folder", "Path", "Help"), colnames(tmp2)), c(25, 60, 150))
      openxlsx::saveWorkbook(wb, locDef2, overwrite = TRUE)
    } else { cat(".\n") }
  }
  #
  tmp2 <- openxlsx::read.xlsx(locDef2)
  W <- which((tmp2$Path == "")|(!dir.exists(tmp2$Path)))
  if (length(W)) {
    for (w in W) {
      tmp2$Path[w] <- rstudioapi::selectDirectory(gsub("( folder)+", " folder",
                                                       paste0("Select ", tmp2$Folder[w], " folder")))
    }
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, "Default_folders")
    openxlsx::writeData(wb, "Default_folders", tmp2)
    openxlsx::addStyle(wb, "Default_folders",
                       openxlsx::createStyle(textDecoration = c("bold", "underline")),
                       1, 1:ncol(tmp2), stack = TRUE)
    openxlsx::addStyle(wb, "Default_folders",
                       openxlsx::createStyle(fontName = "Consolas", textDecoration = "italic"),
                       1:nrow(tmp2) + 1, match("Path", colnames(tmp2)), stack = TRUE)
    openxlsx::addStyle(wb, "Default_folders",
                       openxlsx::createStyle(textDecoration = "italic"),
                       1:nrow(tmp2) + 1, match("Help", colnames(tmp2)), stack = TRUE)
    openxlsx::setColWidths(wb, "Default_folders", match(c("Folder", "Path", "Help"), colnames(tmp2)), c(25, 60, 150))
    openxlsx::saveWorkbook(wb, locDef2, overwrite = TRUE)
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
  # Also copy analysis scripts to home
  scrpts2 <- c("Regulation analysis - master script",
               "Regulation analysis - detailed script",
               "Regulation analysis - detailed script_pepOnly",
               "No replicates analysis - detailed script",
               "Reload_renv_from_lock_file")
  extDr1 <- paste0(proteoPath, "/extdata/R scripts")
  fls <- paste0(extDr1, "/", scrpts2, ".R")
  for (fl in fls) {
    tst <- try(file.copy(fl, homePath, overwrite = TRUE), silent = TRUE)
    print(tst)
  }
  cat("Updated analysis scripts in HOME...")
  #
  # Download ontologies relevant for writing an SDRF file
  #   This sometimes fail with `Error: C stack usage SOMEABSURDLYLARGENUMBER is too close to the limit` when other packages are loaded,
  #   I suspect because of a conflict between internally used functions sharing a same name and not called with package::...
  #   Solution: run this on a cluster without loading those packages!
  #
  # The next step will be to build a shiny app (one more!) into SDRF_4_PRIDE.R to allow selectized choice of valid term(s) for each ontology.
  #
  if (!require(rols)) { pak::pkg_install("rols") }
  tmpCl <- parallel::makeCluster(1, "SOCK")
  get_all_terms <- function(ontology, size = 500) {
    all_terms <- list()
    page <- 0L
    repeat {
      res <- GET(
        sprintf("https://www.ebi.ac.uk/ols4/api/ontologies/%s/terms", ontology),
        query = list(size = size, page = page)
      )
      stop_for_status(res)
      dat <- content(res, as = "parsed", simplifyDataFrame = TRUE)
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
  parallel::clusterExport(tmpCl, list("homePath", "get_all_terms"), envir = environment())
  invisible(parallel::clusterCall(tmpCl, function() {
    library(rols)
    library(httr)
    library(jsonlite)
    return()
  }))
  # - Species: too large, takes forever and not worth it: we already usually have it (taxID is basically it)
  # if ((!file.exists(paste0(homePath, "/Species.csv")))||(updateOntologies)) {
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
  #     write.csv(sp, paste0(homePath, "/Species.csv"), row.names = FALSE)
  #   })
  # }
  # - Tissue
  if ((!file.exists(paste0(homePath, "/Tissues.csv")))||(updateOntologies)) {
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
      write.csv(bto, paste0(homePath, "/Tissues.csv"), row.names = FALSE)
    })
  }
  # - Modifications
  if ((!file.exists(paste0(homePath, "/Modifications.csv")))||(updateOntologies)) {
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
      write.csv(mod, paste0(homePath, "/Modifications.csv"), row.names = FALSE)
    })
  }
  # - Cell type
  if ((!file.exists(paste0(homePath, "/Cell_types.csv")))||(updateOntologies)) {
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
      write.csv(cl, paste0(homePath, "/Cell_types.csv"), row.names = FALSE)
    })
  }
  # - Disease
  if ((!file.exists(paste0(homePath, "/Diseases.csv")))||(updateOntologies)) {
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
      write.csv(dis, paste0(homePath, "/Diseases.csv"), row.names = FALSE)
    })
  }
  # - PRIDE and quant methods
  if ((!file.exists(paste0(homePath, "/PRIDE.csv")))||(!file.exists(paste0(homePath, "/MS_Quant_meth.csv")))||(updateOntologies)) {
    parallel::clusterCall(tmpCl, function() {
      # Get Disease root terms from OLS4
      pride_terms <- get_all_terms("PRIDE", size = 1000)
      pride <- pride_terms[, c("label", "description", "synonyms", "iri", "short_form")]
      w <- which(vapply(colnames(pride), function(x) { "list" %in% class(pride[[x]]) }, TRUE))
      if (length(w)) {
        for (i in w) {
          #pride[[i]] <- vapply(pride[[i]], paste, "", collapse = ";")
          pride[[i]] <- vapply(pride[[i]], function(x) { x <- c(unlist(x), "")[1] }, "") # We don't need all that clutter here
        }
      }
      get_term_children <- function(children_url) {
        res <- httr::GET(children_url)
        stop_for_status(res)
        dat <- jsonlite::fromJSON(httr::content(res, "text"), flatten = TRUE)
        dat$`_embedded`$terms
      }
      quantMeth <- list()
      k <- 0L
      parent <- pride_terms[which(pride_terms$label %in% c("Proteomics data acquisition method", "Quantification method")),]
      repeat {
        if (!"_links.children.href" %in% colnames(parent)) {
          # first loop, dealing with the modified output of get_all_terms()
          children_url <- parent[, grep("^href\\.[0-9]+$", colnames(parent), value = TRUE)]
          wChildr <- grep("children", children_url[1,])
          children_url <- unique(children_url[, wChildr])
        } else {
          # unmodified results
          children_url <- parent$"_links.children.href"
        }
        childr <- lapply(children_url, function(url) { get_term_children(url) })
        childr <- do.call(plyr::rbind.fill, childr)
        childr[, c("label", "obo_id")]
        childr <- childr[grep("[Gg]el-based", childr$label, invert = TRUE),]
        if (!nrow(childr)) break
        k <- k + 1L
        wN <- which(!childr$has_children)
        wY <- which(childr$has_children)
        quantMeth[[k]] <- childr[, c("label", "iri")]
        if (!length(wY)) break
        parent <- childr[wY,]
      }
      quantMeth <- do.call(plyr::rbind.fill, quantMeth)
      if ((!file.exists(paste0(homePath, "/PRIDE.csv")))||(updateOntologies)) {
        write.csv(pride, paste0(homePath, "/PRIDE.csv"), row.names = FALSE)
      }
      if ((!file.exists(paste0(homePath, "/MS_Quant_meth.csv")))||(updateOntologies)) {
        write.csv(quantMeth, paste0(homePath, "/MS_Quant_meth.csv"), row.names = FALSE)
      }
    })
  }
  # - MS instrument ontology
  if ((!file.exists(paste0(homePath, "/MS_models.csv")))||(updateOntologies)) {
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
      write.csv(vendorsInstr, paste0(homePath, "/MS_models.csv"), row.names = FALSE)
      return(0)
    })
  }
  #
  parallel::stopCluster(tmpCl)
  #
  cat("Done\n")
}
