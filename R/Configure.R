#' Configure
#'
#' @description
#' Configuration function which will make sure some files from the package are moved to its subfolder in HOME upon first installation.
#' 
#' @examples
#' proteoCraft::Configure()
#' 
#' @export

Configure <- function() {
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
               "No replicates analysis - detailed script")
  extDr1 <- paste0(proteoPath, "/extdata/R scripts")
  fls <- paste0(extDr1, "/", scrpts2, ".R")
  for (fl in fls) {
    tst <- try(file.copy(fl, homePath, overwrite = TRUE), silent = TRUE)
    print(tst)
  }
  cat("Updated analysis scripts in HOME...")
  #
  # Download MS instrument ontology
  #   This can fail (stack too deep) when some packages are loaded,
  #   probably because of conflict between internally used functions sharing a same name and not called with package::...
  #   Solution I should try: run this on a cluster without loading those packages?
  #
  if (!require(rols)) { pak::pkg_install("rols") }
  ol <- rols::Ontology("ms")
  instrMod <- rols::Term(ol, "MS:1000031")
  vendors <- rols::children(instrMod)  # get all levels
  vendorsLab <- lapply(vendors@x, function(x) {
    x@label
  })
  vendorsInstr <- lapply(1:length(vendorsLab), function(x) { #x <- 2 #x <- 5
    trm <- rols::Term(ol, names(vendorsLab)[[x]])
    models <- rols::children(trm)
    if (!is.null(models)) {
      rs <- vapply(models@x, function(y) { y@label }, "")
      rs2 <- lapply(1:length(rs), function(y) { #y <- g[1]
        try({
          trm2 <- rols::Term(ol, names(rs)[y])
          mods <- rols::children(trm2)
          mods <- vapply(mods@x, function(z) { z@label }, "")
        }, silent = TRUE)
      })
      w <- which(vapply(rs2, function(x) { "try-error" %in% class(x) }, TRUE))
      if (length(w)) { rs2[w] <- rs[w] }
      rs <- unlist(rs2)
      rs <- rs[grep("^MS:[0-9]+$", names(rs))]
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
  write.csv(vendorsInstr, paste0(homePath, "/MS_models.csv"))
  #
  cat("Done\n")
}
