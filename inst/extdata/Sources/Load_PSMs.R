# Load PSMs
source(parSrc)
#
utils::data(modifications, package = "PTMods")
UniMod <- modifications
rm(modifications)
#
searchOutputs %<o% list()
l_inDirs <- length(inDirs)

locDirsFl <- paste0(homePath, "/Default_locations.xlsx")
if (!file.exists(locDirsFl)) { proteoCraft::Configure() }
locDirs <- openxlsx2::read_xlsx(paste0(homePath, "/Default_locations.xlsx"))
archDir <- locDirs$Path[match("Archive folder", locDirs$Folder)]
dbDir <- locDirs$Path[match("Fasta files", locDirs$Folder)]

# Let's load us some PSMs!
# Process sequentially all input directories
psmSrcs <- setNames(paste0(libPath, "/extdata/Sources/Load_PSMs_", c("MQ", "DiaNN", "FP"), ".R"),
                     c("MQ", "DiaNN", "FP"))
for (dir_i in 1L:l_inDirs) { #dir_i <- 1L #dir_i <- 2L
  cat(paste0("\n -> Processing input folder",
             c("", paste0(" #", dir_i))[(l_inDirs > 1L)+1L],
             ":\n\t", inDirs[dir_i], "\n"))
  # Better to run this with a for loop:
  # This will at times ask the user some questions, or will assign some values to the global environment in various ways, etc...
  #
  if (SearchSoft[dir_i] == "MAXQUANT") {
    #rstudioapi::documentOpen(psmSrcs["MQ"])
    source(psmSrcs["MQ"])
  }
  if (SearchSoft[dir_i] == "DIANN") {
    #rstudioapi::documentOpen(psmSrcs["DiaNN"])
    source(psmSrcs["DiaNN"])
  }
  if (SearchSoft[dir_i] == "FRAGPIPE") {
    #rstudioapi::documentOpen(psmSrcs["FP"])
    source(psmSrcs["FP"])
  }
  if (SearchSoft[dir_i] == "PROTEOMEDISCOVERER") { stop("This part has not yet been re-written for Proteome Discoverer!") }
  if (SearchSoft[dir_i] == "ALPHADIA") { stop("This part has not yet been re-written for alphaDIA!") }
  if (SearchSoft[dir_i] == "SPECTRONAUT") { stop("This part has not yet been re-written for Spectronaut!") }
  if (SearchSoft[dir_i] == "PEAKS") { stop("This part has not yet been re-written for Peaks!") }
  if (SearchSoft[dir_i] == "SKYLINE") { stop("This part has not yet been re-written for Skyline (besides, would you really want to use Skyline for this?)") }
  cat("\n\n")
}
#
# Now combine the outputs of the different search directories
#  - Modifications table
Modifs %<o% lapply(searchOutputs, \(x) { x$Mods })
tst <- plyr::rbind.fill(Modifs)
kol1 <- c("Mark", "Mass shift")
kol2 <- c("Mark", "Shifts")
if (!sum(!kol1 %in% colnames(tst))) {
  modsTst1 <- modsTst2 <- do.call(rbind, lapply(Modifs, \(x) { x[, kol1] }))
  modsTst1 <- aggregate(round(modsTst1$`Mass shift`, 4L), list(modsTst1$Mark), unique)
  colnames(modsTst1) <- kol2
  modsTst1$N <- lengths(modsTst1$Shifts)
  w <- which(modsTst1$N > 1L)
  if (length(w)) {
    stop("We need to add code to deal with ambiguous cases where the same mark is used in different searches for different PTMs!!!")
  }
  modsTst2 <- aggregate(modsTst2$Mark, list(round(modsTst2$`Mass shift`, 4L)), unique)
  colnames(modsTst2) <- rev(kol2)
  modsTst2$N <- lengths(modsTst2$Mark)
  w <- which(modsTst2$N > 1L)
  if (length(w)) {
    stop("We need to add code deal with these ambiguous cases where different marks are assigned to the same PTM by different searches!!!")
  }
  modsTst <- do.call(plyr::rbind.fill, Modifs)
  modsTst_a <- aggregate(modsTst$"Full name", list(modsTst$Mark, round(modsTst$"Mass shift", 4L)), \(x) { x[[1L]] })
  colnames(modsTst_a) <- c(kol1, "Full name")
  k <- c("Site", "Site_long", "AA", "UniMod", "Position")
  w <- which(k %in% colnames(modsTst))
  modsTst_b <- aggregate(modsTst[, k[w]], list(modsTst$Mark, round(modsTst$"Mass shift", 4L)), \(x) { unique(unlist(x)) })
  colnames(modsTst_b) <- c(kol1, k[w])
  modsTst_a[, k[w]] <- modsTst_b[, k[w]]
  modsTst_a$Type <- aggregate(modsTst$Type, list(modsTst$Mark, round(modsTst$"Mass shift", 4L)), \(x) {
    x <- unique(x)
    x <- x[!is.na(x)]
    x <- base::setdiff(x, "NA")
    if ((length(x) > 1L) && ("Variable" %in% x)) {
      x <- "Variable"
    }
    return(x)
  })$x
  if ("Mass delta" %in% colnames(modsTst)) {
    modsTst_a$"Mass delta" <- aggregate(1L:nrow(modsTst), list(modsTst$Mark, round(modsTst$"Mass shift", 4L)), \(x) {
      rs <- unique(unlist(modsTst$`Mass delta`[x]))
      rs <- rs[!is.na(rs)]
      if (!length(rs)) {
        rs <- unique(unlist(modsTst$`Mass shift`[x]))
        rs <- rs[!is.na(rs)]
      }
      rs <- if (length(rs)) { mean(rs) } else { NA }
      return(rs)
    })$x
  }
  if ("Max occurences" %in% colnames(modsTst)) {
    modsTst_a$"Max occurences" <- aggregate(modsTst$"Max occurences", list(modsTst$Mark, round(modsTst$"Mass shift", 4L)), \(x) {
      max(x, na.rm = TRUE)
    })$x
  }
  Modifs <- modsTst_a
} else {
  Modifs <- plyr::rbind.fill(Modifs)
}

#  - MatMet text templates
MatMet_Search %<o% setNames(lapply(1L:l_inDirs, \(x) {
  paste0(c("", paste0("Search #", x, " = '", sub(".*/", "", inDirs[x]), "':\n"))[(l_inDirs > 1L)+1L],
         searchOutputs[[x]]$MatMet_txt)
}), inDirs)
# These will need to be updated later
#cat(paste(unlist(MatMet_Search), collapse = "\n"), "\n")
#  - Fastas
fastas_map %<o% list(Original = setNames(lapply(searchOutputs, \(x) { x$Fastas }),
                                         inDirs))
orig_fastas <- fastas <- unique(unlist(fastas_map))
w <- which(!file.exists(fastas))
if (length(w)) {
  dflt <- dbDir
  for (i in w) {
    cat(paste0(" - Fasta file ", fastas[i], " has been relocated since the search",
               c(" was", "es were")[(l_inDirs > 1L)+1L],
               " run and could not be located automatically.\n   Prompting user...\n"))
    tmp <- c()
    while (!length(tmp)) {
      tmp <- rstudioapi::selectFile(paste0("Select missing fasta file ", fastas[i]),
                                    path = dflt,
                                    filter = "fasta file (*.fasta|*.fas|*.fa|*.faa|*.fasta.fas|*.txt)")
    }
    fastas[i] <- tmp
  }
}
fastas_map$Actual <- setNames(lapply(names(fastas_map$Original), \(nm) {
  fastas[match(fastas_map$Original[[nm]], orig_fastas)]
}), names(fastas_map$Original))
fastas %<o% unique(fastas)
if (exists("loadInt") && loadInt && exists("intPrtFst") && file.exists(intPrtFst)) { # Load from input popup at start
  fastas <- base::union(fastas, intPrtFst)
}
if (exists("fastas_reloaded") && file.exists(fastas_reloaded)) { # Re-load from detected local backup
  fastas <- base::union(fastas, fastas_reloaded)
}

#  - FracMap
FracMap %<o% lapply(1L:length(searchOutputs), \(i) {
  res <- searchOutputs[[i]]$FracMap
  res$Search <- i
  return(res)
})
FracMap <- do.call(plyr::rbind.fill, c(FracMap))
if (length(unique(FracMap$Search)) == 1L) { FracMap$Search <- NULL } else {
  FracMap <- FracMap[, c("Search", base::setdiff(colnames(FracMap), "Search"))]
}
tst1 <- "Parent sample" %in% colnames(FracMap)
tst2 <- "MQ.Exp" %in% colnames(FracMap)
if (tst2) {
  if (!tst1) { FracMap$"Parent sample" <- FracMap$`Raw files name` }
  FracMap$MQ.Exp <- FracMap$"Parent sample"
}
if (tst1) {
  if (!tst2) { FracMap$MQ.Exp <- FracMap$`Raw files name` }
  FracMap$"Parent sample" <- FracMap$MQ.Exp
}
w <- which(is.na(FracMap$"Parent sample"))
if (length(w)) { FracMap$"Parent sample"[w] <- FracMap$MQ.Exp[w] }
w <- which(is.na(FracMap$"Parent sample"))
if (length(w)) { FracMap$"Parent sample"[w] <- FracMap$"Raw files name"[w] }
w <- which(is.na(FracMap$MQ.Exp))
if (length(w)) { FracMap$MQ.Exp[w] <- FracMap$"Parent sample"[w] }
#  - LabelType
#    For now we must enforce single values for LabelType!
LabelType %<o% unique(vapply(searchOutputs, \(x) { x$Software_param$LabelType }, ""))
if (length(LabelType) > 1L) {
  stop("It is ok to combine several searches but do it in a responsible manner!\nFor now, it is not possible to combine different searches with different labelling methods.\nWe're not sure that it even makes sense to do it so this is unlikely to change in the future!")
}
#  - MinPepSz
MinPepSz %<o% min(vapply(searchOutputs, \(x) { x$Software_param$MinPepSz }, 1L))
#  - Missed
Missed %<o% max(vapply(searchOutputs, \(x) { x$Software_param$Missed }, 1L))
#  - Raw files
rawFiles %<o% unique(unlist(lapply(searchOutputs, \(x) { x$Raw_files$Full })))
rawFiles2 %<o% unique(unlist(lapply(searchOutputs, \(x) { x$Raw_files$Names })))
#  - PSM files
PSMsFls_2_Engine %<o% setNames(lapply(1L:l_inDirs, \(x) {
  data.frame(File = searchOutputs[[x]]$PSMsFls,
             Search_ID = inDirs[x],
             Search_Engine = paste0(setNames(searchOutputs[[x]]$Software_param$Software, NULL),
                                    "_v", searchOutputs[[x]]$Software_param$Vers))
}), inDirs)
PSMsFls %<o% unlist(lapply(PSMsFls_2_Engine, \(x) { x$File }))
searchEngines %<o% setNames(unlist(lapply(PSMsFls_2_Engine, \(x) { x$Search_Engine })), PSMsFls)
# - isDIA:
isDIA %<o% setNames(vapply(searchOutputs, \(x) { x$isDIA }, TRUE),
                    inDirs)
# - QuantUMS:
QuantUMS %<o% setNames(vapply(searchOutputs, \(x) {
  x <- if ("QuantUMS" %in% names(x$Software_param)) { x$Software_param$QuantUMS } else { FALSE }
  return(x)
}, TRUE), inDirs)
#  - PSMs: 
ev %<o% do.call(plyr::rbind.fill, lapply(searchOutputs, \(x) { x$ev }));ev$id <- 1L:nrow(ev) ### THE 2ND PART IS SO IMPORTANT I AM PUTTING THEM ON 1 ROW
ev$`Raw file` <- sub("\\.[^.]+$", "", sub(".*/", "", ev$`Raw file path`)) # Can be a useful last check
#tstEvs <- do.call(plyr::rbind.fill, lapply(searchOutputs, \(x) { x$ev[1L:10L,] }));View(tstEvs)

if (exists("FracMap_reloaded")) {
  # 1st approach - match by unmodified file path
  srchTst <- sum(c("Search" %in% colnames(FracMap_reloaded),
               "Search" %in% colnames(FracMap)))
  if (srchTst %in% c(0L, 2L)) {
    FracMap_bckp <- FracMap
    mtchKols <- c("Raw file", "Raw files name")
    if (srchTst == 2L) {
      mtchKols <- c("tmp1", "tmp2")
      FracMap_reloaded[[mtchKols[1L]]] <- do.call(paste, c(FracMap_reloaded[, c("Search", "Raw file")], sep = "->"))
      FracMap_reloaded[[mtchKols[2L]]] <- do.call(paste, c(FracMap_reloaded[, c("Search", "Raw files name")], sep = "->"))
      FracMap[[mtchKols[1L]]] <- do.call(paste, c(FracMap[, c("Search", "Raw file")], sep = "->"))
      FracMap[[mtchKols[2L]]] <- do.call(paste, c(FracMap[, c("Search", "Raw files name")], sep = "->"))
    }
    m <- match(FracMap[[mtchKols[1L]]], FracMap_reloaded[[mtchKols[1L]]])
    tst <- sum(is.na(m)) # Are there any NAs?
    gs <- rep(FALSE, 3L) # Should we remove spaces - and if so then how?
    # 2nd approach - match by modified file path
    if (tst) {
      m1 <- match(gsub(" ", "-", FracMap[[mtchKols[1L]]]), gsub(" ", "-", FracMap_reloaded[[mtchKols[1L]]]))
      m2 <- match(gsub(" ", "_", FracMap[[mtchKols[1L]]]), gsub(" ", "_", FracMap_reloaded[[mtchKols[1L]]]))
      m3 <- match(gsub(" ", "", FracMap[[mtchKols[1L]]]), gsub(" ", "", FracMap_reloaded[[mtchKols[1L]]]))
      gs[1L] <- sum(is.na(m1))
      gs[2L] <- sum(is.na(m2))
      gs[3L] <- sum(is.na(m3))
      w <- which(!gs)
      tst <- length(w) == 0L # If we did not find a solution, we need to go on trying
      if (!tst) {
        m <- get(paste0("m", w[1L]))
      }
    }
    # 3rd approach - match by unmodified file name
    if (tst) {
      m <- match(FracMap[[mtchKols[2L]]], FracMap_reloaded[[mtchKols[2L]]])
      tst <- sum(is.na(m))
      gs <- rep(FALSE, 3L)
    }
    # 4th approach - match by modified file name
    if (tst) {
      m1 <- match(gsub(" ", "-", FracMap[[mtchKols[2L]]]), gsub(" ", "-", FracMap_reloaded[[mtchKols[2L]]]))
      m2 <- match(gsub(" ", "_", FracMap[[mtchKols[2L]]]), gsub(" ", "_", FracMap_reloaded[[mtchKols[2L]]]))
      m3 <- match(gsub(" ", "", FracMap[[mtchKols[2L]]]), gsub(" ", "", FracMap_reloaded[[mtchKols[2L]]]))
      gs[1L] <- sum(is.na(m1))
      gs[2L] <- sum(is.na(m2))
      gs[3L] <- sum(is.na(m3))
      w <- which(!gs)
      tst <- length(w) == 0L # If we did not find a solution, we need to go on trying
      if (!tst) {
        m <- get(paste0("m", w[1L]))
      }
    }
    if (!tst) {
      k <- c("Raw file", "Raw files name", "Parent sample", "Fraction", "Use", "PTM-enriched")
      k <- base::intersect(k, colnames(FracMap_reloaded))
      FracMap[, k] <- FracMap_reloaded[m, k]
      if (LabelType == "LFQ") {
        if ("Parent sample" %in% colnames(FracMap_reloaded)) {
          FracMap$"Parent sample" <- FracMap_reloaded$"Parent sample"[m]
        } else {
          if ("MQ.Exp" %in% colnames(FracMap_reloaded)) {
            FracMap$"Parent sample" <- FracMap_reloaded$MQ.Exp[m]
          }
        }
      }
      if (sum(gs)) {
        rp <- c("-", "_", "")[gs][1L]
        m2Ev <- match(gsub(" ", rp, FracMap$`Raw files name`), gsub(" ", rp, ev$`Raw file`))
      } else {
        m2Ev <- match(FracMap$`Raw files name`, ev$`Raw file`)
      }
      tst <- sum(is.na(m2Ev))
      if (!tst) {
        # Fix paths in case we re-mapped them
        FracMap$`Raw file` <- ev$`Raw file path`[m2Ev]
      } else {
        rm(FracMap_reloaded)
        FracMap <- FracMap_bckp
      }
      # Final final test
      rawFiles <- unique(FracMap$`Raw file`)
      rawFiles2 <- unique(FracMap$`Raw files name`)
      stopifnot(sum(!rawFiles %in% ev$`Raw file path`) == 0L,
                sum(!rawFiles2 %in% ev$`Raw file`) == 0L)
      if (srchTst == 2L) {
        if (exists("FracMap_reloaded")) {
          FracMap_reloaded$tmp1 <- NULL
          FracMap_reloaded$tmp2 <- NULL  
        }
        FracMap$tmp1 <- NULL
        FracMap$tmp2 <- NULL
      }
    } else {
      rm(FracMap_reloaded)
    }
  } else {
    rm(FracMap_reloaded)
  }
}
#if (!file.exists(FracMapPath)) {
write.csv(FracMap, file = FracMapPath, row.names = FALSE)
#}
