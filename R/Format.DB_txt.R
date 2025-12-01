#' Format.DB_txt
#'
#' @description 
#' Function to get a table with all entry details from a text version of a UniProtKB proteome.
#' (download all protein entries and choose "text" format)
#'
#' @param txt Text annotations, normally already read into R from your UniProtKB annotations txt file, but it is now also possible to provide a path to a file.
#' @param filter Valid UniProt Accessions. Filter for these to hasten the process.
#' @param GO Set to TRUE (default) to retrieve GO annotations.
#' @param Taxonomy Set to TRUE (default) to retrieve taxonomic annotations.
#' @param InterPro Set to TRUE (default) to retrieve InterPro annotations.
#' @param Pfam Set to TRUE (default) to retrieve Pfam annotations.
#' @param PIRSF Set to TRUE (default) to retrieve PIRSF annotations.
#' @param PROSITE Set to TRUE (default) to retrieve PROSITE annotations.
#' @param EMBL Set to TRUE (default) to retrieve EMBL annotations.
#' @param Ensembl Set to TRUE (default) to retrieve Ensembl annotations.
#' @param MW Set to TRUE (default) to retrieve protein molecular weights.
#' @param Sequence Set to TRUE (default) to retrieve protein sequences.
#' @param PTMs Default = FALSE; if TRUE, will return a map with all PTM sites for which there is evidence.
#' @param PDB Set to TRUE (default) to retrieve PDB IDs.
#' @param TAIR Set to TRUE (default) to retrieve TAIR IDs.
#' @param WormBase Set to TRUE (default) to retrieve WormBase IDs.
#' @param FlyBase Set to TRUE (default) to retrieve FlyBase IDs.
#' @param Features Set to TRUE (default) to retrieve Features.
#' @param Feat_accRgx Normally do not change this, the default regex from https://www.uniprot.org/help/accession_numbers is fine.\cr
#' This is only used if Features = TRUE, and has NO EFFECT on the Accession column itself!\cr
#' Only provide a different regex here if you are parsing a custom file not sourced from UniProtKB where the Accession pattern changes.
#' @param usePar Default = FALSE; if TRUE, uses the next 3 arguments and the "parallel" package for multi-threads parallel processing.
#' @param N.clust A limit on the number of vCPUs to use. If left as NULL (default), uses the number of available clusters - 1, to a minimum of 1.
#' @param N.reserved Default = 1. Number of reserved vCPUs the function is not to use. Note that for obvious reasons it will always use at least one.
#' @param cl Already have a cluster handy? Why spend time making a new one, which on top of that may invalidate the old one. Just pass it along!
#'  
#' @examples

#' data <- Format.DB_txt(txt)
#' 
#' @export

Format.DB_txt <- function(txt,
                          filter = NULL,
                          GO = TRUE,
                          Taxonomy = TRUE,
                          InterPro = TRUE,
                          Pfam = TRUE,
                          PIRSF = TRUE,
                          PROSITE = TRUE,
                          EMBL = TRUE,
                          Ensembl = TRUE,
                          MW = TRUE,
                          Sequence = TRUE,
                          PTMs = FALSE,
                          PDB = TRUE,
                          TAIR = TRUE,
                          WormBase = TRUE,
                          FlyBase = TRUE,
                          Features = FALSE,
                          usePar = TRUE,
                          N.clust,
                          N.reserved = 1,
                          cl,
                          Feat_accRgx = "[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}") {
  TESTING <- FALSE
  #
  #TESTING <- TRUE;proteoCraft::DefArg(proteoCraft::Format.DB_txt);usePar = TRUE;cl = parClust 
  #txt = annot_Fl;Features = TRUE
  #txt = x
  Feat_isoRgx <- paste0("(", Feat_accRgx, ")(-[0-9]+)?")
  #
  #txt <- readLines(AnnotFls)
  if (TESTING) {
    # Note:
    # This is not a perfect alternative to missing but will work in most cases, unless x matches a function imported by a package 
    misFun <- function(x) { return(!exists(deparse(substitute(x)))) }
  } else { misFun <- missing }
  #
  # Create cluster
  tstCl <- stopCl <- misFun(cl)
  if (!misFun(cl)) {
    tstCl <- suppressWarnings(try({
      a <- 1
      parallel::clusterExport(cl, "a", envir = environment())
    }, silent = TRUE))
    tstCl <- !"try-error" %in% class(tstCl)
  }
  if ((misFun(cl))||(!tstCl)) {
    dc <- parallel::detectCores()
    if (misFun(N.reserved)) { N.reserved <- 1 }
    if (misFun(N.clust)) {
      N.clust <- max(c(dc-N.reserved, 1))
    } else {
      if (N.clust > max(c(dc-N.reserved, 1))) {
        warning("More cores specified than allowed, I will ignore the specified number! You should always leave at least one free for other processes, see the \"N.reserved\" argument.")
        N.clust <- max(c(dc-N.reserved, 1))
      }
    }
    cl <- parallel::makeCluster(N.clust, type = "SOCK")
  }
  N.clust <- length(cl)
  #
  # Detect whether txt is actually a path rather than the database already in environment. 
  if ((length(txt) == 1)&&(file.exists(txt))) {
    cat(paste0("File path input detected, reading data at:\n", txt, "\n...\n"))
    txt <- readLines(txt)
  }
  #
  G1 <- grep("^ID ", txt)
  if (!length(G1)) { stop("Invalid file loaded!") }
  lG <- length(G1)
  #
  if (usePar) {
    # Create range in 3 steps to avoid:
    #
    #   Warning message:
    #   In length(txt) * (1:N.clust) : NAs produced by integer overflow
    #
    # This happens for very large files and actually kills the function!!!
    tmp <- round((1:N.clust)/N.clust, 4)
    RG <- length(txt)*tmp
    RG <- round(RG)
    #
    if (N.clust > 1) {
      RG[1:(N.clust-1)] <- vapply(RG[1:(N.clust-1)], function(x) {
        max(G1[which(G1 < x)+1])-1
      }, 1)
    }
    batChes <- setNames(lapply(1:length(RG), function(x) {
      txt[(c(0, RG)[x]+1):(RG[x])]
    }), paste0("Batch", 1:length(RG)))
  } else { batChes <- list(Batch1 = txt) }
  nBatches <- length(batChes)
  currWD <- getwd()
  if (usePar) {
    invisible(lapply(1:nBatches, function(i) {
      readr::write_rds(batChes[[i]], paste0(currWD, "/tmp", i, ".RDS"))
    }))
  }
  F0 <- function(i) { #i <- 1
    if (usePar) {
      btch <- readr::read_rds(paste0(currWD, "/tmp", i, ".RDS"))
    } else {
      btch <- batChes[[i]]
    }
    # UniProtKB IDs
    g1 <- grep("^ID ", btch)
    tbl <- data.frame(ID = btch[g1])
    #should be faster than the previous version (commented line below, in case it causes issues):
    tbl$ID <- gsub(" .+$", "",gsub("^[^ ]+ ", "", gsub(" +", " ", tbl$ID)))
    #tbl$ID <- vapply(strsplit(gsub(" +", " ", tbl$ID), " "), function(x) { unlist(x)[2] }, "")
    # UniProtKB Accessions
    g1 <- c(g1, length(btch)+1)
    tbl$Text <- lapply(1:nrow(tbl), function(x) { btch[g1[x]:(g1[x+1]-1)] })
    tbl$Text2 <- paste0("<->", vapply(tbl$Text, paste, "", collapse = "\n<->"), "\n")
    if (!usePar) {
      cat("-> Parsing txt file...\n")
      cat("   Accessions...\n")
    }
    tbl$Accession <- gsub("<->(?!AC)[^\n]*\n", "", tbl$Text2, perl = TRUE)
    tbl$Accession <- gsub("<->AC +|;\n$", "", gsub(";\n<->AC +", ";", tbl$Accession))
    tbl$Accession <- gsub("; +", ";", tbl$Accession)
    # The whole function takes long, sometimes we can take a shortcut:
    if ((!is.null(filter))&&(!"function" %in% class(filter))) {
      w <- which(vapply(strsplit(tbl$Accession, ";"), function(x) { sum(x %in% filter) > 0 }, TRUE))
      tbl <- tbl[w,]
      btch <- unlist(tbl$Text)
      g1 <- c(grep("^ID ", btch), length(btch)+1)
    }
    # Names
    if (!usePar) { cat("   Names...\n") }
    tbl$Names <- gsub("<->(?!DE +[^;]+Name: Full=)[^\n]*\n", "", tbl$Text2, perl = TRUE)
    tbl$Names <- gsub(";\n", ";", gsub("<->DE +[^;]+Name: Full=| \\{[^\\}]+\\}|;\n$", "", tbl$Names))
    # GO terms
    if (GO) {
      if (!usePar) { cat("   Gene Ontology annotations...\n") }
      tbl$GO <- gsub("<->(?!DR +GO)[^\n]*\n", "", tbl$Text2, perl = TRUE)
      tbl$GO <- gsub(";[^;]+\n", "\n", gsub("<->DR +GO; *", ";", tbl$GO, perl = TRUE))
      tbl$Ontology <- gsub(";$", "", gsub(":[^\n]+\n", ";", gsub(";[^;]+; +", "", tbl$GO)))
      tbl$"GO-ID" <- gsub("^;", "", gsub(";[^;]+\n", "", tbl$GO))
      tbl$GO <- gsub("^;|\n$", "", gsub("\n;GO:[0-9]+; +[CFP]:", ";", paste0("\n", tbl$GO)))
      w <- which(nchar(tbl$GO) > 0)
      tmp <- tbl[w, c("GO", "GO-ID")]
      tmp$GO <- strsplit(tmp$GO, ";")
      tmp$"GO-ID" <- strsplit(tmp$"GO-ID", ";")
      tbl$GO[w] <- apply(tmp, 1, function(x) {
        x1 <- unlist(x[[1]])
        x2 <- unlist(x[[2]])
        l <- length(x1)
        paste(vapply(1:l, function(y) { paste0(x1[y], " [", x2[y], "]") }, ""), collapse = ";")
      })
      tst1 <- unlist(strsplit(tbl$`GO-ID`, ";"))
      tst2 <- unlist(strsplit(tbl$GO, ";"))
      tst3 <- unlist(strsplit(tbl$Ontology, ";"))
      if (sum(c(length(tst1), length(tst2), length(tst3)))) {
        tst4 <- try(aggregate(cbind(tst2, tst3), list(tst1), unique), silent = TRUE)
        if (("try-error" %in% class(tst4))||(length(tst4[[2]]) != length(unique(tst4[[2]])))||(!"character" %in% class(tst4[[2]]))||(!"character" %in% class(tst4[[3]]))) {
          msg <- "Gene Ontology data was not properly formatted in the input txt file!"
          tbl$GO <- gsub("^;|;$", "", gsub(";+", "", tbl$GO))
          tbl$`GO-ID` <- gsub("^;|;$", "", gsub(";+", "", tbl$`GO-ID`))
          tbl$Ontology <- gsub("^;|;$", "", gsub(";+", "", tbl$Ontology))
          tst1 <- unlist(strsplit(tbl$`GO-ID`, ";"))
          tst2 <- unlist(strsplit(tbl$GO, ";"))
          tst3 <- unlist(strsplit(tbl$Ontology, ";"))
          tst4 <- try(aggregate(cbind(tst2, tst3), list(tst1), unique), silent = TRUE)
          if (("try-error" %in% class(tst4))||(length(tst4[[2]]) != length(unique(tst4[[2]])))||(!"character" %in% class(tst4[[2]]))||(!"character" %in% class(tst4[[3]]))) {
            #stop(i)
            stop(msg)
          } else { warning(gsub("\\!$", ", but we managed to make sense of it.", msg)) }
        }
      }
    }
    # Taxonomic information:
    l <- length(g1)-1
    if (Taxonomy) {
      if (!usePar) { cat("   Taxonomic annotations...\n") }
      # Organism
      tbl$Organism <- gsub("<->(?!OS +)[^\n]*\n", "", tbl$Text2, perl = TRUE)
      tbl$Organism <- gsub("^<->OS +|\\.?\n$", "", gsub("\\.?\n<->OS +", ";", tbl$Organism))
      # Taxonomy
      tbl$Taxonomy <- gsub("<->(?!OC +)[^\n]*\n", "", tbl$Text2, perl = TRUE)
      tbl$Taxonomy <- gsub("^<->OC +|\\.?\n$", "", gsub("\n<->OC +", " ", tbl$Taxonomy))
      # TaxID
      tbl$TaxID <- gsub("<->(?!OX +NCBI_TaxID=)[^\n]*\n", "", tbl$Text2, perl = TRUE)
      tbl$TaxID <- gsub("^<->OX +NCBI_TaxID=| *(\\{[^\\}]+\\})?;?\n$", "", gsub("\n<->OX +NCBI_TaxID=", ";", tbl$TaxID))
    }
    # Other annotations:
    for (i in c("InterPro", "Pfam", "PIRSF", "PROSITE")) { #i <- "InterPro" #i <- "PROSITE"
      if (get(i)) {
        if (!usePar) { cat(paste0("   ", i, " annotations...\n")) }
        pat1 <- paste0("<->(?!DR +", i, ")[^\n]*\n")
        pat2 <- paste0("; +[^;\n]+\n<->DR +", i, "; +")
        pat3 <- paste0("^<->DR +", i, "; +|; +[^\n]+\n$")
        tbl[[i]] <- gsub(pat1, "", tbl$Text2, perl = TRUE)
        tbl[[i]] <- gsub(pat3, "", gsub(pat2, ";", tbl[[i]]))
      }
    }
    #View(tbl[, c("InterPro", "Pfam", "PIRSF", "PROSITE")])
    # EMBL
    if (EMBL) {
      if (!usePar) { cat("   EMBL accessions...\n") }
      i <- "EMBL"
      pat1 <- paste0("<->(?!DR +", i, ")[^\n]*\n")
      pat2 <- paste0("; +[^;]+; +[^\\.;]+\\.\n<->DR +", i, "; +")
      pat3 <- paste0("^<->DR +", i, "; +|; +[^;]+; +[^\\.;]+\\.\n$")
      tbl[[i]] <- gsub(pat1, "", tbl$Text2, perl = TRUE)
      tbl[[i]] <- gsub("; +", ";", gsub(pat3, "", gsub(pat2, ";", tbl[[i]])))
    }
    # Ensembl
    if (Ensembl) {
      if (!usePar) { cat("   Ensembl annotations...\n") }
      i <- "Ensembl[A-Z]?[a-z]*"
      pat1 <- paste0("<->(?!DR +", i, ")[^\n]*\n")
      pat2 <- paste0("\n<->DR +", i, "; +")
      pat3 <- paste0("^<->DR +", i, "; +|\n$")
      EnsemblClasses <- c("T", "P", "G")
      enscol <- paste0("Ensembl_", EnsemblClasses)
      tmp <- gsub(pat1, "", tbl$Text2, perl = TRUE)
      tmp <- gsub(pat3, "", gsub(pat2, "<<<>>>", tmp))
      tbl[[enscol[1]]] <- gsub("; +[^;]+; +[^;]+$", "", gsub("; +[^;]+; +[^;<]+<<<>>>", ";", tmp))
      tbl[[enscol[2]]] <- gsub("^[^;]+; +|; +[^;]+$", "", gsub("; +[^;<]+<<<>>>[^;]+; +", ";", tmp))
      tbl[[enscol[3]]] <- gsub("^[^;]+; +[^;]+; +", "", gsub("<<<>>>[^;]+; +[^;]+; +", ";", tmp))
    }
    # MW
    if (MW) {
      if (!usePar) { cat("   molecular weights...\n") }
      tbl$"MW (Da)" <- gsub("<->(?!SQ +SEQUENCE +)[^\n]*\n", "", tbl$Text2, perl = TRUE)
      tbl$"MW (Da)" <- as.integer(gsub("^<->SQ +SEQUENCE +[0-9]+ +AA; +| +MW; +.*", "", tbl$"MW (Da)"))
    }
    # Sequence
    if (Sequence) {
      if (!usePar) { cat("   sequences...\n") }
      # Here we are cheating a bit because we know what to expect
      tbl$Sequence <- gsub(".+<->SQ +SEQUENCE +[^\n]*\n<-> +", "", tbl$Text2)
      tbl$Sequence <- gsub(" |\n<->", "", gsub("\n<->//.*", "", tbl$Sequence))
    }
    # Code after this has not been redone the "new way": it's not that slow
    # Sequence with known PTMs
    if (PTMs) {
      if (!usePar) { cat("   known PTM sites...\n") }
      tbl$"Known PTMs" <- ""
      tbl$tmp <- sapply(tbl$Text, function(x) { grep("^FT +MOD_RES +[0-9]+", x) })
      tmp <- unlist(apply(tbl[, c("tmp", "Text")], 1, function(x) {
        unlist(x[2])[as.numeric(unlist(x[1]))]
      }))
      tst <- vapply(c("^FT +MOD_RES +[0-9]+ +[0-9]+ +.+$", "^FT +MOD_RES +[0-9]+$"), function(x) {
        sum(grepl(x, tmp))
      }, 1) == length(tmp)
      stopifnot(sum(tst) == 1)
      w <- which(vapply(tbl$tmp, length, 1) > 0)
      tbl$"Known PTMs"[w] <- apply(tbl[w, c("Text", "tmp", "Sequence")], 1, function(x) {
        #x <- tbl[w[1], c("Text", "tmp", "Sequence")]
        seq <- paste0("_", unlist(x[[3]]), "_")
        if (which(tst) == 1) {
          tmp <- strsplit(gsub("^FT +MOD_RES +|[\\.,;].*$", "", unlist(x[1])[unlist(x[2])]), " +")
          ptms <- as.data.frame(t(sapply(tmp, function(y) {
            y <- unlist(y)
            stopifnot(y[1] == y[2])
            return(c(y[1], paste(y[3:length(y)], collapse = " ")))
          })))
          colnames(ptms) <- c("Residue", "PTM")
        }
        if (which(tst) == 2) {
          ptms <- data.frame(Residue = gsub("^FT +MOD_RES +", "", unlist(x[1])[unlist(x[2])]),
                             PTM = gsub("^FT +/note=\"|\"$", "", unlist(x[1])[unlist(x[2])+1]))
        }
        ptms <- aggregate(ptms$PTM, list(ptms$Residue), function(y) {
          paste0("([(", paste(y, collapse = "_/_"), ")])")
        })
        ptms <- rbind(c(0, ""), ptms, c(nchar(seq)+1, ""))
        ptms$Group.1 <- as.integer(ptms$Group.1)
        ptms <- ptms[order(ptms$Group.1, decreasing = FALSE),]
        tmp <- vapply(1:(nrow(ptms)-1), function(y) {
          substr(seq, ptms$Group.1[y]+2, ptms$Group.1[y+1])
        }, "")
        tmp2 <- rep("", length(tmp) + nrow(ptms)-1)
        tmp2[2*(1:length(tmp))-1] <- tmp  
        tmp2[2*(1:(nrow(ptms)-2))] <- ptms$x[2:(nrow(ptms)-1)]
        return(paste(tmp2, collapse = ""))
      })
      tbl$tmp <- NULL
      #tst <- unique(unlist(strsplit(tbl$"Known PTMs"[w], "^_?[A-Z]*\\(\\[\\(|\\)\\]\\)[A-Z]*\\(\\[\\(|\\)\\]\\)[A-Z]*_?$|_/_")))
    }
    # PDB and AlphaFold
    if (PDB) {
      pat <- "^DR +PDB; "
      g2a <- grepl(pat, btch)
      if (sum(g2a)) {
        if (!usePar) { cat("   PDB annotations...\n") }
        tbl$PDB <- vapply(tbl$Text, function(x) { #x <- tbl$Text[[1]]
          g <- grep(pat, x, value = TRUE)
          x <- ""
          if (length(g)) { x <- paste(unique(gsub(";.+", "", gsub(pat, "", g))), collapse = ";") }
          return(x)
        }, "")
      }
    }
    if (TAIR) {
      pat <- "^DR +TAIR; "
      g2a <- grepl(pat, btch)
      if (sum(g2a)) {
        if (!usePar) { cat("   TAIR annotations...\n") }
        tbl$TAIR <- vapply(tbl$Text, function(x) { #x <- tbl$Text[[1]]
          g <- grep(pat, x, value = TRUE)
          x <- ""
          if (length(g)) { x <- paste(unique(gsub("locus:[0-9]+; |\\.$", "", gsub(pat, "", g))), collapse = ";") }
          return(x)
        }, "")
      }
    }
    if (WormBase) {
      pat <- "^DR +WormBase; "
      g2a <- grepl(pat, btch)
      if (sum(g2a)) {
        if (!usePar) { cat("   WormBase annotations...\n") }
        tbl[, paste0("WormBase_",
                     c("Primary name", "Sequence", "Gene"))] <- proteoCraft::Isapply(tbl$Text, function(x) { #x <- tbl$Text[[1]]
                       g <- grep(pat, x, value = TRUE)
                       x <- rep("", 3)
                       if (length(g)) {
                         x <- c(paste(gsub(".+; +|\\.$", "", g), collapse = ";"),
                                paste(gsub(";.+", "", gsub(pat, "", g)), collapse = ";"),
                                paste(gsub(".+; +", "", gsub("; [^;]+$", "", g)), collapse = ";"))
                       }
                       return(x)
                     })
      }
    }
    if (FlyBase) {
      pat <- "^DR +FlyBase; "
      g2a <- grepl(pat, btch)
      if (sum(g2a)) {
        if (!usePar) { cat("   FlyBase annotations...\n") }
        tbl[, paste0("FlyBase_",
                     c("Symbol", " FlyBase ID"))] <- proteoCraft::Isapply(tbl$Text, function(x) { #x <- tbl$Text[[1]]
                       g <- grep(pat, x, value = TRUE)
                       x <- rep("", 2)
                       if (length(g)) {
                         x <- c(paste(gsub(".+; +|\\.$", "", g), collapse = ";"),
                                paste(gsub(";.+", "", gsub(pat, "", g)), collapse = ";"))
                       }
                       return(x)
                     })
      }
    }
    # Features
    if (Features) {
      if (!usePar) { cat("   Features...\n") }
      tbl$Features <- gsub("^<->FT +", "", gsub("<->(?!FT +)[^\n]*\n", "", tbl$Text2, perl = TRUE))
      w <- which(tbl$Features != "")
      tbl$Features <- strsplit(tbl$Features, "(\n<->)?FT +")
      tbl$Features[w] <- lapply(w, function(i) { #i <- w[1] #i <- 2 #i <- 116 #i <- 243 #i <- 410 #i <- 475 #i <- 195
        x <- unlist(tbl$Features[i])
        lX <- length(x)
        # The assumption I originally made about the structure was that each item covers 3 rows.
        # Unfortunately, this isn't true...
        pat <- paste0("^[A-Z_]+ +(", Feat_isoRgx, ":)?[\\?<>]?[0-9]*(\\.\\.[\\?<>]?[0-9]*)?$")
        #grep(pat, x, value = TRUE)
        Nms <- grep(pat, x)
        lNms <- length(Nms)
        if (lNms == 1) {
          rgs <- data.frame(first = Nms, last = lX) 
        } else {
          rgs <- data.frame(first = Nms, last = c(Nms[2:lNms]-1, lX))         
        }
        rgs$Name <- gsub(" +.*", "", x[Nms])
        rgs$Start <- suppressWarnings(as.integer(gsub(paste0("^[A-Z_]+ +(", Feat_isoRgx, ":)?[\\?<>]?|\\.\\.[\\?<>]?[0-9]*$"), "", x[Nms])))
        rgs$End <- NA
        g <- grep("\\.\\.\\?$|\\.\\.[\\?<>]?[0-9]+$", x[Nms], invert = TRUE)
        rgs$End[g] <- rgs$Start[g]
        g <- grep("\\.\\.[\\?<>]?[0-9]+$", x[Nms])
        rgs$End[g] <- suppressWarnings(as.integer(gsub(".*\\.\\.[\\?<>]?", "", x[Nms[g]])))
        for (kol in c("Note", "Evidence", "Id")) {
          k <- tolower(kol)
          rgs[[kol]] <- lapply(1:nrow(rgs), function(j) { #j <- 1 #j <- 67 #j <- 2
            y <- rgs[j, c("first", "last")]
            y <- x[(y[[1]]+1):y[[2]]]
            gN <- grep(paste0("/", k, "=\""), y)
            rs <- NA
            lN <- length(gN)
            if (lN) {
              if (lN != 1) { stop(paste(c(i, tbl$Accession[i], kol, j), collapse = " ")) }
              gS <- grep("/", y)
              gS <- gS[which(gS > gN)]
              if (length(gS)) { y <- y[gN:(gS[1]-1)] } else {
                y <- y[gN:length(y)]
              }
              y <- paste(y, collapse = "\n")
              rs <- strsplit(gsub(paste0("^/", k, "=\"|\"$"), "", y), "\n")
            }
            return(rs)
          })
        }
        if (!TESTING) {
          rgs$first <- NULL
          rgs$last <- NULL  
        }
        return(rgs)
      })
    }
    #for (kol in colnames(tbl)) { g <- grep("\n|>>>|<->", tbl[[kol]]) ; if (length(g)) { stop(kol) } } # Check for residuals of our parsing tricks
    rm(list = setdiff(ls(), "tbl"))
    return(tbl)
  }
  if (usePar) {
    environment(F0) <- .GlobalEnv
    exports <- list("GO", "Taxonomy", "InterPro", "Pfam", "PIRSF", "PROSITE", "EMBL", "Ensembl", "MW", "Sequence", "PTMs",
                    "PDB", "TAIR", "WormBase", "FlyBase", "Features", "Feat_isoRgx", "usePar", "filter", "TESTING", "currWD")
    parallel::clusterExport(cl, exports, envir = environment())
    rsTbl <- parallel::parLapply(cl, 1:nBatches, F0)
  } else { rsTbl <- lapply(1:nBatches, F0) }
  if (usePar) {
    invisible(lapply(1:nBatches, function(i) {
      unlink(paste0(currWD, "/tmp", i, ".RDS"))
    }))
  }
  rsTbl <- plyr::rbind.fill(rsTbl)
  if (MW) { rsTbl$"MW (Da)" <- as.numeric(gsub(";.*", "", rsTbl$"MW (Da)")) }
  #Cleanup
  rsTbl$Text <- NULL
  rsTbl$Text2 <- NULL
  # "Melt" table
  tst <- nchar(gsub("[^;]", "", rsTbl$Accession))
  m <- max(tst)
  if (m) {
    kl <- colnames(rsTbl)
    kl <- kl[which(kl != "Accession")]
    w <- which(tst > 0)
    rsTbl$ROW <- 1:nrow(rsTbl)
    tmp <- proteoCraft::listMelt(strsplit(rsTbl$Accession[w], ";"), w, c("Accession", "ROW"))
    tmp[, kl] <- rsTbl[tmp$ROW, kl]
    rsTbl <- rbind(rsTbl[which(tst == 0),],
                   tmp)
    rsTbl <- rsTbl[order(rsTbl$ROW, decreasing = FALSE),]
    rsTbl$ROW <- NULL
  }
  rownames(rsTbl) <- 1:nrow(rsTbl)
  #rsTbl <- rsTbl[match(tstA$Accession, rsTbl$Accession),]
  #
  if (stopCl) { parallel::stopCluster(cl) }
  return(rsTbl)
}
