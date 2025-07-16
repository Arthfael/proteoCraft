#' GO_map
#'
#' @description 
#' A function to map GO terms to proteins and genes.
#' Its main purpose is to make GO_enrich faster.
#' 
#' The functions's code actually contains 2 alternative method to get the same thing done.
#' The 1st, the one actually used by default, is topGO-based. It is slightly faster and used for consistency -
#' since this function's output is later fed to other topGO functions.
#' The 2nd is more DIY (hence the name) and seems to return slightly more mappings (those of which I checked did make sense.)
#' 
#' @param DB The formated protein database, with protein IDs and GO terms.
#' @param db_ID_col Name of the column of Protein IDs (one per row) in the database. Default = "Protein ID"
#' @param db_Gene_col Name of the column of Gene IDs/names in the database. Default = "Gene".
#' @param GO.terms Optional. GO.terms data frame as created by the GO_enrich function.
#' @param method One of "topGO" (default) or "DIY". See description.
#' @param cl Already have a cluster handy? Why spend time making a new one, which on top of that may invalidate the old one. Just pass it along!
#' @param N.clust A limit on the number of vCPUs to use. If left as NULL (default), uses the number of available clusters - 1, to a minimum of 1.
#' @param N.reserved Default = 1. Number of reserved vCPUs the function is not to use. Note that for obvious reasons it will always use at least one.
#' 
#' @return A named list of the following objects
#' \itemize{
#'   \item Mappings A list of up to 2 data frames, "Protein" and optionally "Gene" (if db_Gene_col is a valid column name for DB). Each data frame maps the single GO terms per row in column "GO" to any number of semicolon-separated proteins (resp. genes) in column "Protein" (resp. "Gene").
#'   \item GO.terms - the optional GO.terms data frame if provided as input, with additional "Proteins" and (if db_Gene_col is a valid column name for DB) "Genes" columns.
#' }
#' 
#' @export

GO_map <- function(DB,
                   db_ID_col = "Protein ID",
                   db_Gene_col = "Gene",
                   GO.terms,
                   method = "topGO",
                   cl,
                   N.clust,
                   N.reserved = 1) {
  TESTING <- FALSE
  #proteoCraft::DefArg(proteoCraft::GO_map);DB <- db; TESTING <- TRUE
  if (TESTING) {
    tm1 <<- Sys.time()
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
  method <- toupper(method)
  Ont <- c("BP", "CC", "MF")
  GenTst <- (!is.null(db_Gene_col))&&(db_Gene_col %in% colnames(DB))
  # Note about the code below:
  # For consistency with GO_enrich, we now use topGO's function (which is also slightly faster)
  # However, the code for method "DIY" finds slightly more mappings, the latter which seemed correct when checked.
  GO.mappings <- list()
  if (method == "TOPGO") {
    #a0 <- Sys.time()
    temp <- proteoCraft::listMelt(strsplit(DB$"GO-ID", ";"), DB[[db_ID_col]])
    temp <- temp[which(!is.na(temp$L1)),]
    temp <- proteoCraft::listMelt(strsplit(temp$L1, ";"), temp$value)
    temp <- temp[which(!is.na(temp$L1)),]
    temp <- temp[which(!is.na(temp$value)),]
    temp <- temp[which(temp$L1 != ""),]
    temp <- temp[which(temp$value != ""),]
    Pr <- unique(temp$value)
    Pr2GO <- aggregate(temp$L1, list(temp$value), c)
    Pr2GO <- setNames(Pr2GO$x, Pr2GO$Group.1)
    Pr2GO <- Pr2GO[which(sapply(Pr2GO, length) > 0)]
    Prs2GO <- list(BP = topGO::annFUN.gene2GO("BP", names(Pr2GO), Pr2GO),
                   CC = topGO::annFUN.gene2GO("CC", names(Pr2GO), Pr2GO),
                   MF = topGO::annFUN.gene2GO("MF", names(Pr2GO), Pr2GO))
    # (Alternate way to do it, but results will be identical)
    #GO2Pr <- aggregate(temp$value, list(temp$L1), c)
    #GO2Pr <- setNames(GO2Pr$x, GO2Pr$Group.1)
    #GO2Pr <- GO2Pr[which(sapply(GO2Pr, length) > 0)]
    #GO2Prs <- list(BP = topGO::annFUN.GO2genes("BP", names(Pr2GO), GO2Pr),
    #               CC = topGO::annFUN.GO2genes("CC", names(Pr2GO), GO2Pr),
    #               MF = topGO::annFUN.GO2genes("MF", names(Pr2GO), GO2Pr))
    nms <- unique(unlist(sapply(Ont, function(x) {
      #unique(c(
      names(Prs2GO[[x]])#, names(GO2Prs[[x]])))
    })))
    Prs2GO2 <- setNames(lapply(nms, function(x) {
      sort(unique(c(Prs2GO$BP[[x]], Prs2GO$CC[[x]], Prs2GO$MF[[x]])))
    }), nms) # I am faster if not parallelized!
    #GO2Prs2 <- setNames(lapply(nms, function(x) { sort(unique(c(GO2Prs$BP[[x]], GO2Prs$CC[[x]], GO2Prs$MF[[x]]))) }), nms)
    #tst1 <- sapply(Prs2GO2, paste, collapse = ";")
    #tst2 <- sapply(GO2Prs2, paste, collapse = ";")
    #stopifnot(sum(tst1 != tst2) == 0) # Check that they are identical
    # Check proteins which are not represented: they should all be from entries in the database which are devoid of GO annotation
    misses <- Pr[which(!Pr %in% unlist(Prs2GO2))]
    if (length(misses)) {
      tst <- !DB$"GO-ID"[proteoCraft::grsep(misses, x = DB$`Protein ID`)] %in% c("", NA)
      if (sum(tst)) { warning("Investigate: a few proteins were missed!") }
    }
    # Now check GO terms - some are inexplicably missing: deal with them!
    tst1 <- names(Prs2GO2) # All GO terms mapped
    tst2 <- unique(temp$L1) # All available GO terms
    tst2N <- tst2[which(!tst2 %in% tst1)]
    w <- which(temp$L1 %in% tst2N)
    if (length(w)) {
      temp2 <- aggregate(temp$value[w], list(temp$L1[w]), function(x) { list(unique(x)) })
      temp2 <- setNames(temp2$x, temp2$Group.1)
      Prs2GO2[names(temp2)] <- temp2  
    }
    Prs2GO2 <- Prs2GO2[order(names(Prs2GO2), decreasing = FALSE)]
    temp2 <- data.frame(GO = names(Prs2GO2),
                        Protein = sapply(Prs2GO2, paste, collapse = ";"))
    rownames(temp2) <- NULL
    GO.mappings$Protein <- temp2
    if (GenTst) {
      #ta1 <- Sys.time()
      #temp <- melt(Prs2GO2)
      #temp$Gene <- DB[match(temp$value, DB[[db_ID_col]]), db_Gene_col]
      #temp <- temp[which(!is.na(temp$Gene)),]
      #temp <- temp[which(temp$Gene != ""),]
      #temp <- proteoCraft::listMelt(strsplit(temp$Gene, ";"), temp$L1)
      #temp <- aggregate(temp$value, list(temp$L1), unique)
      #Gns2GO2 <- setNames(temp$x, temp$Group.1)
      #Gns2GO2 <- Gns2GO2[order(names(Gns2GO2), decreasing = FALSE)]
      #ta2 <- Sys.time()
      #print(ta2-ta1)
      #tb1 <- Sys.time()
      #cl <- parClust
      tmp <- DB[, c(db_ID_col, db_Gene_col)]
      parallel::clusterExport(cl, list("tmp", "db_ID_col", "db_Gene_col"), envir = environment())
      f0 <- function(x) { unique(tmp[match(x, tmp[[db_ID_col]]), db_Gene_col]) }
      environment(f0) <- .GlobalEnv
      Gns2GO2a <- parallel::parLapply(cl, Prs2GO2, f0)
      names(Gns2GO2a) <- names(Prs2GO2)
      Gns2GO2a <- Gns2GO2a[order(names(Gns2GO2a), decreasing = FALSE)]
      temp2 <- data.frame(GO = names(Gns2GO2a),
                          Gene = sapply(Gns2GO2a, paste, collapse = ";"))
      rownames(temp2) <- NULL
      GO.mappings$Gene <- temp2
      #tb2 <- Sys.time()
      #print(tb2-tb1)
    }
    #b0 <- Sys.time()
    #b0-a0
  }
  if (method == "DIY") {
    #a1 <- Sys.time()
    for (kl in 1:(1+GenTst)) {
      # 1. Proteins and Genes directly annotated with term in DB
      kol <- get(c("db_ID_col", "db_Gene_col")[kl])
      nm <- c("Protein", "Gene")[kl]
      w <- which(nchar(DB[[kol]]) > 0)
      temp <- magrittr::set_colnames(proteoCraft::listMelt(strsplit(DB[w, kol], ";"), DB$"GO-ID"[w]), c(nm, "GO"))
      temp <- temp[which(!is.na(temp$GO)),]
      temp <- temp[which(nchar(temp$GO) > 0),]
      temp <- proteoCraft::listMelt(strsplit(temp$GO, ";"), temp[[nm]])
      temp <- magrittr::set_colnames(aggregate(temp$L1, list(temp$value), function(x) {
        paste(sort(unique(unlist(x))), collapse = ";")
      }), c("GO", nm))
      #2. Optional: also proteins annotated with an offspring term in DB
      if (OffspringCounts) {
        w <- which(sapply(GO.terms$Offspring, length) > 0)
        temp2 <- proteoCraft::listMelt(GO.terms$Offspring[w], GO.terms$ID[w])
        # Include not just offspring but also term itself:
        temp2 <- rbind(temp2,
                       data.frame(value = GO.terms$ID, L1 = GO.terms$ID))
        temp2[[nm]] <- temp[match(temp2$value, temp$GO), nm]
        temp2 <- temp2[which(!is.na(temp2[[nm]])),]
        temp2[[nm]] <- strsplit(temp2[[nm]], ";")
        temp2 <- magrittr::set_colnames(proteoCraft::listMelt(temp2[[nm]], temp2$L1), c(nm, "GO"))
        temp2 <- magrittr::set_colnames(aggregate(temp2[[nm]], list(temp2$GO), function(x) {
          paste(sort(unique(unlist(x))), collapse = ";")
        }), c("GO", nm))
        temp <- temp2; rm(temp2)
      }
      GO.mappings[[nm]] <- temp
    }
    #b1 <- Sys.time()
    #print(b1-a1)
  }
  Res <- list(Mappings = GO.mappings)
  if (!misFun(GO.terms)) {
    GO.terms$Proteins <- GO.mappings$Protein$Protein[match(GO.terms$ID, GO.mappings$Protein$GO)]
    if (GenTst) { GO.terms$Genes <- GO.mappings$Gene$Gene[match(GO.terms$ID, GO.mappings$Gene$GO)] }
    Res$GO.terms <- GO.terms
  }
  if (TESTING) { tm2 <<- Sys.time() }
  #
  if (stopCl) { parallel::stopCluster(cl) }
  return(Res)
}
