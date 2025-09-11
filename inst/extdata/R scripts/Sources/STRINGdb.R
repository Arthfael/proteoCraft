#### Create STRINGdb graph(s) and Cytoscape networks for regulated proteins
#
# Extend this to modified peptide filters!!! 
#
#
# NB: gets Ratios from the GO enrichment part of the script currently
# If you want to run this, you will thus first need to run the GO enrichment section first
setwd(wd)
# Packages for STRINGdb graphs
packs <- c("rbioapi", "png")
for (pack in packs) { require(pack, character.only = TRUE) }
# Create STRINGdb directory
GraphTypes %<o% c("Functional", "Physical")
dirs <- paste0(wd, "/STRINGdb/", GraphTypes)
for (dir in dirs) { if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) } }
dirlist <- unique(c(dirlist, dir))
# Process filters for STRINGdb
intNets <- list()
Tsts <- c("t-tests", "F-tests", "Localisation", "SAINTexpress")
WhTsts <- which(Tsts %in% names(Reg_filters))
msg <- "Starting STRINGdb network analysis"
ReportCalls <- AddMsg2Report(Space = FALSE)
tmpPG <- listMelt(strsplit(PG$"Leading protein IDs", ";"), PG$id)
tmpPG$L1 <- as.integer(tmpPG$L1)
IDs <- unique(unlist(lapply(names(Reg_filters), function(x) {
  i <- (x == "SAINTexpress")+1
  lapply(Reg_filters[[x]]$`By condition`, function(y) {
    w <- y$Filter
    dat <- get(c("PG", "allSAINTs")[i])
    w <- w[which(dat$"Potential contaminant"[w] != "+")]
    y <- dat[[c("Leading protein IDs", "Protein")[i]]][w]
    if (i == 1) { y <- unlist(strsplit(y, ";")) }
    return(y)
  })
})))
w <- which(db$`Protein ID` %in% IDs)
allTaxIDs <- unique(db$TaxID[w])
tmpPG <- tmpPG[which(tmpPG$value %in% IDs),]
source(parSrc, local = FALSE)
allProteins_mapped <- try(setNames(lapply(allTaxIDs, function(txid) { #txid <- allTaxIDs[1]
  kol <- c("Protein ID", "Common Name", "TAIR")
  kol <- kol[which(kol %in% colnames(db))]
  tmpDB <- db[which((db$`Protein ID` %in% IDs)&(db$TaxID == txid)), kol]
  x <- gsub("^cRAP[0-9]{3}", "", gsub("^CON__", "", tmpDB$`Protein ID`))
  y <- split(x, ceiling(seq_along(x)/20))
  n <- length(y)
  clusterExport(parClust, c("txid", "y", "tmpDB"), envir = environment())
  res <- parLapply(parClust, 1:n, function(i) {
    ids <- y[[i]]
    #a <- try(rbioapi::rba_string_map_ids(ids, txid), silent = TRUE)
    #if ("try-error" %in% class(a)) {
    try({
      rqst <- paste0("https://string-db.org/api/tsv/get_string_ids?identifiers=", paste(ids, collapse = "%0d"),
                     "&species=", txid)
      response <- httr::PUT(rqst, encode = "json")
      if (httr::http_status(response)$category != "Success") { Sys.sleep(2) }
      a <- httr::content(response, "parsed", show_col_types = FALSE)
      a$queryItem <- ids[a$queryIndex+1]
      m <- match(a$queryItem, tmpDB$`Protein ID`)
      a$Name <- tmpDB$`Common Name`[m]
      if ("TAIR" %in% colnames(tmpDB)) { # (For our dear plant-people)
        a$TAIR <- tmpDB$TAIR[m]
      }
    }, silent = TRUE)
    #}
    rs <- list(Outcome = (!"try-error" %in% class(a)))
    if (rs$Outcome) { rs$Result <- a }
    return(rs)
  })
  w <- which(vapply(res, function(rs) { (rs$Outcome)&&("data.frame" %in% class(rs$Result)) }, TRUE))
  res <- res[w]
  if (length(res)) {
    res <- lapply(res, function(rs) { rs$Result })
    res <- plyr::rbind.fill(res)
  } else {
    res <- data.frame(Error = c())
  }
  return(res)
}), paste0("TaxID_", allTaxIDs)), silent = TRUE)
w <- which(!is.null(sapply(allProteins_mapped, nrow)))
allProteins_mapped <- allProteins_mapped[w]
if (4 %in% WhTsts) {
  tmpMap <- listMelt(strsplit(PG$`Leading protein IDs`, ";"), PG$id)
}
# See https://string-db.org/help/api/ for API help
if (length(WhTsts)&&length(allProteins_mapped)) {
  filtersDF <- lapply(WhTsts, function(tt) { #tt <- 1 #tt <- 2 #tt <- 3 #tt <- 4
    Filt <- Reg_filters[[Tsts[tt]]]$"By condition"
    nms <- sort(names(Filt))
    filtDF <- data.frame(Group = Tsts[tt],
                         Name = nms,
                         Type = c(1, 2)[(tt == 4)+1])
    filtDF$Filter <- lapply(filtDF$Name, function(nm) {
      Filt[[nm]]
    })
    filtDF$Reg <- filtDF$TaxID <- c()
    w <- which(vapply(filtDF$Filter, function(flt) {
      sum(c("Ratios", "Filter") %in% names(flt)) == 2
    }, TRUE))
    filtDF <- filtDF[w,]
    nr <- nrow(filtDF)
    if (nr) {
      filtDF$W <- lapply(1:nr, function(x) { #x <- 1
        w <- filtDF$Filter[[x]]$Filter
        typ <- filtDF$Type[x]
        if (length(w)) {
          if (typ == 1) {
            w <- w[which(PG$"Potential contaminant"[w] != "+")]
          }
          if (typ == 2) {
            w <- w[which(allSAINTs$"Potential contaminant"[w] != "+")]
          }
        }
        return(w)
      })
      filtDF$Reg <- lapply(1:nr, function(x) { #x <- 1
        w <- filtDF$W[[x]]
        if (length(w)) {
          typ <- filtDF$Type[x]
          lfc <- filtDF$Filter[[x]]$Ratios[w]
          if (typ == 1) {
            reg <- data.frame("Leading protein IDs" = PG$"Leading protein IDs"[w],
                              "PG id" = PG$id[w],
                              "logFC" = lfc, # Here it would be neat to also have expression in the sample,
                              # but for this I need to rewrite filters!!!!!
                              check.names = FALSE)
            
            reg2 <- set_colnames(listMelt(strsplit(reg$"Leading protein IDs", ";"), 1:nrow(reg)), c("ID", "Row"))
            reg2$logFC <- reg$logFC[reg2$Row]
            reg2$"PG id" <- reg$"PG id"[reg2$Row]
            reg3 <- set_colnames(aggregate(reg2$logFC, list(reg2$ID), mean, na.rm = TRUE), c("ID", "logFC"))
            reg2 <- set_colnames(aggregate(reg2$"PG id", list(reg2$ID), unique), c("ID", "PG id")) # Should be only one value since we are using Leading protein IDs, which are unique
            reg2$logFC <- reg3$logFC[match(reg2$ID, reg3$ID)]
            reg2$Gene <- db$Gene[match(reg2$ID, db$`Protein ID`)]
          }
          if (typ == 2) {
            reg2 <- data.frame("ID" = allSAINTs$Protein[w],
                               "logFC" = lfc, # Here it would be neat to also have expression in the sample,
                               # but for this I need to rewrite filters!!!!!
                               "Gene" = allSAINTs$Gene[w],
                               check.names = FALSE)
            reg2$"PG id" <- tmpMap$L1[match(reg2$ID, tmpMap$value)]
          }
          reg2$TaxID <- db$TaxID[match(reg2$ID, db$"Protein ID")]
          reg2$ID <- gsub("^cRAP[0-9]{3}", "", gsub("^CON__", "", reg2$ID))
          reg2 <- aggregate(1:nrow(reg2), list(reg2$TaxID), function(x) {
            list(reg2[x,])
          })
          reg2 <- setNames(reg2$x, paste0("TaxID_", reg2$Group.1))
        } else {
          reg2 <- NA
        }
        return(reg2)
      })
      filtDF$W <- NULL
      filtDF$Filter <- NULL
      uTx <- unique(unlist(lapply(1:nr, function(x) { names(filtDF$Reg[[x]]) })))
      filtDF <- lapply(uTx, function(tx) { #tx <- uTx[1] #tx <- uTx[2]
        tmp <- filtDF
        nms <- lapply(tmp$Reg, names)
        tmp$Reg <- lapply(1:length(tmp$Reg), function(y) {
          if (tx %in% names(tmp$Reg[[y]])) {
            return(tmp$Reg[[y]][[tx]])
          } else { return() }
        })
        tmp$TaxID <- gsub("^TaxID_", "", tx)
        return(tmp)
      })
      filtDF <- plyr::rbind.fill(filtDF)
      if (!is.null(filtDF$Reg)) {
        filtDF <- filtDF[which(vapply(filtDF$Reg, function(x) { ("data.frame" %in% class(x))&&(nrow(x) > 0) }, TRUE)),]
      }
    }
    return(filtDF)
  })
  filtersDF <- plyr::rbind.fill(filtersDF)
  filtersDF <- filtersDF[which(paste0("TaxID_", filtersDF$TaxID) %in% names(allProteins_mapped)),]
  nr <- nrow(filtersDF)
  if (length(nr)) {
    filtersDF <- rbind(filtersDF, filtersDF)
    filtersDF$GraphType <- GraphTypes[2]
    filtersDF$GraphType[1:nr] <- GraphTypes[1]
    nr <- nr*2
    txidsTst <- (length(unique(filtersDF$TaxID)) > 1)+1
    source(parSrc, local = FALSE)
    clusterExport(parClust, c("wd", "filtersDF", "allProteins_mapped", "tmpPG", "GraphTypes", "Exp", "txidsTst"), envir = environment())
    tstSTRINGs <- parLapply(parClust, 1:nr, function(i) { #i <- 1 #i <- 4 #i <- 10 #i <- 11
      fltNm <- proteoCraft::cleanNms(filtersDF$Name[i], rep = "_")
      regTbl <- filtersDF$Reg[[i]]
      grphType <- filtersDF$GraphType[i]
      tstbee <- filtersDF$Group[i]
      txid <- filtersDF$TaxID[i]
      intNet_I <- NA
      STRINGplot_I <- NA
      img_I <- NA
      grphNm <- paste0(c("", paste0("taxID_", txid, "_"))[txidsTst],
                       gsub("[:\\*\\?<>\\|/| ,]+", "_", tstbee), "_", gsub("[:\\*\\?<>\\|/| ,]+", "_",
                                                                           gsub("\\(|\\)", "",
                                                                                fltNm)))
      imgpath <- paste0(wd, "/STRINGdb/", grphType, "/", grphNm, ".png")
      #cat(grphNm, "\n")
      #typ <- filtersDF$Type[i]
      proteins_mapped <- allProteins_mapped[[paste0("TaxID_", txid)]]
      proteins_mapped <- proteins_mapped[which(proteins_mapped$queryItem %in% regTbl$ID),]
      if (nrow(proteins_mapped) > 1) {
        m <- match(proteins_mapped$queryItem, regTbl$ID)
        proteins_mapped[, c("logFC", "PG id")] <- regTbl[m, c("logFC", "PG id")]
        #if (typ == 2) {
        #  proteins_mapped$PG <- tmpPG$L1[match(proteins_mapped$queryItem, tmpPG$value)]
        #}
        w <- which(nchar(proteins_mapped$Name) > 25)
        proteins_mapped$Name[w] <- paste0(substr(proteins_mapped$Name[w], 1, 22), "...")
        proteins_mapped$Name <- do.call(paste, c(proteins_mapped[, c("queryItem", "Name")],
                                                 sep = "\n"))
        #if (typ == 1) {
        KOL <- "queryItem"
        #}
        #if (typ == 2) { KOL <- "PG" }
        KOL <- c(KOL, "Name")
        if ("TAIR" %in% colnames(proteins_mapped)) { # (For our dear plant-people)
          KOL <- c(KOL, "TAIR")
        }
        if (exists("intNet")) { rm(intNet) }
        # intNet <- try(rbioapi::rba_string_interactions_network(proteins_mapped$stringId,
        #                                                        txid,
        #                                                        500, network_type = tolower(grphType)),
        #               silent = TRUE)
        rqst <- paste0("https://string-db.org/api/tsv/network?identifiers=",
                       paste(proteins_mapped$stringId, collapse = "%0d"),
                       "&species=", txid,
                       "&network_type=", tolower(grphType),
                       "&required_score=500&network_flavor=confidence&show_query_node_labels=1&hide_disconnected_nodes=1")
        try({
          response <- httr::POST(rqst, encode = "json")
          if (httr::http_status(response)$category != "Success") { Sys.sleep(2) }
          intNet <- httr::content(response, "parsed", show_col_types = FALSE)
        })
        if (!"try-error" %in% class(intNet)) { # This may fail if we have too many nodes! The current limit is 2000.
          # They suggest to use their "Cytoscape stringApp", this may be worth looking into.
          if ((!is.null(intNet))&&("data.frame" %in% class(intNet))&&(nrow(intNet))) {
            for (ab in c("A", "B")) {
              m <- match(intNet[[paste0("stringId_", ab)]], proteins_mapped$stringId)
              intNet[, paste0(c("Original_", paste0(KOL, "_")), ab)] <- proteins_mapped[m, c("queryItem", KOL)]
            }
            stopifnot("character" %in% class(intNet$Original_A), "character" %in% class(intNet$Original_B))
            w <- which(proteins_mapped$stringId %in% c(intNet$stringId_A, intNet$stringId_A))
            if (length(w)) {
              intNet_I <- list("Network" = intNet,
                               "Mapping" = proteins_mapped)
              rqst <- paste0("https://string-db.org/api/highres_image/network?identifiers=",
                             paste(proteins_mapped$queryItem[w], collapse = "%0d"),
                             "&species=", txid,
                             "&network_type=", tolower(grphType),
                             "&required_score=500&network_flavor=confidence&show_query_node_labels=1&hide_disconnected_nodes=1")
              try({
                # graph <- rbioapi::rba_string_network_image(proteins_mapped$queryItem[w],
                #                                            "highres_image",
                #                                            imgpath,
                #                                            txid,
                #                                            required_score = 500,
                #                                            network_flavor = "confidence",
                #                                            use_query_labels = TRUE,
                #                                            flat_nodes = TRUE,
                #                                            hide_disconnected_nodes = TRUE # Not good enough on its own, 
                # )
                response <- httr::POST(rqst, encode = "json")
                if (httr::http_status(response)$category != "Success") { Sys.sleep(2) }
                img_I <- httr::content(response, "raw")
                writeBin(img_I, imgpath)
              }, silent = TRUE)
            }
          }
        }
      }
      lst <- list(Name = grphNm,
                  Data = intNet_I,
                  Bin = img_I,
                  STRINGplot = imgpath)
      return(lst)
    })
    STRINGplots %<o% setNames(lapply(GraphTypes, function(grphType) { #grphType <- GraphTypes[1] #grphType <- GraphTypes[2]
      w <- which(filtersDF$GraphType == grphType)
      sapply(w, function(x) { tstSTRINGs[[x]]$STRINGplot })
    }), GraphTypes)
    intNets %<o% setNames(lapply(GraphTypes, function(grphType) {
      w <- which(filtersDF$GraphType == grphType)
      w <- w[which(vapply(w, function(x) { "list" %in% class(tstSTRINGs[[x]]$Data) }, TRUE))]
      setNames(lapply(w, function(x) { tstSTRINGs[[x]]$Data }),
               vapply(w, function(x) { tstSTRINGs[[x]]$Name }, ""))
    }), GraphTypes)
    for (grphType in GraphTypes) { #grphType <- GraphTypes[1]
      if (length(STRINGplots[[grphType]])) {
        fls <- STRINGplots[[grphType]]
        fls <- fls[which(file.exists(fls))]
        if (length(fls) > 1) {
          PNGs <- parLapply(parClust, fls, readPNG)
          g <- parLapply(parClust, PNGs, function(png) { grid::rasterGrob(png, interpolate = TRUE) })
          dim <- ceiling(sqrt(length(g)))
          l <- length(gsub(".*/|\\.png$", "", fls))
          tst <- data.frame(tile = gsub(".*/|\\.png$", "", fls))
          tst$x <- rep(1:dim, dim)[1:l]
          tst$y <- as.numeric(sapply(dim:1, function(x) { rep(x, dim) }))[1:l]
          tst$Grob <- g
          tst$Lab <- gsub("\\)$", "", gsub("\\)_-_\\(", "\nvs\n", gsub("_by_condition_(\\()?", "\n", tst$tile)))
          ttl <- paste0(tolower(grphType), " STRING plots")
          plot <- ggplot(tst) +
            geom_tile(aes(x = x, y = y), width = 0.96, height = 0.96, fill = "white") + coord_equal() + theme_bw() +
            xlim(0.5, dim+0.5) + ylim(0.5, dim+1) +
            theme(axis.line = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(),
                  axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),
                  legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), plot.background = element_blank())
          for (i in 1:nrow(tst)) {
            plot <- plot + annotation_custom(tst$Grob[[i]], xmin = tst$x[i]-0.48, xmax = tst$x[i]+0.48, ymin = tst$y[i]-0.48, ymax = tst$y[i]+0.48)
          }
          plot <- plot +
            geom_text(aes(label = Lab, x = x, y = y + 0.48), cex = 2.5) +
            geom_text(label = paste0(grphType, " interactions"), x = dim/2+0.5, y = dim+0.7, cex = 5)
          poplot(plot)
          ggsave(paste0(wd, "/STRINGdb/", grphType, "/All ", ttl, ".jpeg"), plot,
                 dpi = 120, width = 100, height = 100, units = "in", limitsize = FALSE)
          ggsave(paste0(wd, "/STRINGdb/", grphType, "/All ", ttl, ".pdf"), plot,
                 dpi = 120, width = 100, height = 100, units = "in", limitsize = FALSE)
          ReportCalls <- AddPlot2Report(Space = FALSE, Jpeg  = FALSE)
        }
      }
    }
    ReportCalls <- AddSpace2Report()
    #
    # Cytoscape
    if (CytoScape) {
      wL <- which(vapply(GraphTypes, function(grphType) { length(intNets[[grphType]]) }, 1) > 0)
      if (length(wL)) {
        #
        # Check that Cytoscape is running and available
        Src <- paste0(libPath, "/extdata/R scripts/Sources/Cytoscape_init.R")
        #rstudioapi::documentOpen(Src)
        source(Src, local = FALSE)
        #
        # Create directory for Cytoscape networks
        #dirs <- paste0(wd, "/Cytoscape/", GraphTypes) # No!
        dirs <- paste0(wd, "/STRINGdb/", GraphTypes) # These are Cytoscape versions of STRINGdb graphs,
        # the logical place to save them is in STRINGdb!!!
        for (dir in dirs) { if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) } }
        dirlist <- unique(c(dirlist, dir))
        #
        # Use Cytoscape to represent networks of gated proteins with overlayed logFC
        msg <- "Creating Cytoscape network .cx files..."
        ReportCalls <- AddMsg2Report()
        cat("Creating Cytoscape networks...\n")
        for (grphType in GraphTypes[wL]) { #grphType <- GraphTypes[wL][1]
          cat(" ->", grphType, "\n")
          for (Nm in names(intNets[[grphType]])) { #Nm <- names(intNets[[grphType]])[1]
            cat("   -", Nm, "\n")
            intNet <- intNets[[grphType]][[Nm]]$Network
            mappings <- intNets[[grphType]][[Nm]]$Mapping
            # Legacy code for if doing by protein, not protein groups
            #mappings$name <- vapply(match(mappings$queryItem, db$"Protein ID"), function(m) {
            #  paste(c(paste0("Nm: ", db$"Common Name"[m]),
            #          paste0("Pr: ", db$"Protein ID"[m]),
            #          paste0("Gn: ", db$Gene[m])), collapse = "\n")
            #}, "")
            #if (!grepl("^SAINTexpress_", Nm)) {
            mappings$name <- vapply(match(mappings$queryItem, db$`Protein ID`), function(m) {
              x <- c(paste0("Nm: ", db$`Common Name`[m]),
                     paste0("Pr: ", db$`Protein ID`[m]),
                     paste0("Gn: ", db$Gene[m]))
              if ("TAIR" %in% colnames(mappings)) { x <- c(paste0("TAIR: ", db$TAIR[m], " "), x) }
              return(paste(x, collapse = "\n"))
            }, "")
            #  mappings$PG <- allSAINTs$PG_id[match(mappings$queryItem, allSAINTs$Protein)]
            #} else {
            #   mappings$name <- vapply(match(mappings$PG, PG$id), function(m) {
            #     pr <- unlist(strsplit(PG$"Leading protein IDs"[m], ";"))
            #     if (length(pr) > 1) { pr <- c(pr[1], "...") }
            #     gn <- unlist(strsplit(PG$Genes[m], ";"))
            #     if (length(gn) > 1) { gn <- c(gn[1], "...") }
            #     x <- c(paste0("Nm: ", PG$"Common Name (short)"[m]),
            #            paste0("Pr: ", paste(pr, collapse = ";")),
            #            paste0("Gn: ", paste(gn, collapse = ";")))
            #     if ("TAIR" %in% colnames(mappings)) { x <- c(paste0("TAIR: ", PG$TAIR[m], " "), x) }
            #     return(paste(x, collapse = "\n"))
            #   }, "") 
            #   mappings$PG <- allSAINTs$PG_id[match(mappings$queryItem, allSAINTs$Protein)]
            # }
            mappings$"Av. log10 expression" <- PG$"Av. log10 abundance"[match(mappings$PG, PG$id)]
            #if (!grepl("^SAINTexpress_", Nm)) {
            # intNet$Linkage <- do.call(paste, c(intNet[, paste0("PG_", c("A", "B"))], sep = "_"))
            # intNet <- Isapply(strsplit(unique(intNet$Linkage), "_"), unlist)
            # colnames(intNet) <- paste0("PG_", c("A", "B"))
            # for (ab in c("A", "B")) { #ab <- "A"
            #   intNet[[paste0("Name_", ab)]] <- mappings$name[match(intNet[[paste0("PG_", ab)]], mappings$PG)]
            # }
            #} else {
            intNet$Name_A <- mappings$name[match(intNet$Name_A, mappings$Name)]
            intNet$Name_B <- mappings$name[match(intNet$Name_B, mappings$Name)]
            #}
            kol <- c("Name_A", "Name_B")
            gD <- igraph::simplify(igraph::graph_from_data_frame(intNet[, kol], directed = FALSE))
            seq <- igraph::V(gD)
            logFC <- mappings$logFC[match(names(seq), mappings$name)]
            w <- which(is.na(logFC))
            logFC[w] <- 0
            gD <- igraph::set_vertex_attr(gD, "logFC", igraph::V(gD), logFC)
            Xprss <- mappings$"Av. log10 expression"[match(names(seq), mappings$name)]
            gD <- igraph::set_vertex_attr(gD, "Avg_expression", igraph::V(gD), Xprss)
            #igraph::vcount(gD)
            #igraph::ecount(gD)
            degAll <- igraph::degree(gD, v = igraph::V(gD), mode = "all")
            betAll <- igraph::betweenness(gD, v = igraph::V(gD), directed = FALSE) / (((igraph::vcount(gD) - 1) * (igraph::vcount(gD)-2)) / 2)
            betAll.norm <- (betAll - min(betAll))/(max(betAll) - min(betAll))
            rm(betAll)
            dsAll <- igraph::similarity(gD, vids = igraph::V(gD), mode = "all", method = "dice")
            gD <- igraph::set_vertex_attr(gD, "degree", index = igraph::V(gD), value = degAll)
            gD <- igraph::set_vertex_attr(gD, "betweenness", index = igraph::V(gD), value = betAll.norm)
            #summary(gD)
            F1 <- function(x) {
              data.frame(V4 = dsAll[which(igraph::V(gD)$name == as.character(x$V1)),
                                    which(igraph::V(gD)$name == as.character(x$V2))])
            }
            dataSet.ext <- plyr::ddply(intNet[, kol], .variables = kol, function(x) data.frame(F1(x)))
            gD <- igraph::set_edge_attr(gD, "weight", index = igraph::E(gD), value = 0)
            gD <- igraph::set_edge_attr(gD, "similarity", index = igraph::E(gD), value = 0)
            for (i in 1:nrow(dataSet.ext)) {
              igraph::E(gD)[as.character(dataSet.ext$V1) %--% as.character(dataSet.ext$V2)]$weight <- as.numeric(dataSet.ext$V3)
              igraph::E(gD)[as.character(dataSet.ext$V1) %--% as.character(dataSet.ext$V2)]$similarity <- as.numeric(dataSet.ext$V4)
            }
            rm(dsAll, i, F1)
            #summary(gD)
            exprs <- expression({
              tst <- RCy3::createNetworkFromIgraph(gD, new.title = Nm)
              RCy3::layoutNetwork("force-directed defaultSpringLength=70 defaultSpringCoefficient=0.000003")
              rg <- max(abs(logFC))
              try(RCy3::setNodeColorMapping("logFC", c(-rg, 0, rg),
                                            c("#FF0000", "#999999", "#00FF00"), style.name = "default"), silent = TRUE)
              RCy3::setNodeShapeDefault("ellipse", style.name = "default")
              RCy3::setNodeFontSizeDefault(12, style.name = "default")
              try(RCy3::setNodeSizeMapping("Avg_expression",
                                           style.name = "default"), silent = TRUE)
              RCy3::exportNetwork(paste0(wd, "/Cytoscape/", grphType, "/", Nm), "CX")
              RCy3::deleteAllNetworks()
            })
            tst <- try(eval(exprs), silent = TRUE)
          }
        }
        # Then close Cytoscape:
        try({
          RCy3::closeSession(save.before.closing = FALSE)
          cmd <- paste0("taskkill/im \"", gsub(".+/", "", CytoScExe), "\" /f")
          #cat(cmd)
          shell(cmd)
        }, silent = TRUE)
      }
    }
  }
}
