#### Code chunk - Gene Ontology terms enrichment analysis
if ((Annotate)&&(enrichGO||globalGO)) {
  msg <- "Gene Ontology terms enrichment analysis"
  ReportCalls <- AddMsg2Report(Space = FALSE)
  packs <- c("GO.db", "topGO")
  for (pack in packs) {
    bioc_req <- unique(c(bioc_req, pack))
    biocInstall(pack)
  }
  GO_enrich.dat %<o% list()
  GO_enrich.FCRt %<o% list()
  GO_enrich.tbl %<o% list()
  GO_Plots %<o% list()
  Reg_GO_terms %<o% list()
  GO.enrich.MultiRefs %<o% (("GO.enrichment.Ref.Aggr" %in% colnames(Param))&&(!Param$GO.enrichment.Ref.Aggr %in% c("", "NA", NA)))
  if (GO.enrich.MultiRefs) { parse.Param.aggreg.2("GO.enrichment.Ref.Aggr") }
  #
  if (runClueGO) {
    # Initialize ClueGO
    Src <- paste0(libPath, "/extdata/R scripts/Sources/ClueGO_init.R")
    #rstudioapi::documentOpen(Src)
    source(Src, local = FALSE)
  }
  #
  if (enrichGO) {
    Tsts <- c("t-tests", "F-tests", "Localisation", "SAINTexpress")
    WhTsts <- which(Tsts %in% names(Reg_filters))
    if (!exists("SSD.Root")) { SSD.Root <- "" }
    for (tt in WhTsts) { #tt <- WhTsts[1L] #tt <- WhTsts[2L] #tt <- WhTsts[3L] #tt <- WhTsts[4L]
      tstrt <- Tsts[tt]
      stopifnot(!is.na(tstrt))
      ReportCalls <- AddMsg2Report(Msg = paste0("\n - ", tstrt), Space = FALSE)
      dir <- paste0(wd, "/Reg. analysis/GO enrich/", tstrt)
      if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
      dirlist <- unique(c(dirlist, dir))
      filt <- Reg_filters[[tstrt]]
      #By <- c("By condition", "By reference", "By analysis", "Whole dataset")
      By <- "By condition"
      By <- By[which(By %in% names(filt))]
      if (length(By)) {
        for (bee in By) { #bee <- By[1L]
          flt <- filt[[bee]]
          if (bee == "Whole dataset") { flt <- list("Whole dataset" = flt) }
          tstbee <- paste0(tstrt, "_", tolower(bee))
          if (length(flt)) {
            flt <- if (tt %in% 1L:2L) {
              flt[myContrasts$Contrast[which(myContrasts$Secondary == "")]]
            } else {
              flt[order(names(flt))]
            }
            reg <- setNames(lapply(flt, \(x) { list(x$Columns) }), names(flt))
            reg <- set_colnames(reshape::melt(reg), c("Name", "Bleh", "For"))
            reg$Bleh <- NULL
            fcr <- c(Prot.Rat.Root, Prot.Rat.Root, paste0("Mean ", SSD.Root, "log2(FC) - "))[tt]
            reg$ParentFC <- gsub(paste0(".*Re", c(rep("gulated", 2L), "-localized", "gulated")[tt], " - "), fcr, reg$Name)
            reg$FCname <- paste0(fcr, reg$For)
            tmpdat <- get(c("PG", "F_test_data", "PG", "allSAINTs")[tt])
            UF <- unique(reg$For)
            temp <- as.data.frame(do.call(cbind, lapply(UF, \(x) { #x <- UF[1L]
              x <- reg$ParentFC[which(reg$For == x)]
              x <- if (length(x) > 1L) { apply(tmpdat[, x], 1L, log_ratio_av) } else { tmpdat[[x]] }
              return(x)
            })))
            colnames(temp) <- paste0(fcr, UF)
            for (i in 1L:nrow(reg)) {
              Reg_filters[[tstrt]][[bee]][[reg$For[i]]]$Ratios <- temp[[paste0(fcr, reg$For[i])]]
            }
            if (tt %in% 1L:3L) {
              Kol2Add <- c("Leading protein IDs", "Protein IDs", "id", "Protein names", "No Isoforms", "Names", "Genes",
                           "Common Names", Param$Plot.labels, "GO", "GO-ID", "Potential contaminant")
              Kol2Add <- intersect(Kol2Add, colnames(PG))
              Kol2Add2 <- setdiff(Kol2Add, colnames(tmpdat))
              if (length(Kol2Add2)) {
                if (tt == 3L) { warning("This hasn't been tested for Localisation Analysis yet!") }
                tmpdat[, Kol2Add2] <- PG[match(tmpdat$Protein, PG$`Protein IDs`), Kol2Add2]
              }
              temp[, Kol2Add] <- tmpdat[, Kol2Add]
            } else {
              Kol2Add <- c("No Isoforms", "Gene", "Common Name", "GO", "GO-ID")
              temp$Protein <- allSAINTs$Protein
              temp$PG_id <- allSAINTs$PG_id
              temp[, Kol2Add] <- db[match(allSAINTs$Protein, db$`Protein ID`), Kol2Add]
            }
            GO_enrich.dat[[tstbee]] <- temp
            GO_enrich.FCRt[[tstbee]] <- fcr
            GO_enrich.tbl[[tstbee]] <- reg
            flt <- setNames(lapply(UF, \(x) { flt[[x]]$Filter }), UF)
            #flt <- setNames(lapply(UF, \(x) { flt[[x]]$Filter_down }), UF)
            #flt <- setNames(lapply(UF, \(x) { flt[[x]]$Filter_up }), UF)
            # see function code for defaults)
            ttr <- btr <- ""
            if (length(By) > 1L) { ttr <- btr <- paste0(tolower(bee), "_") }
            Ref.Filt <- tmpFilt <- setNames(lapply(names(flt), \(x) { 1L:nrow(GO_enrich.dat[[tstbee]]) }),
                                            names(flt))
            if (GO.enrich.MultiRefs) {
              Ref.Filt <- try(setNames(lapply(names(flt), \(x) { #x <- names(flt)[1L]
                m <- match(x, myContrasts$Contrast)
                A_ <- myContrasts$A_samples[[m]]
                B_ <- myContrasts$B_samples[[m]]
                y <- Exp.map[which(Exp.map$Ref.Sample.Aggregate %in% c(A_, B_)), GO.enrichment.Ref.Aggr$column]
                z <- Exp.map$Ref.Sample.Aggregate[which(Exp.map[[GO.enrichment.Ref.Aggr$column]] %in% y)]
                w1 <- which(apply(PG[, paste0(Prot.Expr.Root, z)], 1L, \(x) {
                  sum(is.all.good(x, 2L))
                }) > 0L)
                if (tt == 4L) {
                  w1 <- which(temp$PG_id %in% PG$id[w1])
                }
                w2 <- flt[[x]] # Required for if we are imputing missing values
                return(sort(unique(c(w1, w2))))
              }), names(flt)), silent = TRUE)
              if (inherits(Ref.Filt, "try-error")) {
                warning("Invalid \"GO.enrichment.Ref.Aggr\" argument: multiple references for GO enrichment are only feasible if each enrichment filter maps to a single reference! Skipping...")
                GO.enrich.MultiRefs <- FALSE
                Ref.Filt <- tmpFilt
              }
            }
            if ((length(Ref.Filt) > 1L)||(!is.na(Ref.Filt))) {
              flt <- setNames(lapply(names(flt), \(x) { flt[[x]][which(flt[[x]] %in% Ref.Filt[[x]])] }),
                              names(flt))
            }
            # Also save the reference filters 
            nms <- names(Reg_filters[[tstrt]][[bee]])
            for (nm in nms) {
              Reg_filters[[tstrt]][[bee]][[nm]]$Background_filter <- Ref.Filt[[nm]]
            }
            #
            Mode <- "regulated"
            if (tt %in% c(1L, 3L)) { dataType <- "PG" }
            if (tt == 4L) { dataType <- "Prot" }
            if ((!exists("GO_mappings"))&&(file.exists("GO_mappings.RData"))) {
              loadFun("GO_mappings.RData")
            }
            if ((!exists("GO_terms"))&&(file.exists("GO_terms.RData"))) { loadFun("GO_terms.RData") }
            #
            Src <- paste0(libPath, "/extdata/R scripts/Sources/GO_enrich.R")
            #rstudioapi::documentOpen(Src)
            source(Src, local = FALSE)
            #
            if (runClueGO) {
              clueGO_outDir <- dir
              clueGO_type <- "Enrichment (Right-sided hypergeometric test)"
              Src <- paste0(libPath, "/extdata/R scripts/Sources/ClueGO_enrich.R")
              #rstudioapi::documentOpen(Src)
              source(Src, local = FALSE)
            }
            #
            # Cleanup - do it now, not within sources!
            try(rm(list = allArgs), silent = TRUE)
            #
            GO_Plots[[tstbee]] <- goRES
            #
            if (!is.null(names(GO_Plots[[tstbee]]))) {
              if ("All_GO_terms" %in% names(GO_Plots[[tstbee]])) {
                GO_terms <- GO_Plots[[tstbee]]$All_GO_terms
                GO_Plots[[tstbee]]$All_GO_terms <- NULL
              }
              n2 <- names(GO_Plots[[tstbee]]$GO_plots)
              dir2 <- paste0(wd, "/", dir)
              for (ttl in n2) { #ttl <- n2[1L]
                plot <- GO_Plots[[tstbee]]$GO_plots[[ttl]]
                ReportCalls <- AddPlot2Report(Space = FALSE, Jpeg = FALSE)
              }
              if ((create_plotly)&&(!create_plotly_local)) { plot_ly[[paste0("GO plots - Regulated vs Observed - ", tstbee)]] <- GO_Plots[[tstbee]]$GO_plot_ly }
              temp <- GO_Plots[[tstbee]]$GO_terms
              temp$Mapping <- NULL
              if ("Offspring" %in% colnames(temp)) {
                temp$Offspring <- vapply(temp$Offspring, paste, "", collapse = ";")
              }
              temp$`Protein table row(s)` <- NULL
              colnames(temp) <- cleanNms(colnames(temp))
              gn <- grep("^Genes", colnames(temp), value = TRUE)
              pr <- grep("^Proteins", colnames(temp), value = TRUE)
              pg <- grep("^PG IDs", colnames(temp), value = TRUE)
              kn <- grep("^Count", colnames(temp), value = TRUE)
              pv <- grep("^Pvalue", colnames(temp), value = TRUE)
              zs <- grep("^(Z-score|\\(N_Up - N_Down\\)/sqrt\\(Tot\\.\\))", colnames(temp), value = TRUE)
              lf <- grep("^logFC", colnames(temp), value = TRUE)
              si <- grep("^Significance", colnames(temp), value = TRUE)
              #lp <- grep("^Leading protein IDs", colnames(temp), value = TRUE)
              #kl <- colnames(temp)[which(!colnames(temp) %in% c(gn, kn, pv, zs, lf, si, pg, lp))]
              kl <- colnames(temp)[which(!colnames(temp) %in% c(gn, kn, pv, zs, lf, si, pg, pr))]
              #temp <- temp[, c(kl, si, gn, pg, lp, kn, pv, zs, lf)]
              temp <- temp[, c(kl, si, gn, pg, pr, kn, pv, zs, lf)]
              w <- apply(temp[, pv, drop = FALSE], 1L, \(x) { sum(!is.na(x)) }) > 0L
              temp <- temp[w,]
              tst <- apply(temp[, kn, drop = FALSE], 1L, \(x) { sum(x[which(!is.na(x))]) })
              temp <- temp[order(tst, decreasing = TRUE),]
              temp <- temp[order(temp$Ontology, decreasing = FALSE),]
              Reg_GO_terms[[tstbee]] <- temp
              #temp <- Reg_GO_terms[[tstbee]]
              write.csv(temp, file = paste0(dir, "/GO terms - ", tstbee, ".csv"), row.names = FALSE)
              w <- which(vapply(colnames(temp), \(x) { is.character(temp[[x]]) }, TRUE))
              if (length(w)) {
                for (i in w) { #i <- w[1L]
                  w1 <- which(nchar(temp[[colnames(temp)[i]]]) > ExcelMax)
                  if (length(w1)) {
                    temp[[colnames(temp)[i]]][w1] <- paste0(substr(temp[[colnames(temp)[i]]][w1], 1L, ExcelMax-3L), "...")
                  }
                }
              }
              require(openxlsx)
              HdrStl <- createStyle(textDecoration = "bold", halign = "center", valign = "center", wrapText = TRUE,
                                    numFmt = "TEXT", fontSize = 12L)
              wb <- createWorkbook()
              kount <- 0L
              for (ont in names(Ontologies)) { #ont <- names(Ontologies)[1L]
                w <- which(temp$Ontology == Ontologies[ont])
                if (length(w)) {
                  kount <- kount + 1L
                  addWorksheet(wb, ont)
                  writeData(wb, ont, temp[w,])
                  setRowHeights(wb, ont, 1L, 60L)
                  setColWidths(wb, ont, 1L:ncol(temp), 12L)
                  setColWidths(wb, ont, which(colnames(temp) == "Term"), 45L)
                  setColWidths(wb, ont, which(colnames(temp) %in% gn), 20L)
                  setColWidths(wb, ont, which(colnames(temp) %in% pr), 20L)
                  setColWidths(wb, ont, which(colnames(temp) %in% pg), 20L)
                  addStyle(wb, ont, HdrStl, 1L, 1L:ncol(temp))
                  addStyle(wb, ont, createStyle(numFmt = "0"), 2L:(length(w)+1L),
                           which(colnames(temp) %in% kn), gridExpand = TRUE)
                  addStyle(wb, ont, createStyle(numFmt = "0.000"), 2L:(length(w)+1L),
                           which(colnames(temp) %in% c(zs, lf, pv)), gridExpand = TRUE)
                }
              }
              if (kount) {
                saveWorkbook(wb, paste0(dir, "/GO terms - ", tstbee, ".xlsx"), overwrite = TRUE)
                #openXL(paste0(dir, "/GO terms - ", tstbee, ".xlsx"))
                Kol2 <- grep(paste0("^Significance - .+", " ", max(BH.FDR)*100L, "%"),
                             colnames(Reg_GO_terms[[tstbee]]), value = TRUE)
                if (length(Kol2)) {
                  w <- if (length(Kol2) > 1L) {
                    which(apply(Reg_GO_terms[[tstbee]][,Kol2], 1L, \(x) {"+" %in% x}))
                  } else { which(vapply(Reg_GO_terms[[tstbee]][, Kol2], \(x) { "+" %in% x }, TRUE)) }
                  write.csv(Reg_GO_terms[[tstbee]][w,], file = paste0(dir, "/Regulated GO terms - ", tstbee, ".csv"),
                            row.names = FALSE)
                }
                # Summary table and heatmap of number of regulated GO terms
                Kol3 <- grep(paste0("^Significance - [^ ]+ [1-9][0-9]*\\.*[0-9]*%$"),
                             colnames(GO_Plots[[tstbee]]$GO_terms), value = TRUE)
                tst <- as.numeric(gsub(paste0("^Significance - [^ ]+ |%$"), "", Kol3))
                Kol3 <- Kol3[which(tst == max(tst))]
                N <- length(Kol3)
                if (N > 1L) {
                  temp <- as.data.frame(matrix(rep("", (N+1)^2), ncol = N+1L))
                  W <- lapply(Kol3, \(x) { which(GO_Plots[[tstbee]]$GO_terms[[x]] == "+") })
                  names(W) <- gsub(paste0("^Significance - | ", max(BH.FDR)*100L, "%$"), "", Kol3)
                  temp[2L:(N+1L), 1L] <- temp[1L, 2L:(N+1L)] <- names(W)
                  for (i in 2L:(N+1L)) { #i <- 2
                    temp[i, 2L:(N+1L)] <- vapply(Kol3, \(x) {
                      sum((GO_Plots[[tstbee]]$GO_terms[[x]] == "+")&(GO_Plots[[tstbee]]$GO_terms[[Kol3[i]]] == "+"),
                          na.rm = TRUE)
                    }, 1L)
                  }
                  names(W) <- cleanNms(gsub(" [0-9]+%$", "", names(W)))
                  tst <- lengths(strsplit(names(W), " - "))
                  tst <- (min(tst) > 1L)&(length(unique(tst)) == 1L)
                  if (tst) {
                    tst <- as.data.frame(t(sapply(strsplit(names(W), " - "), unlist)))
                    l <- apply(tst, 2L, \(x) { length(unique(x)) })
                    tst <- tst[, which(l > 1L)]
                    names(W) <- do.call(paste, c(tst, sep = " - "))
                  }
                  temp[2L:(N+1L), 1L] <- temp[1L, 2L:(N+1L)] <- names(W)
                  nm <- paste0("N. of co-regulated GO terms\n", tstrt, "\n(", tolower(bee), ")")
                  write.csv(temp, file = paste0(dir, "/", gsub("\n", " - ", gsub("\n\\(", " (", nm)), ".csv"), row.names = FALSE)
                  temp2 <- temp[2L:(N+1L), 2L:(N+1L)]
                  colnames(temp2) <- temp[1L, 2L:(N+1L)]
                  rownames(temp2) <-  temp[2L:(N+1L), 1L]
                  for (i in 1L:nrow(temp2)) { temp2[[i]] <- as.numeric(temp2[[i]]) }
                  if (max(is.all.good(unlist(temp2))) > 0L) {
                    temp2 <- as.matrix(temp2)
                    basic.heatmap(temp2, "N. of co-regulated GO terms", paste0(tstrt, "\n(", tolower(bee), ")"),
                                  save = c("pdf", "jpeg"), folder = dir)
                  }
                } else { warning(paste0(tstrt, " performed for only one condition, skipping.")) }
              }
            }
          } else { warning(paste0("Filter ", tstbee, " has length 0, skipping.")) }
        }
      } else { warning(paste0("No filters available for ", tstrt, ", skipping.")) }
    }
  }
  if (globalGO) {
    msg <- " - Dataset"
    ReportCalls <- AddMsg2Report(Space = FALSE)
    #
    dir <- paste0(wd, "/Reg. analysis/GO enrich/Dataset")
    if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
    dirlist <- unique(c(dirlist, dir))
    #
    Mode <- "dataset"
    dataType <- "PG"
    if ((!exists("GO_mappings"))&&(file.exists("GO_mappings.RData"))) {
      loadFun("GO_mappings.RData")
    }
    if ((!exists("GO_terms"))&&(file.exists("GO_terms.RData"))) {
      loadFun("GO_terms.RData")
    }
    #
    #try(rm(list = allArgs), silent = TRUE)
    #
    Src <- paste0(libPath, "/extdata/R scripts/Sources/GO_enrich.R")
    #rstudioapi::documentOpen(Src)
    source(Src, local = FALSE)
    #
    # Off because it was taking too long and failing sometimes
    # ClueGO really cannot handle too large gene lists, which are the norm here for dataset analysis
    #if (runClueGO) {
    #  clueGO_outDir <- dir
    #  clueGO_type <- "Enrichment/Depletion (Two-sided hypergeometric test)"
    #  Src <- paste0(libPath, "/extdata/R scripts/Sources/ClueGO_enrich.R")
    #  rstudioapi::documentOpen(Src)
    #  source(Src, local = FALSE)
    #}
    #
    # Cleanup - do it now, not within sources!
    try(rm(list = allArgs), silent = TRUE)
    #
    # Quick Fisher exact test on GO terms of interest - looking only at observed data per group
    Src <- paste0(libPath, "/extdata/R scripts/Sources/interestGO_Fisher.R")
    #rstudioapi::documentOpen(Src)
    source(Src, local = FALSE)
    #
    GO_Plots_2 %<o% goRES
    #
    if ((!is.null(GO_Plots_2))&&("All_GO_terms" %in% names(GO_Plots_2))) {
      GO_terms <- GO_Plots_2$All_GO_terms
      GO_Plots_2$All_GO_terms <- NULL
    }
    if ((create_plotly)&&(!create_plotly_local)) { plot_ly$"GO plots - Observed dataset vs Theoretical proteome" <- GO_Plots_2$GO_plot_ly }
  }
  l <- length(DatAnalysisTxt)
  DatAnalysisTxt[l] <- paste0(DatAnalysisTxt[l],
                              " GO terms enrichment analysis was performed, comparing for each test regulated against observed protein groups, using topGO",
                              c("", "and ClueGO")[runClueGO+1L], ".")
}
