#### Code chunk - Gene Ontology terms enrichment analysis
warning("Please, me dude: re-write me, this is getting ugly...")
setwd(wd)
if (Annotate && (enrichGO || globalGO)) {
  if (!exists("GO_mappings_fl")) { GO_mappings_fl %<o% paste0(wd, "/GO_mappings.RDS") }
  if (!exists("GO_terms_fl")) { GO_terms_fl %<o% paste0(wd, "/GO_terms.RDS") }
  #
  cat("Gene Ontology terms enrichment analysis\n")
  packs <- c("GO.db", "topGO")
  for (pack in packs) {
    bioc_req <- unique(c(bioc_req, pack))
    biocInstall(pack)
  }
  GO_Plots <- list()
  Reg_GO_terms <- list()
  GO_enrich.dat <- list()
  GO_enrich.FCRt <- list()
  GO_enrich.tbl <- list()
  GO.enrich.MultiRefs %<o% (("GO.enrichment.Ref.Aggr" %in% colnames(Param)) && (!Param$GO.enrichment.Ref.Aggr %in% c("", "NA", NA)))
  if (GO.enrich.MultiRefs) { parse.Param.aggreg.2("GO.enrichment.Ref.Aggr") }
  #
  if (runClueGO) {
    # Initialize ClueGO
    Src <- paste0(libPath, "/extdata/Sources/ClueGO_init.R")
    #rstudioapi::documentOpen(Src)
    source(Src)
  }
  #
  if (enrichGO) {
    Tsts <- c("t-tests", "F-tests", "Localisation", "SAINTexpress")
    WhTsts <- which(Tsts %in% names(Reg_filters))
    if (!exists("SSD.Root")) { SSD.Root <- "" }
    for (tt in WhTsts) { #tt <- 1L #tt <- 2L #tt <- 3L #tt <- 4L
      tstrt <- Tsts[tt]
      stopifnot(!is.na(tstrt))
      cat("\n - ", tstrt, "\n")
      myDir <- paste0(wd, "/Reg. analysis/GO enrich/", tstrt)
      if (!dir.exists(myDir)) { dir.create(myDir, recursive = TRUE) }
      dirlist <- unique(c(dirlist, myDir))
      filt <- Reg_filters[[tstrt]]
      #By <- c("By condition", "By reference", "By analysis", "Whole dataset")
      By <- "By condition"
      By <- base::intersect(By, names(filt))
      if (length(By)) {
        for (bee in By) { #bee <- By[1L]
          flt <- filt[[bee]]
          if (bee == "Whole dataset") { flt <- list("Whole dataset" = flt) }
          tstbee <- paste0(tstrt, "_", tolower(bee))
          if (length(flt)) {
            flt <- if (tt %in% 1L:2L) {
              flt[base::intersect(myContrasts$Contrast[!myContrasts$isDouble], names(flt))]
            } else {
              if (tt == 4L) {
                flt[myContrasts$Contrast[(!myContrasts$isDouble) & myContrasts$`Up-only`]]
              } else {
                flt[order(names(flt))]                
              }
            }
          }
          if (length(flt)) {
            reg <- lapply(flt, \(x) { x$Columns })
            reg <- listMelt(reg, names(flt), c("Name", "For"))
            fcr <- c(Prot.Rat.Root, Prot.Rat.Root, paste0("Mean ", SSD.Root, "log2(FC) - "), "log2(FC) - ")[tt]
            reg$ParentFC <- sub(paste0(".*Re", c(rep("gulated", 2L), "-localized", "gulated")[tt], " - "), fcr, reg$Name)
            reg$FCname <- paste0(fcr, reg$For)
            tmpdat <- get(c("PG", "F_test_data", "PG", "allSAINTs")[tt])
            UF <- unique(reg$For)
            temp <- as.data.frame(do.call(cbind, lapply(UF, \(x) { #x <- UF[1L]
              x <- reg$ParentFC[reg$For == x]
              x <- if (length(x) > 1L) { apply(tmpdat[, x], 1L, log_ratio_av) } else { tmpdat[[x]] }
              return(x)
            })))
            colnames(temp) <- paste0(fcr, UF)
            for (i in 1L:nrow(reg)) {
              Reg_filters[[tstrt]][[bee]][[reg$For[i]]]$Ratios <- temp[[paste0(fcr, reg$For[i])]]
              flt[[reg$For[i]]]$Ratios <- temp[[paste0(fcr, reg$For[i])]]
            }
            if (tt %in% 1L:3L) {
              Kol2Add <- c("Leading protein IDs", "Protein IDs", "id", "Protein names", "No Isoforms", "Names", "Genes",
                           "Common Names", Param$Plot.labels, "GO", "GO-ID", "Potential contaminant")
              Kol2Add <- base::intersect(Kol2Add, colnames(PG))
              Kol2Add2 <- base::setdiff(Kol2Add, colnames(tmpdat))
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
              #
            }
            GO_enrich.dat[[tstbee]] <- temp
            GO_enrich.FCRt[[tstbee]] <- fcr
            GO_enrich.tbl[[tstbee]] <- reg
            flt <- setNames(lapply(UF, \(x) {
              if (tt == 4L) {
                flt[[x]]$prot_Filter
              } else {
                flt[[x]]$Filter
              }
            }), UF)
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
                y <- Exp.map[Exp.map$Ref.Sample.Aggregate %in% c(A_, B_), GO.enrichment.Ref.Aggr$column]
                z <- Exp.map$Ref.Sample.Aggregate[Exp.map[[GO.enrichment.Ref.Aggr$column]] %in% y]
                w1 <- which(apply(PG[, paste0(Prot.Expr.Root, z)], 1L, \(x) {
                  sum(is.finite(x))
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
            if ((length(Ref.Filt) > 1L) || (!is.na(Ref.Filt))) {
              flt <- setNames(lapply(names(flt), \(x) { base::intersect(flt[[x]], Ref.Filt[[x]]) }),
                              names(flt))
            }
            # Also save the reference filters 
            nms <- names(Reg_filters[[tstrt]][[bee]])
            for (nm in nms) {
              if (tt == 4L) {
                Reg_filters[[tstrt]][[bee]][[nm]]$prot_Background_filter <- Ref.Filt[[nm]]
              } else {
                Reg_filters[[tstrt]][[bee]][[nm]]$Background_filter <- Ref.Filt[[nm]]
              }
            }
            #
            Mode <- "regulated"
            if (tt %in% c(1L, 3L)) { dataType <- "PG" }
            if (tt == 4L) { dataType <- "Prot" }
            if ((!exists("GO_mappings")) && file.exists(GO_mappings_fl)) {
              loadFun(GO_mappings_fl)
            }
            if ((!exists("GO_terms")) && file.exists(GO_terms_fl)) {
              loadFun(GO_terms_fl)
            }
            #
            Src <- paste0(libPath, "/extdata/Sources/GO_enrich.R")
            #rstudioapi::documentOpen(Src)
            source(Src)
            #
            if (runClueGO) {
              clueGO_outDir <- myDir
              clueGO_type <- "Enrichment (Right-sided hypergeometric test)"
              Src <- paste0(libPath, "/extdata/Sources/ClueGO_enrich.R")
              #rstudioapi::documentOpen(Src)
              source(Src)
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
              myDir2 <- paste0(wd, "/", myDir)
              for (ttl in n2) { #ttl <- n2[1L]
                plot <- GO_Plots[[tstbee]]$GO_plots[[ttl]]
              }
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
              #kl <- base::setdiff(colnames(temp), c(gn, kn, pv, zs, lf, si, pg, lp))
              kl <- base::setdiff(colnames(temp), c(gn, kn, pv, zs, lf, si, pg, pr))
              #temp <- temp[, c(kl, si, gn, pg, lp, kn, pv, zs, lf)]
              temp <- temp[, c(kl, si, gn, pg, pr, kn, pv, zs, lf)]
              w <- apply(temp[, pv, drop = FALSE], 1L, sum, na.rm = TRUE) > 0L
              temp <- temp[w,]
              tst <- apply(temp[, kn, drop = FALSE], 1L, sum, na.rm = TRUE)
              temp <- temp[order(tst, decreasing = TRUE),]
              temp <- temp[order(temp$Ontology, decreasing = FALSE),]
              Reg_GO_terms[[tstbee]] <- temp
              #temp <- Reg_GO_terms[[tstbee]]
              write.csv(temp, file = paste0(myDir, "/GO terms - ", tstbee, ".csv"), row.names = FALSE)
              w <- which(vapply(colnames(temp), \(x) { is.character(temp[[x]]) }, TRUE))
              if (length(w)) {
                for (i in w) { #i <- w[1L]
                  w1 <- which(nchar(temp[[colnames(temp)[i]]]) > xl_maxChar)
                  if (length(w1)) {
                    temp[[colnames(temp)[i]]][w1] <- paste0(substr(temp[[colnames(temp)[i]]][w1], 1L, xl_maxChar-3L), "...")
                  }
                }
              }
              #
              wb <- openxlsx2::wb_workbook()
              kount <- 0L
              setW <- \(wb, sheet, cols, width) {
                if (length(cols)) {
                  wb <- openxlsx2::wb_set_col_widths(wb, sheet, cols = cols, widths = width)
                }
                invisible(wb)
              }
              for (ont in names(Ontologies)) { #ont <- names(Ontologies)[1L]
                  wb <- openxlsx2::wb_add_worksheet(wb, ont)
                  wb <- openxlsx2::wb_add_data(wb, ont, temp[w, ])
                  wb <- openxlsx2::wb_set_row_heights(wb, ont, 1L, 60L)
                  wb <- openxlsx2::wb_set_col_widths(wb, ont, seq_len(ncol(temp)), 12L)
                  wb <- setW(wb, ont, match("Term", colnames(temp)), 45L)
                  wb <- setW(wb, ont, which(colnames(temp) %in% gn), 20L)
                  wb <- setW(wb, ont, which(colnames(temp) %in% pr), 20L)
                  wb <- setW(wb, ont, which(colnames(temp) %in% pp), 20L)
                  # header style
                  hdr <- openxlsx2::wb_dims(rows = 1L, cols = seq_len(ncol(temp)))
                  wb <- openxlsx2::wb_add_numfmt(wb, ont, dims = hdr, numfmt = "@")
                  wb <- openxlsx2::wb_add_font(wb, ont, dims = hdr, bold = TRUE, size = 12)
                  wb <- openxlsx2::wb_add_cell_style(wb, ont, dims = hdr, horizontal = "center",
                                                     vertical = "center", wrap_text = TRUE)
                  # number formats for data rows
                  rws <- 2L:(length(w) + 1L)
                  cI  <- which(colnames(temp) %in% kn)
                  cD  <- which(colnames(temp) %in% c(zs, lf, pv))
                  if (length(cI)) {
                    wb <- openxlsx2::wb_add_numfmt(wb, ont, dims = openxlsx2::wb_dims(rows = rws, cols = cI),
                                                   numfmt = "0")
                  }
                  if (length(cD)) {
                    wb <- openxlsx2::wb_add_numfmt(wb, ont, dims = openxlsx2::wb_dims(rows = rws, cols = cD),
                                                   numfmt = "0.000")
                  }
              }
              if (kount) {
                xl_fl <- paste0(myDir, "/GO terms - ", tstbee, ".xlsx")
                openxlsx2::wb_save(wb, xl_fl, overwrite = TRUE)
                #openxlsx2::xl_open(xl_fl)
                Kol2 <- grep(paste0("^Significance - .+", " ", max(BH.FDR)*100L, "%"),
                             colnames(Reg_GO_terms[[tstbee]]), value = TRUE)
                if (length(Kol2)) {
                  w <- if (length(Kol2) > 1L) {
                    which(apply(Reg_GO_terms[[tstbee]][,Kol2], 1L, \(x) {"+" %in% x}))
                  } else { which(vapply(Reg_GO_terms[[tstbee]][, Kol2], \(x) { "+" %in% x }, TRUE)) }
                  write.csv(Reg_GO_terms[[tstbee]][w,], file = paste0(myDir, "/Regulated GO terms - ", tstbee, ".csv"),
                            row.names = FALSE)
                }
                # Summary table and heatmap of number of regulated GO terms
                Kol3 <- grep(paste0("^Significance - [^ ]+ [1-9][0-9]*\\.*[0-9]*%$"),
                             colnames(GO_Plots[[tstbee]]$GO_terms), value = TRUE)
                tst <- as.numeric(gsub(paste0("^Significance - [^ ]+ |%$"), "", Kol3))
                Kol3 <- Kol3[tst == max(tst)]
                N <- length(Kol3)
                if (N > 1L) {
                  temp <- as.data.frame(matrix(rep("", (N+1L)^2L), ncol = N+1L))
                  W <- lapply(Kol3, \(x) { which(GO_Plots[[tstbee]]$GO_terms[[x]] == "+") })
                  names(W) <- gsub(paste0("^Significance - | ", max(BH.FDR)*100L, "%$"), "", Kol3)
                  temp[2L:(N+1L), 1L] <- temp[1L, 2L:(N+1L)] <- names(W)
                  for (i in 2L:(N+1L)) { #i <- 2
                    temp[i, 2L:(N+1L)] <- vapply(Kol3, \(x) {
                      sum((GO_Plots[[tstbee]]$GO_terms[[x]] == "+") & (GO_Plots[[tstbee]]$GO_terms[[Kol3[i]]] == "+"),
                          na.rm = TRUE)
                    }, 1L)
                  }
                  names(W) <- cleanNms(sub(" [0-9]+%$", "", names(W)))
                  tst <- lengths(strsplit(names(W), " - "))
                  tst <- (min(tst) > 1L) & (length(unique(tst)) == 1L)
                  if (tst) {
                    tst <- as.data.frame(t(sapply(strsplit(names(W), " - "), unlist)))
                    l <- apply(tst, 2L, \(x) { length(unique(x)) })
                    tst <- tst[, which(l > 1L)]
                    names(W) <- do.call(paste, c(tst, sep = " - "))
                  }
                  temp[2L:(N+1L), 1L] <- temp[1L, 2L:(N+1L)] <- names(W)
                  nm <- paste0("N. of co-regulated GO terms\n", tstrt, "\n(", tolower(bee), ")")
                  write.csv(temp, file = paste0(myDir, "/", gsub("\n", " - ", gsub("\n\\(", " (", nm)), ".csv"), row.names = FALSE)
                  temp2 <- temp[2L:(N+1L), 2L:(N+1L)]
                  colnames(temp2) <- temp[1L, 2L:(N+1L)]
                  rownames(temp2) <-  temp[2L:(N+1L), 1L]
                  for (i in 1L:nrow(temp2)) { temp2[[i]] <- as.numeric(temp2[[i]]) }
                  temp2ul <- unlist(temp2)
                  if (max(temp2ul[is.finite(temp2ul)]) > 0L) {
                    temp2 <- as.matrix(temp2)
                    basic.heatmap(temp2, "N. of co-regulated GO terms", paste0(tstrt, "\n(", tolower(bee), ")"),
                                  "svg", folder = myDir)
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
    cat("\n - Dataset\n")
    #
    myDir <- paste0(wd, "/Reg. analysis/GO enrich/Dataset")
    if (!dir.exists(myDir)) { dir.create(myDir, recursive = TRUE) }
    dirlist <- unique(c(dirlist, myDir))
    #
    Mode <- "dataset"
    dataType <- "PG"
    if ((!exists("GO_mappings")) && file.exists(GO_mappings_fl)) {
      loadFun(GO_mappings_fl)
    }
    if ((!exists("GO_terms")) && file.exists(GO_terms_fl)) {
      loadFun(GO_terms_fl)
    }
    #
    #try(rm(list = allArgs), silent = TRUE)
    #
    Src <- paste0(libPath, "/extdata/Sources/GO_enrich.R")
    #rstudioapi::documentOpen(Src)
    source(Src)
    #
    # Off because it was taking too long and failing sometimes
    # ClueGO really cannot handle too large gene lists, which are the norm here for dataset analysis
    #if (runClueGO) {
    #  clueGO_outDir <- myDir
    #  clueGO_type <- "Enrichment/Depletion (Two-sided hypergeometric test)"
    #  Src <- paste0(libPath, "/extdata/Sources/ClueGO_enrich.R")
    #  rstudioapi::documentOpen(Src)
    #  source(Src)
    #}
    #
    # Cleanup - do it now, not within sources!
    try(rm(list = allArgs), silent = TRUE)
    #
    # Quick Fisher exact test on GO terms of interest - looking only at observed data per group
    Src <- paste0(libPath, "/extdata/Sources/interestGO_Fisher.R")
    #rstudioapi::documentOpen(Src)
    source(Src)
    #
    if ((!is.null(goRES)) && ("All_GO_terms" %in% names(goRES))) {
      GO_terms <- goRES$All_GO_terms
    }
  }
  if (dataType == "PG") {
    DatAnalysisTxt$GO <- paste0("GO terms enrichment analysis was performed with topGO",
                                c("", " and ClueGO (Cytoscape)")[runClueGO+1L], ", comparing for each contrast ",
                                c("differentially expressed",
                                  "up-regulated")[(sum(myContrasts$`Up-only`) == nrow(myContrasts)) + 1L],
                                " against observed protein groups.")
  }
}
