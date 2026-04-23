#### Code chunk - Gene Ontology terms enrichment analysis
if (globalGO) {
  # Global dataset GO enrichment - expression per sample vs total proteome
  ref <- rev(PG.int.cols[which(PG.int.cols != paste0("Imput. ", PG.int.cols["Original"]))])[1L]
  xprsFilt <- setNames(lapply(Exp, \(e) { which(PG[[paste0(ref, e)]] > 0L) }), Exp)
  dir <- paste0(wd, "/GO enrichment analysis/Sample vs total proteome")
  if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
  setwd(dir)
  #temp <- PG
  #temp$`Leading protein IDs` <- gsub("CON__[^;]+", "", temp$`Leading protein IDs`) # Remove contaminants
  #
  Mode <- "dataset"
  dataType <- "PG"
  #
  Src <- paste0(libPath, "/extdata/R scripts/Sources/GO_enrich.R")
  #rstudioapi::documentOpen(Src)
  source(Src, local = FALSE)
  #
  if (runClueGO) {
    clueGO_outDir <- dir
    clueGO_type <- "Enrichment/Depletion (Two-sided hypergeometric test)"
    Src <- paste0(libPath, "/extdata/R scripts/Sources/ClueGO_enrich.R")
    #rstudioapi::documentOpen(Src)
    source(Src, local = FALSE)
  }
  #
  # Cleanup - do it now, not within sources!
  suppressWarnings(try(rm(list = allArgs), silent = TRUE))
  #
  temp <- goRES
  #
  setwd(wd)
  if ((is.list(goRES))&&("GO_terms" %in% names(goRES))) {
    temp <- goRES$GO_terms
    temp$"Protein group IDs" <- vapply(temp$`Protein table row(s)`, \(x) {
      paste(PG$id[as.numeric(x)], collapse = ";")
    }, "")
    temp$Mapping <- NULL
    if ("Offspring" %in% colnames(temp)) {
      temp$Offspring <- vapply(temp$Offspring, paste, "", collapse = ";")
    }
    temp$`Protein table row(s)` <- NULL
    gn <- grep("^Genes", colnames(temp), value = TRUE)
    pr <- grep("^Proteins", colnames(temp), value = TRUE)
    pg <- grep("^Protein group IDs", colnames(temp), value = TRUE)
    kn <- grep("^Count", colnames(temp), value = TRUE)
    pv <- grep("^Pvalue", colnames(temp), value = TRUE)
    zs <- grep("^(Z-score|\\(N_Up - N_Down\\)/sqrt\\(Tot\\.\\))", colnames(temp), value = TRUE)
    lf <- grep("^logFC", colnames(temp), value = TRUE)
    si <- grep("^Significance", colnames(temp), value = TRUE)
    kl <- colnames(temp)[which(!colnames(temp) %in% c(gn, kn, pv, zs, lf, si, pg, pr))]
    temp <- temp[, c(kl, si, gn, pg, pr, kn, pv, zs, lf)]
    w <- apply(temp[, pv, drop = FALSE], 1L, \(x) { sum(!is.na(x)) }) > 0L
    temp <- temp[w,]
    tst <- apply(temp[, kn, drop = FALSE], 1L, \(x) { sum(x[which(!is.na(x))]) })
    temp <- temp[order(tst, decreasing = TRUE),]
    temp <- temp[order(temp$Ontology, decreasing = FALSE),]
    w <- which(vapply(colnames(temp), \(x) { is.character(temp[[x]]) }, TRUE))
    if (length(w)) {
      for (i in w) { #i <- w[1L]
        w1 <- which(nchar(temp[[colnames(temp)[i]]]) > ExcelMax)
        if (length(w1)) {
          temp[[colnames(temp)[i]]][w1] <- paste0(substr(temp[[colnames(temp)[i]]][w1], 1L, ExcelMax-3L), "...")
        }
      }
    }
    wb <- wb_workbook()
    myFont <- create_font(sz = "12", color = wb_color("#FFFFFF"))
    wb$styles_mgr$add(myFont, "Header_font")
    HdrStl <- create_cell_style(font_id = 1L,
                                num_fmt_id = "General",
                                horizontal = "justify",
                                vertical = "top",
                                text_rotation = 60,
                                wrap_text = TRUE)
    wb$styles_mgr$add(HdrStl, "Header_style")
    nTabs <- 0L
    for (ont in names(Ontologies)) { #ont <- names(Ontologies)[1L]
      w <- which(temp$Ontology == Ontologies[ont])
      if (length(w)) {
        sheetNm <- ont
        if (sheetNm %in% wb_get_sheet_names(wb)) {
          wb <- wb_remove_worksheet(wb, sheetNm)
          nTabs <- nTabs + 1L
        }
        wb <- wb_add_worksheet(wb, sheetNm)
        nTabs <- nTabs + 1L
        wb <- wb_add_data_table(wb, sheetNm, temp[w,], wb_dims(2L, 2L))
        wb <- wb_set_row_heights(wb, sheetNm, 2L, 120L)
        wb <- wb_set_col_widths(wb, sheetNm, 1L, 1L)
        wb <- wb_set_col_widths(wb, sheetNm, 1L+1L:ncol(temp), 12L)
        wb <- wb_set_col_widths(wb, sheetNm, 1L+which(colnames(temp) == "Term"), 45L)
        wb <- wb_set_col_widths(wb, sheetNm, 1L+which(colnames(temp) %in% c(gn, pr)), 20L)
        dms <- wb_dims(2L, 1L+1L:ncol(temp))
        wb <- wb_set_cell_style(wb, sheetNm, dms, wb$styles_mgr$get_xf_id("Header_style"))
        dms <- wb_dims(1L+seq_along(w), 1L+which(colnames(temp) %in% kn))
        wb <- wb_add_numfmt(wb, sheetNm, dms, "0")
        dms <- wb_dims(1L+seq_along(w), 1L+which(colnames(temp) %in% c(zs, lf, pv)))
        wb <- wb_add_numfmt(wb, sheetNm, dms, "0.000")
      }
    }
    if (nTabs) {
      goXLfl <- paste0(dir, "/GO terms - Samples vs Proteome.xlsx")
      wb_save(wb, goXLfl)
      #xl_open(goXLfl)
    }
  }
  if (enrichGO) {
    for (grp in rat.grps) { #grp <- rat.grps[1L]
      wb <- wb_workbook()
      myFont <- create_font(sz = "12", color = wb_color("#FFFFFF"))
      wb$styles_mgr$add(myFont, "Header_font")
      HdrStl <- create_cell_style(font_id = 1L,
                                  num_fmt_id = "General",
                                  horizontal = "justify",
                                  vertical = "top",
                                  text_rotation = 60,
                                  wrap_text = TRUE)
      wb$styles_mgr$add(HdrStl, "Header_style")
      nTabs <- 0L
      SmplMp <- SamplesMap[which(SamplesMap$`Ratios group` == grp),]
      tmp <- paste0("Comparisons to ", SmplMp$Experiment[which(SmplMp$Reference)])
      if (length(rat.grps) > 1L) { tmp <- paste0("Group", grp, " (", tmp, ")") }
      dir <- paste0(wd, "/GO enrichment analysis/", tmp)
      if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
      setwd(dir)
      db_obs %<o% db[which(db$"Protein ID" %in% unlist(strsplit(PG$"Leading protein IDs", ";"))),]
      fcFilt <- FC_filt[which(names(FC_filt) %in% SmplMp$Experiment)]
      nms <- names(fcFilt)
      filt2 <- setNames(lapply(nms, \(x) {
        unique(unlist(xprsFilt[SmplMp$Experiment]))
      }), nms)
      #
      Mode <- "regulated"
      dataType <- "PG"
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
      suppressWarnings(try(rm(list = allArgs), silent = TRUE))
      #
      setwd(wd)
      if ((is.list(goRES))&&("GO_terms" %in% names(goRES))) {
        temp <- goRES$GO_terms
        temp$"Protein group IDs" <- vapply(temp$`Protein table row(s)`, \(x) {
          paste(PG$id[as.numeric(x)], collapse = ";")
        }, "")
        temp$Mapping <- NULL
        if ("Offspring" %in% colnames(temp)) {
          temp$Offspring <- vapply(temp$Offspring, paste, "", collapse = ";")
        }
        temp$`Protein table row(s)` <- NULL
        gn <- grep("^Genes", colnames(temp), value = TRUE)
        pr <- grep("^Proteins", colnames(temp), value = TRUE)
        pg <- grep("^Protein group IDs", colnames(temp), value = TRUE)
        kn <- grep("^Count", colnames(temp), value = TRUE)
        pv <- grep("^Pvalue", colnames(temp), value = TRUE)
        zs <- grep("^(Z-score|\\(N_Up - N_Down\\)/sqrt\\(Tot\\.\\))", colnames(temp), value = TRUE)
        lf <- grep("^logFC", colnames(temp), value = TRUE)
        si <- grep("^Significance", colnames(temp), value = TRUE)
        kl <- colnames(temp)[which(!colnames(temp) %in% c(gn, kn, pv, zs, lf, si, pg, pr))]
        temp <- temp[, c(kl, si, gn, pg, pr, kn, pv, zs, lf)]
        w <- apply(temp[, pv, drop = FALSE], 1L, \(x) { sum(!is.na(x)) }) > 0L
        temp <- temp[w,]
        tst <- apply(temp[, kn, drop = FALSE], 1L, \(x) { sum(x[which(!is.na(x))]) })
        temp <- temp[order(tst, decreasing = TRUE),]
        temp <- temp[order(temp$Ontology, decreasing = FALSE),]
        w <- which(vapply(colnames(temp), \(x) { is.character(temp[[x]]) }, TRUE))
        if (length(w)) {
          for (i in w) { #i <- w[1L]
            w1 <- which(nchar(temp[[colnames(temp)[i]]]) > ExcelMax)
            if (length(w1)) {
              temp[[colnames(temp)[i]]][w1] <- paste0(substr(temp[[colnames(temp)[i]]][w1], 1L, ExcelMax-3L), "...")
            }
          }
        }
        for (ont in names(Ontologies)) { #ont <- names(Ontologies)[1L]
          w <- which(temp$Ontology == Ontologies[ont])
          if (length(w)) {
            sheetNm <- ont
            if (sheetNm %in% wb_get_sheet_names(wb)) {
              wb <- wb_remove_worksheet(wb, sheetNm)
              nTabs <- nTabs-1L
            }
            wb <- wb_add_worksheet(wb, sheetNm)
            nTabs <- nTabs + 1L
            wb <- wb_add_data_table(wb, sheetNm, temp[w,], wb_dims(2L, 2L))
            wb <- wb_set_row_heights(wb, sheetNm, 2L, 120L)
            wb <- wb_set_col_widths(wb, sheetNm, 1L, 1L)
            wb <- wb_set_col_widths(wb, sheetNm, 1L+1L:ncol(temp), 12L)
            wb <- wb_set_col_widths(wb, sheetNm, 1L+which(colnames(temp) == "Term"), 45L)
            wb <- wb_set_col_widths(wb, sheetNm, 1L+which(colnames(temp) %in% c(gn, pr)), 20L)
            dms <- wb_dims(2L, 1L+1L:ncol(temp))
            wb <- wb_set_cell_style(wb, sheetNm, dms, wb$styles_mgr$get_xf_id("Header_style"))
            dms <- wb_dims(1L+seq_along(w), 1L+which(colnames(temp) %in% kn))
            wb <- wb_add_numfmt(wb, sheetNm, dms, "0")
            dms <- wb_dims(1L+seq_along(w), 1L+which(colnames(temp) %in% c(zs, lf, pv)))
            wb <- wb_add_numfmt(wb, sheetNm, dms, "0.000")
          }
        }
        if (nTabs) {
          goXLfl <- paste0(dir, "/GO terms - Comp group ", grp, ".xlsx")
          wb_save(wb, goXLfl)
          #xl_open(goXLfl)
        }
      }
    }
    if (PTMriched) {
      for (ptm in EnrichedPTMs) { #ptm <- EnrichedPTMs[1L]
        ptmpep <- PTMs_pep[[ptm]]
        for (grp in rat.grps) { #grp <- rat.grps[1L]
          SmplMp <- SamplesMap[which(SamplesMap$`Ratios group` == grp),]
          tmp <- paste0("Comparisons to ", SmplMp$Experiment[which(SmplMp$Reference)])
          if (length(rat.grps) > 1L) { tmp <- paste0("Group", grp, " (", tmp, ")") }
          nms <- names(PTMs_FC_filt[[ptm]])
          filt3 <- setNames(lapply(nms, \(x) {
            # Filter: any expressed in the same group, i.e. at least with 1 non null value in the group
            which(rowSums(ptmpep[, paste0(PTMs_intRf[length(PTMs_intRf) # Not log!
            ], " - ", SmplMp$Experiment)], na.rm = TRUE) > 0L)
          }), nms)
          dir <- paste0(wd, "/GO enrichment analysis/", tmp)
          if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
          setwd(dir)
          #
          Mode <- "regulated"
          dataType <- "modPeptides"
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
          suppressWarnings(try(rm(list = allArgs), silent = TRUE))
          #
          setwd(wd)
          if ((is.list(goRES))&&("GO_terms" %in% names(goRES))) {
            temp <- goRES$GO_terms
            temp$"Protein group IDs" <- vapply(temp$`Protein table row(s)`, \(x) {
              paste(PG$id[as.numeric(x)], collapse = ";")
            }, "")
            temp$Mapping <- NULL
            if ("Offspring" %in% colnames(temp)) {
              temp$Offspring <- vapply(temp$Offspring, paste, "", collapse = ";")
            }
            temp$`Protein table row(s)` <- NULL
            gn <- grep("^Genes", colnames(temp), value = TRUE)
            pr <- grep("^Proteins", colnames(temp), value = TRUE)
            pg <- grep("^Protein group IDs", colnames(temp), value = TRUE)
            kn <- grep("^Count", colnames(temp), value = TRUE)
            pv <- grep("^Pvalue", colnames(temp), value = TRUE)
            zs <- grep("^(Z-score|\\(N_Up - N_Down\\)/sqrt\\(Tot\\.\\))", colnames(temp), value = TRUE)
            lf <- grep("^logFC", colnames(temp), value = TRUE)
            si <- grep("^Significance", colnames(temp), value = TRUE)
            kl <- colnames(temp)[which(!colnames(temp) %in% c(gn, kn, pv, zs, lf, si, pg, pr))]
            temp <- temp[, c(kl, si, gn, pg, pr, kn, pv, zs, lf)]
            w <- apply(temp[, pv, drop = FALSE], 1L, \(x) { sum(!is.na(x)) }) > 0L
            temp <- temp[w,]
            tst <- apply(temp[, kn, drop = FALSE], 1L, \(x) { sum(x[which(!is.na(x))]) })
            temp <- temp[order(tst, decreasing = TRUE),]
            temp <- temp[order(temp$Ontology, decreasing = FALSE),]
            w <- which(vapply(colnames(temp), \(x) { is.character(temp[[x]]) }, TRUE))
            if (length(w)) {
              for (i in w) { #i <- w[1L]
                w1 <- which(nchar(temp[[colnames(temp)[i]]]) > ExcelMax)
                if (length(w1)) {
                  temp[[colnames(temp)[i]]][w1] <- paste0(substr(temp[[colnames(temp)[i]]][w1], 1L, ExcelMax-3L), "...")
                }
              }
            }
            for (ont in names(Ontologies)) { #ont <- names(Ontologies)[1L]
              w <- which(temp$Ontology == Ontologies[ont])
              if (length(w)) {
                sheetNm <- ont
                if (sheetNm %in% wb_get_sheet_names(wb)) {
                  wb <- wb_remove_worksheet(wb, sheetNm)
                  nTabs <- nTabs-1L
                }
                wb <- wb_add_worksheet(wb, sheetNm)
                nTabs <- nTabs + 1L
                wb <- wb_add_data_table(wb, sheetNm, temp[w,], wb_dims(2L, 2L))
                wb <- wb_set_row_heights(wb, sheetNm, 2L, 120L)
                wb <- wb_set_col_widths(wb, sheetNm, 1L, 1L)
                wb <- wb_set_col_widths(wb, sheetNm, 1L+1L:ncol(temp), 12L)
                wb <- wb_set_col_widths(wb, sheetNm, 1L+which(colnames(temp) == "Term"), 45L)
                wb <- wb_set_col_widths(wb, sheetNm, 1L+which(colnames(temp) %in% c(gn, pr)), 20L)
                dms <- wb_dims(2L, 1L+1L:ncol(temp))
                wb <- wb_set_cell_style(wb, sheetNm, dms, wb$styles_mgr$get_xf_id("Header_style"))
                dms <- wb_dims(1L+seq_along(w), 1L+which(colnames(temp) %in% kn))
                wb <- wb_add_numfmt(wb, sheetNm, dms, "0")
                dms <- wb_dims(1L+seq_along(w), 1L+which(colnames(temp) %in% c(zs, lf, pv)))
                wb <- wb_add_numfmt(wb, sheetNm, dms, "0.000")
              }
            }
            if (nTabs) {
              goXLfl <- paste0(dir, "/GO terms - ", ptm, " - Comp group ", grp, ".xlsx")
              wb_save(wb, goXLfl)
              #xl_open(goXLfl)
            }
          }
        }
      }
    }
  }
}
