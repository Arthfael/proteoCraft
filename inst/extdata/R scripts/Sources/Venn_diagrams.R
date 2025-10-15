#### Venn diagrams
# cleanNms2 function specifically designed to clean names for Venn diagrams
library(venn)
library(ggplot2)
library(openxlsx2)
cleanNms2 %<o% function(names, groups = VPAL, sep = "\n/vs/\n", simplify) {
  isList <- ("list" %in% class(names))
  if (missing(simplify)) { simplify <- !isList }
  if (!isList) {
    names <- sapply(strsplit(gsub("^\\(|\\)$", "", names), "\\) - \\("), unlist)
  }
  uNms <- unique(unlist(names))
  nuNms <- as.data.frame(t(as.data.frame(strsplit(uNms, "___"))))
  colnames(nuNms) <- groups$names
  nuNms$Full <- uNms
  w <- which(vapply(groups$names, function(x) { length(unique(nuNms[[x]])) }, 1) > 1)
  nuNms$New <- do.call(paste, c(nuNms[, groups$names[w], drop = FALSE], sep = ""))
  nuNames <- lapply(names, function(x) { nuNms$New[match(x, nuNms$Full) ]})
  if (simplify) { nuNames <- vapply(nuNames, paste, "", collapse = sep) }
  return(nuNames)
}
#
ReportCalls <- AddSpace2Report()
msg <- "Venn diagrams"
ReportCalls <- AddMsg2Report()
HdrStlVenn <- openxlsx2::create_cell_style(num_fmt_id = "General",
                                           horizontal = "left",
                                           vertical = "bottom",
                                           text_rotation = 60,
                                           wrap_text = TRUE)
Mod2Venn <- c()
scrptMtch <- match(scrptTypeFull,
                   c("withReps_PTMs_only",
                     "withReps_PG_and_PTMs"))
if (("Venn.Groups" %in% colnames(Param))&&(Param$Venn.Groups != "")) {
  VennGrp <- Param$Venn.Groups
} else { VennGrp <- Param$Ratios.Groups }
if (toupper(VennGrp) == "GLOBAL") { VennGrp <- "GLOBAL" }
#if (grepl("^[A-Z][a-z]{2}(;[A-Z][a-z]{2})*$", VennGrp)) {
VennGrp2 <- parse.Param.aggreg(Param_filter(VennGrp, "Rep"))
#}
VennMx <- 7
#
II <- setNames(1,
               c("All peptidoforms", "Protein groups")[scrptMtch])
if ((exists("PTMs_pep"))&&(length(PTMs_pep))) {
  Mod2Venn <- names(PTMs_pep)
  II[paste0(Mod2Venn, "-mod. pept.")] <- 1+(seq_along(length(Mod2Venn)))
}
if (scrptMtch == 1) { II <- II[1:length(II)] } # For now
for (ii in II) { #ii <- II[1] #ii <- II[2] #ii <- II[3]
  if ((scrptMtch == 2)&&(ii == 1)) { # For now, we may extend to withReps_PTMs_only later
    dir <- paste0(wd, "/Venn diagrams")
    ttest_Filt <- Reg_filters$"t-tests"$"By condition"
    vennRoot <- ""
    myData <- PG
    myRef <- paste0("Mean ", prtRfRoot)
    infoKol <- "Genes"
    idKol <- "Leading protein IDs"
    # Draft alt. code for pep
    #myData <- pep
    #myRef <- grep("^Mean [^ ]+ log10\\(", colnames(myData), value = TRUE)
    #myRef <- paste0(unique(gsub(" - .*", "", myRef)), " - ")
    #myRef <- myRef[length(myRef)]
    #idKol <- "Modified sequence"
    #infoKol <- "Proteins"
    if (F.test) {
      Ftest_Filt <- Reg_filters$"F-tests"$"By condition"
      myFData <- F_test_data
    }
  } else {
    Ptm <- Mod2Venn[ii-1]
    dir <- paste0(wd, "/Reg. analysis/", Ptm, "/Venn diagrams")
    ttest_Filt <- PTMs_Reg_filters[[Ptm]]$"t-tests"$"By condition"
    vennRoot <- paste0(Ptm, " ")
    myData <- PTMs_pep[[Ptm]]
    myRef <- grep("^Mean ([^ ]+ )?log10\\(", colnames(myData), value = TRUE)
    myRef <- paste0(unique(gsub(" - .*", "", myRef)), " - ")
    myRef <- myRef[length(myRef)]
    infoKol <- "Proteins"
    idKol <- "Modified sequence"
    if (F.test) {
      Ftest_Filt <- PTMs_Reg_filters[[Ptm]]$"F-tests"$"By condition"
      myFData <- PTMs_F_test_data[[Ptm]]
    }
  }
  topTitle <- paste0("Venn diagram - ", names(II)[ii])
  infoKol <- infoKol[which(infoKol %in% colnames(myData))]
  if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
  dirlist <- unique(c(dirlist, dir))
  #
  setwd(dir)
  #View(myData[, paste0(myRef, VPAL$values)])
  comp_list <- setNames(lapply(VPAL$values, function(grp) { #grp <- VPAL$values[1]
    x <- myData[[paste0(myRef, grp)]]
    which(is.all.good(as.numeric(x), 2))
  }), cleanNms2(VPAL$values))
  w <- which(vapply(comp_list, length, 1) > 0)
  VennExp <- names(comp_list)[w]
  lV <- length(VennExp)
  OK <- lV > 1
  ttl <- paste0(vennRoot, "LFQ Venn diagram, all")
  msg <- paste0(" - ", names(II)[ii], " - observations per sample groups")
  ReportCalls <- AddMsg2Report(Space = FALSE)
  if (lV > VennMx) {
    msg <- paste0("Too many groups, select at least 2 and up to ", VennMx,
                  " to include in Venn diagram ", ttl, ", or cancel to skip.")
    opt <- vapply(VennExp, function(x) { paste(c(x, rep(" ", 200-nchar(x))), collapse = "") }, "")
    VennExp <- VennExp[match(dlg_list(opt, opt[1:VennMx], TRUE, msg)$res, opt)]
    if (length(VennExp) <= 1) {
      warning("Skipping per-sample-group observations Venn diagrams: you should have selected at least 2 groups!")
      OK <- FALSE
    }
    if ((length(VennExp) > VennMx)&&(length(VennExp) < 1)) {
      msg <- paste0("Skipping per-sample-group observations Venn diagrams: you should have selected ",
                    VennMx, " groups at most!")
      warning(msg)
      OK <- FALSE
    }
  }
  if (OK) {
    comp_list <- comp_list[VennExp]
    plot <- venn::venn(comp_list, ilabels = "counts", ellipse = TRUE, zcolor = "style", ggplot = TRUE)
    plot <- plot + ggtitle(topTitle, subtitle = "Global, LFQ") +
      theme(plot.title = element_text(size = 15),
            plot.subtitle = element_text(size = 10))
    poplot(plot)
    suppressMessages({
      ggsave(paste0(dir, "/", ttl, ".jpg"), plot, dpi = 150)
      ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 150)
    })
    #system(paste0("open \"", dir, "/", ttl, ".jpg", "\""))
    ReportCalls$Calls <- append(ReportCalls$Calls, paste0("body_add_img(Report, \"", dir, "/", ttl, ".jpg\", height = 6, width = 6)"))
    #
    wb <- openxlsx2::wb_workbook()
    wb$styles_mgr$add(HdrStlVenn, "Header_style")
    SheetNm <- "Sample groups composition"
    if (SheetNm %in% wb_get_sheet_names(wb)) { wb <- openxlsx2::wb_remove_worksheet(wb, SheetNm) }
    wb <- openxlsx2::wb_add_worksheet(wb, SheetNm)
    l <- length(comp_list)
    tmp <- sapply(names(comp_list), function(grp) {
      res <- rep("", nrow(myData))
      res[comp_list[[grp]]] <- "+"
      return(res)
    })
    tmp <- as.data.frame(tmp)
    tmp[[idKol]] <- myData[[idKol]]
    tmp <- tmp[, c(idKol, names(comp_list))]
    wb <- openxlsx2::wb_add_data(wb, SheetNm, tmp, openxlsx2::wb_dims(2, 2))
    wb <- openxlsx2::wb_set_row_heights(wb, SheetNm, 2, 120)
    wb <- openxlsx2::wb_set_col_widths(wb, SheetNm, 1, 3)
    wb <- openxlsx2::wb_set_col_widths(wb, SheetNm, 1, 10)
    hdrDms <- openxlsx2::wb_dims(2, 1:(l+1) + 1)
    wb <- openxlsx2::wb_add_font(wb, SheetNm, hdrDms, size = 12, bold = TRUE)
    wb <- openxlsx2::wb_set_cell_style(wb, SheetNm, hdrDms, wb$styles_mgr$get_xf_id("Header_style"))
    vennXLfl <- paste0(dir, "/Venn diagrams.xlsx")
    openxlsx2::wb_save(wb, vennXLfl)
    #openxlsx2::xl_open(vennXLfl)
  } else {
    msg <- paste0("     ", ttl, ": not enough groups to compare!")
    ReportCalls <- AddMsg2Report(Space = FALSE)
  }
  ReportCalls <- AddSpace2Report()
  #
  dir2 <- paste0(dir, "/t-tests")
  if (!dir.exists(dir2)) { dir.create(dir2, recursive = TRUE) }
  dirlist <- unique(c(dirlist, dir2))
  setwd(dir2)
  wb <- openxlsx2::wb_workbook()
  wb$styles_mgr$add(HdrStlVenn, "Header_style")
  wbKount <- 0
  vennlev <- c("up", "down")[1:(1+TwoSided)]
  flt_nmz <- names(ttest_Filt)
  msg <- " - t-tests"
  ReportCalls <- AddMsg2Report(Space = FALSE)
  if (length(flt_nmz) > 1) {
    for (r in vennlev) { #r <- "up"
      m <- (length(vennlev) > 1)+1
      levTxt0 <- c("", r)[m]
      levTxt1 <- c("", paste0(" ", r))[m]
      levTxt2 <- c("", paste0(", ", r))[m]
      ReportCalls <- AddMsg2Report(Msg = paste0("   -> ", levTxt0), Offset = TRUE, Space = FALSE)
      for (grp in VennGrp2$values) { #grp <- VennGrp2$values[1]
        n <- (length(VennGrp2$values) > 1)+1
        grpTxt <- cleanNms(grp)
        grpTxt0 <- c("", grpTxt)[m]
        grpTxt1 <- c("", paste0(" ", grpTxt))[n]
        grpTxt2 <- c("", paste0(", ", grpTxt))[n]
        if (n == 2) {
          ReportCalls <- AddMsg2Report(Msg = paste0("     + ", grpTxt), Offset = TRUE, Space = FALSE)
        }
        em <- Exp.map[which(Exp.map[[VennGrp2$column]] == grp),]
        nmz <- flt_nmz[which(flt_nmz %in% em[[VPAL$column]])]
        comp_list <- setNames(lapply(nmz, function(x) {
          ttest_Filt[[x]][[paste0("Filter_", r)]]
        }), cleanNms2(nmz))
        w <- which(vapply(comp_list, length, 1) > 0)
        VennExp <- names(comp_list)[w]
        lV <- length(VennExp)
        OK <- lV > 1
        ttl <- gsub("^g", "G", paste0(vennRoot, "t-tests Venn diagram", grpTxt2, levTxt2))
        if (lV > VennMx) {
          msg <- paste0("Too many groups, select at least 2 and up to ", VennMx,
                        " to include in Venn diagram ", ttl)
          opt <- vapply(VennExp, function(x) { paste(c(x, rep(" ", 200-nchar(x))), collapse = "") }, "")
          VennExp <- VennExp[match(dlg_list(opt, opt[1:VennMx], TRUE, msg)$res, opt)]
          if (length(VennExp) <= 1) {
            warning("       ", ttl,
                    " - skipping: you should have selected at least 2 groups!")
            OK <- FALSE
          }
          if ((length(VennExp) > VennMx)&&(length(VennExp) < 1)) {
            msg <- paste0("       ", ttl, " - skipping: you should have selected ",
                          VennMx, " groups at most!")
            warning(msg)
            OK <- FALSE
          }
        }
        if (OK) {
          comp_list <- comp_list[VennExp]
          kr <- paste0("Regulated - ", nmz)
          tmp <- myData[unique(unlist(comp_list)), c(idKol, infoKol, "id", kr)]
          if (r == "up") { good <- grep("^up|^Specific", unique(unlist(tmp[, kr])), value = TRUE) }
          if (r == "down") { good <- grep("^up|^Specific", unique(unlist(tmp[, kr])), value = TRUE) }
          tst <- apply(tmp[, kr], 1, function(x) { length(w[which(x %in% good)]) })
          tmp <- tmp[order(tst, decreasing = TRUE),]
          write.csv(tmp, paste0(dir2, "/", ttl, " - table.csv"), row.names = FALSE)
          plot <- venn(comp_list, ilabels = "counts", ellipse = TRUE, zcolor = "style", ggplot = TRUE)
          plot <- plot +
            ggtitle(topTitle, subtitle = paste0("t-tests, ", r)) +
            theme(plot.title = element_text(size = 15),
                  plot.subtitle = element_text(size = 10))
          poplot(plot)
          suppressMessages({
            ggsave(paste0(dir2, "/", ttl, ".jpg"), plot, dpi = 150)
            ggsave(paste0(dir2, "/", ttl, ".pdf"), plot, dpi = 150)
          })
          #system(paste0("open \"", dir2, "/", ttl, ".jpg", "\""))
          ReportCalls$Calls <- append(ReportCalls$Calls, paste0("body_add_img(Report, \"", dir2, "/", ttl,
                                                                ".jpg\", height = 6, width = 6)"))
          #
          wbKount <- wbKount+1
          SheetNm <- paste0("t-test", grpTxt1, levTxt1)
          if (SheetNm %in% wb_get_sheet_names(wb)) { wb <- openxlsx2::wb_remove_worksheet(wb, SheetNm) }
          wb <- openxlsx2::wb_add_worksheet(wb, SheetNm)
          l <- length(comp_list)
          tmp <- sapply(names(comp_list), function(grp) {
            res <- rep("", nrow(myData))
            res[comp_list[[grp]]] <- "+"
            return(res)
          })
          tmp <- as.data.frame(tmp)
          tmp[[idKol]] <- myData[[idKol]]
          tmp <- tmp[, c(idKol, names(comp_list))]
          wb <- openxlsx2::wb_add_data(wb, SheetNm, tmp, openxlsx2::wb_dims(2, 2))
          wb <- openxlsx2::wb_set_row_heights(wb, SheetNm, 2, 120)
          wb <- openxlsx2::wb_set_col_widths(wb, SheetNm, 1, 3)
          wb <- openxlsx2::wb_set_col_widths(wb, SheetNm, 1, 10)
          hdrDms <- openxlsx2::wb_dims(2, 1:(l+1) + 1)
          wb <- openxlsx2::wb_add_font(wb, SheetNm, hdrDms, size = 12, bold = TRUE)
          wb <- openxlsx2::wb_set_cell_style(wb, SheetNm, hdrDms, wb$styles_mgr$get_xf_id("Header_style"))
        } else {
          msg <- paste0("       ", ttl, ": not enough groups with regulated proteins to compare!")
          ReportCalls <- AddMsg2Report(Space = FALSE)
        }
      }
    }
  }
  if (wbKount) {
    vennXLfl <- paste0(dir2, "/Venn diagrams.xlsx")
    openxlsx2::wb_save(wb, vennXLfl)
    #openxlsx2::xl_open(vennXLfl)
  }
  ReportCalls <- AddSpace2Report()
  #
  if (F.test) {
    vennlev <- c("up", "down") # For F-test we always have both 
    dir2 <- paste0(dir, "/F-test")
    if (!dir.exists(dir2)) { dir.create(dir2, recursive = TRUE) }
    dirlist <- unique(c(dirlist, dir2))
    setwd(dir2)
    wb <- openxlsx2::wb_workbook()
    wb$styles_mgr$add(HdrStlVenn, "Header_style")
    wbKount <- 0
    #
    nmz <- names(Ftest_Filt)
    nmz <- nmz[which(!grepl(" VS ", nmz))]
    OK <- length(nmz) > 1
    if (length(nmz) > VennMx) {
      msg <- paste0("Too many groups, select at least 2 and up to ", VennMx,
                    " to include in Venn diagram ", ttl, ", or cancel to skip.")
      opt <- vapply(nmz, function(x) { paste(c(x, rep(" ", 200-nchar(x))), collapse = "") }, "")
      nmz <- nmz[match(dlg_list(opt, opt[1:VennMx], TRUE, msg)$res, opt)]
      if (length(nmz) <= 1) {
        warning("   Skipping F-test Venn diagrams: you should have selected at least 2 groups!")
        OK <- FALSE
      }
      if ((length(nmz) > VennMx)&&(length(nmz) < 1)) {
        msg <- paste0("   Skipping F-test Venn diagrams: you should have selected ",
                      VennMx, " groups at most!")
        warning(msg)
        OK <- FALSE
      }
    }
    if (OK) {
      msg <- " - F-tests"
      ReportCalls <- AddMsg2Report(Space = FALSE)
      for (r in vennlev) { #r <- "up"
        m <- (length(vennlev) > 1)+1
        levTxt0 <- c("", r)[m]
        levTxt1 <- c("", paste0(" ", r))[m]
        levTxt2 <- c("", paste0(", ", r))[m]
        ReportCalls <- AddMsg2Report(Msg = paste0("   -> ", levTxt0), Offset = TRUE, Space = FALSE)
        ttl <- paste0(vennRoot, "Venn diagram - F-test, global -", levTxt1)
        nmz2 <- cleanNms2(nmz)
        comp_list <- setNames(lapply(nmz, function(x) {
          Ftest_Filt[[x]][[paste0("Filter_", r)]]
        }), nmz2)
        tst <- vapply(comp_list, length, 1)
        w <- which(tst > 0)
        if (length(w) > 1) {
          if (length(w) <= VennMx) {
            comp_list <- comp_list[w]
            w <- unique(unlist(comp_list))
            wN <- which(!infoKol %in% colnames(myFData))
            if (length(wN)) {
              myFData[, infoKol[wN]] <- myData[match(myFData[[idKol]], myData[[idKol]]), infoKol]
            }
            tmp <- myFData[w, c(idKol, infoKol)]
            tmp$id <- myData$id[w]
            kr <- paste0("mod. F-test Regulated - ", nmz)
            tmp[, kr] <- myFData[w, kr]
            if (r == "up") { good <- grep("^up|^Specific", unique(unlist(tmp[, kr])), value = TRUE) }
            if (r == "down") { good <- grep("^up|^Specific", unique(unlist(tmp[, kr])), value = TRUE) }
            tst <- apply(tmp[, kr], 1, function(x) { length(w[which(x %in% good)]) })
            tmp <- tmp[order(tst, decreasing = TRUE),]
            write.csv(tmp, paste0(dir2, "/", ttl, " - table.csv"), row.names = FALSE)
            names(comp_list) <- cleanNms(names(comp_list))
            plot <- venn(comp_list, ilabels = "counts", ellipse = TRUE, zcolor = "style", ggplot = TRUE)
            plot <- plot +
              ggtitle(topTitle, subtitle = paste0(vennRoot, "F-test", levTxt2)) +
              theme(plot.title = element_text(size = 15),
                    plot.subtitle = element_text(size = 10))
            poplot(plot)
            suppressMessages({
              ggsave(paste0(dir2, "/", ttl, ".jpg"), plot, dpi = 150)
              ggsave(paste0(dir2, "/", ttl, ".pdf"), plot, dpi = 150)
            })
            #system(paste0("open \"", dir2, "/", ttl, ".jpg", "\""))
            ReportCalls$Calls <- append(ReportCalls$Calls, paste0("body_add_img(Report, \"", dir2, "/",
                                                                  ttl, ".jpg\", height = 6, width = 6)"))
            wbKount <- wbKount+1
            SheetNm <- paste0("F-test", levTxt1)
            if (SheetNm %in% wb_get_sheet_names(wb)) { wb <- openxlsx2::wb_remove_worksheet(wb, SheetNm) }
            wb <- openxlsx2::wb_add_worksheet(wb, SheetNm)
            l <- length(comp_list)
            tmp <- sapply(names(comp_list), function(grp) {
              res <- rep("", nrow(myData))
              res[comp_list[[grp]]] <- "+"
              return(res)
            })
            tmp <- as.data.frame(tmp)
            tmp[[idKol]] <- myData[[idKol]]
            tmp <- tmp[, c(idKol, names(comp_list))]
            wb <- openxlsx2::wb_add_data(wb, SheetNm, tmp, openxlsx2::wb_dims(2, 2))
            wb <- openxlsx2::wb_set_row_heights(wb, SheetNm, 2, 120)
            wb <- openxlsx2::wb_set_col_widths(wb, SheetNm, 1, 3)
            wb <- openxlsx2::wb_set_col_widths(wb, SheetNm, 1, 10)
            hdrDms <- openxlsx2::wb_dims(2, 1:(l+1) + 1)
            wb <- openxlsx2::wb_add_font(wb, SheetNm, hdrDms, size = 12, bold = TRUE)
            wb <- openxlsx2::wb_set_cell_style(wb, SheetNm, hdrDms, wb$styles_mgr$get_xf_id("Header_style"))
          } else {
            msg <- paste0("     Could not draw global F-tests Venn diagram: more than ", VennMx,
                          " groups to compare!")
            ReportCalls <- AddMsg2Report(Space = FALSE)
          }
          ReportCalls <- AddSpace2Report()
        } else {
          msg <- paste0("   Could not draw global F-tests Venn diagram", levTxt1,
                        ": not enough groups with regulated proteins to compare!")
          ReportCalls <- AddMsg2Report()
        }
      }
    } else {
      msg <- "   Could not draw global F-tests Venn diagram: not enough groups to compare!"
      ReportCalls <- AddMsg2Report()
    }
    ReportCalls <- AddSpace2Report()
    if (wbKount) {
      vennXLfl <- paste0(dir2, "/Venn diagrams.xlsx")
      openxlsx2::wb_save(wb, vennXLfl)
      #openxlsx2::xl_open(vennXLfl)
    }
    #
  }
  setwd(wd)
  ReportCalls <- AddSpace2Report()
}
