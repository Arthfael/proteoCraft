#### Code chunk - Modified peptides analysis
if (scrptTypeFull == "withReps_PG_and_PTMs") {
  myIDcol <- "Leading proteins"
}
if (scrptTypeFull == "withReps_PTMs_only") {
  myIDcol <- "Proteins"
  PTM_normalize %<o% FALSE
}
#
PepLabKol %<o% setNames(c("Mod. sequence", "Proteins", "Common protein names", "PEP"),
                        c("Mod. sequence", "Protein(s)", "Protein name(s)", "PEP"))
pep$"Mod. sequence" <- gsub("^_|_$", "", pep$"Modified sequence")
if (("Phospho.analysis" %in% colnames(Param))&&(Param$Phospho.analysis)) {
  warning("This argument is deprecated - use the \"PTM.analysis\" argument now.\nI understand you want to perform phosphopeptides analysis though so will be doing the editing for you... this time!\nCarry on please...")
  if ("PTM.analysis" %in% colnames(Param)) {
    a <- toupper(unlist(strsplit(Param$PTM.analysis, ";")))
    if (!"PHOSPHO" %in% a) { Param$PTM.analysis <- paste(c(a, "PHOSPHO"), collapse = ";") }
  } else { Param$PTM.analysis <- "PHOSPHO" }
}
if ("PTM.analysis" %in% colnames(Param)) {
  PTMs %<o% unlist(strsplit(Param$PTM.analysis, ";"))
  if (length(PTMs)) {
    msg <- "Modified peptides analysis"
    ReportCalls <- AddMsg2Report(Space = FALSE)
    PTMs_ref.ratios %<o% list()
    PTMs_FDR.thresholds %<o% list()
    PTMs_pep %<o% list()
    PTMs_Reg_filters %<o% list()
    PTMs_int.ref %<o% list()
    PTMs_rat.ref %<o% list()
    PTM_normalize_NAsReplaceMethod %<o% list()
    PTMs_SAM_thresh %<o% list()
    if (F.test) {
      PTMs_F_test_data %<o% list()
      #PTMs_F_test_ref_ratios %<o% list() # Not needed
    }
    if (("PTM.analysis_Norm" %in% colnames(Param))&&(Param$PTM.analysis_Norm != "")) {
      tmp <- as.logical(unlist(strsplit(as.character(Param$PTM.analysis_Norm), ";")))
      if ((!length(tmp) %in% c(1, length(PTMs)))||(NA %in% tmp)) {
        warning("I cannot make sense of parameter PTM.analysis_Norm, defaulting to TRUE")
        tmp <- TRUE
      }
      if ((length(tmp) == 1)&&(length(PTMs) > 1)) {
        warning(paste0("Recycling PTM.analysis_Norm (value = ", as.character(tmp), ") over PTMs..."))
        tmp <- rep(tmp, length(PTMs))
      }
      PTM_normalize %<o% setNames(tmp, PTMs)
    } else { PTM_normalize %<o% setNames(rep(TRUE, length(PTMs)), PTMs) }
    if (enrichGO) {
      PTMs_GO_Plots %<o% list()
      PTMs_Reg_GO_terms %<o% list()
      PTMs_GO_enrich.dat %<o% list()
      PTMs_GO_enrich.FCRt %<o% list()
      PTMs_GO_enrich.tbl %<o% list()
      PTMs_GO_Plots %<o% list()
      PTMs_Reg_GO_terms %<o% list()
      ### Check that Cytoscape is installed and can run, then launch it.
      Src <- paste0(libPath, "/extdata/R scripts/Sources/Cytoscape_init.R")
      #rstudioapi::documentOpen(Src)
      source(Src, local = FALSE)
      # Initialize ClueGO
      Src <- paste0(libPath, "/extdata/R scripts/Sources/ClueGO_init.R")
      #rstudioapi::documentOpen(Src)
      source(Src, local = FALSE)
      #
    }
    tmp <- PTMs; l <- length(tmp)
    if (l > 1) { tmp <- paste0(paste(tmp[seq_len(l-1)], collapse = ", "), " and ", tmp[l], "") }
    tmp2 <- setNames(sapply(PTMs, function(Ptm) { PTM_normalize[[Ptm]] }), PTMs); tmp2u <- unique(tmp2); l2u <- length(tmp2u)
    if ((l2u == 1)&&(tmp2u)) {
      tmp2 <- ", re-normalizing values to account for average parent protein group(s) fold change."
    } else {
      if (l2u > 1) {
        tmp2 <- names(tmp2)[which(tmp2)]; l2 <- length(tmp2)
        if (l2 > 1) { tmp2 <- paste0(paste(tmp2[seq_len(l-1)], collapse = ", "), " and ", tmp2[l], "") }
        tmp2 <- ", normalizing values to correct for parent protein group(s) fold change for "
      } else { tmp2 <- "." }
    }
    DatAnalysisTxt <- paste0(DatAnalysisTxt, " Statistical analysis was performed for ", c("", "PTMs ")[(l > 1)+1], tmp,
                             "-modified peptides as for protein groups", tmp2)
    for (ptm in PTMs) { #ptm <- PTMs[1]
      ptms.ref %<o% pep.ref[length(pep.ref)]
      ptms.ratios.ref %<o% setNames(pep.ratios.ref[length(pep.ratios.ref)],
                                    "Original")
      a <- unlist(strsplit(gsub("\\)$", "", ptm), "\\("))
      if (length(a) > 1) {
        Ptm <- paste0(toupper(substr(a[1], 1, 1)), substr(a[1], 2, nchar(ptm)), "(", a[2], ")")
      } else {
        Ptm <- paste0(toupper(substr(a[1], 1, 1)), substr(a[1], 2, nchar(ptm)))
      }
      w <- which(Modifs$"Full name" == Ptm)
      if (length(w) == 1) {
        p <- Modifs$Mark[w]
        ppat <- paste0("\\(", p, "\\)|\\(", p, ",|,", p, "\\)|,", p, ",") # Pattern to catch all instances of the mod
        pep[[Ptm]] <- grepl(ppat, pep$"Modified sequence")
        g <- which(pep[[Ptm]])
        if (length(g)) {
          source(parSrc, local = FALSE)
          ReportCalls <- AddMsg2Report(Msg = paste0(" - ", ptm), Space = FALSE)
          dir <- c("", "/t-tests")
          if (F.test) { dir <- c(dir, "/F-tests") }
          dir <-  paste0(wd, "/Reg. analysis/", ptm, dir)
          for (d in dir) { if (!dir.exists(d)) { dir.create(d, recursive = TRUE) }}
          dirlist <- unique(c(dirlist, dir))
          pep[[paste0(Ptm, " ID")]] <- ""
          ptmpep <- pep[g,]
          pep[g, paste0(Ptm, " ID")] <- ptmpep$ModPep_ID <- seq_len(nrow(ptmpep))
          ptmsh <- substr(p, 1, 1)
          temp <- ptmpep[, c("Modified sequence", myIDcol)]
          temp[[myIDcol]] <- strsplit(temp[[myIDcol]], ";")
          temp$"Modified sequence" <- gsub(paste0("[^A-Z", ptmsh, "]"), "",
                                           gsub(ppat, ptmsh, temp$"Modified sequence"))
          #ptmpep[, c("Match(es)", paste0(Ptm, "-site(s)"))] <- ""
          dbsmall <- db[which(db$"Protein ID" %in% unique(unlist(temp[[myIDcol]]))), c("Protein ID", "Sequence")]
          # On I/L ambiguity remaining even with newer DIA methods taking into account RT, IM and fragments intensity, see https://github.com/vdemichev/DiaNN/discussions/1631
          dbsmall$"Seq*" <- gsub("I", "L", dbsmall$Sequence)
          temp$"ModSeq*" <- gsub("I", "L", temp$"Modified sequence")
          #
          kol <- c(myIDcol, "ModSeq*")
          temp$`ModSeq*` <- strsplit(temp$`ModSeq*`, "")
          temp2 <- temp[, kol]
          clusterExport(parClust, list("temp2", "dbsmall", "ptmsh"), envir = environment())
          temp3 <- parApply(parClust, temp2, 1, function(x) {
            #x <- temp[1, kol]
            #A <- sapply(seq_len(nrow(temp)), function(i) {
            #print(i)
            #x <- temp[i, kol]
            m <- unlist(x[[2]])
            m <- data.frame(Seq = m, Mod.seq = m, Test = FALSE)
            w1 <- which(m$Seq == ptmsh)
            w2 <- which(m$Seq != ptmsh)
            m$Mod.seq[w1-1] <- paste0(ptmsh, m$Mod.seq[w1-1])
            m$Test[w1-1] <- TRUE
            m <- m[w2,]
            l <- nrow(m)
            m$Offset <- 0:(l-1)
            q <- unlist(x[[1]])
            mtch <- match(q, dbsmall$"Protein ID")
            wN <- which(!is.na(mtch))
            mtch <- mtch[wN]
            q <- dbsmall$"Protein ID"[mtch[wN]]
            if (length(mtch)) {
              seq <- strsplit(dbsmall$"Seq*"[mtch], "")
              matches <- lapply(seq, function(S) { #S <- seq[1]
                S <- unlist(S)
                lS <- length(S)
                m1 <- m
                m1$Match <- apply(m1[, c("Seq", "Offset")], 1, function(y) { which(S == y[1]) - as.numeric(y[2]) })
                M <- unlist(m1$Match)
                M <- M[which(M > 0)]
                M <- aggregate(M, list(M), length)
                M <- M[order(-M$x),]
                M <- M$Group.1[which(M$x == l)]
                # Check that peptides are tryptic:
                #test <- sapply(M, function(y) {
                #  # r1: on the N-terminal end, is the peptide preceded by K, R or (if starting at position 2, M)?
                #  if (y > 1) {
                #    if (y == 2) { r1 <- S[y-1] %in% c("K", "R", "M") } else { r1 <- S[y-1] %in% c("K", "R") }
                #  } else { r1 <- TRUE }
                #  # r2: on the C-terminal end, is this a tryptic peptide or the last peptide in the protein?
                #  r2 <- (m1$Seq[l] %in% c("K", "R"))|(y+l-1 == lS)
                #  return(r1+r2 == 2)
                #})
                #M <- M[which(test)]
                return(M)
              })
              names(matches) <- q
              matches <- matches[which(sapply(matches, length) > 0)]
              if (length(matches)) {
                matches <- proteoCraft::listMelt(matches, ColNames = c("Match", "Protein"))
                matches <- aggregate(matches$Protein, list(matches$Match), paste, collapse = ";")
                colnames(matches) <- c("Match", "Proteins")
                w <- which(m$Test)
                matches$Sites <- sapply(matches$Match, function(y) {
                  y <- paste(sapply(w, function(z) {paste0(m$Mod.seq[z], y+m$Offset[z])}), collapse = "-")
                })
                matches$Match <- apply(matches[, c("Match", "Proteins")], 1, paste, collapse = " ")
                matches$Sites <- apply(matches[, c("Sites", "Proteins")], 1, paste, collapse = " ")
                matches <- apply(matches[, c("Match", "Sites")], 2, paste, collapse = "/")
              } else { matches <- c(NA, NA) }
            } else { matches <- c(NA, NA) }
            return(matches)
            #})
          })
          ptmpep[, c("Match(es)", paste0(Ptm, "-site(s)"))] <- as.data.frame(t(temp3))
          ptmpep[[paste0(Ptm, "-site")]] <- gsub(" .+", "", ptmpep[[paste0(Ptm, "-site(s)")]])
          ptmpep <- ptmpep[which(!is.na(ptmpep$`Match(es)`)),]
          ptmpep$tmp1 <- gsub("^_|_$", "", ptmpep$"Modified sequence")
          ptmpep$tmp2 <- ptmpep[[paste0(Ptm, "-site(s)")]]
          nc <- nchar(ptmpep$tmp2)
          w <- which(nc > 25)
          ptmpep$tmp2[w] <- paste0(substr(ptmpep$tmp2[w], 1, 22), "...")
          ptmpep$Code <- apply(ptmpep[, paste0("tmp", 1:2)], 1, paste, collapse = "\n")
          ptmpep$tmp1 <- NULL
          ptmpep$tmp2 <- NULL
          ptmpep$Name <- ""
          w <- which(ptmpep[[myIDcol]] != "")
          temp2 <- as.data.frame(t(sapply(strsplit(gsub("[/,;].+$", "", ptmpep[w, paste0(Ptm, "-site(s)")]), " "), unlist)))
          colnames(temp2) <- c("Site", "Protein")
          temp2$Protein <- db$"Common Name"[match(temp2$Protein, db$"Protein ID")]
          temp2$ModSeq <- gsub("^_|_$", "", ptmpep$`Modified sequence`[w])
          ptmpep$Name[w] <- do.call(paste, c(temp2[, c("Site", "ModSeq", "Protein")], sep = "\n"))
          ptmpep$Name[which(ptmpep$Name == "")] <- paste0("Unknown source ", ptm, "-modified peptide #", seq_along(which(ptmpep$Name == "")))
          #View(ptmpep[,c("Match(es)", "Modified sequence", "Code", paste0(Ptm, "-site(s)"))])
          if (grepl("^[P,p]hospho( \\([A-Z]+\\))?$", ptm)) {
            p_col <- paste0(gsub(" |\\(|\\)", ".", ptm), ".Probabilities")
            scd_col <- paste0(gsub(" |\\(|\\)", ".", ptm), ".Score.Diffs")
            if (sum(c(p_col, scd_col) %in% colnames(ptmpep)) == 2) {
              temp <- try(phos_QC(ptmpep, P_col = p_col, ScD_col = scd_col), silent = TRUE)
              if ("try-error" %in% class(temp)) {
                warning("No phospho QC performed, check colnames or the phos_QC function!")
              } else { ptmpep$High_Quality_ptmpep <- temp }
            }
          }
          # Convert intensities to log10
          # I've had it! I should've done this way earlier
          ref <- ptms.ref[1]
          ref2 <- paste0("log10(", gsub(" - $", ") - ", gsub("Evidence intensities", "PSMs. int.", ref)))
          kols <- grep(topattern(ref), colnames(ptmpep), value = TRUE)
          for (kol in kols) {
            smpl <- gsub(topattern(ref), "", kol)
            ptmpep[[paste0(ref2, smpl)]] <- log10(ptmpep[[kol]])
            ptmpep[[kol]] <- NULL
          }
          ptms.ref <- c("Original" = ref2)
          #
          # Optional: normalize to parent protein group
          if ((scrptTypeFull == "withReps_PG_and_PTMs")&&(PTM_normalize[[Ptm]])) {
            # Step 1: normalize ratios:
            # We essentially want to correct the fold change of each modified peptide by that of the parent protein
            #
            #
            #
            ################################
            #         IMPORTANT!!!         #
            ################################
            #
            #  This can only work if any protein groups-level re-normalisation is propagated back onto peptides!!!
            #  This is currently the case. Make sure it stays so!
            #
            #
            #
            #
            kolRPp <- grep(topattern(ptms.ratios.ref["Original"]), colnames(ptmpep), value = TRUE)
            grps <- gsub(topattern(ptms.ratios.ref["Original"]), "", kolRPp)
            kolRPr <- paste0(Prot.Rat.Root, grps)
            stopifnot(sum(!kolRPr %in% colnames(PG)) == 0)
            kolXprPr <- grep(topattern(Prot.Expr.Root), colnames(PG), value = TRUE)
            grps <- gsub("_REF.to.REF_[0-9]+$", "", grps)
            w1 <- which(grps %in% RSA$values)
            w2 <- which(grps %in% RRG$values)
            grps[w1] <- Exp.map[match(grps[w1], Exp.map[[RSA$column]]),
                                RG$column]
            grps[w2] <- Exp.map[match(grps[w2], Exp.map[[RRG$column]]),
                                RG$column]
            tmpPG <- PG[, c("id", kolXprPr, kolRPr)]
            clusterExport(parClust, list("tmpPG", "is.all.good", "kolXprPr", "kolRPr"), envir = environment())
            #
            uPG <- unique(ptmpep$"Protein group IDs")
            tmpPGXpr <- parSapply(parClust, strsplit(uPG, ";"), function(x) {
              #x <- strsplit(uPG, ";")[1]
              x <- unlist(x)
              y <- tmpPG[match(x, tmpPG$id), kolXprPr, drop = FALSE]
              if (length(x) > 1) { y <- apply(y, 2, function(x) { mean(is.all.good(x)) }) }
              return(unlist(y))
            })
            tmpPGXpr <- as.data.frame(t(tmpPGXpr))
            tmpPGRat <- parSapply(parClust, strsplit(uPG, ";"), function(x) {
              #x <- strsplit(uPG, ";")[1]
              x <- unlist(x)
              y <- tmpPG[match(x, tmpPG$id), kolRPr, drop = FALSE]
              if (length(x) > 1) { y <- apply(y, 2, function(x) { mean(is.all.good(x)) }) }
              return(unlist(y))
            })
            tmpPGRat <- as.data.frame(t(tmpPGRat))
            #
            # We have some missing values here, which will mean loss of values for some peptides.
            # We will need to replace them... with what?
            # In most experiments most proteins do not change.
            #  - For now we replace missing values with the median of the column (NAsReplMeth == 2)
            #  - Alternatively, there is also code for imputation (NAsReplMeth == 1),
            #    but in this context this means adding random variation and does not feel like a good idea.
            #
            # (In the future, it would be possible to have this under the control of a parameter.)
            #
            # Note that for now the correction is sample-specific, but it could also be sample group-specific
            # (use average value). There are pros and cons to this.
            #
            #PTM_normalize_NAsReplaceMethod[[Ptm]] <- 0
            PTM_normalize_NAsReplaceMethod[[Ptm]] <- NAsReplMeth
            #PTM_normalize_NAsReplaceMethod[[Ptm]] <- 2
            if (PTM_normalize_NAsReplaceMethod[[Ptm]] == 1) {
              for (grp in RG$values) { #grp <- RG$values[1] # Won't work with RRG because of the nested case (one sample per group breaks the imputation)
                em <- Exp.map[which(Exp.map[[RG$column]] == grp),]
                kolEPG0 <- paste0(Prot.Expr.Root, em$Ref.Sample.Aggregate[which(!em$Reference)])
                kolEPG1 <- paste0(Prot.Expr.Root, em$Ref.Sample.Aggregate[which(em$Reference)])
                kolRtPG0 <- paste0(Prot.Rat.Root, em$Ref.Sample.Aggregate[which(!em$Reference)])
                w0 <- which(kolRtPG0 %in% colnames(tmpPGRat))
                kolRtPG1 <- paste0(Prot.Rat.Root, em$Ref.Sample.Aggregate[which(em$Reference)])
                tmp <- Data_Impute2(tmpPGXpr[, c(kolEPG0, kolEPG1)],
                                    c(rep(0, length(kolEPG0)), rep(1, length(kolEPG1))))
                tmpPGXpr[, c(kolEPG0, kolEPG1)] <- tmp$Imputed_data
                kolEPGRf <- paste0(Prot.Expr.Root, grp, ".REF")
                tmpPGXpr[[kolEPGRf]] <- apply(tmpPGXpr[, kolEPG0, drop = FALSE], 1, mean, na.rm = TRUE)
                tmpPGRat[, kolRtPG1] <- sweep(tmpPGXpr[, kolEPG1], 1, tmpPGXpr[[kolEPGRf]], "-")/log10(2)
                if (length(w0)) {
                  # Ratios for individual controls
                  tmpPGRat[, kolRtPG0] <- sweep(tmpPGXpr[, kolEPG0], 1, tmpPGXpr[[kolEPGRf]], "-")/log10(2)
                } # Else no ratios for individual controls!
                #tmpPGRat[, colnames(temp2)] <- temp2
              }
            }
            if (PTM_normalize_NAsReplaceMethod[[Ptm]] == 2) {
              for (k in kolRPr) { #k <- kolRPr[1]
                w <- which(is.na(tmpPGRat[[k]]))
                tmpPGRat[w, k] <- median(PG[[k]], na.rm = TRUE)
              }
              for (k in kolXprPr) { #k <- kolXprPr[1]
                w <- which(is.na(tmpPGXpr[[k]]))
                tmpPGXpr[w, k] <- median(PG[[k]], na.rm = TRUE)
              }
            }
            #sum(is.na(unlist(tmpPGRat)))
            #sum(is.na(unlist(tmpPGXpr)))
            # Ref-to-Ref ratios:
            temp <- make_RefRat(data = tmpPGXpr,
                                int.root = Prot.Expr.Root,
                                rat.root = Prot.Rat.Root,
                                logInt = 10)
            #sum(is.na(unlist(temp)))
            tmpPGRat[, colnames(temp)] <- temp
            #
            tmpPGRat <- tmpPGRat[match(ptmpep$"Protein group IDs", uPG), ]
            tmpPGXpr <- tmpPGXpr[match(ptmpep$"Protein group IDs", uPG), ]
            #sum(is.na(unlist(tmpPGRat)))
            #sum(is.na(unlist(tmpPGXpr)))
            #
            #View(tmpPGRat[grsep2(prot.list[1], ptmpep$Proteins),])
            #
            #tst1 <- apply(ptmpep[, kol0], 2, summary);View(tst1)
            ptmpep[, paste0("ReNorm. ", kolRPp)] <- ptmpep[, kolRPp] - tmpPGRat[, kolRPp]
            #sum(is.na(unlist(ptmpep[, paste0("ReNorm. ", kolRPp)]))) == sum(is.na(unlist(ptmpep[, kolRPp])))
            #
            #View(ptmpep[, paste0("ReNorm. ", kolRPp)])
            #tst2 <- apply(ptmpep[, paste0("ReNorm. ", kol0)], 2, summary);View(tst2)
            #tstPlot <- FALSE
            #if (tstPlot) {
            df1 <- ptmpep[, kolRPp]
            df2 <- ptmpep[, paste0("ReNorm. ", kolRPp)]
            # Check that there was no data loss at this stage
            tst <- length(is.all.good(unlist(df1))) == length(is.all.good(unlist(df2)))
            if (!tst) {
              warning(paste0("Are you expecting re-normalisation of ", ptm,
                             "-modified peptides to results in loss of valid ratio values? If not, investigate, because it happened!"))
            }
            pepPlotFun(df1,
                       df2,
                       "Re-normalisation (ratios)",
                       dir[1])
            # The distributions should remain barely changed
            ptms.ref["ReNorm."] <- paste0("ReNorm. ", ptms.ref["Original"])
            ptms.ratios.ref["ReNorm."] <- paste0("ReNorm. ", ptms.ratios.ref["Original"])
            # Step 2: adjust expression values to reflect the re-normalized ratios:
            # There would have been several ways to do this.
            # However, ultimately the easiest is to use the values we just calculated
            # and subtract them from the peptide quant values
            # (so not keeping the row sums constant)
            # This ensures the correction is applied consistently.
            #
            kolRPp1 <- grep("_REF\\.to\\.REF_", kolRPp, value = TRUE, invert = TRUE) # Pep. orig. rat. per sample
            kolEPp1 <- gsub(topattern(ptms.ratios.ref["Original"]),
                            ptms.ref["Original"], kolRPp1) # Pep. orig. int. per sample
            kolRNrmEPp1 <- paste0("ReNorm. ", kolEPp1) # Pep. renorm. int. per sample
            ptmpep[, kolRNrmEPp1] <- ptmpep[, kolEPp1] - tmpPGRat[, kolRPp1]/log2(10)
            # Nested case:
            # (deal with samples without ratios, i.e. ref. samples in nested case)
            kolEPp0 <- paste0(ptms.ref["Original"], Exp.map$Ref.Sample.Aggregate)
            kolEPp0 <- kolEPp0[which(!kolEPp0 %in% kolEPp1)]
            kolRNrmEPp <- kolRNrmEPp1
            if (length(kolEPp0)) {
              kolRNrmEPp0 <- paste0("ReNorm. ", kolEPp0)
              ptmpep[, kolRNrmEPp0] <- ptmpep[, kolEPp0] # These do not change after renorm (parent PG logFC = 0)
              kolRNrmEPp <- c(kolRNrmEPp0, kolRNrmEPp)
            }
            #
            kolEPp <- c(kolEPp0, kolEPp1)
            temp1 <- ptmpep[, kolEPp]
            temp2 <- ptmpep[, kolRNrmEPp]
            w1 <- which(is.infinite(as.matrix(temp1)), arr.ind = TRUE)
            w2 <- which(is.infinite(as.matrix(temp2)), arr.ind = TRUE)
            if (nrow(w1)) { temp1[w1] <- NA }
            if (nrow(w2)) { temp2[w2] <- NA }
            #sum(is.na(unlist(temp1))) == sum(is.na(unlist(temp1)))
            beforSum <- rowSums(temp1, na.rm = TRUE)
            afterSum <- rowSums(temp2, na.rm = TRUE)
            #test <- all(is.all.good(beforSum, 2) == is.all.good(afterSum, 2)); print(test)
            #View(data.frame(Before = beforSum, After = afterSum))
            ptmpep[, kolRNrmEPp] <- sweep(ptmpep[, kolRNrmEPp], 1, beforSum-afterSum, "-")
            #sum(is.na(unlist(ptmpep[, kolRNrmEPp])))
            # What about ref columns?
            for (grp in RRG$values) { #grp <- RRG$values[1]
              em <- Exp.map[which(Exp.map[[RRG$column]] == grp),]
              w0 <- which(em$Reference)
              kolEPpRf <- paste0(ptms.ref["Original"], grp, ".REF")
              kolRNrmEPpRf <- paste0(ptms.ref["ReNorm."], grp, ".REF")
              kolRNrmE0 <- paste0(ptms.ref["ReNorm."], em$Ref.Sample.Aggregate[w0])
              if (length(w0) == 1) {
                # Nested
                ptmpep[[kolRNrmEPpRf]] <- ptmpep[[kolRNrmE0]]
              } else {
                # Non-nested
                ptmpep[[kolRNrmEPpRf]] <- parApply(parClust, ptmpep[, kolRNrmE0], 1, log_ratio_av)
              }
            }
            #if (tstPlot) {
            k1 <- grep(topattern(ptms.ref["Original"]), colnames(ptmpep), value = TRUE)
            k2 <- gsub(topattern(ptms.ref["Original"]), ptms.ref["ReNorm."], k1)
            k2a <- grep(topattern(ptms.ref["ReNorm."]), colnames(ptmpep), value = TRUE)
            sum(!k2 %in% k2a) == 0
            sum(!k2a %in% k2) == 0
            df1 <- ptmpep[, k1]
            df2 <- ptmpep[, k2]
            w <- which(is.na(df1), arr.ind = TRUE)
            if (nrow(w)) { df2[w] <- df1[w] }
            #View(df1); View(df2)
            # Check that there was no data loss at this stage
            tst <- length(is.all.good(unlist(df1))) == length(is.all.good(unlist(df2)))
            if (!tst) {
              warning(paste0("Are you expecting re-normalisation of ", ptm,
                             "-modified peptides to results in losses of valid intensity values? If not, investigate, because it happened!"))
            }
            pepPlotFun(df1,
                       df2,
                       "Re-normalisation (intensities)",
                       dir[1])
            #}
            #kol <- grep("ReNorm. log2", colnames(ptmpep), value = TRUE)
            #tst <- apply(ptmpep[, kol], 2, function(x) { summary(is.all.good(x)) })
            if (length(prot.list_pep)) {
              source(parSrc, local = FALSE)
              pepHtmp(prot.list_pep,
                      ptmpep,
                      ptms.ref["Original"],
                      paste0(dir[1], "/Heatmaps"),
                      paste0(ptm, "-mod. pept. log2 heatmap, original"),
                      is.log = TRUE,
                      cl = parClust)
              pepHtmp(prot.list_pep,
                      ptmpep,
                      ptms.ref["ReNorm."],
                      paste0(dir[1], "/Heatmaps"),
                      paste0(ptm, "-mod. pept. log2 heatmap, re-normalised"),
                      is.log = TRUE,
                      cl = parClust)
            }
          }
          #
          pepRf <- ptms.ref[length(ptms.ref)]
          pepRatRf <- ptms.ratios.ref[length(ptms.ratios.ref)]
          #View(ptmpep[, grep(topattern(pepRf), colnames(ptmpep), value = TRUE)])
          #View(ptmpep[, grep(topattern(pepRatRf), colnames(ptmpep), value = TRUE)])
          #
          # Calculate average intensities and ratios, as well as Welch's t-test and moderated P-values;
          # For unpaired replicates a permutations t-test is also performed.
          samDir <- paste0(dir[2], "/SAM")
          ebamDir <- paste0(dir[2], "/EBAM")
          dataType <- "modPeptides"
          #
          Src <- paste0(libPath, "/extdata/R scripts/Sources/Av_and_Stat_tests.R")
          #rstudioapi::documentOpen(Src)
          source(Src, local = FALSE)
          #
          #kol <- grep(topattern(paste0("Mean ", pepRf)), colnames(ptmpep), value = TRUE)
          #tst <- apply(ptmpep[, kol], 2, function(x) { summary(is.all.good(log10(x))) })
          #kol2 <- grep(topattern(paste0("Mean ", ptms.ratios.ref[length(ptms.ratios.ref)])), colnames(ptmpep), value = TRUE)
          #tst2 <- apply(ptmpep[, kol2], 2, function(x) { summary(is.all.good(x)) })
          #
          # Also mean expression over whole dataset
          a <- grep(topattern(pepRf), colnames(ptmpep), value = TRUE)
          ptmpep$"Mean Expr." <- apply(ptmpep[, a], 1, function(x) { mean(is.all.good(unlist(x))) })
          # Create list of control ratio values for the purpose of identifying thresholds for plots:
          #
          if (Param$Ratios.Thresholds == "Absolute log2 FC threshold") {
            ref.rat <- NULL
          }
          if (Param$Ratios.Thresholds == threshMsg) {
            PTMs_ref.ratios[[Ptm]] <- ref.rat <- setNames(lapply(VPAL$values, function(x) { #x <- VPAL$values[1]
              if (RatConGrps == "Ratio groups") {
                x1 <- unique(Exp.map[which(Exp.map[[VPAL$column]] == x), RG$column])
              }
              if (RatConGrps == "Experiments") {
                x1 <- unique(Exp.map$Experiment[which(Exp.map[[VPAL$column]] == x)])
                x1 <- unique(Exp.map[which(Exp.map$Experiment == x1), RG$column])
              }
              if (RatConGrps == "Whole dataset") {
                x1 <- unique(Exp.map[[RG$column]])
              }
              x <- grep(paste0(topattern(paste0(ptms.ratios.ref[length(ptms.ratios.ref)], x1, "_REF.to.REF_")), "[0-9]+"),
                        colnames(ptmpep), value = TRUE)
              if (length(x)) { x <- is.all.good(as.numeric(unlist(ptmpep[, x]))) } else { x <- NULL }
              return(x)
            }), VPAL$values)
          }
          #
          # Estimate P-value significance for a set of accepted FDRs:
          ## NB: For graphical reasons (volcano plots), there is only support for 4 different FDR values. This should suffice anyway.
          temp_thrsh <- c()
          A <- VPAL$values
          test <- sapply(A, function(x) { #x <- A[6]
            x <- paste0(pvalue.col[which(pvalue.use)], x)
            r <- x %in% colnames(ptmpep)
            if (r) { r <- length(is.all.good(as.numeric(ptmpep[[x]]))) > 0 }
            return(r)
          })
          A <- A[which(test)]
          for (a in A) { #a <- A[1]
            temp <- FDR(data = ptmpep,
                        aggregate = a,
                        pvalue_root = pvalue.col[which(pvalue.use)],
                        fdr = BH.FDR, returns = c(TRUE, TRUE), method = "BH")
            ptmpep[, colnames(temp$`Significance vector`)] <- temp$`Significance vector`
            temp_thrsh <- c(temp_thrsh, temp$Thresholds)
            PTMs_FDR.thresholds[[Ptm]] <- temp_thrsh
          }
          ptmpep$"1-PEP" <- 1 - ptmpep$PEP
          ptmpep$"log10(1-PEP)" <- log10(ptmpep$"1-PEP")
          a <- grep(topattern(pepRf), colnames(ptmpep), value = TRUE)
          a <- a[which(!grepl("\\.REF$", a))]
          ptmpep$"Av. log10 abundance" <- apply(ptmpep[, a], 1, function(x) { mean(is.all.good(unlist(x))) })
          ptmpep$"Rel. av. log10 abundance" <- ptmpep$"Av. log10 abundance"/max(is.all.good(ptmpep$"Av. log10 abundance"))
          # Arbitrary thresholds
          P <- Param
          P$Plot.labels <- "Name"
          P$Plot.metrics <- paste0("X:", paste0(ptms.ratios.ref[length(ptms.ratios.ref)], "Mean."),
                                   ";Y:", pvalue.col[which(pvalue.use)])
          ReportCalls <- AddMsg2Report(Msg = " - t-tests", Space = FALSE)
          #k1 <- grep(topattern(paste0("Mean ", ptms.ratios.ref[length(ptms.ratios.ref)])), colnames(ptmpep), value = TRUE)
          #df1 <- ptmpep[, k1]
          #subDr <- gsub(topattern(wd), "", dir[2])
          subDr <- dir[2]
          if (useSAM) {
            # In this case, we bypass the original decision and base it off SAM even though we plot Student's P-values
            for (i in names(PTMs_SAM_thresh[[Ptm]])) { #i <- names(PTMs_SAM_thresh[[Ptm]])[1]
              dec <- PTMs_SAM_thresh[[Ptm]][[i]]$decision
              mKol <- rev(colnames(dec))[1]
              FCkol <- paste0("Mean ", ptms.ratios.ref[length(ptms.ratios.ref)], i)
              stopifnot(FCkol %in% names(ptmpep))
              regKol <- paste0("Regulated - ", i)
              ptmpep[[regKol]] <- "non significant"
              fdrs <- as.numeric(gsub("FDR$", "", colnames(dec)[which(colnames(dec) != mKol)]))
              fdrs <- sort(fdrs, decreasing = TRUE)
              for (f in fdrs) { #f <- fdrs[1]
                w <- which(ptmpep[[mKol]] %in% dec[which(dec[[paste0(f, "FDR")]] == "+"), mKol])
                if (length(w)) {
                  ptmpep[which(ptmpep[w, FCkol] > 0), regKol] <- paste0("up, FDR = ", f*100, "%")
                  if (TwoSided) {
                    ptmpep[which(ptmpep[w, FCkol] < 0), regKol] <- paste0("down, FDR = ", f*100, "%")
                  }
                }
              }
            }
          }
          stopCluster(parClust)
          source(parSrc)
          tempVPptm <- Volcano.plot(ptmpep, "custom",
                                    paste0("Mean ", ptms.ratios.ref[length(ptms.ratios.ref)]),
                                    pvalue.col[which(pvalue.use)],
                                    experiments.map = Exp.map,
                                    aggregate.map = Aggregate.map,
                                    aggregate.list = Aggregate.list,
                                    aggregate.name = VPAL$aggregate,
                                    parameters = P,
                                    save = c("jpeg", "pdf"),
                                    labels = c("FDR", "both")[useSAM+1],
                                    Ref.Ratio.values = ref.rat,
                                    Ref.Ratio.method = paste0("obs", RefRat_Mode),
                                    ratios.FDR = as.numeric(Param$Ratios.Contamination.Rates),
                                    FDR.thresh = PTMs_FDR.thresholds[[Ptm]], arbitrary.lines = arbitrary.thr,
                                    proteins = prot.list, IDs.col = "Code", Proteins.col = "Proteins",
                                    proteins_split = protsplit, return = TRUE, return.plot = TRUE,
                                    title = paste0(Ptm, " volcano plot "), subfolder = subDr,
                                    subfolderpertype = FALSE, Symmetrical = TwoSided,
                                    Size = "Rel. av. log10 abundance", Size.max = 2,
                                    plotly = create_plotly, plotly_local = create_plotly_local,
                                    plotly_labels = c(PepLabKol, paste0(Ptm, "-site")),
                                    cl = parClust,
                                    SAM = useSAM, curved_Thresh = PTMs_SAM_thresh[[Ptm]]
                                    #, saveData = TRUE
          )
          #k2 <- grep(topattern(paste0("Mean ", ptms.ratios.ref[length(ptms.ratios.ref)])), colnames(ptmpep), value = TRUE)
          #df2 <- ptmpep[, k2]
          #pepPlotFun(df1, df2, "Before VS after volc. plot", dir[1], FALSE)
          stopCluster(parClust)
          source(parSrc)
          #
          # Save plotly plots
          dr <- subDr
          myPlotLys <- tempVPptm
          Src <- paste0(libPath, "/extdata/R scripts/Sources/save_Volcano_plotlys.R")
          #rstudioapi::documentOpen(Src)
          source(Src, local = FALSE)
          #
          VP_list <- tempVPptm
          insrt <- ""
          Src <- paste0(libPath, "/extdata/R scripts/Sources/thresholds_Excel.R")
          #rstudioapi::documentOpen(Src)
          source(Src, local = FALSE)
          #
          #
          ptmpep <- tempVPptm$Protein_groups_file
          volcano.plots[[Ptm]] <- tempVPptm$Plots
          n2 <- names(volcano.plots[[Ptm]]$Labelled)
          for (ttl in n2) {
            plot <- volcano.plots[[Ptm]]$Labelled[[ttl]]
            ReportCalls <- AddPlot2Report(Space = FALSE, Jpeg = FALSE)
          }
          if ((create_plotly)&&(!create_plotly_local)) {
            plot_ly[[paste0(Ptm, "_Volcano plots (t-tests)")]] <- tempVPptm$"Plotly plots"
          }
          # Specificity mark for untested proteins
          # We can assume that proteins with PSMs only in the specific pull-down samples are actually specifically enriched!
          # -> Label them as such!
          # NB: This used to be for pull-downs only but I am extending it to the whole workflow for now.
          #     This could also be parameter-controlled.
          # NB: For the F-test, this is done within the function
          #if (IsPullDown) {
          for (a in RG$values) { #a <- RG$values[1]
            e <- Exp.map[which(Exp.map[[RG$column]] == a),]
            grp1 <- unique(e[which(!e$Reference %in% c(TRUE, "TRUE")), VPAL$column])
            kolTR <- paste0("Regulated - ", grp1)
            w <- which(kolTR %in% colnames(ptmpep))
            if (length(w)) {
              grp1 <- grp1[w]
              kolTR <- kolTR[w]
              grp0 <- unique(e[which(e$Reference %in% c(TRUE, "TRUE")), VPAL$column])
              e0 <- unique(e$Ref.Sample.Aggregate[which(e[[VPAL$column]] == grp0)])
              kole0 <- paste0("Evidence IDs - ", e0)
              tst0 <- rowSums(ptmpep[, kole0] == "") == length(kole0)
              for (i in seq_along(grp1)) { #i <- 1
                e1i <- unique(e$Ref.Sample.Aggregate[which(e[[VPAL$column]] == grp1[i])])
                kole1i <- paste0("Evidence IDs - ", e1i)
                tst1i <- rowSums(ptmpep[, kole1i] != "") == length(kole1i)
                w1i <- which(tst0 & tst1i)
                if (length(w1i)) {
                  evcount <- apply(ptmpep[w1i, kole1i, drop = FALSE], 1, function(x) { length(unlist(strsplit(x, ";"))) })
                  txtup <- paste0("Specific: at least ", evcount, " PSMs/sample")
                  ptmpep[w1i, kolTR[i]] <- txtup
                }
              }
            } else { warning("I would expect \"Regulated ...\" columns in the modified peptides table by this stage!") }
          }
          #}
          # Create t-test filters:
          ## These can then be used for further steps down the line, such as volcano plots, etc...
          g <- grep("^Regulated - ", colnames(ptmpep), value = TRUE)
          g1 <- gsub("^Regulated - ", "", g)
          up <- grep("^up|^Specific", unique(unlist(ptmpep[, g])), value = TRUE)
          down <- grep("^down|^Anti-specific", unique(unlist(ptmpep[, g])), value = TRUE) # The "anti-specific" part will only become relevant if in future I disconnect symmetry and/or adding specific tags from pull-down experiments
          PTMs_Reg_filters[[Ptm]] <- list()
          PTMs_Reg_filters[[Ptm]]$"t-tests" <- list()
          if ("con" %in% filter_types) {
            PTMs_Reg_filters[[Ptm]]$"t-tests"$"By condition" <- list()
            for (i in seq_along(g)) {
              PTMs_Reg_filters[[Ptm]]$"t-tests"$"By condition"[[g1[i]]] <- list(Columns = g[i],
                                                                                Filter_up = sort(which(ptmpep[[g[i]]] %in% up)),
                                                                                Filter_down = sort(which(ptmpep[[g[i]]] %in% down)),
                                                                                Filter = sort(which(ptmpep[[g[i]]] %in% c(up, down))))
            }
          }
          if ("ref" %in% filter_types) {
            PTMs_Reg_filters[[Ptm]]$"t-tests"$"By reference" <- list()
            g2 <- as.data.frame(t(as.data.frame(strsplit(g1, "___"))))
            colnames(g2) <- VPAL$names
            tst <- apply(Exp.map[, RRG$names, drop = FALSE], 1, paste, collapse = "___")
            tmp <- apply(g2[,RRG$names, drop = FALSE], 1, paste, collapse = "___")
            g2$Ref <- sapply(tmp, function(x) { y <- unique(tst[which((Exp.map$Reference)&(tst == x))]) })
            for (i in unique(g2$Ref)) {
              w <- which(g2$Ref == i)
              PTMs_Reg_filters[[Ptm]]$"t-tests"$"By reference"[[i]] <- list(Columns = g[w],
                                                                            Filter_up = sort(which(apply(ptmpep[, g[w], drop = FALSE], 1, function(x) {
                                                                              length(which(x %in% up))
                                                                            }) > 0)),
                                                                            Filter_down = sort(which(apply(ptmpep[, g[w], drop = FALSE], 1, function(x) {
                                                                              length(which(x %in% down))
                                                                            }) > 0)),
                                                                            Filter = sort(which(apply(ptmpep[, g[w], drop = FALSE], 1, function(x) {
                                                                              length(which(x %in% c(up, down)))
                                                                            }) > 0)))
            }
          }
          if (sum(c("dat", "dat2") %in% filter_types)) {
            PTMs_Reg_filters[[Ptm]]$"t-tests"$"Whole dataset" <- list(Columns = g,
                                                                      Filter_up = sort(which(apply(ptmpep[, g, drop = FALSE], 1, function(x) {
                                                                        length(which(x %in% up))
                                                                      }) > 0)),
                                                                      Filter_down = sort(which(apply(ptmpep[, g, drop = FALSE], 1, function(x) {
                                                                        length(which(x %in% down))
                                                                      }) > 0)),
                                                                      Filter = sort(which(apply(ptmpep[, g, drop = FALSE], 1, function(x) {
                                                                        length(which(x %in% c(up, down)))
                                                                      }) > 0)))
          }
          if (F.test) {
            #kol <- grep(topattern(pepRf), colnames(ptmpep), value = TRUE)
            #View(ptmpep[, kol])
            ReportCalls <- AddMsg2Report(Msg = " - F-tests", Space = FALSE)
            # NB: id.col below should be unique!
            stopifnot(length(unique(ptmpep$Name)) == nrow(ptmpep))
            #
            expMap <- Exp.map[match(rownames(designMatr), Exp.map$Ref.Sample.Aggregate),]
            Group_ <- do.call(paste, c(expMap[, Coefficients, drop = FALSE], sep = "_"))
            Group_ <- as.factor(Group_)
            expMap$Group_ <- Group_
            dataType <- "modPeptides"
            #
            FSrc %<o% paste0(libPath, "/extdata/R scripts/Sources/run_F_test.R")
            #rstudioapi::documentOpen(FSrc)
            tstFtst <- try(source(FSrc, local = FALSE), silent = TRUE)
            #
            if (!"try-error" %in% class(tstFtst)) {
              #F_test_ref_ratios %<o% F_volc$`Reference ratios` # Not needed
              volcano.plots[[Ptm]]$"F-tests_Unlabelled" <- F_volc$Plots$"Unlabelled"
              volcano.plots[[Ptm]]$"F-tests_Labelled" <- F_volc$Plots$"Labelled"
              n2 <- names(volcano.plots[[Ptm]]$"F-tests_Labelled")
              dir <- paste0(dir, "/Reg. analysis/F-tests")
              for (ttl in n2) {
                plot <- volcano.plots[[Ptm]]$"F-tests_Labelled"[[ttl]]
                ReportCalls <- AddPlot2Report(Space = FALSE, Jpeg = FALSE)
              }
              # Legacy code for web-hosted plotly plots:
              if ((create_plotly)&&(!create_plotly_local)) { plot_ly[[paste0(Ptm, "_Volcano plots (F-tests)")]] <- F_volc$"Plotly plots" }
              # Also a posteriori F-test P- or Q- values histogram:
              nbin <- 20
              bd <- (0:nbin)/nbin
              if (F_Root %in% colnames(ptmpep)) {
                temp <- data.frame(value = is.all.good(10^(-ptmpep[[F_Root]])))
                ttl <- paste0("Histogram: F-test moderated Pvalue (", ptm, ")")
                plot <- ggplot(temp, aes(x = value)) +
                  geom_histogram(bins = nbin, colour = "black", alpha = 0.25, fill = "green") +
                  guides(fill = "none") + theme_bw() + ggtitle(ttl)
                poplot(plot)
                dir <- paste0(wd, "/Workflow control/Protein groups/P-values")
                dirlist <- unique(c(dirlist, dir))
                if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
                ttla <- gsub(": *", " - ", ttl)
                suppressMessages({
                  ggsave(paste0(dir, "/", ttla, ".jpeg"), plot, dpi = 300)
                  ggsave(paste0(dir, "/", ttla, ".pdf"), plot, dpi = 300)
                })
                ReportCalls <- AddPlot2Report(Title = ttla)
              }
              # Create F-test filters:
              g <- grep("^mod\\. F-test Regulated - ", colnames(PTMs_F_test_data[[Ptm]]), value = TRUE)
              g1 <- gsub("^mod\\. F-test Regulated - ", "", g)
              PTMs_Reg_filters[[Ptm]]$"F-tests" <- list()
              if ("con" %in% filter_types) {
                PTMs_Reg_filters[[Ptm]]$"F-tests"$"By condition" <- list()
                for (i in seq_along(g)) {
                  PTMs_Reg_filters[[Ptm]]$"F-tests"$"By condition"[[g1[i]]] <- list(Columns = g[i],
                                                                                    Filter_up = sort(which(PTMs_F_test_data[[Ptm]][[g[i]]] %in% up)),
                                                                                    Filter_down = sort(which(PTMs_F_test_data[[Ptm]][[g[i]]] %in% down)),
                                                                                    Filter = sort(which(PTMs_F_test_data[[Ptm]][[g[i]]] %in% c(up, down))))
                }
              }
              if ("ref" %in% filter_types) {
                PTMs_Reg_filters[[Ptm]]$"F-tests"$"By reference" <- list()
                g2 <- as.data.frame(t(as.data.frame(strsplit(g1, " - "))))
                row.names(g2) <- NULL
                colnames(g2) <- c("Group", "Analysis", "Condition")
                g2$Ref <- ""
                w <- grep(" VS ", g2$Condition)
                if (length(w)) {
                  g2$tmp[w] <- sapply(strsplit(gsub("\\)$", "", g2$Condition[w]), split = " VS "), function(x) {x[[2]]})
                  g2$Ref[w] <- apply(g2[w, c("Group", "Analysis", "tmp")], 1, function(x) {
                    paste0(x[[1]], " - ", x[[2]], ", ref: ", x[[3]])
                  })
                }
                w <- which((grepl(" vs ", g2$Condition))&(!grepl(" VS ", g2$Condition)))
                if (length(w)) {
                  g2$tmp[w] <- sapply(strsplit(g2$Condition[w], split = " vs "), function(x) {x[[2]]})
                  g2$Ref[w] <- apply(g2[w, c("Group", "Analysis", "tmp")], 1, function(x) {
                    paste0(x[[1]], ", ref: ", x[[2]])
                  })
                }
                w <- which(!grepl(" vs ", g2$Condition))
                if (length(w)) { g2$Ref[w] <- g1[w] }
                for (i in unique(g2$Ref)) {
                  w <- which(g2$Ref == i)
                  PTMs_Reg_filters[[Ptm]]$"F-tests"$"By reference"[[i]] <- list(Columns = g[w],
                                                                                Filter_up = sort(which(apply(PTMs_F_test_data[[Ptm]][, g[w], drop = FALSE], 1, function(x) {
                                                                                  length(which(x %in% up))
                                                                                }) > 0)),
                                                                                Filter_down = sort(which(apply(PTMs_F_test_data[[Ptm]][, g[w], drop = FALSE], 1, function(x) {
                                                                                  length(which(x %in% down))
                                                                                }) > 0)),
                                                                                Filter = sort(which(apply(PTMs_F_test_data[[Ptm]][, g[w], drop = FALSE], 1, function(x) {
                                                                                  length(which(x %in% c(up, down)))
                                                                                }) > 0)))
                }
              }
              if ("dat" %in% filter_types) {
                PTMs_Reg_filters[[Ptm]]$"F-tests"$"By analysis" <- list()
                g2 <- as.data.frame(t(as.data.frame(strsplit(g1, " - "))))
                colnames(g2) <- c("Group", "Analysis", "Condition")
                g2$Group_Analysis <- apply(g2[,c("Group", "Analysis")], 1, paste, collapse = "_")
                for (i in unique(g2$Group_Analysis)) {
                  w <- which(g2$Group_Analysis == i)
                  PTMs_Reg_filters[[Ptm]]$"F-tests"$"By analysis"[[i]] <- list(Columns = g[w],
                                                                               Filter_up = sort(which(apply(PTMs_F_test_data[[Ptm]][, g[w], drop = FALSE], 1, function(x) {
                                                                                 length(which(x %in% up))
                                                                               }) > 0)),
                                                                               Filter_down = sort(which(apply(PTMs_F_test_data[[Ptm]][, g[w], drop = FALSE], 1, function(x) {
                                                                                 length(which(x %in% down))
                                                                               }) > 0)),
                                                                               Filter = sort(which(apply(PTMs_F_test_data[[Ptm]][, g[w], drop = FALSE], 1, function(x) {
                                                                                 length(which(x %in% c(up, down)))
                                                                               }) > 0)))
                }
              }
              if ("dat2" %in% filter_types) {
                PTMs_Reg_filters[[Ptm]]$"F-tests"$"Whole dataset" <- list(Columns = g,
                                                                          Filter_up = sort(which(apply(PTMs_F_test_data[[Ptm]][, g, drop = FALSE], 1, function(x) {
                                                                            length(which(x %in% up))
                                                                          }) > 0)),
                                                                          Filter_down = sort(which(apply(PTMs_F_test_data[[Ptm]][, g, drop = FALSE], 1, function(x) {
                                                                            length(which(x %in% down))
                                                                          }) > 0)),
                                                                          Filter = sort(which(apply(PTMs_F_test_data[[Ptm]][, g, drop = FALSE], 1, function(x) {
                                                                            length(which(x %in% c(up, down)))
                                                                          }) > 0)))
              }
            } else { warning("The F-test did not generate any plots! Investigate!") }
          }
          #
          # Gene-Set Enrichment Analysis (GSEA)
          if (runGSEA) {
            dataType <- "modPeptides"
            GSEAmode <- "standard"
            Src <- paste0(libPath, "/extdata/R scripts/Sources/GSEA.R")
            #rstudioapi::documentOpen(Src)
            source(Src, local = FALSE)
          }
          #
          # Heatmap
          g <- paste0(pepRf, RSA$values)
          grps <- Exp.map[match(RSA$values, Exp.map$Ref.Sample.Aggregate), VPAL$column]
          w <- which(g %in% colnames(ptmpep))
          g <- g[w]; grps <- grps[w]
          temp <- ptmpep[, g]
          rownames(temp) <- gsub("\n", " ", ptmpep$Name)
          kol <- gsub(topattern(pepRf), "", g)
          colnames(temp) <- gsub(topattern(pepRf), "", colnames(temp))
          kol <- cleanNms(kol)
          colnames(temp) <- cleanNms(colnames(temp))
          w <- which(!is.finite(as.matrix(temp)), arr.ind = TRUE)
          temp[w] <- NA
          rwMns <- rowMeans(temp[, kol], na.rm = TRUE)
          w <- which(!is.na(rwMns))
          temp <- temp[w,]
          rwMns <- rwMns[w]
          temp[, kol] <- sweep(temp[, kol], 1, rwMns, "-")
          temp2 <- as.matrix(Data_Impute2(temp, grps)$Imputed_data)
          hcluster <- hclust(dist(t(temp2)))
          hdendro <- as.dendrogram(hcluster)
          hord <- order.dendrogram(hdendro)
          hord <- colnames(temp)[hord]
          vcluster <- hclust(dist(temp2))
          vdendro <- as.dendrogram(vcluster)
          vord <- order.dendrogram(vdendro)
          vord <- rownames(temp)[vord]
          require(ggdendro)
          hdendro.plot <- ggdendrogram(data = hdendro) + theme(axis.text.y = element_text(size = 0.1),
                                                               plot.margin = margin(0, 0, 0, 0, "cm"))
          vdendro.plot <- ggdendrogram(data = vdendro, rotate = TRUE, labels = FALSE, leaf_labels = FALSE) +
            theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
                  panel.background = element_rect(fill = "transparent", colour = NA), 
                  plot.background = element_rect(fill = "transparent", colour = NA),
                  plot.margin = margin(0, 0, 0, 0, "cm"))
          # Data wrangling
          temp2 <- set_colnames(reshape::melt.data.frame(temp), c("Sample", "value"))
          temp2$Label <- rownames(temp)
          temp2$Sample <- as.character(temp2$Sample)
          temp2$Colour <- "grey"
          temp2$Size <- 1
          # Extract the order of the tips in the dendrograms
          # Order the levels according to their position in the clusters
          temp2$Xmin <- match(temp2$Sample, hord)-1
          temp2$Ymin <- match(temp2$Label, vord)-1
          Xscale <- length(unique(temp2$Sample))
          temp2$Label2 <- temp2$Label
          w <- which(nchar(temp2$Label2) > 25)
          temp2$Label2[w] <- paste0(substr(temp2$Label2[w], 1, 22), "...")
          # Create heatmap plot
          w1 <- which(temp2$Colour == "green")
          w2 <- which((temp2$Xmin == max(temp2$Xmin))&(temp2$Colour == "green"))
          nm <- paste0("Heatmap\n", Ptm, "-modified peptides")
          heatmap.plot <- ggplot(temp2) +
            geom_rect(aes(xmin = Xmin, xmax = Xmin+1, ymin = Ymin, ymax = Ymin+1, fill = value)) +
            geom_text(data = temp2, aes(y = Ymin+0.5, label = Label2),
                      x = length(unique(temp2$Sample)), colour = "grey", hjust = 0, vjust = 0.5, size = 2) +
            theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
                  axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
                  panel.background = element_rect(fill = "transparent", colour = NA),
                  plot.margin = margin(0, 0, 0, 0, "cm")) +
            scale_fill_gradient2(low = "darkblue", mid = "lightgrey", high = "darkred") +
            xlab(NULL) + ylab(NULL)
          leg <- get_legend(heatmap.plot)
          heatmap.plot <- heatmap.plot + theme(legend.position = "none")
          htmp <- arrangeGrob(grobs = list(hdendro.plot, leg, heatmap.plot, vdendro.plot),
                              widths = c(1, 2*Xscale, 1, 3), heights = c(3, 5),
                              layout_matrix = rbind(c(NA, 1, NA, 2), c(3, 3, 3, 4)),
                              padding = 5)
          nm <- gsub("\n", " - ", nm)
          suppressMessages({
            ggsave(paste0(dir[1], "/", nm, ".jpeg"), htmp, width = 20, height = 20, units = "in", dpi = 600)
            ggsave(paste0(dir[1], "/", nm, ".pdf"), htmp, width = 20, height = 20, units = "in", dpi = 600)
          })
          ReportCalls <- AddPlot2Report(Title = nm, Dir = dir[1])
          #system(paste0("open \"", dir[1], "/", nm, ".jpeg", "\""))
          #system(paste0("open \"", dir[1], "/", nm, ".pdf", "\""))
          #
          # Gene Ontology terms enrichment analysis
          if (enrichGO) {
            if ((!exists("GO_mappings"))&&(file.exists("GO_mappings.RData"))) {
              loadFun("GO_mappings.RData")
            }
            if ((!exists("GO_terms"))&&(file.exists("GO_terms.RData"))) {
              loadFun("GO_terms.RData")
            }
            p <- strsplit(ptmpep$Proteins, ";")
            test <- sapply(annot.col, function(x) { x %in% colnames(db) })
            annot.col2 <- annot.col[which(test)]
            tmpDB <- db[which(db$Observed), c("Protein ID", annot.col2)]
            clusterExport(parClust, c("tmpDB", "annot.col2"), envir = environment())
            temp <- parLapply(parClust, p, function(x) {
              m <- match(x, tmpDB$"Protein ID")
              y <- tmpDB[m, annot.col2]
              if (length(m) > 1) {
                y <- apply(y, 2, function(z) {
                  z <- z[which(!is.na(z))]
                  paste(sort(unique(unlist(strsplit(as.character(z), ";")))), collapse = ";")
                })
              }
              y <- setNames(as.character(unlist(y)), annot.col2)
              return(y)
            })
            temp <- as.data.frame((do.call(rbind, temp)))
            for (i in annot.col2) {
              w <- which(is.na(temp[[i]]))
              if (length(w)) { temp[w, i] <- "" }
            }
            ptmpep[, annot.col2] <- temp
            Kol <- unique(c("Modified sequence", "Evidence IDs", "Sequence",
                            myIDcol, "Proteins", "Leading razor proteins",
                            "PEP", "id", "Protein group IDs", "Razor protein group ID", 
                            "Leading razor protein", "Gene names", "Protein names", "Gene names (all)",
                            "Protein names (all)", "Weights", "Match(es)", "Code", "Name", annot.col2,
                            "Potential contaminant"))
            Kol <- Kol[which(Kol %in% colnames(ptmpep))]
            PTMs_GO_enrich.dat[[Ptm]] <- list()
            PTMs_GO_enrich.FCRt[[Ptm]] <- list()
            PTMs_GO_enrich.tbl[[Ptm]] <- list()
            PTMs_GO_Plots[[Ptm]] <- list()
            PTMs_Reg_GO_terms[[Ptm]] <- list()
            Tsts <- c("t-tests", "F-tests")
            WhTsts <- which(Tsts %in% names(PTMs_Reg_filters[[Ptm]]))
            for (tt in WhTsts) { #tt <- 1 #tt <- 2
              tstrt <- Tsts[tt]
              stopifnot(!is.na(tstrt))
              dir <- paste0(wd, "/Reg. analysis/", ptm, "/GO enrich/", tstrt)
              if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
              dirlist <- unique(c(dirlist, dir))
              filt <- PTMs_Reg_filters[[Ptm]][[tstrt]]
              #By <- c("By condition", "By reference", "By analysis", "Whole dataset")
              By <- "By condition"
              By <- By[which(By %in% names(filt))]
              if (length(By)) {
                for (bee in By) { #bee <- By[1]
                  flt <- filt[[bee]]
                  if (bee == "Whole dataset") { flt <- list("Whole dataset" = flt) }
                  tstbee <- paste0(tstrt, "_", tolower(bee))
                  if (length(flt)) {
                    if (tt == 1) { tmpdat <- ptmpep }
                    if (tt == 2) { tmpdat <- PTMs_F_test_data[[Ptm]] }
                    for (kol in Kol) { tmpdat[, kol] <- ptmpep[, kol] }
                    flt <- flt[order(names(flt))]
                    reg <- setNames(lapply(flt, function(x) { list(x$Columns) }), names(flt))
                    reg <- set_colnames(reshape2::melt(reg), c("Name", "Bleh", "For"))
                    reg$Bleh <- NULL
                    # PTMs_GO_enrich.FCRt[[Ptm]][[tstbee]] <- paste0("Mean ", c(ptms.ratios.ref[length(ptms.ratios.ref)],
                    #                                                           "log2(Ratio) - ")[tt])
                    PTMs_GO_enrich.FCRt[[Ptm]][[tstbee]] <- gsub("(ReNorm\\. )+", "ReNorm. ",
                                                                 paste0("Mean ",
                                                                        c("", "ReNorm. ")[PTM_normalize[[Ptm]]+1],
                                                                        c(ptms.ratios.ref[length(ptms.ratios.ref)], "log2(Ratio) - ")[tt]))
                    reg$ParentFC <- gsub(".*Regulated - ", PTMs_GO_enrich.FCRt[[Ptm]][[tstbee]], reg$Name)
                    reg$FCname <- paste0(PTMs_GO_enrich.FCRt[[Ptm]][[tstbee]], reg$For)
                    tmpdat[, Kol] <- ptmpep[, Kol]
                    #tmpdat <- get(c("ptmpep", "PTMs_F_test_data[[Ptm]]", "ptmpep", "PTMs_allSAINTs")[tt]) # PTMs_allSAINTs doesn't exist
                    UF <- unique(reg$For)
                    temPTM <- as.data.frame(sapply(UF, function(x) { #x <- UF[1]
                      x <- reg$ParentFC[which(reg$For == x)]
                      if (length(x) > 1) {
                        x <- apply(tmpdat[, x], 1, log_ratio_av)
                      } else { x <- tmpdat[[x]] }
                      return(x)
                    }))
                    colnames(temPTM) <- paste0(PTMs_GO_enrich.FCRt[[Ptm]][[tstbee]], UF)
                    temPTM[, Kol] <- ptmpep[, Kol]
                    #temPTM$"First protein" <- sapply(strsplit(temPTM[[myIDcol]], ";"), function(x) { unlist(x)[1] })
                    if (scrptTypeFull == "withReps_PG_and_PTMs") {
                      temp <- listMelt(strsplit(temPTM$`Protein group IDs`, ";"), temPTM$id)
                      temp$Genes <- PG$Genes[match(temp$value, PG$id)]
                      temp <- aggregate(temp$Genes, list(temp$L1), function(x) {
                        paste(sort(unique(unlist(strsplit(x, ";")))), collapse = ";")
                      })
                    }
                    if (scrptTypeFull == "withReps_PTMs_only") {
                      temp <- listMelt(strsplit(temPTM$Proteins, ";"), temPTM$id)
                      temp$Genes <- db$Gene[match(temp$value, db$`Protein ID`)]
                      temp <- aggregate(temp$Genes, list(temp$L1), function(x) {
                        paste(sort(unique(unlist(strsplit(x, ";")))), collapse = ";")
                      })
                    }
                    temPTM$Genes <- temp$x[match(temPTM$id, temp$Group.1)]
                    PTMs_GO_enrich.dat[[Ptm]][[tstbee]] <- temPTM
                    PTMs_GO_enrich.tbl[[Ptm]][[tstbee]] <- reg
                    flt <- setNames(lapply(UF, function(x) { flt[[x]]$Filter }), UF)
                    ttr <- btr <- ""
                    if (length(By) > 1) { ttr <- btr <- paste0(tolower(bee), "_") }
                    Pep.Ref.Filt <- tmpFilt <- setNames(lapply(names(flt), function(x) { seq_len(nrow(PTMs_GO_enrich.dat[[Ptm]][[tstbee]])) }), names(flt))
                    if (tt == 1) {
                      if (GO.enrich.MultiRefs) {
                        Pep.Ref.Filt <- try(setNames(lapply(names(flt), function(x) { #x <- names(flt)[1]
                          m <- Exp.map[which(Exp.map[[VPAL$column]] == x), , drop = FALSE]
                          y <- unique(m[[GO.enrichment.Ref.Aggr$column]])
                          stopifnot(length(y) == 1)
                          w <- which(Exp.map[[GO.enrichment.Ref.Aggr$column]] == y)
                          w1 <- which(apply(ptmpep[, paste0(ptms.ref, Exp.map$Ref.Sample.Aggregate[w])], 1, function(x) {
                            sum(is.all.good(x, 2))
                          }) > 0)
                          w2 <- flt[[x]] # Required for if we are imputing missing values
                          return(sort(unique(c(w1, w2))))
                        }), names(flt)), silent = TRUE)
                        if ("try-error" %in% class(Pep.Ref.Filt)) {
                          warning("Invalid \"GO.enrichment.Ref.Aggr\" argument: multiple references for GO enrichment are only feasible if each enrichment filter maps to a single reference! Skipping...")
                          GO.enrich.MultiRefs <- FALSE
                          Pep.Ref.Filt <- tmpFilt
                        }
                      }
                    }
                    if (tt == 2) {
                      cM <- as.data.frame(contrMatr_F)
                      Pep.Ref.Filt <- setNames(lapply(names(flt), function(x) { #x <- names(flt)[1]
                        x <- expMap$Ref.Sample.Aggregate[which(expMap$Group_ %in% unlist(expContrasts_F$All[match(x, expContrasts_F$name)]))]
                        kol <- paste0(ptms.ref[1], x)
                        which(parApply(parClust, ptmpep[, kol], 1, function(x) { length(proteoCraft::is.all.good(x)) }) > 0)
                      }), names(flt))
                    }
                    if ((length(Pep.Ref.Filt) > 1)||(!is.na(Pep.Ref.Filt))) {
                      flt <- setNames(lapply(names(flt), function(x) { flt[[x]][which(flt[[x]] %in% Pep.Ref.Filt[[x]])] }), names(flt))
                    }
                    # Also save the reference filters 
                    nms <- names(PTMs_Reg_filters[[Ptm]][[tstrt]][[bee]])
                    for (nm in nms) {
                      PTMs_Reg_filters[[Ptm]][[tstrt]][[bee]][[nm]]$Background_filter <- Pep.Ref.Filt[[nm]]
                    }
                    #
                    ReportCalls <- AddMsg2Report(Msg = " - GO terms enrichment analysis", Space = FALSE)
                    #
                    Mode <- "regulated"
                    dataType <- "modPeptides"
                    if ((!exists("GO_mappings"))&&(file.exists("GO_mappings.RData"))) {
                      loadFun("GO_mappings.RData")
                    }
                    if ((!exists("GO_terms"))&&(file.exists("GO_terms.RData"))) {
                      loadFun("GO_terms.RData")
                    }
                    #
                    Src <- paste0(libPath, "/extdata/R scripts/Sources/GO_enrich.R")
                    #rstudioapi::documentOpen(Src)
                    source(Src, local = FALSE)
                    #
                    clueGO_outDir <- dir
                    clueGO_type <- "Enrichment (Right-sided hypergeometric test)"
                    Src <- paste0(libPath, "/extdata/R scripts/Sources/ClueGO_enrich.R")
                    #rstudioapi::documentOpen(Src)
                    source(Src, local = FALSE)
                    #
                    # Cleanup - do it now, not within sources!
                    try(rm(list = allArgs), silent = TRUE)
                    #
                    PTMs_GO_Plots[[Ptm]][[tstbee]] <- goRES
                    #
                    if (!is.null(names(PTMs_GO_Plots[[Ptm]][[tstbee]]))) {
                      if ("All_GO_terms" %in% names(PTMs_GO_Plots[[Ptm]][[tstbee]])) {
                        GO_terms <- PTMs_GO_Plots[[Ptm]][[tstbee]]$All_GO_terms
                        PTMs_GO_Plots[[Ptm]][[tstbee]]$All_GO_terms <- NULL  
                      }
                      n2 <- names(PTMs_GO_Plots[[Ptm]][[tstbee]]$GO_plots)
                      for (ttl in n2) { #ttl <- n2[1]
                        plot <- PTMs_GO_Plots[[Ptm]][[tstbee]]$GO_plots[[ttl]]
                        ReportCalls <- AddPlot2Report(Space = FALSE, Jpeg = FALSE)
                      }
                      if ((create_plotly)&&(!create_plotly_local)) {
                        plot_ly[[paste0(Ptm, "_GO plots - Regulated vs Observed - ", tstbee)]] <- PTMs_GO_Plots[[Ptm]][[tstbee]]$GO_plot_ly
                      }
                      temp <- PTMs_GO_Plots[[Ptm]][[tstbee]]$GO_terms
                      temp$Mapping <- NULL
                      if ("Offspring" %in% colnames(temp)) {
                        temp$Offspring <- parSapply(parClust, temp$Offspring, paste, collapse = ";")
                      }
                      temp$`Protein table row(s)` <- NULL
                      colnames(temp) <- cleanNms(colnames(temp), start = FALSE)
                      gn <- grep("^Genes", colnames(temp), value = TRUE)
                      pr <- grep("^Proteins", colnames(temp), value = TRUE)
                      pp <- grep("^Peptide IDs", colnames(temp), value = TRUE)
                      kn <- grep("^Count", colnames(temp), value = TRUE)
                      pv <- grep("^Pvalue", colnames(temp), value = TRUE)
                      zs <- grep("^(Z-score|\\(N_Up - N_Down\\)/sqrt\\(Tot\\.\\))", colnames(temp), value = TRUE)
                      lf <- grep("^logFC", colnames(temp), value = TRUE)
                      si <- grep("^Significance", colnames(temp), value = TRUE)
                      #lp <- grep("^Leading protein IDs", colnames(temp), value = TRUE)
                      #kl <- colnames(temp)[which(!colnames(temp) %in% c(gn, kn, pv, zs, lf, si, pp, lp))]
                      kl <- colnames(temp)[which(!colnames(temp) %in% c(gn, kn, pv, zs, lf, si, pp, pr))]
                      #temp <- temp[, c(kl, si, gn, pp, lp, kn, pv, zs, lf)]
                      temp <- temp[, c(kl, si, gn, pp, pr, kn, pv, zs, lf)]
                      w <- apply(temp[, pv, drop = FALSE], 1, function(x) { sum(!is.na(x)) }) > 0
                      temp <- temp[w,]
                      tst <- apply(temp[, kn, drop = FALSE], 1, function(x) { sum(x[which(!is.na(x))]) })
                      temp <- temp[order(tst, decreasing = TRUE),]
                      temp <- temp[order(temp$Ontology, decreasing = FALSE),]
                      PTMs_Reg_GO_terms[[Ptm]][[tstbee]] <- temp
                      write.csv(temp, file = paste0(dir, "/", Ptm, " GO terms - ", tstbee, ".csv"), row.names = FALSE)
                      w <- which(sapply(colnames(temp), function(x) { class(temp[[x]]) }) == "character")
                      if (length(w)) {
                        for (i in w) { #i <- w[1]
                          w1 <- which(nchar(temp[[colnames(temp)[i]]]) > ExcelMax)
                          if (length(w1)) {
                            temp[[colnames(temp)[i]]][w1] <- paste0(substr(temp[[colnames(temp)[i]]][w1], 1, ExcelMax-3), "...")
                          }
                        }
                      }
                      require(openxlsx)
                      HdrStl <- createStyle(textDecoration = "bold", halign = "center", valign = "center",
                                            wrapText = TRUE, numFmt = "TEXT", fontSize = 12)
                      wb <- createWorkbook()
                      kount <- 0
                      for (ont in names(Ontologies)) { #ont <- names(Ontologies)[1]
                        w <- which(temp$Ontology == Ontologies[ont])
                        if (length(w)) {
                          kount <- kount + 1
                          addWorksheet(wb, ont)
                          writeData(wb, ont, temp[w,])
                          setRowHeights(wb, ont, 1, 60)
                          setColWidths(wb, ont, seq_len(ncol(temp)), 12)
                          setColWidths(wb, ont, which(colnames(temp) == "Term"), 45)
                          setColWidths(wb, ont, which(colnames(temp) %in% gn), 20)
                          setColWidths(wb, ont, which(colnames(temp) %in% pr), 20)
                          setColWidths(wb, ont, which(colnames(temp) %in% pp), 20)
                          addStyle(wb, ont, HdrStl, 1, seq_len(ncol(temp)))
                          addStyle(wb, ont, createStyle(numFmt = "0"), 2:(length(w)+1),
                                   which(colnames(temp) %in% kn), gridExpand = TRUE)
                          addStyle(wb, ont, createStyle(numFmt = "0.000"), 2:(length(w)+1),
                                   which(colnames(temp) %in% c(zs, lf, pv)), gridExpand = TRUE)
                        }
                      }
                      if (kount) {
                        saveWorkbook(wb, paste0(dir, "/", Ptm, " GO terms - ", tstbee, ".xlsx"), overwrite = TRUE)
                        if (tt == 1) {
                          Kol2 <- paste0("Significance - ", cleanNms(VPAL$values), " ", max(BH.FDR)*100, "%")
                          Kol2 <- Kol2[which(Kol2 %in% colnames(PTMs_Reg_GO_terms[[Ptm]][[tstbee]]))]
                        }
                        if (tt == 2) {
                          Kol2 <- paste0("Significance - ", names(PTMs_Reg_filters[[Ptm]][[tstrt]][[bee]]), " ", max(BH.FDR)*100, "%")
                          Kol2 <- Kol2[which(Kol2 %in% colnames(PTMs_Reg_GO_terms[[Ptm]][[tstbee]]))]
                        }
                        if (length(Kol2)) {
                          if (length(Kol2) > 1) {
                            w <- which(apply(PTMs_Reg_GO_terms[[Ptm]][[tstbee]][, Kol2], 1, function(x) {"+" %in% x}))
                          } else { w <- which(sapply(PTMs_Reg_GO_terms[[Ptm]][[tstbee]][,Kol2], function(x) {"+" %in% x})) }
                          write.csv(PTMs_Reg_GO_terms[[Ptm]][[tstbee]][w,],
                                    file = paste0(dir, "/", Ptm, " regulated GO terms - ", tstbee, ".csv"),
                                    row.names = FALSE)
                        }
                        # Summary table and heatmap of number of co-regulated GO terms
                        Kol3 <- Kol2
                        N <- length(Kol3)
                        if (N > 1) {
                          temp <- as.data.frame(matrix(rep("", (N+1)^2), ncol = N+1))
                          W <- lapply(Kol3, function(x) { which(PTMs_GO_Plots[[Ptm]][[tstbee]]$GO_terms[[x]] == "+") })
                          names(W) <- gsub(paste0("^Significance - | ", max(BH.FDR)*100, "%$"), "", Kol3)
                          temp[2:(N+1), 1] <- temp[1, 2:(N+1)] <- names(W)
                          for (i in 2:(N+1)) { #i <- 2
                            temp[i, 2:(N+1)] <- sapply(Kol3, function(x) {
                              sum((PTMs_GO_Plots[[Ptm]][[tstbee]]$GO_terms[[x]] == "+")&(PTMs_GO_Plots[[Ptm]][[tstbee]]$GO_terms[[Kol3[i]]] == "+"),
                                  na.rm = TRUE)
                            })
                          }
                          names(W) <- cleanNms(gsub(" [0-9]+%$", "", names(W)))
                          tst <- sapply(strsplit(names(W), " - "), length)
                          tst <- (min(tst) > 1)&(length(unique(tst)) == 1)
                          if (tst) {
                            tst <- as.data.frame(t(sapply(strsplit(names(W), " - "), unlist)))
                            l <- apply(tst, 2, function(x) { length(unique(x)) })
                            tst <- tst[, which(l > 1), drop = FALSE]
                            names(W) <- apply(tst, 1, paste, collapse = " - ")
                          }
                          temp[2:(N+1), 1] <- temp[1, 2:(N+1)] <- names(W)
                          nm <- paste0("N. of co-regulated GO terms\n", tstrt, "\n(", tolower(bee), ")")
                          write.csv(temp, file = paste0(dir, "/", gsub("\n", " - ", gsub("\n\\(", " (", nm)), ".csv"), row.names = FALSE)
                          temp2 <- temp[2:(N+1), 2:(N+1)]
                          colnames(temp2) <- temp[1, 2:(N+1)]
                          rownames(temp2) <-  temp[2:(N+1), 1]
                          for (i in seq_len(nrow(temp2))) { temp2[[i]] <- as.numeric(temp2[[i]]) }
                          if (max(is.all.good(unlist(temp2)))) {
                            temp2 <- as.matrix(temp2)
                            basic.heatmap(temp2, "N. of co-regulated GO terms", paste0(tstrt, "\n(", tolower(bee), ")"),
                                          save = c("pdf", "jpeg"), folder = dir)
                          } else { warning(paste0(tstrt, " co-regulated GO terms heatmap: nothing to plot")) }
                        } else { warning(paste0(tstrt, " performed for only one condition, skipping.")) }
                      }
                    }
                  } else { warning(paste0("Filter ", tstbee, " has length 0, skipping.")) }
                }
              } else { warning(paste0("No filters available for ", tstrt, ", skipping.")) }
            }
          }
          PTMs_pep[[Ptm]] <- ptmpep
          ReportCalls <- AddSpace2Report()
          PTMs_int.ref[[Ptm]] <- ptms.ref
          PTMs_rat.ref[[Ptm]] <- ptms.ratios.ref
          #write.csv(ptmpep, paste0("Tables/", Ptm, "-modified peptides.csv"), row.names = FALSE)
          #w <- grsep2(prot.list, ptmpep$Proteins)
          #if (length(w)) {
          #  write.csv(ptmpep[w,], paste0("Tables/", Ptm, "-modified peptide - proteins in list.csv"), row.names = FALSE)
          #}
        } else {
          warning(paste0("There are no ", ptm, "-modified peptides to report! Check your parameters!"))
        }
      } else {
        warning(paste0("PTM \"", ptm, "\" could not be recognized, please use names consistent with those specified in the search."))
      }
    }
  } else {
    if (scrptTypeFull == "withReps_PTMs_only") {
      stop("Really? There is no PTM-modified class of peptides to analyze? Why did you run this workflow then?")
    }
  }
} else {
  if (scrptTypeFull == "withReps_PTMs_only") {
    stop("Really? There is no PTM-modified class of peptides to analyze? Why did you run this workflow then?")
  }
}

# Once we are done with using GO stuff, backup the final versions
if (Annotate) {
  # (Testing for existence in case we are rerunning part of the script)
  if (exists("GO_terms")) {
    saveFun(GO_terms, file = "GO_terms.RData")
    rm(GO_terms)
  }
  if (exists("GO_mappings")) {
    saveFun(GO_mappings, file = "GO_mappings.RData")
    rm(GO_mappings)
  }
  # otherwise they are also unused further down.
  .obj <- .obj[which(!.obj %in% c("GO_terms", "GO_mappings"))]
}
if (CytoScape) {
  # Then close Cytoscape:
  tst <- try({
    invisible({
      RCy3::closeSession(save.before.closing = FALSE)
      cmd <- paste0("taskkill/im \"", gsub(".+/", "", CytoScVrs), "\" /f")
      #cat(cmd)
      shell(cmd)
    })
  })
  if ("try-error" %in% class(tst)) {
    cat("(This error can be ignored as it does not interrupt script execution.\nIt looks like you are trying to close Cytoscape although it is already open (maybe running this script by bits after an interruption?)\n\n")
  }
}
