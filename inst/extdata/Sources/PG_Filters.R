#### Optional - Apply True-Discovery or negative filter
# Let's assume that our data is contaminated (e.g. impure fractions),
# but that we have another source to assess the validity of protein groups discoveries.
# Here we can load a table of known valid protein groups, with one column for each value of RG
# Any protein group not matching that table, and its peptides, will have all quantitative values set to NA from here on.
if (DiscFilt) {
  # Create filter
  DiscFiltFilt %<o% strsplit(PG$"Leading protein IDs", ";")
  DiscFiltFilt <- listMelt(DiscFiltFilt, PG$id)
  colnames(DiscFiltFilt) <- c("Leading protein ID", "PG ID")
  w <- which(DiscFiltFilt$"Leading protein ID" %in% DiscFiltTbl$"Protein ID")
  if (DiscFiltMode %in% DiscFiltModes[1L:2L]) {
    # Apply TRUE/FALSE from loaded filter
    for (grp in RG$values) { #grp <- RG$values[1L]
      DiscFiltFilt[[grp]] <- c(FALSE, TRUE)[match(DiscFiltMode, DiscFiltModes)]
      tmp <- DiscFiltTbl[match(DiscFiltFilt$"Leading protein ID"[w], DiscFiltTbl$`Protein ID`), grp]
      if (DiscFiltMode == DiscFiltModes[2L]) { tmp <- !tmp }
      DiscFiltFilt[w, grp] <- tmp
    }
    if (DiscFiltMode == DiscFiltModes[1L]) {
      # We only remove a PG if no leading protein is TRUE in the filter
      DiscFiltFilt <- aggregate(DiscFiltFilt[, RG$values], list(DiscFiltFilt$"PG ID"), function(x) { as.logical(max(x)) })
    }
    if (DiscFiltMode == DiscFiltModes[2L]) {
      # We remove a PG if any leading protein is FALSE in the filter
      DiscFiltFilt <- aggregate(DiscFiltFilt[, RG$values], list(DiscFiltFilt$"PG ID"), function(x) { as.logical(min(x)) })
    }
    colnames(DiscFiltFilt) <- c("PG ID", RG$values)
    DiscFiltFilt <- DiscFiltFilt[match(PG$id, DiscFiltFilt$"PG ID"),] # Re-order
    # Apply filter to quantitative data
    for (grp in RG$values) { #grp <- RG$values[1L]
      w <- which(!DiscFiltFilt[[grp]])
      em <- Exp.map[which(Exp.map[[RG$column]] == grp),]
      kol <- c(paste0(Prot.Expr.Root, em$Ref.Sample.Aggregate),
               grep(topattern(paste0(Prot.Expr.Root, grp, ".REF")), colnames(quantData), value = TRUE),
               paste0(Prot.Rat.Root, em$Ref.Sample.Aggregate),
               grep(topattern(paste0(Prot.Rat.Root, grp, "_REF.to.REF_")), colnames(quantData), value = TRUE))
      kol <- kol[which(kol %in% colnames(quantData))]
      quantData[w, kol] <- NA
    }
    l <- length(DatAnalysisTxt)
    DatAnalysisTxt[l] <- paste0(DatAnalysisTxt[l],
                                " Data was filtered to only include proteins identified in the provided true-discovery filter.")
    ReportCalls <- AddTxt2Report("Removing proteins not included in the TRUE-Discovery filter!")
  }
  if (DiscFiltMode == DiscFiltModes[3L]) {
    if (length(unique(RG$values)) == 1L) {
      DiscFiltCols <- DiscFiltCol
    } else { DiscFiltCols <- paste0(DiscFiltCol, " - ", RG$values) }
    DiscFiltFilt[, RG$values] <- ""
    for (grp in RG$values) { DiscFiltTbl[[grp]] <- c("", "+")[DiscFiltTbl[[grp]]+1L] }
    DiscFiltFilt[w, RG$values] <- DiscFiltTbl[match(DiscFiltFilt$`Leading protein ID`[w], DiscFiltTbl$`Protein ID`), RG$values]
    DiscFiltFilt <- aggregate(DiscFiltFilt[, RG$values], list(DiscFiltFilt$`PG ID`), function(x) {
      c("", "+")[("+" %in% unlist(x))+1L]
    })
    colnames(DiscFiltFilt) <- c("id", DiscFiltCols)
    PG[, DiscFiltCols] <- ""
    w <- which(PG$id %in% DiscFiltFilt$id)
    PG[w, DiscFiltCols] <- DiscFiltFilt[match(PG$id[w], DiscFiltFilt$id), DiscFiltCols]
  }
}
