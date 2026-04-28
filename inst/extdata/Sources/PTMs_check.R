
# In case we have enriched for a PTM, it helps to check how good the enrichment was:
# This should be before any PSMs are filtered out - we want to look at the data "straight out of the MS"
tstEnrich <- unique(FracMap$`PTM-enriched`)
tstEnrich <- tstEnrich[which((!is.na(tstEnrich))&(tstEnrich != "NA"))]
if (length(tstEnrich)) {
  dir <- paste0(wd, "/Summary plots")
  dirlist <- unique(c(dirlist, dir))
  if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
  for (Mod in tstEnrich) { #Mod <- tstEnrich[1L]
    pat <- paste0("\\(", Modifs$Mark[match(Mod, Modifs$`Full name`)], "\\)")
    ev[[Mod]] <- grepl(pat, ev$"Modified sequence")
    for (Type in c("PSMs", "Pep")) { #Type <- "PSMs"
      if (Type == "PSMs") {
        tst <- data.table(Mod = ev[[Mod]],
                          MS_file = ev$"Raw file path",
                          temp = 1L)
        tst <- tst[, list(Count = sum(temp)), by = list(MS_file = MS_file, Mod = Mod)]
        tst <- as.data.frame(tst)
        Root <- "PSM"
      }
      if (Type == "Pep") {
        tst <- data.table(Mod = ev[[Mod]],
                          MS_file = ev$"Raw file path",
                          ModSeq = ev$`Modified sequence`)
        tst <- tst[, list(Count = length(unique(ModSeq))), by = list(MS_file = MS_file, Mod = Mod)]
        tst <- as.data.frame(tst)
        Root <- "peptidoform"
      }
      colnames(tst)[2L] <- Mod
      m <- match(tst$MS_file, FracMap$"Raw file")
      tst$Sample <- FracMap$MQ.Exp[m]
      tst <- tst[which(!is.na(tst$Sample)),]
      tst <- tst[which(tst$Sample %in% unlist(Exp.map$MQ.Exp)),]
      tst <- tst[which(tst$Sample %in% unlist(Exp.map$MQ.Exp)),]
      #which(vapply(Exp.map$MQ.Exp, \(y) { x %in% unlist(y) }, TRUE))
      tst2 <- reshape2::melt(tst)
      colnames(tst2) <- gsub("\\(|\\)|\\[|\\]", "", gsub(" ", "_", colnames(tst2)))
      frml <- as.formula(paste0("MS_file ~ `", gsub("\\(|\\)|\\[|\\]", "", gsub(" ", "_", Mod)), "`"))
      tst2 <- cast(tst2, frml, fun.aggregate = sum)
      kN <- paste0("Non-", Mod, "-modified")
      kY <- paste0(Mod, "-modified")
      colnames(tst2)[which(colnames(tst2) == "FALSE")] <- kN
      colnames(tst2)[which(colnames(tst2) == "TRUE")] <- kY
      tst2[[paste0(Mod, " [%]")]] <- signif(100*tst2[[kY]]/(tst2[[kY]]+tst2[[kN]]), 3L)
      tst2$Sample <- tst$Sample[match(tst2$MS_file, tst$MS_file)]
      tst2 <- tst2[, c("Sample", "MS_file", kN, kY, paste0(Mod, " [%]"))]
      write.csv(tst2, paste0(dir, "/", Mod, "-", Root, "s per MS file.csv"), row.names = FALSE)
      rw <- unique(tst$MS_file)
      for (i in c(TRUE, FALSE)) {
        w <- which(tst[[Mod]] == i)
        w2 <- which(!rw %in% tst$MS_file[w])
        if (length(w2)) {
          tmp <- tst[which((tst$MS_file %in% rw[w2])&(tst[[Mod]] == !i)),]
          tmp[[Mod]] <- i
          tmp$Count <- 0L
          tst <- rbind(tst, tmp)
        }
      }
      tst[[Mod]] <- factor(c("-", "+")[tst[[Mod]]+1L], levels = c("-", "+"))
      ttl <- paste0(Mod, "-", Root, "s per MS file")
      plot <- ggplot(tst) + geom_bar(stat = "identity", position = "dodge",
                                     aes(x = MS_file, y = Count, fill = .data[[Mod]])) +
        theme_bw() + scale_fill_viridis(discrete = TRUE) + ggtitle(ttl) +
        theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1, size = 5L),
              plot.margin = unit(c(0L, 0L, 0L, 3L), "in"))
      #poplot(plot, 12, 20)
      ggsave(paste0(dir, "/", ttl, ".jpg"), plot, dpi = 300L)
      ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300L)
      ReportCalls <- AddPlot2Report()
    }
  }
}

# Write PTMs table
temp <- Modifs
w <- which(vapply(colnames(Modifs), \(x) { inherits(Modifs[[x]], "list") }, TRUE))
for (i in w) { temp[[i]] <- vapply(temp[[i]], paste, "", collapse = ", ") }
dir <- paste0(wd, "/Workflow control")
dirlist <- unique(c(dirlist, dir))
if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
write.csv(temp, paste0(dir, "/Modifications.csv"), row.names = FALSE)
ReportCalls$Calls <- AddTxt2Report("PTMs table:")
ReportCalls$Objects$AABiases <- temp
ReportCalls$Calls <- AddTbl2Report(ReportCalls$Calls, "AABiases")
ReportCalls <- AddSpace2Report()

#
if (LabelType == "LFQ") {
  aggrCol <- "RSA"
  tstMQXp2 <- setNames(vapply(MQ.Exp, \(x) {
    unique(Exp.map$Ref.Sample.Aggregate[tstMQXp[[x]]])
  }, ""), MQ.Exp) # Here for LFQ we should not have multiples, this should fail if it is the case!!!
  ev[[aggrCol]] <- tstMQXp2[ev$MQ.Exp]
  tmp <- ev[, c("Proteins", "Intensity", aggrCol)]
}
if (LabelType == "Isobaric") {
  # Isobaric case - subtle difference
  aggrCol <- "Iso"
  tstMQXp2 <- setNames(vapply(MQ.Exp, \(x) {
    unique(Exp.map$Isobaric.set[tstMQXp[[x]]])
  }, 1L), MQ.Exp)
  tmp <- ev[, c("Proteins", "Intensity")]
  tmp[[aggrCol]] <- tstMQXp2[ev$MQ.Exp]
}
# Plot of contamination levels per sample
dir <- paste0(wd, "/Summary plots")
dirlist <- unique(c(dirlist, dir))
if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
tmp2 <- listMelt(strsplit(tmp$Proteins, ";"), 1L:nrow(tmp), c("Protein", "Row"))
m <- match(tmp2$Protein, db$`Protein ID`)
w <- which(!is.na(m))
tmp2 <- tmp2[w,]; m <- m[w]
# For our purpose here we must match contaminant proteins.
tmp2$Cont <- db$`Potential contaminant`[m]
if (tstOrg) {
  tmp2$Organism <- db[m, dbOrgKol]
  tmp2 <- as.data.table(tmp2)
  f0 <- \(x) { c("Target", "Contaminant")[("Contaminant" %in% x)+1L] }
  tmp2 <- tmp2[, list(x = f0(Organism)), by = list(Group.1 = Row)]
  tmp2 <- as.data.frame(tmp2)
} else {
  tmp2 <- as.data.table(tmp2)
  f0 <- \(x) { c("Target", "Contaminant")[("+" %in% x)+1L] }
  tmp2 <- tmp2[, list(x = f0(Cont)), by = list(Group.1 = Row)]
  tmp2 <- as.data.frame(tmp2)
}
tmp$Organism <- tmp2$x[match(1L:nrow(tmp), tmp2$Group.1)]
tmp <- tmp[which(!is.na(tmp$Organism)),]
tmp$Organism <- factor(tmp$Organism, levels = c("Contaminant", "Target"))
tmp$Intensity <- as.numeric(tmp$Intensity)
tmp <- aggregate(tmp$Intensity, list(tmp[[aggrCol]], tmp$Organism), sum, na.rm = TRUE)
if (LabelType == "LFQ") {
  k <- "Sample"
  colnames(tmp) <- c(k, "Organism", "Total intensity")
  tmp[[k]] <- factor(cleanNms(tmp[[k]]), levels = cleanNms(unique(Exp.map$Ref.Sample.Aggregate)))
}
if (LabelType == "Isobaric") {
  k <- "Isobaric.set"
  colnames(tmp) <- c(k, "Organism", "Total intensity")
  tmp[[k]] <- factor(cleanNms(tmp[[k]]), levels = Iso)
}
ttl <- "Contributions to TIC"
plot <- ggplot(tmp) +
  geom_bar(stat = "identity", aes(x = .data[[k]], y = `Total intensity`, fill = Organism)) +
  theme_bw() + scale_fill_viridis(discrete = TRUE, begin = 0.8, end = 0.2) +
  ggtitle(ttl, subtitle = "Summed TIC for each class of identified peptides") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
print(plot) # This type of QC plot does not need to pop up, the side panel is fine
ggsave(paste0(wd, "/Summary plots/", ttl, ".jpeg"), plot, dpi = 150L, width = 10L, height = 10L, units = "in")
ggsave(paste0(wd, "/Summary plots/", ttl, ".pdf"), plot, dpi = 150L, width = 10L, height = 10L, units = "in")

# Time points
if (exists("Tim")) {
  if (("Time.Points" %in% colnames(Param))&&(Param$Time.Points != "")) {
    tp %<o% as.character(sort(as.numeric(unlist(strsplit(Param$Time.Points, ";")))))
    if (("Time.Point.Names" %in% colnames(Param))&&(Param$Time.Point.Names != "")) {
      names(tp) <- gsub(" ", ".", unlist(strsplit(Param$Time.Point.Names, ";")))
    } else {
      names(tp) <- tp
    }
    if (sum(tp != Tim)) {
      stop("Review your time points!")
    } else {
      Tim <- tp
    }
  } else {
    Tim <- as.character(sort(as.numeric(Tim)))
  }
}

rm(list = ls()[which(!ls() %in% .obj)])
Script <- readLines(ScriptPath)
saveImgFun(BckUpFl)
#loadFun(BckUpFl)
source(parSrc, local = FALSE)
