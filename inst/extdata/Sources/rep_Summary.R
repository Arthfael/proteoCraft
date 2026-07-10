#### Summary table and QC plots
mods <- setNames(Modifs$Mark[which(Modifs$Type == "Variable")],
                 nm = Modifs$"Full name"[which(Modifs$Type == "Variable")])
tmp <- aggregate(Frac.map$"Raw file", list(Frac.map$MQ.Exp), length)
tmp <- round(mean(tmp$x)) # Size of a full fraction set, rounding for cases where we removed some fractions
defSc <- 60L # (Non-strict) default maximum number of files to look at per plot
if (tmp > defSc) {
  # If one set of fractions is larger than defSc
  sc <- tmp
} else {
  # What is closest to default: n or n+1 full sets of fractions?
  tst <- (defSc %% tmp) >= defSc/2L
  # Identify N = fixed number of files from full fraction sets we can fit in one plot
  N <- c(floor(defSc/tmp), ceiling(defSc/tmp))[tst+1L]*tmp
  # If we divide the total number of files by that number, how many plots do we have?
  Nplts <- ceiling(nrow(Frac.map)/N)
  # So that makes how many files per plot:
  sc <- ceiling(nrow(Frac.map)/Nplts)
}
sc <- max(c(sc, 1L))
source(parSrc, local = FALSE)
tmp <- MQ.summary(ev = ev, pg = PG, wd = wd, mods = mods,
                  raw.files = rawFiles, sc = sc, cl = parClust,
                  MQtxt = inDirs[which(SearchSoft == "MAXQUANT")])
Exp_summary %<o% tmp$table
if (!exists("QC_plotLys")) { QC_plotLys %<o% list() }
QC_plotLys[names(tmp$plotLy)] <- tmp$plotLy
Exp_summary$"Biological sample" <- ""
if ("Parent sample" %in% colnames(Frac.map)) {
  Exp_summary$"Biological sample" <- Frac.map$"Parent sample"[match(Exp_summary$Sample, Frac.map$"Raw file")]
}
if ("MQ.Exp" %in% colnames(Frac.map)) {
  Exp_summary$"Biological sample" <- Frac.map$MQ.Exp[match(Exp_summary$Sample, Frac.map$"Raw file")]
}
Exp_summary$"Biological sample"[1L] <- "All samples"
Exp_summary <- Exp_summary[, c("Sample", "Biological sample",
                               colnames(Exp_summary)[which(!colnames(Exp_summary) %in% c("Sample", "Biological sample"))])]
write.csv(Exp_summary, paste0(wd, "/Workflow control/Summary.csv"), row.names = FALSE)
#Exp_summary <- read.csv(paste0(wd, "/Workflow control/Summary.csv"), check.names = FALSE)
