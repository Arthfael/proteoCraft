# Prepare for statistical testing
# At least, Welch's t-test and moderated t-test; for unpaired replicates a permutations t-test is also performed.
# By default in the "t.test" function var.equal is set to FALSE, which performs a Welch's t.test.
# This is as well, since Welch's t.test apparently never performs worse than Student's.
Av_SE_fun %<o% function(vect) {
  res <- proteoCraft::is.all.good(as.numeric(vect))
  if (length(res)) { res <- c(mean(res), sd(res)/sqrt(length(res))) } else {
    res <- unique(vect[which((!is.nan(vect))&(!is.na(vect)))])
    if (length(res) == 1) { return(c(res, NA)) } else { return(c(NA, NA)) }
    return(res)
  }
}
# 
# Internal function nicked from https://github.com/cran/miRtest/blob/master/R/miRtest.R
# Barely modified
limma.one.sided %<o% function(myFit, lower) {
  #se.coef <- sqrt(myFit$s2.post) * myFit$stdev.unscaled
  df.total <- myFit$df.prior + myFit$df.residual
  rs <- pt(myFit$t, df = df.total, lower.tail = lower)
  return(rs[, colnames(myFit$p.value), drop = FALSE])
}
ReportCalls <- AddTxt2Report("Calculating average intensities and ratios and performing statistical tests...")
# Alternative hypothesis for tests:
#AltHyp <- c("two.sided", "greater")[IsPullDown+1]
AltHyp %<o% c(c("greater", "lower")[Mirror.Ratios+1], "two.sided")[TwoSided+1]
# Why use two-sided in most cases?
# For pull-downs, we can still learn something on the left, e.g.:
# - what is enriched specifically in ctrl?
# - see how much variability we see on the left to compare with right
# In addition, it is difficult to implement for the other tests.
#
# Define design matrix and contrasts (limma) 
#  From Exp.map to design matrix
#Coefficients %<o% Factors[which(!Factors %in% c("Experiment", "Replicate"))]
Coefficients %<o% VPAL$names
w <- which(vapply(Coefficients, function(x) { length(unique(FactorsLevels[[x]])) < nrow(Exp.map) }, TRUE))
Coefficients <- Coefficients[w]
w <- which(vapply(Coefficients, function(x) { length(unique(Exp.map[[x]])) > 1 }, TRUE))
Coefficients <- Coefficients[w]
expMap %<o% Exp.map
expMap <- expMap[order(#expMap[[RRG$column]], # Do not use RRG here
  expMap[[RG$column]], # Cf. below: safe because each RG contains at least one ref samples group
  expMap$Reference, # Very important: if any level is dropped from the design matrix, it must be a reference!
  # (Otherwise every gets confusing and  my head starts hurting...)
  expMap[[VPAL$column]],
  expMap$Replicate),]
#
# Replace hyphens by dots to avoid issues with evaluating contrasts
for (Coeff in Coefficients) {
  l <- length(grep("-", Exp.map[[Coeff]])) # Not expMap in case we are re-running a small chunk
  if (l) {
    nuCoeff <- paste0(Coeff, "___")
    stopifnot(!nuCoeff %in% colnames(Exp.map)) # Not expMap in case we are re-running a small chunk
    Coefficients[which(Coefficients == Coeff)] <- nuCoeff
    expMap[[nuCoeff]] <- gsub("-", ".", expMap[[Coeff]])
  }
}
#
Group_ <- do.call(paste, c(expMap[, Coefficients, drop = FALSE], sep = "_"))
Group_ref <- unique(Group_[which(expMap$Reference)])
Group_rst <- unique(Group_[which(!expMap$Reference)])
Group_ <- factor(Group_,
                 levels = c(Group_ref, Group_rst))
expMap$Group_ <- Group_
if (Nested) {
  expMap$Replicate_ <- Replicate_ <- as.factor(expMap$Replicate)
  designMatr %<o% model.matrix(~0 + Replicate_ + Group_)
} else {
  designMatr %<o% model.matrix(~0 + Group_)
}
rownames(designMatr) <- expMap$Ref.Sample.Aggregate
#
# Define contrasts
expContrasts %<o% list()
for (ratGrp in RG$values) { #ratGrp <- RG$values[1]
  em <- expMap[which(expMap[[RG$column]] == ratGrp),]
  grp1 <- unique(em$Group_[which(!em$Reference)])
  vpal1 <- em[match(grp1, em$Group_), VPAL$column]
  grp0 <- unique(em$Group_[which(em$Reference)])
  expContrasts[[ratGrp]] <- plyr::rbind.fill(lapply(grp0, function(g0) {
    data.frame(x1 = paste0("Group_", grp1),
               x0 = paste0("Group_", g0),
               name = vpal1)
  }))
}
expContrasts <- plyr::rbind.fill(expContrasts)
expContrasts$All <- lapply(1:nrow(expContrasts), function(x) { gsub("^Group_", "", expContrasts[x, c("x1", "x0")]) })
expContrasts$Contrasts <- apply(expContrasts[, c("x1", "x0")], 1, function(x) {
  x <- x[which(x %in% colnames(designMatr))]
  paste(x, collapse = " - ")
})
contrCall <- paste0("contrMatr %<o% makeContrasts(",
                    paste(expContrasts$Contrasts, collapse = ", "),
                    ", levels = designMatr)")
#cat(contrCall, "\n")
eval(parse(text = contrCall), envir = .GlobalEnv)
# (NB: Contrasts could be renamed to something shorter, e.g. makeContrasts(Comp1 = A - B, Comp2 = A - C)
#

bhFDRs %<o% sort(BH.FDR, decreasing = FALSE)
samRoots %<o% c(samRoot,
                paste0("SAM regulated-FDR=", paste(100*bhFDRs, collapse = "/"), "% FDR - "))
samSubDir %<o% "Reg. analysis/SAM"
ebamSubDir %<o% "Reg. analysis/EBAM"
ebamRoot %<o% paste0("EBAM regulated-FDR=", paste(100*bhFDRs, collapse = "/"), "% FDR - ")
PolySTestRoot %<o% paste0("EBAM regulated-FDR=", paste(100*bhFDRs, collapse = "/"), "% FDR - ")
#
samDir <- paste0(wd, "/", samSubDir)
ebamDir <- paste0(wd, "/", ebamSubDir)
