# Create list of control ratio values for the purpose of identifying vertical thresholds for plots:
## This should be, for each volcano plot, a list of ratios with the same name.
# Note: It's taken me forever to get this right, but I think this should be ok now!!!
if (Param$Ratios.Thresholds == threshMsg) {
  warning(paste0("Parameter Ratios.Thresholds = \"", threshMsg, "\" is deprecated!"))
}
if (Param$Ratios.Thresholds == "Absolute log2 FC threshold") {
  plotMetr <- as.data.frame(strsplit(unlist(strsplit(Param$Plot.threshold.metrics, ";")), ":"))
  plotMetr <- as.data.frame(t(plotMetr)) 
  rownames(plotMetr) <- NULL
  colnames(plotMetr) <- c("Levels", "Axis")
  a2 <- set_colnames(as.data.frame(t(sapply(strsplit(unlist(strsplit(Param$Plot.threshold.values, split = "; *")), split = ": *"), unlist))),
                     c("Direction", "Text.value"))
  plotMetr$Text.value <- a2$Text.value[match(plotMetr$Levels, a2$Direction)]
  w <- which(plotMetr$Axis == "X")
  m <- w[match(c("down", "up"), plotMetr$Levels[w])]
  plotMetr$Text.value[w] <- as.character(c(-Param$Ratios.Contamination.Rates, Param$Ratios.Contamination.Rates))
  Param$Plot.threshold.values <- do.call(paste, c(plotMetr[, c("Levels", "Text.value")], sep = ": ", collapse = ";"))
  Ref.Ratios %<o% NULL
}
if (Param$Ratios.Thresholds == threshMsg) {
  stop("This option is deprecated!")
  # Ref.Ratios %<o% setNames(lapply(VPAL$values, \(x) { #x <- VPAL$values[1L]
  #   if (RatConGrps == "Ratio groups") {
  #     x1 <- unique(Exp.map[which(Exp.map[[VPAL$column]] == x), RG$column])
  #   }
  #   if (RatConGrps == "Experiments") {
  #     x1 <- unique(Exp.map$Experiment[which(Exp.map[[VPAL$column]] == x)])
  #     x1 <- unique(Exp.map[which(Exp.map$Experiment == x1), RG$column])
  #   }
  #   if (RatConGrps == "Whole dataset") {
  #     x1 <- unique(Exp.map[[RG$column]])
  #   }
  #   x <- unique(Exp.map[which(Exp.map[[VPAL$column]] == x), RG$column])
  #   x <- grep(paste0(topattern(paste0(Prot.Rat.Root, x1, "_REF.to.REF_")), "[0-9]+"), colnames(quantData), value = TRUE)
  #   x <- if (length(x)) { is.all.good(as.numeric(unlist(quantData[, x]))) } else { NULL }
  #   return(x)
  # }), VPAL$values)
}
