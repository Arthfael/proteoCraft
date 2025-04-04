#' cleanNms3
#'
#' @description
#' Meant as an internal function to clean up "___"-aggregated names for nicer plot annotations.
#' (NB: should work for other separators, we are just using "___" here).
#' For instance, it will turn "Exp1___Treated___KO___1" and "Exp1___Treated___KO" into the much user friendlier "Exp1 Treated KO 1" and "Exp1 Treated KO",
#' or, if - as is usually the case - there is only one Experiment, into "Treated KO 1" and "Treated KO".
#' 
#' Should be an improvement over cleanNms, which has the same purpose but a different code.
#' With this one it should be possible to feed a vector of data, example the colnames from a data.frame, and to fix all these names at once.
#' 
#' This is called cleanNms3 because cleanNms2 is a different function with a similar purpose used in the workflows for Venn Diagrams.
#' 
#' @param names The vector of sample names, or of column names containing sample names. Each name should be made of "___"-aggregated factor levels.
#' @param rep Replacement, user-friendlier separator. Default = " "
#' @param experiments.map The experiments map. Default = Exp.map
#' @param aggregate.map The aggregate map. Default = Aggregate.map
#' @param aggregate.list The named list of aggregates. Default = Aggregate.list
#' @param RSA List object describing individual samples and how they relate to the experimental factors. Default = Ref.Sample.Aggregate
#' 
#' @export

cleanNms3 <- function(names,
                      rep = " ",
                      experiments.map = Exp.map,
                      aggregate.map = Aggregate.map,
                      aggregate.list = Aggregate.list,
                      RSA = Ref.Sample.Aggregate) {
  warning("I CANNOT DEAL YET WITH DOUBLE CONTRAST NAMES CONTAINING A VS! HELP ME DEAL WITH THOSE!!!")
  # E.g.:
  #names <- colnames(PG)
  namesBck <- names
  tstA <- aggregate(1:length(namesBck), list(namesBck), c)
  g <-  grep("___", names)
  rts <- gsub(" [^ ]+$", "", names[g])
  names2Fix <- gsub(".* ", "", names[g])
  gR1 <- grep("\\.REF$", names2Fix)
  gR2 <- grep("_REF\\.to\\.REF_[0-9]+$", names2Fix)
  names2Fix[gR1] <- gsub("\\.REF$", "", names2Fix[gR1])
  num2 <- gsub(".*_REF\\.to\\.REF_", "", names2Fix[gR2])
  names2Fix[gR2] <- gsub("_REF\\.to\\.REF_[0-9]+$", "", names2Fix[gR2])
  kol <- grep("^([A-Z][a-z]{2})+$", colnames(experiments.map), value = TRUE)
  #
  Aggr <- unique(unlist(aggregate.map$Characteristics))
  w1 <- which(sapply(Aggr, function(x) {
    length(unique(experiments.map[[x]]))
  }) == 1)
  if (!length(w1)) {
    names2Fix <- gsub("___", rep, names2Fix)
  } else {
    names2Fix <- sapply(names2Fix, function(x) { #x <- names2Fix[1]
      tst <- setNames(lapply(kol, function(y) { which(experiments.map[[y]] == x) }), kol)
      tst <- tst[which(sapply(tst, length) > 0)[1]]
      w <- tst[[1]]
      k <- names(tst)
      k <- aggregate.map$Characteristics[[match(k, aggregate.map$Aggregate.Name)]]
      k <- k[which(!k %in% Aggr[w1])]
      k <- RSA$names[which(RSA$names %in% k)] # for order
      return(paste(experiments.map[w, k], collapse = rep))
    })
  }
  names2Fix[gR1] <- paste0(names2Fix[gR1], ".REF")
  names2Fix[gR2] <- paste0(names2Fix[gR2], "_REF.to.REF", num2)
  names[g] <- paste0(rts, " ", names2Fix)
  # Check for unicity
  tstB <- aggregate(1:length(names), list(names), c)
  tstB$Orig <- namesBck[sapply(tstB$x, function(x) { x[[1]] })]
  tstB <- tstB[match(tstA$Group.1, tstB$Orig),]
  if (!identical(tstA$x, tstB$x)) {
    stop("Values corruption! Simplification is creating identity where there was none in original input!")
  }
  return(names)
}
