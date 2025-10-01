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
  #proteoCraft::DefArg(cleanNms3)
  #
  # For instance:
  #names <- colnames(PG)
  #names <- filtersDF$Name
  #
  nmsBckp <- names
  tstA <- aggregate(1:length(nmsBckp), list(nmsBckp), c)
  g <-  grep("___", names) # Names to fix are the ones where the triad "___" occurs
  # We need to protect the ") - (" from double contrasts, as it is the only allowed occurrence of a space in the names we want to correct:
  # But unfortunately not all occurrences of ") - (" are double contrasts, only ever the last
  # For instance, there are two in "Mean log2(Ratio) - (Exp1___Treatment1) - (Exp1___Treatment2)"
  # However, we know that if it is present once, then we have a double contrast!
  g2c <- grep("\\) - \\(", names[g])
  if (length(g2c)) {
    names[g][g2c] <- vapply(strsplit(names[g][g2c], "\\) - \\("), function(x) { #x <- strsplit(names[g][g2c], "\\) - \\(")[1]
      x <- unlist(x)
      l <- length(x)
      rg <- (l-1):l
      x[l-1] <- paste(x[rg], collapse = "<-VERSUS->")
      x <- x[1:(l-1)]
      return(paste(x, collapse = ") - (")) 
    }, "")
  }
  gSp <- grep(" ", names[g])
  rts <- gsub(" [^ ]+$", "", names[g]) # Roots are defined as anything before the last occurrence of a space
  names2Fix <- gsub(".* ", "", names[g]) # Names are after that space
  #
  # We need to treat separately these stupid .REF and _REF.to.REF_ columns
  gR1 <- grep("\\.REF$", names2Fix)
  gR2 <- grep("_REF\\.to\\.REF_[0-9]+$", names2Fix)
  names2Fix[gR1] <- gsub("\\.REF$", "", names2Fix[gR1])
  num2 <- gsub(".*_REF\\.to\\.REF_", "", names2Fix[gR2])
  names2Fix[gR2] <- gsub("_REF\\.to\\.REF_[0-9]+$", "", names2Fix[gR2])
  #
  # Also deal with double contrast names
  names2Fix <- strsplit(gsub("^\\(|\\)$", "", names2Fix), "<-VERSUS->")
  iNms <- 1:length(names2Fix)
  names2Fix_map <- listMelt(names2Fix, iNms)
  names2Fix2 <- unique(names2Fix_map$value)
  #
  kol <- grep("^([A-Z][a-z]{2})+$", colnames(experiments.map), value = TRUE)
  Aggr <- unique(unlist(aggregate.map$Characteristics))
  w1 <- which(vapply(Aggr, function(x) {
    length(unique(experiments.map[[x]]))
  }, 1) == 1)
  if (!length(w1)) {
    names2Fix2new <- gsub("___", rep, names2Fix2)
  } else {
    names2Fix2new <- vapply(names2Fix2, function(x) { #x <- names2Fix2[1] #x <- rev(names2Fix2)[1]
      tst <- setNames(lapply(kol, function(y) { which(experiments.map[[y]] == x) }),
                      kol)
      tst <- tst[which(vapply(tst, length, 1) > 0)[1]]
      w <- tst[[1]]
      k <- names(tst)
      k <- aggregate.map$Characteristics[[match(k, aggregate.map$Aggregate.Name)]]
      k <- k[which(!k %in% Aggr[w1])]
      k <- RSA$names[which(RSA$names %in% k)] # for order
      return(unique(do.call(paste, c(experiments.map[w, k, drop = FALSE], sep = rep))))
    }, "")
  }
  names2Fix_map$New <- names2Fix2new[match(names2Fix_map$value, names2Fix2)]
  names2Fix_map <- aggregate(names2Fix_map$New, list(names2Fix_map$L1), function(x) {
    x <- unlist(x)
    if (length(x) == 2) { x <- paste0("(", x[1], ") - (", x[2], ")") }
    return(x)
  })
  names2Fix_map <- names2Fix_map[order(names2Fix_map$Group.1),]
  w <- which(iNms %in% names2Fix_map$Group.1)
  names2Fix[w] <- names2Fix_map$x[match(iNms[w], names2Fix_map$Group.1)]
  #
  names2Fix[gR1] <- paste0(names2Fix[gR1], ".REF")
  names2Fix[gR2] <- paste0(names2Fix[gR2], "_REF.to.REF", num2)
  names2Fix <- unlist(names2Fix)
  #
  names[g] <- names2Fix
  names[g][gSp] <- paste0(rts, " ", names[g][gSp])
  # Check for unicity
  tstB <- aggregate(1:length(names), list(names), c)
  tstB$Orig <- nmsBckp[vapply(tstB$x, function(x) { x[[1]] }, 1)]
  tstB <- tstB[match(tstA$Group.1, tstB$Orig),]
  if (!identical(tstA$x, tstB$x)) {
    stop("Values corruption! Simplification is creating identity where there was none in original input!")
  }
  return(names)
}
