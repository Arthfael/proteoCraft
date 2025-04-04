#' Format_Mod_Pep
#'
#' @description 
#' A function to correct the evidences table's Modified sequence file, because since MaxQuant 1.6.?.? a new format is being used which prevents my previous scripts from working.
#' Will output the corrected modified sequences as well as a modifications table for inspection.
#' This is a duplicate, currently unused, instead I use the more current cor_mod_seq: I appear to have created 2 functions for the same purpose.
#' Candidate for deletion.
#'
#' @param x The Modified sequences from the evidence file.
#' 
#' @examples
#' test <- Format_Mod_Pep(ev$Modified.sequence)
#' View(test$"Modifications table")
#' ev$Modified.sequence <- test$"Corrected sequences"
#'
#' @export

Format_Mod_Pep <- function(x) {
  x1 <- unique(x)
  x2 <- gsub("^_|_$", "", x1)
  for (a in AA) { x2 <- gsub(a, paste0("-", a, "-"), x2) }
  x2 <- strsplit(x2, "-+")
  y <- unique(unlist(x2))
  y <- y[which(!y %in% AA)]
  y <- unique(unlist(strsplit(gsub("\\(", "-(", gsub("\\)", ")-", y)), "-")))
  y <- y[which(y != "")]
  g <- grepl("^\\([^A-Z]{2}\\)$", y)
  n <- nchar(y)
  w <- which(nchar(y[which(!g)]) > 1)
  if (length(w) > 0) {
    warning("Check your sequences, there are abnormal series of characters!\nWe expect to only see:\n- left- and right- terminal \"_\"\n- capital letters for amino acids\n- modification \"Mod\" indicated by either \"(mo)\" or \"m\" - the aim of this function being to correct the former to the latter!")
  }
  mod <- data.frame(Full = y[which(g)])
  mod$Short <- substr(mod$Full, 2, 2)
  mod$Short[which(! mod$Short %in% y)] <- NA
  if (length(which(!y %in% unlist(mod))) > 0) { warning("Something went awry!") }
  mod$Specificity <- NA
  W <- which(!is.na(mod$Short))
  xtemp <- x1
  for (w in W) {
    #w <- W[1]
    g <- grep(paste(paste0(mod$Short[w], AA), collapse = "|"), xtemp)
    x3 <- strsplit(gsub(mod$Short[w], paste0("-", mod$Short[w]), xtemp[g]), "-")
    x3 <- proteoCraft::Isapply(x3, function(z) {
      z <- unlist(z)
      l <- length(z)
      res <- sapply(z[2:l], function(s) {
        n <- nchar(s)
        s1 <- substr(s, 2, 2)
        if (n > 2) { s2 <- substr(s, 3, n) } else { s2 <- "" }
        s <- paste0(s1, mod$Full[w], s2)
        return(list(c(s, s1)))
      })
      res1 <- paste(c(z[1], sapply(res, function(r) { r[[1]] })), collapse = "")
      res2 <- list(unique(sapply(res, function(r) { r[[2]] })))
      return(c(res1, res2))
    })
    x3a <- unlist(x3[,1])
    if (length(x3a) != nrow(x3)) { stop("Something went awry!") } else {
      xtemp[g] <- x3a
    }
    mod$Specificity[w] <- paste(sort(unique(unlist(x3[,2]))), collapse = ";")
  }
  x <- xtemp[match(x, x1)]
  return(list("Corrected sequences" = x, "Modifications table" = mod))
}
