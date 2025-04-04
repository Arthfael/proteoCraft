#' Param.load
#'
#' @description 
#' A function to quickly import analysis Parameters into the environment.
#' 
#' @param file The parameters file. Must be a txt or csv file.
#' @param filter.deprecated Should we remove deprecated parameters? Default = TRUE
#' @param WD_detect If set to TRUE (default), will either set work directory to that contained in the parameter file (if there), or will set the parameter value "WD" to the current one.
#' 
#' @examples
#' Param <- Param.load()
#' 
#' @export

Param.load <- function(file = "Parameters.csv", filter.deprecated = TRUE, WD_detect = TRUE) {
  #proteoCraft::DefArg(proteoCraft::Param.load)
  ext <- unlist(strsplit(file, split = "\\."))
  ext <- ext[length(ext)]
  if (ext == "csv") { Param <- read.csv(file, header = FALSE) } else {
    if (ext == "txt") { Param <- read.delim(file, header = FALSE) } else { stop("I only accept txt or csv files") }
  }
  # The commented lines are from older code, which now causes more harm than good
  #if ((ncol(Param) > 2)&&(nrow(Param) > 2)) {
  #  if ((ncol(Param) == 3)&&(nrow(Param) > 3)) {
  #    Param <- Param[,1:2]
  #  } else {
  #    if ((nrow(Param) == 3)&&(ncol(Param) > 3)) {
  #      Param <- Param[1:2,]
  #    } else {
  #      stop("There was a mistake processing this parameters file!")
  #    }
  #  }
  #}
  #if (ncol(Param) == 2) {
  #  a <- Param[,1]
  #  Param <- as.data.frame(t(Param[,2]))
  #  colnames(Param) <- a
  #}
  #if (nrow(Param) == 2) {
  #  a <- Param[1,]
  #  Param <- as.data.frame(Param[2,])
  #  colnames(Param) <- a
  #}
  prm <- as.data.frame(matrix(Param[[2]], nrow = 1))
  Param <- magrittr::set_colnames(prm, Param[[1]])
  if (filter.deprecated) { Param <- Param[,which(!grepl("^deprecated", Param[1,], ignore.case = TRUE)), drop = FALSE] }
  Param[,which(Param[1,] == "T")] <- TRUE
  Param[,which(Param[1,] == "F")] <- FALSE
  if (!is.null(Param$WD)) {
    Param$WD <- gsub("\\\\","/",Param$WD)
    Param$WD <- gsub("^ | $","",Param$WD)
  }
  if (WD_detect) {
    if ((!is.null(Param$WD))&&(dir.exists(Param$WD))) {
      setwd(Param$WD)
    } else {
      tmpwd <- dirname(file)
      if (!grepl("^[A-Z]:", tmpwd)) { tmpwd <- paste0(getwd(), "/", gsub("^\\.$", "", dirname(file))) }
      if (grepl("^[A-Z]:", tmpwd)) { Param$WD <- tmpwd }
    }
  }
  #if (("Scatter" %in% colnames(Param))&&(Param$Scatter)) {
  #  Param.Scatter <- read.csv("Param.Scatter.csv")
  #}
  if ("MQ.Experiments.to.discard" %in% colnames(Param)) { Param$MQ.Experiments.to.discard <- gsub("-| ", ".", Param$MQ.Experiments.to.discard) }
  Param[,which(is.na(Param[1,]))] <- ""
  if ("Type" %in% colnames(Param)) { Param$Type <- toupper(gsub(" ", "", Param$Type)) }
  # Convert characters with numeric values to numerics:
  w <- which(!Param[1,] %in% c(NA, "NA", NaN, "NaN"))
  if (length(w)) {
    w2 <- which(sapply(Param[1, w], function(x) { "character" %in%  class(x) }))
    test1 <- suppressWarnings(as.numeric(Param[1,w[w2]]))
    w3 <- which(!is.na(test1))
    if (length(w3)) {
      w4 <- which(suppressWarnings(as.character(as.numeric(Param[1,w[w2[w3]]]))) == Param[1,w[w2[w3]]])
      if (length(w4)) { Param[,w[w2[w3[w4]]]] <- as.numeric(Param[,w[w2[w3[w4]]]]) }
    }
  }
  # Convert characters to logicals where relevant:
  w <- which(sapply(colnames(Param), function(x) { "character" %in% class(x) })&(toupper(Param[1,]) %in% c("FALSE", "TRUE")))
  if (length(w)) { Param[, w] <- as.logical(Param[, w]) }
  # Test and normalize paths if necessary
  w <- which(sapply(colnames(Param), function(x) { "character" %in% class(Param[[x]]) }))
  w <- w[which(!w %in% match("MQ.Experiments", colnames(Param)))]
  g <- grep("\\\\", Param[1, w])
  if (length(g)) {
    test <- lapply(strsplit(as.character(Param[1, w[g]]), ";"), function(x) {
      #x <- strsplit(as.character(Param[1, w[g]])[1], ";")
      x <- unlist(x)
      tst <- sum(sapply(x, function(y) {
        suppressWarnings(!"try-error" %in% class(try(normalizePath(y, winslash = "/"), silent = TRUE)))
      }))
      if (tst == length(x)) {
        x <- suppressWarnings(normalizePath(x, winslash = "/"))
      }
      return(paste(x, collapse = ";"))
    })
    Param[,w[g]] <- test
  }
  #
  return(Param)
}
