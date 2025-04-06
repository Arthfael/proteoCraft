#' TopN
#'
#' @description 
#' A function to calculate protein Intensity-based TopN abundance estimator from individual peptide intensities.
#' 
#' @param N How many peptides should be averaged?
#' @param Prot The protein/protein groups file.
#' @param Peptide.IDs The name of the Protein/Protein groups file's peptide IDs column. Default = "Peptide IDs"
#' @param Pep The peptides file.
#' @param id The name of the Peptides file IDs column. Default = "id"
#' @param Pep.Intens.Nms Name of the peptides intensity column(s)
#' @param log.Pep.Intens Set to 0 or FALSE if input peptide intensities are linear scale. If the data is already log scale, set to the relevant base. Default = 10
#' @param Prot.TopN.to.log Should the output protein TopN values be log-scale or not? If set to a number, this is assumed to be the base of the output log scale. Default = 10
#' @param Mods Which modifications should be included? If set to FALSE, will not filter any modifications.
#' @param Mod.Nms Default = "Modified sequence". The name of the column containing the Modified sequence in the peptides file. Can be set to a non-modified sequence if Mods = FALSE
#' @param Min.Pep.Nb How many peptides should at least be present? Should be at the very least 1 (set do default). If greater than N, will be ignored.
#' @param In.Norm Should the Input be re-normalised? Default = FALSE
#' @param Out.Norm Should the Output be re-normalised? Default = TRUE
#' @param Method Which averaging method should be used. The default, "mean", is an arithmetic mean on log values, equivalent to a geometric mean for linear values. Can be instead set to "median", or the name of any averaging function of choice which R knows.
#' @param corr Should cases where there are less than N peptides be corrected? Values can be "local" (default), "global" or FALSE. Local will do the correction for factors infered for each column. Global will average for all columns. FALSE will not do it. Ignored if N = 1
#' 
#' @examples
#' test <- TopN(N = 3, Prot = prot, Peptide.IDs = "Peptide IDs", Pep = pep, id = "id", Pep.Intens.Nms = Pep.Intensity.Names,
#'              log.Pep.Intens = 10, Prot.TopN.to.log = 10,
#'              Mods = Mod4Quant, Mod.Nms = "Modified sequence", Min.Pep.Nb = 1, In.Norm = TRUE, Out.Norm = TRUE, Method = "mean", corr = "local")
#' colnames(test) <- gsub("\\.Intensity\\.", ".", colnames(test))
#' prot[, colnames(test)] <- test
#' 
#' @export

TopN <- function(N = 3,
                 Prot,
                 Peptide.IDs = "Peptide IDs",
                 Pep,
                 id = "id",
                 Pep.Intens.Nms,
                 log.Pep.Intens = 10,
                 Prot.TopN.to.log = 10,
                 Mods,
                 Mod.Nms = "Modified sequence",
                 Min.Pep.Nb = 1,
                 In.Norm = FALSE,
                 Out.Norm = TRUE,
                 Method = "mean",
                 corr = "local") {
  #proteoCraft::DefArg(proteoCraft::TopN)
  #N = 3
  #N = 1
  #Prot = PG;Pep = pep; Pep.Intens.Nms = paste0(pep.ref[length(pep.ref)], Ref.Sample.Aggregate$values); log.Pep.Intens = FALSE; Mods = Mod4Quant; Out.Norm = FALSE; corr = "global"
  if (Min.Pep.Nb > N) { Min.Pep.Nb <- N }
  Method <- get(Method)
  test <- rep(0, nrow(Pep))
  if (sum(Mods != FALSE) != 0) {
    a <- paste(c(paste(c(proteoCraft::AA, "_"), collapse = "|"), sapply(Mods, function(x) {
      paste0("\\(", x, "\\)")
    })), collapse = "|")
    test <- sapply(gsub(a, "", Pep[[Mod.Nms]]), function(x) {nchar(x)})
  }
  pep1 <- Pep[which(!test),]
  log.Pep.Intens <- as.integer(log.Pep.Intens)
  if (!log.Pep.Intens) {
    pep1[, Pep.Intens.Nms] <- suppressWarnings(log10(pep1[, Pep.Intens.Nms]))
    base <- 10
  } else { base <- log.Pep.Intens }
  if (In.Norm) {
    for (i in Pep.Intens.Nms) {
      m <- median(proteoCraft::is.all.good(pep1[[i]]))
      pep1[[i]] <- pep1[[i]] - m
    }
  }
  temp.ids <- strsplit(Prot[[Peptide.IDs]], ";")
  res <- sapply(temp.ids, function(x) { #x <- temp.ids[1]
    x <- which(pep1[[id]] %in% unlist(x))
    if (length(x)) {
      if (length(Pep.Intens.Nms) > 1) {
        temp <- as.data.frame(t(apply(pep1[x, Pep.Intens.Nms], 2, function(y) {
          y <- proteoCraft::is.all.good(y)
          if (length(y) >= N) {
            return(list(sapply(1:N, function(z) { Method(sort(y, decreasing = TRUE)[1:min(length(y), z)]) })))
          } else {
            if (length(y) < Min.Pep.Nb) { return(list(rep(NA, N))) } else {
              return(list(c(rep(paste0("L", length(y)), N - 1), Method(y))))
            }
          }
        })))
      } else {
        temp <- proteoCraft::is.all.good(pep1[x, Pep.Intens.Nms])
        if (length(temp) >= N) {
          temp <- sapply(1:N, function(z) {Method(sort(temp, decreasing = TRUE)[1:min(length(temp), z)])})
        } else {
          if (length(temp) < Min.Pep.Nb) { temp <- list(rep(NA, N)) } else {
            temp <- list(c(rep(paste0("L", length(temp)), N - 1), Method(temp)))
          }
        }
      }
    } else {
      if (length(Pep.Intens.Nms) > 1) {
        temp <- as.data.frame(t(apply(pep1[x, Pep.Intens.Nms], 2, function(y) {return(list(rep(NA, N)))})))
      } else { temp <- list(rep(NA, N)) }
    }
    return(temp)
  })
  if (length(Pep.Intens.Nms) > 1) { res <- as.data.frame(t(res)) } else { res <- as.data.frame(t(t(res))) }
  colnames(res) <- Pep.Intens.Nms
  nk <- max(c(1, N-1))
  temp <- as.data.frame(matrix(rep(0, length(Pep.Intens.Nms)*nk), nrow = length(Pep.Intens.Nms)))
  rownames(temp) <- Pep.Intens.Nms
  colnames(temp) <- 1:nk
  if (N > 1) {
    if (corr == "local") {
      for (o in Pep.Intens.Nms) {
        y <- res[[o]]
        for (n in 1:(N-1)) {
          y1 <- sapply(y, function(x) {unlist(x)[n]})
          y1 <- as.data.frame(y1)
          y1$N <- sapply(y, function(x) {unlist(x)[N]})
          options(warn = -1)
          y2 <- apply(y1, 1, function(x) {
            x <- proteoCraft::is.all.good(as.numeric(x))
            if (length(x) == 2) {return(T)} else {return(F)}
          })
          options(warn = 0)
          y1 <- y1[which(y2 == TRUE),]
          temp[o,n] <- median(as.numeric(y1[,1]) - as.numeric(y1[,2]))
        }
      }
    }
    if (corr == "global") {
      temp2 <- res[,1]
      if (length(Pep.Intens.Nms) > 1) {
        for (i in 2:length(Pep.Intens.Nms)) { temp2 <- c(temp2, res[,i]) }
      }
      y <- temp2
      for (n in 1:(N-1)) {
        y1 <- sapply(y, function(x) {unlist(x)[n]})
        y1 <- as.data.frame(y1)
        y1$N <- sapply(y, function(x) {unlist(x)[N]})
        options(warn = -1)
        y2 <- apply(y1, 1, function(x) { length(proteoCraft::is.all.good(as.numeric(x))) }) == 2
        options(warn = 0)
        y1 <- y1[which(y2),]
        temp[,n] <- median(as.numeric(y1[,1]) - as.numeric(y1[,2]))
      }
    }
  }
  res2 <- sapply(1:ncol(res), function(x) {
    x1 <- temp[which(rownames(temp) == Pep.Intens.Nms[x]),]
    x <- res[,x]
    x <- sapply(1:length(x), function(y) {
      y <- unlist(x[y])
      options(warn = -1)
      y1 <- proteoCraft::is.all.good(as.numeric(y))
      options(warn = 0)
      if (length(y1) == N) { ret <- y[length(y)] } else {
        if (length(which(is.na(y))) == length(y)) { ret <- NA } else {
          options(warn = -1)
          y2 <- unique(y[which(!proteoCraft::is.all.good(as.numeric(y), 2))])
          options(warn = 0)
          y2 <- as.numeric(unlist(strsplit(y2, split = "L"))[2])
          y2 <- x1[y2]
          ret <- y1+y2
        }
      }
      return(ret)
    })
    return(unlist(x))
  })
  if ("numeric" %in% class(res2)) { res2 <- t(as.matrix(res2)) }
  res2 <- as.data.frame(res2)
  if (!Prot.TopN.to.log) { colnames(res2) <- paste0("Top", N, " - ", Pep.Intens.Nms) } else {
    colnames(res2) <- paste0("Top", N, " log", base, " - ", Pep.Intens.Nms)
  }
  if (Out.Norm) {
    for (i in colnames(res2)) {
      a <- median(proteoCraft::is.all.good(unlist(res2[[i]])))
      res2[[i]] <- res2[[i]] - a
    }
  }
  if (!Prot.TopN.to.log) {
    res2 <- base^res2
  } else {
    if (Prot.TopN.to.log != base) {
      res2 <- res2/log(Prot.TopN.to.log, base)
    }
  }
  return(res2)
}
