#' pepHtmp
#' 
#' @description 
#' Peptides heatmap function - useful for checking normalisations.
#' 
#' @param intProt Vector of protein accessions of interest.
#' @param Pep Peptides data.frame, default = pep
#' @param ref Roof of intensity column names, default = pep.ref[length(pep.ref)]
#' @param dstDir Destination directory.
#' @param ttlRoot Root of the plot title, default = "Peptides log2 heatmap"
#' @param DB Parsed fasta data.base, default = db
#' @param Experiment = Exp
#' @param RSA Ref.Sample.Aggregate object, default = Ref.Sample.Aggregate
#' @param is.log Logical, is the data log-transformed? (default = FALSE)
#' @param cl Already have a cluster handy? Why spend time making a new one, which on top of that may invalidate the old one. Just pass it along!
#' @param N.clust A limit on the number of vCPUs to use. If left as NULL (default), uses the number of available clusters - 1, to a minimum of 1.
#' @param N.reserved Default = 1. Number of reserved vCPUs the function is not to use. Note that for obvious reasons it will always use at least one.
#' 
#' @returns
#' This function does not return anything.
#' 
#' @export

pepHtmp <- function(intProt = prot.list_pep,
                    Pep = pep,
                    ref = pep.ref[length(pep.ref)],
                    dstDir,
                    ttlRoot = "Peptides log2 heatmap",
                    DB = db,
                    Experiment = Exp,
                    RSA = Ref.Sample.Aggregate,
                    is.log = FALSE,
                    cl,
                    N.clust,
                    N.reserved = 1) {
  TESTING <- FALSE
  #TESTING <- TRUE
  if (TESTING) {
    # Note:
    # This is not a perfect alternative to missing but will work in most cases, unless x matches a function imported by a package 
    misFun <- function(x) { return(!exists(deparse(substitute(x)))) }
  } else { misFun <- missing }
  stopifnot(!is.null(intProt), length(intProt) > 0,
            !is.null(Pep), length(Pep) > 0)
  if ((is.null(ttlRoot))||(!"character" %in% class(ttlRoot))||(!nchar(ttlRoot))) {
    warning("Argument \"ttlRoot\" cannot be NULL!")
    ttlRoot <- "Peptides log2 heatmap"
  }
  if (!dir.exists(dstDir)) { dir.create(dstDir, recursive = TRUE) }
  if (exists("dirlist")) { assign("dirlist", unique(c(dirlist, dstDir)), envir = .GlobalEnv) }
  #
  g <- paste0(ref, RSA$values)
  g <- g[which(g %in% colnames(Pep))]
  g1 <- as.data.frame(t(as.data.frame(strsplit(sub(topattern(ref), "", g), "___"))))
  colnames(g1) <- RSA$names
  g1$Col <- g
  test <- RSA$names[which(RSA$names != "Replicate")]
  for (i in rev(test)) { g1 <- g1[order(g1[[i]]),] }
  g <- g1$Col
  g1 <- sub(topattern(ref), "", g)
  g1 <- cleanNms(g1, Experiment = Experiment)
  DB <- DB[match(intProt, DB$"Protein ID"),
           c("Common Name", "Protein ID", "Sequence")]
  tmpPep <- Pep[grsep2(intProt, Pep$Proteins),
                c("Proteins", "Sequence", "Modified sequence", g)]
  cleanUp <- FALSE
  if (nrow(tmpPep)) {
    # Create cluster (some steps are slow otherwise)
    if (misFun(cl)) {
      dc <- parallel::detectCores()
      if (misFun(N.reserved)) { N.reserved <- 1 }
      if (misFun(N.clust)) {
        N.clust <- max(c(dc-N.reserved, 1))
      } else {
        if (N.clust > max(c(dc-N.reserved, 1))) {
          warning("More cores specified than allowed, I will ignore the specified number! You should always leave at least one free for other processes, see the \"N.reserved\" argument.")
          N.clust <- max(c(dc-N.reserved, 1))
        }
      }
      cl <- parallel::makeCluster(N.clust, type = "SOCK")
      cleanUp <- TRUE
    }
    WD <- getwd()
    tmpFl1 <- paste0(wd, "/tmp1.rds")
    tmpFl2 <- paste0(wd, "/tmp2.rds")
    readr::write_rds(DB, tmpFl1)
    readr::write_rds(tmpPep, tmpFl2)
    parallel::clusterExport(cl, list("tmpFl1", "tmpFl2", "dir", "g", "g1", "ttlRoot", "is.log", "dstDir"), envir = environment())
    invisible(parallel::clusterCall(cl, function(x) {
      DB <<- readr::write_rds(tmpFl1)
      tmpPep <<- readr::write_rds(tmpFl2)
      return()
    }))
    l <- length(intProt)
    tempDat <- data.frame(Protein = rep(intProt, 2),
                          Method = c(rep("Mean", l), rep("ZSc", l)))
    f0 <- function(x) { #x <- tempDat[1,]
      plp <- unlist(x[[1]])
      meth <- x[[2]]
      Plp <- paste(DB[which(DB$"Protein ID" == plp),
                      c("Common Name", "Protein ID")], collapse = " - ")
      grs <- grsep2(plp, tmpPep$Proteins)
      if (length(grs)) {
        Seq <- DB$Sequence[match(plp, DB$"Protein ID")]
        temp <- tmpPep[grs, c("Sequence", "Modified sequence", g)]
        temp$tst1 <- vapply(temp$Sequence, function(x) { nchar(unlist(strsplit(Seq, x))[1]) }, 1)
        temp$tst2 <- vapply(temp$Sequence, nchar, 1)
        temp$tst3 <- vapply(temp$"Modified sequence", nchar, 1)
        temp <- temp[order(temp$tst1, temp$tst2, temp$tst3),]
        w <- which(rowSums(temp[, g], na.rm = TRUE) == 0)
        if (length(w)) {
          warning("Peptide(s) found with all invalid intensity values!")
          w <- which(rowSums(temp[, g], na.rm = TRUE) > 0)
          temp <- temp[w,]
        }
        if (nrow(temp)) {
          colnames(temp)[match(g, colnames(temp))] <- g1
          # Create heatmap
          temp <- temp[, c("Modified sequence", g1)]
          if (!is.log) {
            warning(" Converting data to log2...")
            temp[, g1] <- suppressWarnings(log2(temp[, g1]))  
          }
          w <- which(!is.finite(as.matrix(temp[, g1])), arr.ind = TRUE)
          temp[, g1][w] <- NA
          M <- rowMeans(temp[, g1], na.rm = TRUE)
          if (meth == "ZSc") { SD <- apply(temp[, g1], 1, sd, na.rm = TRUE) }
          temp[, g1] <- sweep(temp[, g1], 1, M, "-")
          if (meth == "ZSc") { temp[, g1] <- sweep(temp[, g1], 1, SD, "/") }
          temp2 <- dfMelt(temp[, g1], c("Sample", "value"))
          temp2$"Modified sequence" <- temp$"Modified sequence"
          temp2$Sample <- as.character(temp2$Sample)
          temp2$Xmin <- match(temp2$Sample, colnames(temp))-1
          temp2$Xmax <- temp2$Xmin+1
          temp2$Ymax <- nrow(temp):1
          temp2$Ymin <- temp2$Ymax-1
          hlab <- aggregate(temp2[, c("Xmin", "Xmax")], list(temp2$Sample), unique)
          colnames(hlab)[1] <- "Sample"
          hlab$X <- rowMeans(hlab[, c("Xmin", "Xmax")])
          vlab <- aggregate(temp2[, c("Ymin", "Ymax")], list(temp2$"Modified sequence"), unique)
          colnames(vlab)[1] <- "Modified sequence"
          vlab$Y <- rowMeans(vlab[, c("Ymin", "Ymax")])
          Xscale <- max(temp2$Xmax)
          Yscale <- max(temp2$Ymax)
          # Create heatmap plot
          Xlim <- c(-1, Xscale+15)
          Ylim <- c(-6, Yscale+20)
          temp2a <- temp2[, c("Xmin", "Ymin", "value")]
          Splits <- 20
          XScale2 <- Xscale*0.1/Splits
          temp2b <- data.frame(Xmin = Xscale/2 - XScale2*((-Splits/2):(Splits/2)),
                               Ymin = -5)
          temp2b$Xmax <- temp2b$Xmin + XScale2
          temp2b$value <- min(temp2a$value, na.rm = TRUE) + (Splits:0)*(max(temp2a$value, na.rm = TRUE)-min(temp2a$value, na.rm = TRUE))/Splits
          temp2b$Label <- ""
          temp2b$Label[c(1, Splits/2+1, Splits+1)] <- round(temp2b$value[c(1, Splits/2+1, Splits+1)], 1)
          # Create graph
          ttl <- paste0(ttlRoot, " - ", gsub("/", "-", Plp), c("", " (Z-scored)")[(meth == "ZSc")+1])
          heatmap.plot <- ggplot2::ggplot(temp2a) +
            ggplot2::geom_rect(ggplot2::aes(xmin = Xmin, xmax = Xmin+1, ymin = Ymin, ymax = Ymin+1, fill = value)) +
            ggplot2::geom_rect(data = temp2b, ggplot2::aes(xmin = Xmin, xmax = Xmax, ymin = Ymin, ymax = Ymin+1, fill = value)) +
            ggplot2::geom_text(data = hlab, ggplot2::aes(x = X, label = Sample), y = Yscale+1, angle = 60, hjust = 0, size = 2.5) +
            ggplot2::geom_text(data = vlab, ggplot2::aes(y = Y, label = `Modified sequence`), x = Xscale+0.5, hjust = 0, size = 1.7) +
            ggplot2::geom_text(data = temp2b, ggplot2::aes(x = Xmin, label = Label), y = -3.5, hjust = 0.5, size = 2) +
            ggplot2::ggtitle(paste0(ttlRoot, ", ", c("mean-normalized", "Z-scored")[(meth == "ZSc")+1]),
                             subtitle = Plp) +
            ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(), 
                           axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(), 
                           panel.background = ggplot2::element_rect(fill = "transparent", color = NA),
                           plot.margin = ggplot2::margin(0, 0, 0, 0, "cm")) +
            ggplot2::coord_fixed(0.5) +
            ggplot2::scale_fill_gradient2(low = "red", mid = "black", high = "green", na.value = "lightblue") +
            ggplot2::xlab(NULL) + ggplot2::ylab(NULL) + ggplot2::theme(legend.position = "none") +
            ggplot2::xlim(Xlim[1], Xlim[2]) + ggplot2::ylim(Ylim[1], Ylim[2])
          #poplot(heatmap.plot)
          suppressMessages({
            ggplot2::ggsave(paste0(dstDir, "/", ttl, ".jpeg"), heatmap.plot,
                            dpi = 600, width = 20, height = 12, units = "in")
            ggplot2::ggsave(paste0(dstDir, "/", ttl, ".pdf"), heatmap.plot,
                            dpi = 600, width = 20, height = 12, units = "in")
          })
          #system(paste0("open \"", dstDir, "/", ttl, ".jpeg", "\""))
          #system(paste0("open \"", dstDir, "/", ttl, ".pdf", "\""))
        }
      }
      return()
    }
    environment(f0) <- .GlobalEnv
    invisible(parallel::parApply(cl, tempDat, 1, f0))
  }
  if (cleanUp) { parallel::stopCluster(cl) }
  return()
}
