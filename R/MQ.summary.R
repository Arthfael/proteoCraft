#' MQ.summary
#'
#' @description 
#' A function to summarize the data from an MS experiment. Originally limited to MaxQuant searches, hence the name.
#' Now also works partially for the output of DIA-NN if converted to MaxQuant-like format.
#' 
#' @param wd Working directory. Automatically set to local if left empty.
#' @param ev The evidence file. Will try to load from the local directory if value is missing.
#' @param pg The protein groups file. Will try to load from the local directory if value is missing.
#' @param filter Should we filter the files for contaminants, reverse hits, NA values, etc... Default = FALSE
#' @param mods (Ideally named) vector of 2-letter modification codes (as used by older MaxQuant versions for Modified Sequences) of interest. The default is the default MQ PTMs list which we use in most projects.
#' @param raw.files Vector of full raw file paths, sorted in the order you will want to use.
#' @param plot Should we create plot?
#' @param save Set this to a vector of acceptable file extensions to save the graph to the corresponding file format, or to FALSE if you do not want to save it. Default = "jpeg" 
#' @param sc Scale, the max number of raw files to be plotted together in one plot. Default = 60. 
#' @param N.clust A limit on the number of vCPUs to use. If left as NULL (default), uses the number of available clusters - 1, to a minimum of 1.
#' @param N.reserved Default = 1. Number of reserved vCPUs the function is not to use. Note that for obvious reasons it will always use at least one.
#' @param cl Already have a cluster handy? Why spend time making a new one, which on top of that may invalidate the old one. Just pass it along!
#' @param MQtxt For MaxQuant only: the directory in which to search for MaxQuant evidence.txt, msms.txt and msmsScans.txt files.
#' 
#' @examples
#' ## To load data from local directory and filter it:
#' MQ.summary(filter = TRUE)
#' ## To summarize processed data, already loaded in the environment, and look at phospho-peptides
#' MQ.summary(ev = ev, pg = PG, mods = "ph")
#' 
#' @export

MQ.summary <- function(wd, ev, pg, filter = FALSE,
                       mods = setNames(c("ox", "ac", "de", "gl", "ph"),
                                       c("Oxidation (M)", "Acetyl (Protein N-term)", "Deamidation (NQ)", "Gln->pyro-Glu", "Phospho (STY)")),
                       raw.files,
                       subfolder = "Summary plots",
                       plot = TRUE, save = "pdf", sc = 60,
                       N.clust,
                       N.reserved = 1,
                       cl,
                       MQtxt = indir) {
  TESTING <- FALSE
  #proteoCraft::DefArg(proteoCraft::MQ.summary)
  #pg = PG; mods = setNames(Modifs$Mark, Modifs$"Full name"); raw.files = rawFiles; sc = sc
  #pg = PG; mods = setNames(Modifs$Mark, Modifs$"Full name"); raw.files = rawFiles; sc = max(c(20, round(length(rawFiles2)/length(Exp)))); save = c("jpeg", "pdf")
  #TESTING <- TRUE
  if (TESTING) {
    # Note:
    # This is not a perfect alternative to missing but will work in most cases, unless x matches a function imported by a package 
    misFun <- function(x) { return(!exists(deparse(substitute(x)))) }
  } else { misFun <- missing }
  if ((is.null(MQtxt))||(!dir.exists(MQtxt))) { MQtxt <- getwd() }
  wd0 <- getwd()
  if (misFun(wd)) { WD <- wd0 } else { WD <- wd }
  WD <- gsub("/+$", "", normalizePath(WD, "/"))
  if (!is.null(subfolder)) {
    subfolder <- gsub("^/+|/+$", "", gsub("\\\\", "/", subfolder)) # Here do not use normalizePath... try, you will see why ^^
    WD <- paste0(WD, "/", subfolder, "/")
  }
  if (!dir.exists(WD)) { dir.create(WD) }
  sv <- (length(save) > 1)||((length(save) == 1)&(!toupper(as.character(save)) %in% c("FALSE", "NaN", "NA")))
  if (sv) { save <- unique(gsub("^jpg$", "jpeg", gsub("^\\.", "", tolower(save)))) }
  if (misFun(ev)) {
    ev <- proteoCraft::MQ.load(return = TRUE, assign = FALSE, pep = FALSE, prot = FALSE)$evidences
  }
  ev$"Protein group IDs" <- as.character(ev$"Protein group IDs")
  if (misFun(pg)) {
    pg <- proteoCraft::MQ.load(return = TRUE, assign = FALSE, ev = FALSE, pep = FALSE)$protein.groups
  }
  if (filter) {
    W <- colnames(ev)[which(toupper(colnames(ev)) %in% c("CONTAMINANT", "POTENTIAL.CONTAMINANT"))]
    for (w in W) { ev <- ev[which((ev[[w]] != "+")|(is.na(ev[[w]]))),] }
    ev <- ev[which((ev$Reverse != "+")|(is.na(ev$Reverse))),]
    ev <- ev[which((proteoCraft::is.all.good(ev$Intensity, 2))&(ev$Intensity > 0)),]
    W <- colnames(pg)[which(toupper(colnames(pg)) %in% c("CONTAMINANT", "POTENTIAL CONTAMINANT"))]
    for (w in W) { pg <- pg[which((pg[[w]] != "+")|(is.na(pg[[w]]))),] }
    pg <- pg[which(pg$id %in% unique(unlist(strsplit(ev$"Protein group IDs", ";")))),]
  }
  Res <- data.frame(Sample = "Whole dataset",
                    Evidences = nrow(ev),
                    Peptides = length(unique(ev$"Modified sequence")),
                    "Protein groups" = nrow(pg),
                    check.names = FALSE)
  if ("Leading Protein IDs" %in% colnames(pg)) { Res$Proteins <- length(unique(unlist(strsplit(pg$"Leading Protein IDs", ";")))) }
  UseMods <- ((!is.null(mods))&&(length(mods)))
  if (UseMods) {
    if (is.null(names(mods))) { names(mods) <- mods }
    for (i in 1:length(mods)) {
      pat <- paste0("\\(", mods[i], "\\)")
      temp <- grep(pat, ev$"Modified sequence", value = TRUE)
      l <- length(temp); lu <- length(unique(temp))
      Res[[paste0(names(mods)[i], " - evidences")]] <- l
      Res[[paste0(names(mods)[i], " - peptides")]] <- lu
      Res[[paste0(names(mods)[i], " - % ev.")]] <- round(100*l/Res$Evidences[1], 2)
      Res[[paste0(names(mods)[i], " - % pep.")]] <- round(100*lu/Res$Peptides[1], 2)
    }
  }
  if ("Raw file" %in% colnames(ev)) {
    Raw <- data.frame("Raw file" = gsub("\\.((raw)|(mzX?ML)|(d)|(dia))$", "", basename(raw.files)),
                      "Directory" = dirname(raw.files),
                      "Path" = raw.files,
                      "Extension" = sapply(strsplit(raw.files, "\\."), function(x) { rev(unlist(x))[1] }),
                      "Exists" = file.exists(raw.files),
                      check.names = FALSE)
    if (length(which(!Raw$Exists))) {
      warning(paste0("Could not find the following raw files:\n",
                     paste(paste0(" - ", Raw$Path[which(!Raw$Exists)]), collapse = "\n")))
      #Raw$Path <- NULL
    }
    # Use of file name, with or without extension, or full path, is inconsistent. remediate this:
    tstPth <- sapply(1:3, function(x) {
      kol <- c("Raw file", "Path", "Path")[x]
      raws <- Raw[[kol]]
      if (x == 3) { raws <- paste0(raws, ".", Raw$Extension) }
      return(sum(raws %in% ev$`Raw file`))
    })
    L <- length(unique(ev$`Raw file`))
    if (!L %in% tstPth) { warning("Couldn't find raw paths in evidences file.") }
    w <- which(tstPth == L)[1]
    if (!is.na(w)) {
      rawFls <- Raw[[c("Raw file", "Path", "Path")[w]]]
      if (w == 3) { rawFls <- paste0(rawFls, ".", Raw$Extension) }
    } else {
      raw <- unique(ev$"Raw file")
      Raw <- data.frame("Raw file" = gsub("\\.((raw)|(mzX?ML)|(d)|(dia))$", "", basename(raw)),
                        "Directory" = dirname(raw),
                        "Path" = raw,
                        "Extension" = sapply(strsplit(raw, "\\."), function(x) { rev(unlist(x))[1] }),
                        "Exists" = file.exists(raw),
                        check.names = FALSE)
    }
    rawFls <- factor(rawFls, levels = unique(rawFls))
    for (r in as.character(rawFls)) {# r <- as.character(rawFls)[1]
      Res <- rbind(Res, rep(NA, ncol(Res)))
      e <- ev[which(ev$"Raw file" == r),]
      p <- pg[which(pg$id %in% unlist(strsplit(e$"Protein group IDs", ";"))),]
      n <- nrow(Res)
      Res$Sample[n] <- r
      Res$Evidences[n] <- nrow(e)
      Res$Peptides[n] <- length(unique(e$"Modified sequence"))
      if ("Leading Protein IDs" %in% colnames(pg)) {
        Res$Proteins[n] <- length(unique(unlist(strsplit(p$"Leading Protein IDs", ";"))))
      }
      Res$"Protein groups"[n] <- nrow(p)
      if (UseMods) {
        for (i in 1:length(mods)) {
          pat <- paste0("\\(", mods[i], "\\)")
          temp <- grep(pat, e$"Modified sequence", value = TRUE)
          l <- length(temp); lu <- length(unique(temp))
          Res[n, paste0(names(mods)[i], " - evidences")] <- l
          Res[n, paste0(names(mods)[i], " - peptides")] <- lu
          Res[n, paste0(names(mods)[i], " - % ev.")] <- round(100*l/Res$Evidences[n], 2)
          Res[n, paste0(names(mods)[i], " - % pep.")] <- round(100*lu/Res$Peptides[n], 2)
        }
      }
    }
  } else { stop("Expecting a \"Raw file\" column in the evidence file!!!") }
  if (plot) {
    tstVir <- require(viridis)
    # Plot - Peptides composition
    temp <- Res[which(!Res$Sample %in% c("Whole dataset", "")), which(!grepl("%", colnames(Res)))]
    if (nrow(temp)) {
      temp$Proteins <- NULL; temp$"Protein groups" <- NULL
      temp <- reshape::melt.data.frame(temp, id.vars = "Sample")
      if (UseMods) {
        temp$Modification <- sapply(as.character(temp$variable), function(x) {
          x <- unlist(strsplit(x, " - "))
          if (length(x) == 1) { x <- "All" } else { x <- x[1] }
          return(x)
        })
        temp$Modification <- factor(temp$Modification, levels = c("All", names(mods)))
      }
      temp$variable <- apply(cbind(grepl("evidences", temp$variable, ignore.case = TRUE),
                                   grepl("peptides", temp$variable, ignore.case = TRUE)),
                             1, function(x) {c("Evidences", "Peptides")[which(unlist(x))]})
      temp$Sample <- factor(temp$Sample, levels = levels(rawFls))
    }
    if (length(levels(rawFls)) > sc) {
      Raw$iter <- ceiling((1:nrow(Raw))/sc)
      if (nrow(temp)) { temp$iter <- Raw$iter[match(temp$Sample, Raw$"Raw file")] }
    } else {
      Raw$iter <- 1
      if (nrow(temp)) { temp$iter <- 1 }
    }
    iters <- unique(temp$iter)
    for (i in iters) { #i <- 1
      if (length(iters) == 1) { ttl <- "Peptides composition" } else { ttl <- paste0("Peptides composition - ", i) }
      tmp <- temp[which(temp$iter == i),]
      if (UseMods) {
        plot <- ggplot2::ggplot(tmp) +
          ggplot2::geom_col(ggplot2::aes(x = Sample, y = value, fill = Modification),
                            colour = NA) +
          ggplot2::facet_grid(variable~Modification)
      } else {
        plot <- ggplot2::ggplot(tmp) +
          ggplot2::geom_col(ggplot2::aes(x = Sample, y = value, fill = variable),
                            colour = NA) +
          ggplot2::facet_grid(variable~.)
      }
      plot <- plot + ggplot2::ggtitle(ttl) + ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 66, hjust = 1, size = 7))
      if (tstVir) {
        plot <- plot + viridis::scale_fill_viridis(begin = 0.2, discrete = TRUE, option = "D")
      }
      proteoCraft::poplot(plot, 12, 22)
      if (sv) { for (s in save) {
        setwd(WD)
        if (s %in% c("jpeg", "tiff", "png", "bmp")) {
          ggplot2::ggsave(paste0("Summary plots - ", ttl, ".", s), plot,
                 dpi = 300, width = 10, height = 10, units = "in")
        } else {
          ggplot2::ggsave(paste0("Summary plots - ", ttl, ".", s), plot)
        }
        setwd(wd0)
      } }
    } 
    # Plot - Protein coverage
    g <- grep("^Sequence coverage \\[%\\]( - )?", colnames(pg), value = TRUE)
    if (length(g)) {
      temp <- reshape::melt.data.frame(pg[, g])
      temp$Sample <- gsub("^Sequence coverage \\[%\\]( - )?", "", temp$variable)
      temp$Sample[which(temp$Sample == "")] <- "Global"
      temp$Sample <- gsub("___", " ", temp$Sample)
      temp <- aggregate(temp$value, list(round(temp$value), temp$Sample), length)
      colnames(temp) <- c("Sequence coverage [%]", "Sample", "Count")
      temp <- temp[which(!is.na(temp$"Sequence coverage [%]")),]
      temp <- temp[which(temp$"Sequence coverage [%]" > 0),]
      meds <- aggregate(temp$"Sequence coverage [%]", list(temp$Sample), median)
      colnames(meds) <- c("Sample", "Median")
      meds$Label <- paste0("Median = ", signif(meds$Median, 3), "%")
      ttl <- "Sequence coverage"
      sttl <- "(1st protein in group)"
      lev <- unique(temp$Sample)
      lev <- c(lev[which(lev == "Global")], lev[which(lev != "Global")])
      temp$Sample <- factor(temp$Sample, levels = lev)
      meds$Sample <- factor(meds$Sample, levels = lev)
      plot <- ggplot2::ggplot(temp) +
        ggplot2::geom_bar(stat = "Identity",
                          ggplot2::aes(x = `Sequence coverage [%]`, y = Count,
                                       fill = `Sequence coverage [%]`)) +
        ggplot2::geom_text(data = meds,
                           ggplot2::aes(label = Label),
                           x = max(temp$`Sequence coverage [%]`)*0.9,
                           y = max(temp$Count)*0.9, hjust = 1, size = 3) +
        ggplot2::facet_wrap(~Sample) +
        ggplot2::ggtitle(paste0(ttl, " (%)"), sttl) + ggplot2::theme_bw()
      if (tstVir) {
        plot <- plot + viridis::scale_fill_viridis(begin = 0.2, discrete = FALSE, option = "B")
      }
      proteoCraft::poplot(plot, 12, 22)
      if (sv) { for (s in save) {
        setwd(WD)
        if (s %in% c("jpeg", "tiff", "png", "bmp")) {
          ggplot2::ggsave(paste0("Summary plots - ", ttl, ".", s), plot,
                 dpi = 300, width = 10, height = 10, units = "in")
        } else {
          ggplot2::ggsave(paste0("Summary plots - ", ttl, ".", s), plot)
        }
        setwd(wd0)
      } }
    }
    # Plot - Missed cleavages
    if (length(levels(rawFls)) > sc) { ev$iter <- Raw$iter[match(ev$"Raw file", Raw$"Raw file")] } else { ev$iter <- 1 }
    for (i in iters) { #i <- 1
      basettl <- ttl <- "Missed cleavages"
      if (length(iters) > 1) { ttl <- paste0(basettl, " - ", i) }
      w <- which(ev$iter == i)
      e <- ev[w,]
      temp <- data.frame(`Raw file` = rawFls, check.names = FALSE)
      for (n in 0:max(e$"Missed cleavages")) {
        w2 <- which(e$"Missed cleavages" == n)
        temp2 <- aggregate(w2, list(e$"Raw file"[w2]), length)
        temp[[as.character(n)]] <- temp2$x[match(as.character(temp$`Raw file`), temp2$Group.1)]
      }
      temp <- reshape::melt.data.frame(temp, id.vars = "Raw file")
      temp$`Raw file` <- factor(temp$`Raw file`, levels = rawFls)
      colnames(temp)[which(colnames(temp) == "variable")] <- "Missed cleavages"
      colnames(temp)[which(colnames(temp) == "value")] <- "Number of identifications"
      temp <- temp[which(!is.na(temp$"Number of identifications")),]
      plot <- ggplot2::ggplot(temp) +
        ggplot2::geom_bar(stat = "identity",
                          ggplot2::aes(x = `Missed cleavages`, y = `Number of identifications`,
                                       fill = `Missed cleavages`, group = `Missed cleavages`)) +
        ggplot2::facet_wrap(~`Raw file`) + ggplot2::ggtitle(ttl) + ggplot2::theme_bw()
      if (tstVir) {
        plot <- plot + viridis::scale_fill_viridis(begin = 0.2, discrete = TRUE, option = "E")
      }
      proteoCraft::poplot(plot, 12, 22)
      if (sv) { for (s in save) {
        setwd(WD)
        if (s %in% c("jpeg", "tiff", "png", "bmp")) {
          ggplot2::ggsave(paste0("Summary plots - ", ttl, ".", s), plot,
                 dpi = 300, width = 10, height = 10, units = "in")
        } else {
          ggplot2::ggsave(paste0("Summary plots - ", ttl, ".", s), plot)
        }
        setwd(wd0)
      } }
    }
    # Plot - Peptide length
    for (i in iters) { #i <- 1
      basettl <- ttl <- "Peptide length"
      if (length(iters) > 1) { ttl <- paste0(basettl, " - ", i) }
      w <- which(ev$iter == i)
      e <- ev[w,]
      temp <- data.frame("Raw file" = as.character(rawFls[which(rawFls %in% e$`Raw file`)]), check.names = FALSE)
      for (n in min(e$Length):max(e$Length)) {
        w2 <- which(e$Length == n)
        if (length(w2)) {
          temp2 <- aggregate(w2, list(e$"Raw file"[w2]), length)
          temp[[as.character(n)]] <- temp2$x[match(temp$"Raw file", temp2$Group.1)]
        } else { temp[[as.character(n)]] <- 0 }
      }
      temp <- reshape::melt.data.frame(temp, id.vars = "Raw file")
      temp$`Raw file` <- factor(temp$`Raw file`, levels = rawFls)
      colnames(temp)[which(colnames(temp) == "variable")] <- "Length"
      colnames(temp)[which(colnames(temp) == "value")] <- "Number of peptides"
      temp <- temp[which(!is.na(temp$"Number of peptides")),]
      plot <- ggplot2::ggplot(temp) +
        ggplot2::geom_bar(stat = "identity",
                          ggplot2::aes(x = Length, y = `Number of peptides`, fill = Length,
                                       group = Length)) +
        ggplot2::facet_wrap(~`Raw file`) + ggplot2::ggtitle(ttl) + ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(hjust = 1, size = 5))
      if (tstVir) {
        plot <- plot + viridis::scale_fill_viridis(begin = 0.2, discrete = TRUE, option = "D")
      }
      proteoCraft::poplot(plot, 12, 22)
      if (sv) { for (s in save) {
        setwd(WD)
        if (s %in% c("jpeg", "tiff", "png", "bmp")) {
          ggplot2::ggsave(paste0("Summary plots - ", ttl, ".", s), plot,
                 dpi = 300, width = 10, height = 10, units = "in")
        } else {
          ggplot2::ggsave(paste0("Summary plots - ", ttl, ".", s), plot)
        }
        setwd(wd0)
      } }
    }
    # Plot - TIC and Base peak
    if ("Path" %in% colnames(Raw)) {
      we <- which(file.exists(gsub("\\.[^\\.]+$", ".raw", Raw$Path)))
      if (length(we)) {
        #
        # Create cluster
        tstCl <- stopCl <- misFun(cl)
        if (!misFun(cl)) {
          tstCl <- suppressWarnings(try({
            a <- 1
            clusterExport(cl, "a", envir = environment())
          }, silent = TRUE))
          tstCl <- !"try-error" %in% class(tstCl)
        }
        if ((misFun(cl))||(!tstCl)) {
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
        }
        N.clust <- length(cl)
        #
        if (!suppressWarnings(require(rawrr))) {
          pak::pkg_install("rawrr", ask = FALSE)
          #install.packages('http://fgcz-ms.uzh.ch/~cpanse/rawrr_0.2.1.tar.gz', repo = NULL)
          #rawrr::installRawFileReaderDLLs(sourceUrl = rawrr::.thermofisherlsmsUrl()) # Deprecated
          rawrr::installRawrrExe()
        }
        chromtypes <- setNames(c("tic", "bpc"), c("TIC", "Base peak"))
        exports <- list("Raw", "rawFls", "we")
        clusterExport(cl, exports, envir = environment())
        #clusterCall(cl, function() { library(rawrr); return(0) })
        f0 <- function(x) {
          x2 <- rawrr::readChromatogram(Raw$Path[x], type = "tic")
          return(data.frame("Raw file" = Raw$`Raw file`[x],
                            "Raw file path" = rawFls[x],
                            "Retention time" = as.numeric(x2$times),
                            "Intensity" = as.numeric(x2$intensities),
                            check.names = FALSE))
        }
        f1 <- function(x) {
          x2 <- rawrr::readChromatogram(Raw$Path[x], type = "bpc")
          return(data.frame("Raw file" = Raw$`Raw file`[x],
                            "Raw file path" = rawFls[x],
                            "Retention time" = as.numeric(x2$times),
                            "Intensity" = as.numeric(x2$intensities),
                            check.names = FALSE))
        }
        environment(f0) <- .GlobalEnv
        environment(f1) <- .GlobalEnv
        tic <- try(parallel::parLapply(cl, we, f0), silent = TRUE)
        bpc <- try(parallel::parLapply(cl, we, f1), silent = TRUE)
        if (!"try-error" %in% c(class(tic), class(bpc))) {
          for (chrmtp in names(chromtypes)) { #chrmtp <- names(chromtypes)[1]
            temp <- get(chromtypes[chrmtp])
            temp <- plyr::rbind.fill(temp)
            temp$"Raw file path" <- factor(temp$"Raw file path", levels = levels(rawFls))
            temp$Label <- signif(temp$"Retention time", 3)
            winsz <- 0.25
            temp$tst <- FALSE
            temp2 <- temp[, c("Raw file path", "Retention time", "Intensity")]
            clusterExport(cl, list("temp2", "winsz"), envir = environment())
            f2 <- function(x) { #x <- rawFls[1]
              w <- data.frame(Wh = which(temp2$"Raw file path" == x))
              w$tst <- sapply(1:length(w$Wh), function(x) {
                m <- max(temp2$"Retention time"[w$Wh])
                rng <- c(max(c(0, temp2$"Retention time"[w$Wh[x]]-winsz/2)),
                         min(c(temp2$"Retention time"[w$Wh[x]]+winsz/2, m)))
                w2 <- which((temp2$"Retention time"[w$Wh] >= rng[1])&(temp2$"Retention time"[w$Wh] <= rng[2]))
                return(temp2$Intensity[w$Wh[x]] == max(temp2$Intensity[w$Wh][w2]))
              })
              return(w)
            }
            tmp <- parLapply(cl, rawFls, f2)
            tmp <- do.call(rbind, tmp)
            temp$tst <- tmp$tst[match(1:nrow(temp), tmp$Wh)]
            temp$tst <- (temp$tst)&(temp$Intensity > max(temp$Intensity)/10)
            yscl <- max(temp$Intensity)
            if (length(levels(rawFls)) > sc) {
              temp$iter <- Raw$iter[match(temp$"Raw file path", Raw$Path)]
              w <- which(is.na(temp$iter))
              temp$iter[w] <- Raw$iter[match(temp$"Raw file path"[w], Raw$"Raw file")]
            } else { temp$iter <- 1 }
            #aggregate(temp$iter, list(temp$"Raw file"), unique)
            for (i in iters) { #i <- 1
              basettl <- ttl <- chrmtp
              if (length(iters) > 1) { ttl <- paste0(basettl, " - ", i) }
              w <- which(temp$iter == i)
              w2 <- which(temp$tst[w])
              #length(unique(temp$"Raw file"[w]))
              Xmax <- max(temp$`Retention time`)
              Ymax <- max(temp$Intensity)
              Rat <- Xmax/(2*Ymax)
              plot <- ggplot2::ggplot(temp[w,]) +
                ggplot2::geom_segment(ggplot2::aes(x = `Retention time`,
                                                   xend = `Retention time`, yend = Intensity),
                                      y = 0, colour = "grey", size = 0.1) +
                ggplot2::geom_line(ggplot2::aes(x = `Retention time`, y = Intensity),
                                   linewidth = 0.1) +
                ggrepel::geom_text_repel(data = temp[w[w2],],
                                         ggplot2::aes(x = `Retention time`, y = Intensity,
                                                      label = Label),
                                angle = 45, hjust = 0, cex = 2.8, force = 0.0005,
                                direction = "y", ylim = c(1, yscl), min.segment.length = 0) +
                ggplot2::facet_wrap(~`Raw file path`) + ggplot2::ggtitle(ttl) +
                ggplot2::ylab("Intensity") + ggplot2::coord_fixed(Rat) + ggplot2::theme_bw() # +
              #ggplot2::scale_x_discrete(breaks = 1:floor(max(temp$"Retention time"))) # NB: Here use the global values for consistency between iterations
              proteoCraft::poplot(plot, 12, 22)
              if (sv) { for (s in save) {
                setwd(WD)
                if (s %in% c("jpeg", "tiff", "png", "bmp")) {
                  ggplot2::ggsave(paste0("Summary plots - ", ttl, ".", s), plot,
                         dpi = 300, width = 10, height = 10, units = "in")
                } else {
                  ggplot2::ggsave(paste0("Summary plots - ", ttl, ".", s), plot)
                }
                setwd(wd0)
              } }
            }
          }
        }
      }
    }
    # Number of Retention time bins
    nb <- 200
    # MSMS file
    if ("msms.txt" %in% list.files(MQtxt)) {
      msms <- try(data.table::fread(paste0(MQtxt, "/msms.txt"), integer64 = "numeric", check.names = FALSE, data.table = FALSE), silent = TRUE)
      if (class(msms) != "try-error") {
        msms$"Raw file" <- factor(msms$"Raw file", levels = levels(rawFls))
        msms <- msms[which(!is.na(msms$"Raw file")),]
        msms <- msms[order(msms$"Retention time", decreasing = FALSE),]
        msms <- msms[order(msms$"Raw file"),]
        msms$Unique_scan <- gsub("___ +", "___", do.call(paste, c(msms[, c("Raw file", "Scan number")], sep = "___")))
        msms$"Retention time bin" <- as.factor(floor(msms$"Retention time"))
        if (length(levels(rawFls)) > sc) {
          msms$iter <- Raw$iter[match(msms$"Raw file", Raw$"Raw file")]
          w <- which(is.na(msms$iter))
          if ((length(w))&&("Raw file path" %in% names(msms))) {
            msms$iter[w] <- Raw$iter[match(msms$"Raw file path"[w], Raw$"Raw file")]
          }
        } else { msms$iter <- 1 }
        # Plot - Mass error distribution
        for (i in iters) { #i <- 1
          basettl <- ttl <- "Mass error distribution (Da)"
          if (length(iters) > 1) { ttl <- paste0(basettl, " - ", i) }
          wi <- which(msms$iter == i)
          w1 <- which(!is.na(msms$"Mass error [Da]"[wi]))
          plot <- ggplot2::ggplot(msms[wi[w1],]) +
            ggplot2::geom_boxplot(ggplot2::aes(x = `Retention time bin`, y = `Mass error [Da]`),
                                  fill = "orange", outlier.size = 0.1) +
            ggplot2::facet_wrap(~`Raw file`) + ggplot2::ggtitle(ttl) + ggplot2::theme_bw() +
            ggplot2::xlab("Retention time") +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 5))
          if (nb > 20) {
            lev <- levels(msms$"Retention time bin")
            l <- length(lev)
            lev <- lev[which(((1:l) %% 5) == 0)]
            plot <- plot + ggplot2::scale_x_discrete(breaks = lev)
          }
          proteoCraft::poplot(plot, 12, 22)
          if (sv) { for (s in save) {
            setwd(WD)
            if (s %in% c("jpeg", "tiff", "png", "bmp")) {
              ggplot2::ggsave(paste0("Summary plots - ", ttl, ".", s), plot,
                     dpi = 300, width = 10, height = 10, units = "in")
            } else {
              ggplot2::ggsave(paste0("Summary plots - ", ttl, ".", s), plot)
            }
            setwd(wd0)
          } }
        }
      } else { warning("NB: \"msms.txt\" could not be loaded and may be corrupted...") }
    } else { message("NB: \"msms.txt\" could not be found.") }
    #
    # MSMS scans file
    if ("msmsScans.txt" %in% list.files(MQtxt)) {
      msmsScans <- try(data.table::fread(paste0(MQtxt, "/msmsScans.txt"), integer64 = "numeric", check.names = FALSE, data.table = FALSE), silent = TRUE)
      if (class(msmsScans) != "try-error") {
        msmsScans <- msmsScans[order(msmsScans$"Retention time", decreasing = FALSE),]
        msmsScans <- msmsScans[order(msmsScans$"Raw file"),]
        msmsScans$Unique_scan <- gsub("___ +", "___", do.call(paste, c(msmsScans[, c("Raw file", "Scan number")], sep = "___")))
        if ("msms.txt" %in% list.files(MQtxt)) {
          msmsScans$Identified <- c("-", "+")[msmsScans$Unique_scan %in% msms$Unique_scan+1]
        }
        msmsScans$Identified <- factor(msmsScans$Identified, levels = c("-", "+"))
        msmsScans$"Raw file" <- factor(msmsScans$"Raw file", levels = levels(rawFls))
        msmsScans <- msmsScans[which(!is.na(msmsScans$`Raw file`)),]
        msmsScans$"Retention time bin" <- as.factor(floor(msmsScans$"Retention time"))
        bw <- (max(msmsScans$"Retention time")-min(msmsScans$"Retention time"))/nb
        if (length(levels(rawFls)) > sc) {
          msmsScans$iter <- Raw$iter[match(msmsScans$"Raw file", Raw$"Raw file")]
          w <- which(is.na(msmsScans$iter))
          if ((length(w))&&("Raw file path" %in% names(msmsScans))) {
            msmsScans$iter[w] <- Raw$iter[match(msmsScans$"Raw file path"[w], Raw$"Raw file")]
          }
        } else { msmsScans$iter <- 1 }
        # Plot - MSMS identification success
        for (i in iters) { #i <- 1
          basettl <- ttl <- "MSMS identification success"
          if (length(iters) > 1) { ttl <- paste0(basettl, " - ", i) }
          myColors <- setNames(c("grey", "red"), c("-", "+"))
          fillScale <- ggplot2::scale_fill_manual(name = "Identified", values = myColors)
          plot <- ggplot2::ggplot(msmsScans[which(msmsScans$iter == i),]) +
            ggplot2::geom_histogram(ggplot2::aes(x = `Retention time`, fill = Identified),
                                    binwidth = bw) +
            ggplot2::facet_wrap(~`Raw file`) + ggplot2::ggtitle(ttl) + ggplot2::theme_bw() +
            fillScale
          proteoCraft::poplot(plot, 12, 22)
          if (sv) { for (s in save) {
            setwd(WD)
            if (s %in% c("jpeg", "tiff", "png", "bmp")) {
              ggplot2::ggsave(paste0("Summary plots - ", ttl, ".", s), plot,
                     dpi = 300, width = 10, height = 10, units = "in")
            } else {
              ggplot2::ggsave(paste0("Summary plots - ", ttl, ".", s), plot)
            }
            setwd(wd0)
          } }
        }
        # Plot - Precursor apex intensities
        temp <- data.frame(Retention.time = c(min(msmsScans$"Retention time")+bw*c(1:(nb-1)), max(msmsScans$"Retention time")))
        # (NB on above: I have to do it that way to avoid an issue with rounding numbers.)
        for (r in sort(unique(msmsScans$"Raw file"))) { #r <- sort(unique(msmsScans$"Raw file"))[1]
          m <- msmsScans[which(msmsScans$"Raw file" == r),]
          w1 <- which(m$Identified == "-")
          w2 <- which(m$Identified == "+")
          temp[[paste0(r, "___Intensity_Not.identified")]] <- NA
          temp[[paste0(r, "___Intensity_Identified")]] <- NA
          if (length(w1)) {
            tst1 <- sapply(m$"Retention time"[w1], function(x) { temp$Retention.time[which(temp$Retention.time >= x)[1]] })
            temp1 <- aggregate(m$"Precursor intensity"[w1]/m$"Precursor apex fraction"[w1], list(tst1),
                               function(x) { sum(proteoCraft::is.all.good(x)) })
            w1 <- which(temp$Retention.time %in% temp1$Group.1)
            temp[w1, paste0(r, "___Intensity_Not.identified")] <- temp1$x[match(temp$Retention.time[w1], temp1$Group.1)]
          }
          if (length(w2)) {
            tst2 <- sapply(m$"Retention time"[w2], function(x) { temp$Retention.time[which(temp$Retention.time >= x)[1]] })
            temp2 <- aggregate(m$"Precursor intensity"[w2]/m$"Precursor apex fraction"[w2], list(tst2),
                               function(x) { sum(proteoCraft::is.all.good(x)) })
            w2 <- which(temp$Retention.time %in% temp2$Group.1)
            temp[w2, paste0(r, "___Intensity_Identified")] <- temp2$x[match(temp$Retention.time[w2], temp2$Group.1)]
          }
        }
        temp <- reshape::melt.data.frame(temp, id.vars = "Retention.time")
        temp$"Raw file" <- factor(gsub("___Intensity_.+", "", temp$variable), levels = levels(rawFls))
        temp$variable <- gsub(".+___", "", temp$variable) 
        temp$Identified <- c("-", "+")[match(temp$variable, c("Intensity_Not.identified", "Intensity_Identified"))]
        temp$Identified <- factor(temp$Identified, levels = c("-", "+"))
        colnames(temp)[which(colnames(temp) == "value")] <- "Intensity"
        colnames(temp) <- gsub("\\.", " ", colnames(temp))
        if (length(levels(rawFls)) > sc) {
          temp$iter <- Raw$iter[match(temp$"Raw file", Raw$"Raw file")]
          w <- which(is.na(temp$iter))
          temp$iter[w] <- Raw$iter[match(temp$"Raw file path"[w], Raw$"Raw file")]
        } else { temp$iter <- 1 }
        for (i in iters) { #i <- 1
          basettl <- ttl <- "Summed Precursor apex intensities"
          if (length(iters) > 1) { ttl <- paste0(basettl, " - ", i) }
          wi <- which(temp$iter == i)
          w1 <- which(!is.na(temp$Intensity[wi]))
          plot <- ggplot2::ggplot(temp[wi[w1],]) +
            suppressWarnings(ggplot2::geom_histogram(stat = "identity",
                                                     ggplot2::aes(x = `Retention time`,
                                                                  y = Intensity,
                                                                  fill = Identified))) +
            ggplot2::facet_wrap(~`Raw file`) + ggplot2::ggtitle(ttl) + ggplot2::theme_bw() +
            fillScale
          proteoCraft::poplot(plot, 12, 22)
          if (sv) { for (s in save) {
            setwd(WD)
            if (s %in% c("jpeg", "tiff", "png", "bmp")) {
              ggplot2::ggsave(paste0("Summary plots - ", ttl, ".", s), plot,
                     dpi = 300, width = 10, height = 10, units = "in")
            } else {
              ggplot2::ggsave(paste0("Summary plots - ", ttl, ".", s), plot)
            }
            setwd(wd0)
          } }
        }
        # Plot - Number of MSMS per duty cycle
        tst <- as.numeric(msmsScans$"Scan event number") > c(as.numeric(msmsScans$"Scan event number"[2:nrow(msmsScans)]), 1)
        msmsScans$"Scan event number" <- factor(msmsScans$"Scan event number", levels = sort(unique(msmsScans$"Scan event number")))
        msmsScans$`Raw file` <- factor(msmsScans$`Raw file`, levels = rawFls)
        for (i in iters) { #i <- 1
          basettl <- ttl <- "Number of MSMS per duty cycle"
          if (length(iters) > 1) { ttl <- paste0(basettl, " - ", i) }
          plot <- ggplot2::ggplot(msmsScans[which((msmsScans$iter == i)&(tst)),]) +
            ggplot2::geom_bar(ggplot2::aes(x = `Scan event number`,
                                           fill = `Scan event number`)) +
            ggplot2::facet_wrap(~`Raw file`) + ggplot2::ggtitle(ttl) + ggplot2::theme_bw()
          if (tstVir) {
            plot <- plot + viridis::scale_fill_viridis(begin = 0.2, discrete = TRUE, option = "D")
          }
          proteoCraft::poplot(plot, 12, 22)
          if (sv) { for (s in save) {
            setwd(WD)
            if (s %in% c("jpeg", "tiff", "png", "bmp")) {
              ggplot2::ggsave(paste0("Summary plots - ", ttl, ".", s), plot,
                     dpi = 300, width = 10, height = 10, units = "in")
            } else {
              ggplot2::ggsave(paste0("Summary plots - ", ttl, ".", s), plot)
            }
            setwd(wd0)
          } }
        }
        # Plot - AGC fill
        if ("AGC Fill" %in% colnames(msmsScans)) {
          w <- which(!is.na(msmsScans$`AGC Fill`))
          if (length(w)) {
            for (i in iters) { #i <- 1
              basettl <- ttl <- "AGC fill"
              if (length(iters) > 1) { ttl <- paste0(basettl, " - ", i) }
              plot <- ggplot2::ggplot(msmsScans[which(msmsScans$iter == i),]) +
                ggplot2::geom_boxplot(ggplot2::aes(x = `Retention time bin`, y = `AGC Fill`),
                                      fill = "purple", outlier.size = 0.1) +
                ggplot2::ylab(basettl) + ggplot2::facet_wrap(~`Raw file`) +
                ggplot2::ggtitle(ttl) + ggplot2::theme_bw() + ggplot2::xlab("Retention time") +
                ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                                                                   hjust = 1,
                                                                   size = 5)) +
                ggplot2::guides(color = "none")
              if (nb > 20) {
                lev <- levels(msmsScans$"Retention time bin")
                l <- length(lev)
                lev <- lev[which(((1:l) %% 5) == 0)]
                plot <- plot + ggplot2::scale_x_discrete(breaks = lev)
              }
              proteoCraft::poplot(plot, 12, 22)
              if (sv) { for (s in save) {
                setwd(WD)
                if (s %in% c("jpeg", "tiff", "png", "bmp")) {
                  ggplot2::ggsave(paste0("Summary plots - ", ttl, ".", s), plot,
                         dpi = 300, width = 10, height = 10, units = "in")
                } else {
                  ggplot2::ggsave(paste0("Summary plots - ", ttl, ".", s), plot)
                }
                setwd(wd0)
              } }
            }
          }
        }
      } else { warning("NB: \"msmsScans.txt\" could not be loaded and may be corrupted...") }
    } else { message("NB: \"msmsScans.txt\" could not be found.") }
    # Plot - original evidences contamination level
    if ("evidence.txt" %in% list.files(MQtxt)) {
      ev2 <- try(data.table::fread(paste0(MQtxt, "/evidence.txt"), integer64 = "numeric", check.names = FALSE, data.table = FALSE), silent = TRUE)
      if (!"try-error" %in% class(ev2)) {
        ev2$"Raw file" <- factor(ev2$"Raw file", levels = levels(rawFls))
        ev2$Class <- "Sample"
        w1 <- which(ev2$Reverse == "+")
        w2 <- which(ev2$`Potential contaminant` == "+")
        ev2$Class[w1] <- "Reverse"
        ev2$Class[w2] <- "Potential\ncontaminant"
        # (I assume that if an evidence is both reverse and contaminant, the latter takes priority.)
        ev2$Class <- factor(ev2$Class, levels = c("Sample", "Reverse", "Potential\ncontaminant"))
        temp <- aggregate(1:nrow(ev2), list(ev2$`Raw file`, ev2$Class), length)
        colnames(temp) <- c("Raw file", "Class", "Count")
        ttl <- "Evidences QC"
        plot <- ggplot2::ggplot(temp) +
          ggplot2::geom_bar(stat = "identity",
                            ggplot2::aes(x = Class, y = Count, fill = Class)) +
          ggplot2::ggtitle(ttl) + ggplot2::facet_wrap(~ `Raw file`) + ggplot2::theme_bw() +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 66, hjust = 1, size = 7))
        if (tstVir) {
          plot <- plot + viridis::scale_fill_viridis(begin = 0.2, discrete = TRUE, option = "D")
        }
        proteoCraft::poplot(plot, 12, 22)
        if (sv) { for (s in save) {
          setwd(WD)
          if (s %in% c("jpeg", "tiff", "png", "bmp")) {
            ggplot2::ggsave(paste0("Summary plots - ", ttl, ".", s), plot,
                            dpi = 300, width = 10, height = 10, units = "in")
          } else {
            ggplot2::ggsave(paste0("Summary plots - ", ttl, ".", s), plot)
          }
          setwd(wd0)
        } }
      }
    }
    #if ("allPeptides.txt" %in% list.files(MQtxt)) {
    #  require(fGarch)
    #  allpep <- read.delim("allPeptides.txt", check.names = FALSE)
    #  allpep$Retention.times <- apply(allpep[, c("Min scan number", "Max scan number")], 1, function(x) { 
    #    m <- match(x[[1]]:x[[2]], ms$`Scan number`)
    #    m <- m[which(!is.na(m))]
    #    paste(ms$`Retention time`[m], collapse = ";")
    #  })
    #  tst1 <- sapply(strsplit(allpep$Intensities, ";"), length)
    #  tst2 <- sapply(strsplit(allpep$Retention.times, ";"), length)
    #  allpep <- allpep[which(tst1 == tst2),]
    #  psnorm_diff <- function(rt, mean, sd, xi, sampling) {
    #    psnorm(rt+0.5*sampling, mean, sd, xi) - psnorm(rt-0.5*sampling, mean, sd, xi)
    #  }
    #  kount <- 0
    #  for (i in 1:nrow(allpep)) {
    #    temp <- data.frame(RT = as.numeric(unlist(strsplit(allpep$Retention.times[i], ";"))),
    #                       Intensities = as.numeric(unlist(strsplit(allpep$Intensities[i], ";"))))
    #    if (length(which(temp$RT > 0)) > 10) {
    #      Sampling <- mean(temp$RT[2:nrow(temp)] - temp$RT[1:(nrow(temp)-1)])
    #      S <- sum(temp$Intensities)
    #      temp$Rel.intens <- temp$Intensities/S
    #      Mx <- max(temp$Intensities)
    #      M <- mean(temp$RT[which(temp$Intensities == Mx)])
    #      SD <- sd(temp$Rel.intens)
    #      test <- try(nls(Rel.intens ~ psnorm_diff(RT, mean, sd, xi, Sampling), data = temp,
    #                      start = list(mean = M, sd = SD, xi = 1)), silent = TRUE)
    #      if (class(test) == "nls") {
    #        if (exists("fwhm")) { rm(fwhm) }
    #        res <- as.data.frame(summary(test)$parameters)
    #        res <- setNames(res$Estimate, rownames(res))
    #        Mx2 <- psnorm_diff(res[["mean"]], res[["mean"]], res[["sd"]], res[["xi"]], Sampling)*S
    #        extr <- c(min(temp$RT), max(temp$RT))
    #        rg <- data.frame(RT = ((0:100)/100)*(extr[2]-extr[1])+extr[1])
    #        rg$value <- psnorm_diff(rg$RT, res[["mean"]], res[["sd"]], res[["xi"]], Sampling)*S
    #        limz <- suppressWarnings(c(L1 = max(which((rg$value < Mx2/2)&(rg$RT < res[["mean"]]))),
    #                                   L2 = min(which((rg$value > Mx2/2)&(rg$RT < res[["mean"]]))),
    #                                   R1 = max(which((rg$value > Mx2/2)&(rg$RT > res[["mean"]]))),
    #                                        R2 = min(which((rg$value < Mx2/2)&(rg$RT > res[["mean"]])))))
    #        if (length(proteoCraft::is.all.good(limz)) == 4) {
    #          diff <- c(psnorm_diff(mean(rg$RT[limz[1:2]]), res[["mean"]], res[["sd"]], res[["xi"]], Sampling)*S,
    #                    psnorm_diff(mean(rg$RT[limz[3:4]]), res[["mean"]], res[["sd"]], res[["xi"]], Sampling)*S)
    #          if (max(abs((diff-Mx2/2)*2/Mx2)) > 0.01) {
    #            while (max(abs((diff-Mx2/2)*2/Mx2)) > 0.01) {
    #              rg <- data.frame(RT = c(((0:100)/100)*(rg$RT[limz[2]]-rg$RT[limz[1]])+rg$RT[limz[1]],
    #                                      ((0:100)/100)*(rg$RT[limz[4]]-rg$RT[limz[3]])+rg$RT[limz[3]]))
    #                rg$value <- psnorm_diff(rg$RT, res[["mean"]], res[["sd"]], res[["xi"]], Sampling)*S
    #              limz <- suppressWarnings(c(L1 = max(which((rg$value < Mx2/2)&(rg$RT < res[["mean"]]))),
    #                                         L2 = min(which((rg$value > Mx2/2)&(rg$RT < res[["mean"]]))),
    #                                         R1 = max(which((rg$value > Mx2/2)&(rg$RT > res[["mean"]]))),
    #                                         R2 = min(which((rg$value < Mx2/2)&(rg$RT > res[["mean"]])))))
    #              diff <- c(psnorm_diff(mean(rg$RT[limz[1:2]]), res[["mean"]], res[["sd"]], res[["xi"]], Sampling)*S,
    #                        psnorm_diff(mean(rg$RT[limz[3:4]]), res[["mean"]], res[["sd"]], res[["xi"]], Sampling)*S)
    #            }
    #          }
    #          fwhm <- mean(rg$RT[limz[3:4]])-mean(rg$RT[limz[1:2]])
    #          kount <- kount+1
    #          # Neat plot to check (for one iteration only) that the fit is good:
    #          #plot <- ggplot2::ggplot(temp) +
    #          #  ggplot2::geom_line(ggplot2::aes(x = RT, y = Rel.intens)) +
    #          #  ggplot2::geom_vline(xintercept = allpep$`Retention time`[i], colour = "red") +
    #          #  ggplot2::stat_function(fun = function(rt) {
    #          #    psnorm_diff(rt, res[["mean"]], res[["sd"]], res[["xi"]], Sampling)
    #          #  }, colour = "blue") +
    #          #  ggplot2::geom_segment(x = mean(rg$RT[limz[1:2]]), y = Mx2*0.5/S,
    #          #                        xend = mean(rg$RT[limz[3:4]]), yend = Mx2*0.5/S,
    #          #                        xcolour = "green") +
    #          #  ggplot2::geom_text(x = mean(rg$RT[limz[1:2]]) + fwhm/2, y = Mx2*0.525/S,
    #          #                     label = "approx. FWHM", colour = "green") +
    #          #  ggplot2::theme_bw()
    #          #proteoCraft::poplot(plot, 12, 22)
    #          #Add title + save chunk!
    #          if (kount == 1) {
    #            Res <- data.frame(fwhm = fwhm, skewness = res[["xi"]], corr.RT = res[['mean']])
    #          } else { Res[kount,] <- c(fwhm = fwhm, skewness = res[["xi"]], corr.RT = res[['mean']]) }
    #        }
    #      }
    #    }
    #  }
    #}
  }
  setwd(wd0)
  #
  if (stopCl) { parallel::stopCluster(cl) }
  return(Res)
}
