#' Prot.Ruler
#' 
#' @description 
#' In house implementation of the Proteomic ruler (https://www.mcponline.org/content/13/12/3497.long), with minor modifications.
#' Values output are estimates of copy number per cell.
#'
#' @param Prot # The Protein Groups table.
#' @param DB # The processed Fasta proteome database.
#' @param Expr.roots # The roots of the expression columns which should be processed. Default = c("Top3.log10.Intensity.", "Mean.Top3.log10.Intensity.")
#' @param Proteins.col Default = "Protein.IDs"
#' @param NuclL Nucleosome length + estimate of linker region length; default = 196 (146 + 50)
#' @param log.Expr Set to 0 or FALSE if input expression values are linear scale. If the data is already log scale, set to the relevant scale's base. Default = 10
#' @param norm Are expression columns normalised? If TRUE (default), then we will use as estimates the data for histones averaged over all columns. If FALSE, each expression column will be treated independently.
#' @param norm_filter If TRUE (default), and all core histones are missing from a sample, then that sample will be excluded from the results.
#' @param MaxOrg Sometimes, because of contaminants, there are a lot of organisms in an experiment. This is the maximum number (default = 3) of organisms that the function will try to calculate a value for. If there are more, it will prompt the user for which to include.
#' @param getLatest If TRUE, will attempt to get genome sizes from the latest version of the online file at "ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/overview.txt"; if FALSE (default), will use the version re-distributed with the package.
#' 
#' @examples
#' PG <- Prot.Ruler(PG, db)
#' 
#' @export

Prot.Ruler <- function(Prot,
                       DB,
                       Expr.roots = c("Top3 log10(Intensity) ", "Mean Top3 log10(Intensity) "),
                       Proteins.col = "Protein IDs",
                       NuclL = 196,
                       log.Expr = 10,
                       norm = TRUE,
                       norm_filter = TRUE,
                       MaxOrg = 3,
                       getLatest = FALSE) {
  #proteoCraft::DefArg(proteoCraft::Prot.Ruler)
  #Prot <- PG; DB <- db; Expr.roots = c(Prot.Expr.Root, paste0("Mean ", Prot.Expr.Root))
  #Prot <- PG; DB <- db; Expr.roots = ref
  #Prot <- temp; DB <- db; Expr.roots = exprsRt; NuclL = ProtRulNuclL
  PrRulerRoot <- "log10(est. copies/cell) - "
  if (is.logical(log.Expr)) {
    if (!log.Exp) { PrRulerRoot <- "est. copies/cell - " }
    log.Expr.base <- 10
  } else { log.Expr.base <- log.Expr }
  for (i in 1:length(Expr.roots)) { #i <- 1
    kol <- grep(proteoCraft::topattern(Expr.roots[i]), colnames(Prot), value = TRUE)
    kol <- kol[which(!grepl("\\.REF$", kol))]
    # Log-transform now if data wasn't, we will de-log later
    if ((is.logical(log.Expr))&&(!log.Exp)) { Prot[, kol]  <- log(Prot[, kol], log.Expr.base) }
    PRcol <- gsub(proteoCraft::topattern(Expr.roots[i]), PrRulerRoot, kol)
    Samples <- gsub(proteoCraft::topattern(Expr.roots[i]), "", kol)
    Prot[, PRcol] <- NA
    temp <- data.frame(Root = rep(Expr.roots[i], length(kol)), Expr.cols = kol, Pr.Ruler.cols = PRcol, Sample = Samples)
    if (i == 1) { Expr.cols <- temp } else { Expr.cols <- rbind(Expr.cols, temp) }
  }
  if (getLatest) {
    url <- "ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/overview.txt"
    orgmap <- try(read.delim(url), silent = TRUE)
    if ("try-error" %in% class(orgmap)) {
      orgmap <- try(curl::curl_download(url, "tmpOrg.txt"), silent = TRUE)
      if ("try-error" %in% class(orgmap)) {
        warning("Failed to get the latest version of Genome Reports from NCBI, defaulting to re-distributed version...")
        getLatest <- FALSE
      } else {
        orgmap <- read.delim("tmpOrg.txt")
        unlink("tmpOrg.txt")
      }
    }
  }
  if (!getLatest) {
    RPath <- as.data.frame(library()$results)
    RPath <- normalizePath(RPath$LibPath[match("proteoCraft", RPath$Package)], winslash = "/")
    libPath <- paste0(RPath, "/proteoCraft")
    orgmap <- read.delim(paste0(libPath, "/extdata/NCBI_GENOME_REPORTS_overview.txt"))
  }
  orgmap$Linnean_name <- gsub("\\'.*$", "", gsub("^\\'", "", orgmap$X.Organism.Name))
  if ("Organism_Full" %in% colnames(DB)) {
    Orgs <- unique(unlist(strsplit(Prot[[Proteins.col]], ";")))
    Orgs <- Orgs[which(!is.na(Orgs))]
    Orgs <- DB$Organism_Full[match(Orgs, DB$`Protein ID`)]
    Orgs <- aggregate(Orgs, list(Orgs), length)
    Orgs <- data.frame(Organism = Orgs$Group.1[order(Orgs$x, decreasing = TRUE)])
    if ("TaxID" %in% colnames(DB)) {
      Orgs$TaxID <- sapply(Orgs$Organism, function(o) { unique(DB$TaxID[match(o, DB$Organism_Full)]) })
    }
    Orgs[, c("Kingdom", "Genome size")] <- orgmap[match(Orgs$Organism, orgmap$Linnean_name),
                                                  c("Kingdom", "Size..Mb.")]
    Orgs$"Genome size" <- as.numeric(Orgs$"Genome size")*10^6
    Orgs <- Orgs[which(!is.na(Orgs$"Genome size")),]
    w <- which(Orgs$Kingdom != "Eukaryota")
    if (length(w)) {
      if (length(w) == 1) {
        msg <- paste0("Skipping non-eukaryotic organism \"", Orgs$Organism[w], "\"!")
      } else {
        msg <- paste0("Skipping non-eukaryotic organisms ", paste(paste0("\"", Orgs$Organism[w[1:(length(w)-1)]], "\""), collapse = ", "),
                      " and \"", Orgs$Organism[w[length(w)]], "\"!")
      }
      warning(msg)
      Orgs <- Orgs[which(Orgs$Kingdom == "Eukaryota"),]
    }
    if (nrow(Orgs) > MaxOrg) {
      opt <- vapply(Orgs$Organism, function(x) { paste(c(x, rep(" ", max(c(1, 250-nchar(x))))), collapse = "") }, "")
      orgChc <- svDialogs::dlg_list(opt, opt[1], TRUE, "Choose organism(s) for which you want to calculate a Proteomic Ruler value")$res
      orgChc <- Orgs$Organism[match(orgChc, opt)]
      Orgs <- Orgs[which(Orgs$Organism %in% orgChc),]
    }
  } else {
    Orgs <- dlgInput("Enter the Linnaean name of the relevant (Eukaryotic) organism:", "Homo sapiens")$res
    m <- match(Orgs, orgmap$Linnean_name)
    Kngdm <- orgmap$Kingdom[m]
    while (Kngdm != "Eukaryota") {
      Orgs <- dlgInput("Enter the Linnaean name of the relevant organism:\nThis must be a Eukaryote:\n - Bacteria do not have histones.\n - Archaea do, but inter-nucleosome distance does not seem to be as consistent than in Eukaryotes.\n", "Homo sapiens")$res
      m <- match(Orgs, orgmap$Linnean_name)
      Kngdm <- orgmap$Kingdom[m]
    }
    Orgs <- data.frame(Organism = Orgs, Kingdom = Kngdm)
    Orgs$TaxID <- sapply(Orgs$Organism, function(o) { unique(DB$TaxID[grep(o, DB$Header)]) })
    test <- sapply(Orgs$TaxID, length) != 1
    if (sum(test)) {
      warning("Some of the provided Linnaean names could not be found in the database or have duplicate TaxIDs!")
      Orgs <- Orgs[which(!test), , drop = FALSE]
    }
    DB$Organism_Full <- Orgs$Organism[match(DB$TaxID, Orgs$TaxID)]
    Orgs$"Genome size" <- as.numeric(orgmap$Size..Mb.[m])*10^6
    Orgs <- Orgs[which(!is.na(Orgs$"Genome size")),]
  }
  if (nrow(Orgs)) {
    for (i in 1:nrow(Orgs)) { #i <- 1
      GSize <- Orgs$"Genome size"[i]
      NNucl <- round(GSize/NuclL, 0)
      # Identify affected protein groups (usually, all):
      temp <- proteoCraft::listMelt(strsplit(Prot[[Proteins.col]], ";"), Prot$id)
      temp$Organism <- DB$Organism_Full[match(temp$value, DB$"Protein ID")]
      Pr.Ruler.cols <- Expr.cols$Pr.Ruler.cols
      if (nrow(Orgs) > 1) { Pr.Ruler.cols <- paste0(Orgs$Organism[i], ": ", Pr.Ruler.cols) }
      wpg <- which(Prot$id %in% unique(temp$L1[which(temp$Organism == Orgs$Organism[i])]))
      if (length(wpg)) {
        # Histone IDs:
        # NB: We will not include H1, as its presence on nucleosomes is optional!
        histones <- data.frame(Type = paste0("H", c("2A", "2B", "3", "4")))
        histones$IDs <- lapply(histones$Type, function(x) {
          DB$"Protein ID"[grep(paste0("^Histone ", x, "\\.?[0-9]*$"), DB$"Common Name")]
        })
        temp <- proteoCraft::listMelt(strsplit(PG[[Proteins.col]], ";"), PG$id)
        histones$"Protein Group IDs" <- lapply(histones$IDs, function(x) {
          unique(temp$L1[which(temp$value %in% x)])
        })
        ltest <- sapply(histones$"Protein Group IDs", length)
        if (sum(ltest)) {
          histones <- histones[which(ltest > 0),]
          exprsdata <- Prot[, Expr.cols$Expr.cols, drop = FALSE]
          R <- sapply(Expr.cols$Expr.cols, function(r) { #r <- Expr.cols$Expr.cols[1]
            sapply(histones$"Protein Group IDs", function(h) { mean(proteoCraft::is.all.good(as.numeric(exprsdata[match(h, Prot$id), r]))) })
          })
          if (norm) {
            if (norm_filter) {
              wR <- which(apply(R, 2, function(x) { length(proteoCraft::is.all.good(x)) }) > 0)
              exprsdata <- exprsdata[, wR, drop = FALSE]
            } else { wR <- 1:length(Pr.Ruler.cols) }
            R <- mean(proteoCraft::is.all.good(as.numeric(unlist(R))))
          } else {
            R <- apply(R, 2, function(x) { mean(proteoCraft::is.all.good(x)) })
            wR <- which(proteoCraft::is.all.good(R, 2))
          }
          Pr.Ruler.cols <- Pr.Ruler.cols[wR]
          Prot[, Pr.Ruler.cols] <- NA
          if (length(wR)) {
            if (norm) {
              Prot[wpg, Pr.Ruler.cols] <- exprsdata[wpg, Expr.cols$Expr.cols[wR], drop = FALSE] + log(NNucl, log.Expr.base)-R+1
            } else {
              Prot[wpg, Pr.Ruler.cols] <- sweep(exprsdata[wpg, Expr.cols$Expr.cols[wR], drop = FALSE], 2,
                                                log(NNucl, log.Expr.base)-R+1, "+")
            }
            if ((is.logical(log.Expr))&&(!log.Expr)) {
              for (k in Pr.Ruler.cols) { Prot[wpg, k] <- log.Expr.base^Prot[wpg, k] }
            }
          } else {
            warning(paste0("None of organism ", Orgs$Organism[i], "'s core histones have valid quantitative data for any sample => skipping!"))
          }
        } else { warning(paste0("None of organism ", Orgs$Organism[i], "'s core histones were identified in the dataset => skipping!")) }
      }
    }
    if (norm_filter) { for (k in Pr.Ruler.cols) { if (sum(is.na(Prot[[k]])) == nrow(Prot)) { Prot[[k]] <- NULL} } }
    Res <- list(Protein.groups = Prot, Database = DB)
  } else {
    warning("Not a single valid organisms, returning NA!")
    Res <- NA
  }
  return(Res)
}
