#' annot_from_GTF
#'
#' @description 
#' Function to parse a GTF (or GFF) file to get annotations (for now only GO-terms).
#'
#' @param file GFF or GTF (GFF v2) file.
#' @param CDS_only Logical: keep only CDS entries? This is a proteomics package so this is TRUE by default: this is usually all we care about.
#' @param mode One of "GFF" or "GTF". If missing, we will attempt to detect it from the file extension.
#' 
#' @examples
#' data <- annot_from_GTF(file)
#' 
#' @export

annot_from_GTF <- function(file,
                           CDS_only = TRUE,
                           mode) {
  TESTING <- FALSE
  #
  #TESTING <- TRUE;proteoCraft::DefArg(proteoCraft::Format.DB_GTF) 
  #file <- paste0(annotDir, "/genomic.gtf")
  #file <- paste0(annotDir, "/genomic.gff")
  allModes <- c("GTF", "GFF")
  if (!missing("mode")) { mode <- toupper(mode) }
  if ((missing("mode"))||(mode %in% allModes)) {
    tst <- setNames(c(grepl("\\.gtf$", file), grepl("\\.gff$", file)),
                    allModes)
    stopifnot(sum(tst) == 1)
    mode <- names(tst)[which(tst)]
  }
  txt <- readLines(file)
  txt <- grep("^#", txt, invert = TRUE, value = TRUE)
  dat <- as.data.frame(t(as.data.frame(strsplit(txt, "\t"))))
  colnames(dat) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
  rownames(dat) <- NULL 
  if (CDS_only) {
    dat <- dat[which(dat$feature == "CDS"),]
  }
  if (mode == "GTF") {
    tmp <- strsplit(gsub("\"$", "", dat$attribute), "\"; *")
    tmp <- lapply(tmp, function(x) { #x <- tmp[[1]]
      x1 <- gsub(" *\".*", "", x)
      x2 <- gsub("^[^ ]+ *\"", "", x)
      x2 <- data.frame(t(x2))
      colnames(x2) <- x1
      return(x2)
    })
  }
  if (mode == "GFF") {
    tmp <- strsplit(dat$attribute, ";")
    tmp <- lapply(tmp, function(x) { #x <- tmp[[1]]
      x1 <- gsub("=.*", "", x)
      x2 <- gsub(".*=", "", x)
      x2 <- data.frame(t(x2))
      colnames(x2) <- x1
      return(x2)
    })
  }
  tmp <- plyr::rbind.fill(tmp)
  if (!"protein_id" %in% colnames(tmp)) { return() }
  tmp$Accession <- tmp$protein_id
  tmp$protein_id <- NULL
  goKol1 <- c("go_component", "go_function", "go_process")
  goKol1 <- goKol1[which(goKol1 %in% colnames(tmp))]
  kol <- c("product", #"Ontology_term",
           goKol1)
  kol <- kol[which(kol %in% colnames(tmp))]
  # Sometimes we have nothing to work with...
  if (!length(kol)) { return() }
  # ... but usually we should've something
  tmp <- tmp[, c("Accession", kol)]
  if (!length(goKol1)) { return() }
  for (k in goKol1) {
    stopifnot(sum(grepl(";", tmp[[k]])) == 0) # Would indicate that our assumption that there is only one entry per row in incorrect
    tmp[[k]] <- gsub("\\|", " [GO:", gsub("\\|\\|.*", "]", tmp[[k]]))
  }
  tmp2 <- reshape::melt(tmp[, c("Accession", goKol1)], id.vars = "Accession")
  tmp2 <- tmp2[which(!is.na(tmp2$value)),]
  tmp2$"GO-ID" <- gsub(".* \\[", "", gsub("\\]$", "", tmp2$value))
  tst <- unique(gsub("^GO:[0-9]{7}$", "", tmp2$"GO-ID"))
  stopifnot(length(tst) == 1,
            tst == "") # Again, would indicate false parsing assumptions
  tmp2 <- aggregate(tmp2[, c("value", "GO-ID")], list(tmp2$Accession), function(x) {
    paste(unique(x), collapse = ";")
  })
  colnames(tmp2) <- c("Accession", "GO", "GO-ID")
  if ("product" %in% colnames(tmp)) {
    tmp2$"Common Name" <- tmp$product[match(tmp2$Accession, tmp$Accession)]
  }
  return(tmp2)
}
