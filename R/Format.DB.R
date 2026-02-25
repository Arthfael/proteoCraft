#' Format.DB
#'
#' @description 
#' A function to load Fasta databases from the work directory or the environment and format them into a data frame.
#' Depending on the mode you select, specific regex rules are applied to the header to extract information.
#' Additionally, three rules can be custom specified by the arguments.
#' 
#' @param file Name (as character) of the Fasta file.
#' @param mode One of "Uniprot", "Ensembl", "RefSeq-RNA", "RefSeq-Protein", "RefSeq-CDS", "NCBI", "TAIR" or "custom".
#' @param in.env If TRUE, the function will attempt processing the character string in one of two ways: a) if it recognizes the character string as a Fasta file, it will directly process it; else, it will treat it as the name of a Fasta file already in the environment. If set to FALSE (default), it will attempt to load the file from the work directory.
#' @param Full.ID.rule Regex rule used to parse the header and extract the full ID. Only used if operating in custom mode. Default is "^>([^ ]*)" (UniProt).
#' @param Protein.ID.rule Regex rule used to parse the header and extract protein ID. Only used if operating in custom mode. Default is "^>[a-z]{2}\\|([^\\|]*)" (UniProt).
#' @param Name.rule Regex rule used to parse the header and extract the protein name. Only used if operating in custom mode. Default is "^>[a-z]{2}\\|[^\\|]*\\|([^ ]*)" (UniProt).
#' @param Unique If TRUE (default), will filter fasta to provide only one accession per sequence. (This may mean that only 1 gene will be retained where several encode the exact same protein, but reducing redundancy is good for FDR.)
#' @param IDs_only Allows extracting just the IDs (default = FALSE). Used internally only. Ignores Unique!
#' @param parallel Default = FALSE; if TRUE, uses parallel processing and the next few arguments.
#' @param N.clust A limit on the number of vCPUs to use. If left as NULL (default), uses the number of available clusters - 1, to a minimum of 1.
#' @param N.reserved Default = 1. Number of reserved vCPUs the function is not to use. Note that for obvious reasons it will always use at least one.
#' @param cl Already have a cluster handy? Why spend time making a new one, which on top of that may invalidate the old one. Just pass it along!
#' @param trimName Logical (TRUE by default). Remove anything after the last underscore from the name, e.g. "ACTIN_HUMAN" become ACTIN. A more complex but safer way was used previously (up to 3.1.4).
#' 
#' @returns
#' A data frame of parsed annotations.
#' 
#' @examples
#' db <- Format.DB(file = "Search database.FASTA", in.env = TRUE, species = "HUMAN")
#' db <- Format.DB(file = "Search database2.FASTA", in.env = FALSE, species = "RAT;MOUSE")
#' CDS <- Format.DB(list.files()[grep(".*\\.CDS\\..*\\.FASTA$", toupper(list.files()))], mode = "ensembl")
#' 
#' @export

Format.DB <- function(file,
                      mode = "UniProt",
                      in.env = FALSE,
                      Full.ID.rule = "^>([^ ]*)",
                      Protein.ID.rule = "^>[a-z]{2}\\|([^\\|]*)",
                      Name.rule = "^>[a-z]{2}\\|[^\\|]*\\|([^ ]*)",
                      Unique = TRUE,
                      IDs_only = FALSE,
                      parallel = FALSE,
                      N.clust,
                      N.reserved = 1,
                      cl,
                      trimName = TRUE) {
  TESTING <- FALSE
  #TESTING <- TRUE;DefArg(Format.DB)
  #file = unlist(fastasTbl$Data[[i]]); in.env = TRUE; mode = fastasTbl$Type[i]; parallel = TRUE; cl = parClust
  if (TESTING) {
    # Note:
    # This is not a perfect alternative to missing but will work in most cases, unless x matches a function imported by a package 
    misFun <- function(x) { return(!exists(deparse(substitute(x)))) }
  } else { misFun <- missing }
  #
  if (!misFun(cl)) { # Override "parallel" argument if a cluster was provided
    parallel <- TRUE
  }
  #
  if ((is.logical(trimName))||(length(trimName) != 1)||(is.na(trimName))) {
    trimName <- TRUE
  }
  #
  mode2 <- gsub("-|_|\\.| ", "", toupper(mode))
  # Check mode
  if (!mode2 %in% c("UNIPROT", "UNIPROTKB", "SWISSPROT", "TREMBL", "ENSEMBL", "REFSEQRNA", "REFSEQPROTEIN",
                    "REFSEQCDS", "NCBI", "TAIR", "CUSTOM")) {
    warning(paste0("The value provided for \"mode\" (\"", mode, "\") could not be recognized, defaulting to \"custom\"!"))
    mode <- "CUSTOM"
  } else { mode <- mode2 }
  if (mode %in% c("UNIPROT", "UNIPROTKB", "SWISSPROT", "TREMBL")) { mode <- "UNIPROTKB" }
  # Load file
  if (in.env) {if (grepl(">", file[1])) { DB <- file } else { DB <- get(file) }
  } else {
    DB <- readLines(file[1])
    if (length(file) > 1) { for (i in 2:length(file)) { DB <- c(DB, readLines(file[i])) } }
  }
  # Rules
  Roolz <- list(`Full ID` = Full.ID.rule,
                `Protein ID` = Protein.ID.rule,
                Name = Name.rule)
  nms.list <- c("Name", "Protein ID", "Full ID")
  # Ensembl-specific behaviour
  if (mode == "ENSEMBL") {
    Roolz$"Full ID" <- "^>([^ ]*)"
    Roolz$"Protein ID" <- NULL
    Roolz$Name <- NULL
    Roolz$Chromosome <- "^>.*chromosome:([^ ]*)"
    Roolz$Gene <- "^>.*gene:([^ ]*)"
    Roolz$Gene_biotype <- "^>.*gene_biotype:([^ ]*)"
    Roolz$Transcript_biotype <- "^>.*transcript_biotype:([^ ]*)"
    Roolz$Gene_symbol <- "^>.*gene_symbol:([^ ]*)"
    nms.list <- c("Name", "GenBank", "Gene ID", "Protein ID", "No Isoforms")
  }
  # RefSeq-RNA specific behaviour
  if (mode %in% c("REFSEQRNA", "REFSEQPROTEIN")) {
    Roolz$"Full ID" <- "^>([^ ]*)"
    Roolz$"Protein ID" <- Roolz$Name <- NULL
    nms.list <- c("Full Name", "Full ID")
  }
  # RefSeq-Protein specific behaviour
  if (mode == "REFSEQPROTEIN") {
    Roolz$"Full ID" <- "^>([^ ]*)"
    Roolz$"Protein ID" <- "^>gi\\|[0-9]+\\|[a-z]+\\|([^ ]+(?![ \\|])[^ \\|])"
    Roolz$Name <- "^>gi\\|[0-9]+\\|[a-z]+\\|[^ ]+ ([^\\[]+(?! \\[))"
    Roolz$Organism <- "^>gi\\|[0-9]+\\|[a-z]+\\|[^\\[]+\\[([^\\]]+)\\]"
    nms.list <- c("Full Name", "Protein ID", "No Isoforms", "Full ID")
  }
  # RefSeq-CDS specific behaviour
  if (mode == "REFSEQCDS") {
    Roolz$"Full ID" <- "^>([^ ]*)"
    Roolz$Name <- NULL
    Roolz$Protein <- "^>.*\\[protein=([^\\]]+)\\]"
    Roolz$"Protein ID" <- "^>.*\\[protein_id=([^\\]]+)\\]"
    Roolz$Gene <- "^>.*\\[gene=([^\\]]+)\\]"
    nms.list <- c("Name", "Full ID", "Protein", "Protein ID", "No Isoforms", "Gene", "Gene ID")
  }
  # NCBI-specific behaviour
  if (mode == "NCBI") {
    Roolz$"Full ID" <- "^>([^ ]*)"
    Roolz$"Protein ID" <- "^>([^ ]*)"
    Roolz$Name <- "^>[^ ]+ *(.+) \\["
    Roolz$GenBank <- "^>.*(gb\\|[^\\|]+)"
    Roolz$"Gene ID" <- "^>(gi\\|[^\\|]+)"
    Roolz$Organism <- "^>.*\\[([^\\]]+)\\]"
    nms.list <- c("Name", "GenBank", "Gene ID", "Protein ID", "No Isoforms")
  }
  headRs <- grep("^>", DB)
  if ((parallel)&&(length(headRs) < 500)) {
    parallel <- FALSE
  }
  # Create cluster
  if (parallel) {
    tstCl <- stopCl <- misFun(cl)
    if (!misFun(cl)) {
      tstCl <- suppressWarnings(try({
        a <- 1
        parallel::clusterExport(cl, "a", envir = environment())
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
  }
  #
  lH <- length(headRs)
  if (parallel) {
    n <- min(c(length(headRs), N.clust))
    RG <- unique(round(length(DB)*(1:n)/n))
    RG[1:(n-1)] <- vapply(RG[1:(n-1)], function(x) {
      max(headRs[which(headRs < x)+1])-1
    }, 1)
    batChes <- setNames(lapply(1:length(RG), function(x) {
      DB[(c(0, RG)[x]+1):(RG[x])]
    }), paste0("Batch", 1:length(RG)))
  } else { batChes <- list(Batch1 = DB) }
  if (parallel) {
    res <- parallel::clusterApplyLB(cl,
                                    batChes,
                                    .Format.DB_worker,
                                    Roolz = Roolz,
                                    mode = mode,
                                    IDs_only = IDs_only,
                                    trimName = trimName,
                                    file = file)
  } else {
    res <- lapply(batChes,
                  .Format.DB_worker,
                  Roolz = Roolz,
                  mode = mode,
                  IDs_only = IDs_only,
                  trimName = trimName,
                  file = file)
  }
  #
  res <- plyr::rbind.fill(res)
  #
  # Backup code in case we do not have an Organism column yet
  if (!"Organism" %in% colnames(res)) {
    # Not clean but a last resort cheat - this is likely to go wrong often
    res$Organism <- gsub("_", " ", unlist(strsplit(gsub(".*/", "", file), "\\."))[1])
  }
  #
  # Check unicity of sequences
  if (Unique) {
    seq <- unique(res$Sequence)
    if (length(seq) < nrow(res)) {
      wY <- match(seq, res$Sequence)
      wN <- setdiff(1:nrow(res), wY)
      lN <- length(wN)
      cat(paste0(lN, " duplicate sequence", c(" was", "s were")[(lN > 1)+1], " removed.\n"))
      res <- res[wY, , drop = FALSE]
    }
  }
  #
  # Check unicity of IDs!
  IDsCheck <- aggregate(res$Sequence, list(res$"Protein ID"), function(x) { length(unique(x)) })
  if (max(IDsCheck$x) > 1) {
    stop("Database is corrupt: multiple entries were found with distinct sequences and the same protein accession!")
  }
  #
  if ((parallel)&&(stopCl)) {
    parallel::stopCluster(cl)
  }
  return(res)
}
