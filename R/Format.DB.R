#' Format.DB
#'
#' @description 
#' A function to load Fasta databases from the work directory or the environment and format them into a data frame.
#' Depending on the mode you select, specific regex rules are applied to the header to extract information.
#' Additionally, three rules can be custom specified by the arguments.
#' 
#' @param file Name (as character) of the Fasta file.
#' @param mode One of "Uniprot", "Ensembl", "RefSeq-RNA", "RefSeq-Protein", "RefSeq-CDS", "NCBI", "TAIR" or "custom" (synonym: "unknown").
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
#' @param revString Character, default = "rev_". If provided, entries with the prefix will be removed by default.
#' @param revDrop Logical, default = TRUE. Assuming that revString is provided, should we remove matching entries? Default = TRUE.
#' @param verbose Logical (FALSE by default).
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
                      #Protein.ID.rule = "^>[a-z]{2}\\|([^\\|]+)",
                      Protein.ID.rule = "^>(?:[a-z]{2}\\|)?([^\\| ]+)",
                      #Name.rule = "^>[a-z]{2}\\|[^\\|]+\\|([^ ]*)", 
                      #Name.rule = "^>[^ ]+ ([^ ]*)", 
                      Name.rule = "^>(?:[a-z]{2}\\|[^\\|]+\\|)?([^ ]*)", 
                      Unique = TRUE,
                      IDs_only = FALSE,
                      parallel = FALSE,
                      N.clust, 
                      N.reserved = 1L,
                      cl,
                      trimName = TRUE,
                      revString = "rev_",
                      revDrop = TRUE,
                      verbose = FALSE) {
  TESTING <- FALSE
  #DefArg(Format.DB); TESTING <- TRUE
  misFun <- if (TESTING) {
    function(x) { return(!exists(deparse(substitute(x)))) }
  } else { missing }
  if (!misFun(cl)) { parallel <- TRUE }
  #
  if ((is.logical(trimName))||(length(trimName) != 1L)||(is.na(trimName))) {
    trimName <- TRUE
  }
  MODE <- gsub("-|_|\\.| ", "", toupper(mode))
  if (MODE == "UNKNOWN") { MODE <- "CUSTOM" }
  if (!MODE %in% c("UNIPROT", "UNIPROTKB", "SWISSPROT", "TREMBL", 
                    "ENSEMBL", "REFSEQRNA", "REFSEQPROTEIN", "REFSEQCDS", 
                    "NCBI", "TAIR", "CUSTOM")) {
    warning(paste0("The value provided for \"mode\" (\"", mode, "\") could not be recognized, defaulting to \"custom\"!"))
    MODE <- "CUSTOM"
  }
  if (MODE %in% c("UNIPROT", "UNIPROTKB", "SWISSPROT", "TREMBL", "UNIPARC")) {
    # For now we assume we can use global rules which can deal together with UniParc and UniProtKB.
    # However the headers are not structured exactly the same way and UniParc may have to be treated separately in future versions...
    MODE <- "UNIPROTKB"
  }
  if (in.env) {
    DB <- if (grepl(">", file[1L])) { file } else { get(file) }
  } else {
    DB <- readLines(file[1L])
    if (length(file) > 1L) {
      for (i in 2L:length(file)) { DB <- c(DB, readLines(file[i])) }
    }
  }
  revTest <- (is.character(revString))&&(length(revString) == 1L)&&(nchar(revString))
  if (revTest) {
    Protein.ID.rule <- gsub("^\\^>", paste0("^>(?:", revString, ")?"), Protein.ID.rule)
    Name.rule <- gsub("^\\^>", paste0("^>(?:", revString, ")?"), Name.rule)
  }
  # print(revString)
  # print(revTest)
  # print(Full.ID.rule)
  # print(Protein.ID.rule)
  # print(Name.rule)
  Roolz <- list(`Full ID` = Full.ID.rule,
                `Protein ID` = Protein.ID.rule, 
                Name = Name.rule)
  nms.list <- c("Name", "Protein ID", "Full ID")
  if (MODE == "ENSEMBL") {
    Roolz$"Full ID" <- "^>([^ ]*)"
    Roolz$"Protein ID" <- NULL
    Roolz$Name <- NULL
    Roolz$Chromosome <- "^>.*chromosome:([^ ]*)"
    Roolz$Gene <- "^>.*gene:([^ ]*)"
    Roolz$Gene_biotype <- "^>.*gene_biotype:([^ ]*)"
    Roolz$Transcript_biotype <- "^>.*transcript_biotype:([^ ]*)"
    Roolz$Gene_symbol <- "^>.*gene_symbol:([^ ]*)"
    nms.list <- c("Name", "GenBank", "Gene ID", "Protein ID", 
                  "No Isoforms")
  }
  if (MODE %in% c("REFSEQRNA", "REFSEQPROTEIN")) {
    Roolz$"Full ID" <- "^>([^ ]*)"
    Roolz$"Protein ID" <- Roolz$Name <- NULL
    nms.list <- c("Full Name", "Full ID")
  }
  if (MODE == "REFSEQPROTEIN") {
    Roolz$"Full ID" <- "^>([^ ]*)"
    Roolz$"Protein ID" <- "^>gi\\|[0-9]+\\|[a-z]+\\|([^ ]+(?![ \\|])[^ \\|])"
    Roolz$Name <- "^>gi\\|[0-9]+\\|[a-z]+\\|[^ ]+ ([^\\[]+(?! \\[))"
    Roolz$Organism <- "^>gi\\|[0-9]+\\|[a-z]+\\|[^\\[]+\\[([^\\]]+)\\]"
    nms.list <- c("Full Name", "Protein ID", "No Isoforms", 
                  "Full ID")
  }
  if (MODE == "REFSEQCDS") {
    Roolz$"Full ID" <- "^>([^ ]*)"
    Roolz$Name <- NULL
    Roolz$Protein <- "^>.*\\[protein=([^\\]]+)\\]"
    Roolz$"Protein ID" <- "^>.*\\[protein_id=([^\\]]+)\\]"
    Roolz$Gene <- "^>.*\\[gene=([^\\]]+)\\]"
    nms.list <- c("Name", "Full ID", "Protein", "Protein ID", 
                  "No Isoforms", "Gene", "Gene ID")
  }
  if (MODE == "NCBI") {
    Roolz$"Full ID" <- "^>([^ ]*)"
    Roolz$"Protein ID" <- "^>([^ ]*)"
    Roolz$Name <- "^>[^ ]+ *(.+) \\["
    Roolz$GenBank <- "^>.*(gb\\|[^\\|]+)"
    Roolz$"Gene ID" <- "^>(gi\\|[^\\|]+)"
    Roolz$Organism <- "^>.*\\[([^\\]]+)\\]"
    nms.list <- c("Name", "GenBank", "Gene ID", "Protein ID", 
                  "No Isoforms")
  }
  headRs <- grep("^>", DB)
  if ((parallel) && (length(headRs) < 500L)) {
    parallel <- FALSE
  }
  if (parallel) {
    stopCl <- FALSE
    if ((is.null(cl))||(!inherits(cl, "cluster"))) {
      dc <- parallel::detectCores()
      if (misFun(N.reserved)) { N.reserved <- 1L }
      nMax <- max(c(dc - N.reserved, 1L))
      if (misFun(N.clust)) { N.clust <- nMax } else {
        if (N.clust > nMax) {
          warning("More cores specified than allowed, I will ignore the specified number! You should always leave at least one free for other processes, see the \"N.reserved\" argument.")
          N.clust <- nMax
        }
      }
      cl <- parallel::makeCluster(N.clust, type = "SOCK")
      stopCl <- TRUE
    }
    N.clust <- length(cl)
  }
  lH <- length(headRs)
  if (parallel) {
    n <- min(c(length(headRs), N.clust))
    RG <- unique(round(length(DB) * (1L:n)/n))
    RG[1L:(n - 1L)] <- vapply(RG[1L:(n - 1L)], \(x) {
      max(headRs[which(headRs < x) + 1L]) - 1L
    }, 1L)
    batChes <- setNames(lapply(1L:length(RG), \(x) {
      DB[(c(0L, RG)[x] + 1L):(RG[x])]
    }), paste0("Batch", 1L:length(RG)))
  } else {
    batChes <- list(Batch1 = DB)
  }
  if (parallel) {
    xports <- list("Roolz", "IDs_only", "MODE", "Unique", "nms.list", "trimName", "Isapply")
    parallel::clusterExport(cl, xports, envir = environment())
    if (verbose) { cat(" Running parallel...\n") }
    f0 <- .bind_worker(.Format.DB_worker,
                       list(Roolz = Roolz,
                            IDs_only = IDs_only,
                            MODE = MODE,
                            Unique = Unique,
                            nms.list = nms.list,
                            trimName = trimName))
    if (verbose) { cat("    (cluster size = ", length(cl), ")\n") }
    res <- parallel::clusterApplyLB(cl,
                                    batChes,
                                    f0,
                                    Roolz = Roolz,
                                    mode = MODE,
                                    IDs_only = IDs_only,
                                    trimName = trimName)
  } else {
    if (verbose) { cat(" Running serial...\n") }
    res <- lapply(batChes,
                  .Format.DB_worker,
                  Roolz = Roolz,
                  mode = MODE,
                  IDs_only = IDs_only,
                  trimName = trimName)
  }
  res <- plyr::rbind.fill(res)
  #
  # Backup code in case we do not have an Organism column yet
  if (!"Organism" %in% colnames(res)) {
    # Not clean but a last resort cheat - this is likely to go wrong often
    res$Organism <- gsub("_", " ", unlist(strsplit(gsub(".*/", "", file), "\\."))[1L])
  }
  #
  # Check unicity of sequences
  if (Unique) {
    seq <- unique(res$Sequence)
    if (length(seq) < nrow(res)) {
      wY <- match(seq, res$Sequence)
      wN <- setdiff(1L:nrow(res), wY)
      lN <- length(wN)
      cat(paste0("Removing ", lN, " duplicate sequence", c("", "s")[(lN > 1L)+1L], "\n"))
      res <- res[wY, , drop = FALSE]
    }
  }
  #
  # Deal with rev
  if (revTest) {
    if ((!is.logical(revDrop))||(length(revDrop) != 1L)||(is.na(revDrop))) { revDrop <- TRUE }
    if (revDrop) { # Either remove rev entries...
      gN <- grep(paste0("^>", revString), res$Header, invert = TRUE)
      res <- res[gN,]
    } else { # ... or mark them as such
      gY <- grep(paste0("^>", revString), res$Header)
      res$`Protein ID`[gY] <- paste0(revString, res$`Protein ID`[gY])
      for (k in c("Name", "Full Name", "Common Name")) {
        if (k %in% colnames(res)) {
          res[gY, k] <- paste0(res$Name[gY, k], " (reversed)")
        }
      }
    }
  }
  #res <- res[order(res$Sequence),]
  #View(res)
  IDsCheck <- aggregate(res$Sequence, list(res$"Protein ID"), \(x) { length(unique(x)) })
  #i <- 1L
  #i <- i+1L
  #IDsCheck$Group.1[i]
  #res$Sequence[which(res$`Protein ID` == IDsCheck$Group.1[i])]
  IDsCheck <- IDsCheck[order(IDsCheck$x, decreasing = TRUE),]
  # print(IDsCheck[1:5,])
  if (max(IDsCheck$x) > 1L) {
    stop("Incorrect Protein.ID.rule or corrupt database: multiple entries with distinct sequences have the same protein accession!")
  }
  if ((parallel) && (stopCl)) {
    parallel::stopCluster(cl)
  }
  return(res)
}
