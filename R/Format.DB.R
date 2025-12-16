#' Format.DB
#'
#' @description 
#' A function to load Fasta databases from the work directory or the environment and format them into a data frame.
#' Depending on the mode you select, specific regex rules are applied to the header to extract information.
#' Additionally, three rules can be custom specified by the arguments.
#' Now parallelized!!!
#' 
#' @param file Name (as character) of the Fasta file.
#' @param mode One of "Uniprot", "Ensembl", "RefSeq-RNA", "RefSeq-Protein", "RefSeq-CDS", "NCBI", "TAIR" or "custom".
#' @param in.env If TRUE, the function will attempt processing the character string in one of two ways: a) if it recognizes the character string as a Fasta file, it will directly process it; else, it will treat it as the name of a Fasta file already in the environment. If set to FALSE (default), it will attempt to load the file from the work directory.
#' @param Full.ID.rule Regex rule used to parse the header and extract the full ID. Only used if operating in custom mode. Default is "^>([^ ]*)" (UniProt).
#' @param Protein.ID.rule Regex rule used to parse the header and extract protein ID. Only used if operating in custom mode. Default is "^>[a-z]{2}\\|([^\\|]*)" (UniProt).
#' @param Name.rule Regex rule used to parse the header and extract the protein name. Only used if operating in custom mode. Default is "^>[a-z]{2}\\|[^\\|]*\\|([^ ]*)" (UniProt).
#' @param Unique If TRUE (default), will filter fasta to provide only one accession per sequence. (This may mean that only 1 gene will be retained where several encode the exact same protein, but reducing redundancy is good for FDR.)
#' @param IDs_only Allows extracting just the IDs (default = FALSE). Used internally only.
#' @param parallel Default = FALSE; if TRUE, uses parallel processing and the next few arguments.
#' @param N.clust A limit on the number of vCPUs to use. If left as NULL (default), uses the number of available clusters - 1, to a minimum of 1.
#' @param N.reserved Default = 1. Number of reserved vCPUs the function is not to use. Note that for obvious reasons it will always use at least one.
#' @param cl Already have a cluster handy? Why spend time making a new one, which on top of that may invalidate the old one. Just pass it along!
#' @param trimName Logical (TRUE by default). Remove anything after the last underscore from the name, e.g. "ACTIN_HUMAN" become ACTIN. A more complex but safer way was used previously (up to 3.1.4).
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
  #TESTING <- TRUE;proteoCraft::DefArg(proteoCraft::Format.DB)
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
  F0 <- function(btch) { #btch <- batChes[[1]]
    hdrs <- grep("^>", btch)
    temp1 <- data.frame(Header = btch[hdrs])
    # Batches
    #temp1 <- temp1[1:10000, , drop = FALSE]
    if (IDs_only) { Roolz <- Roolz["Protein ID"] }
    if (mode != "TAIR") {
      for (i in 1:length(Roolz)) {
        temp1[[names(Roolz)[i]]] <- vapply(temp1$Header, function(x) {
          y <- regexpr(Roolz[[i]], x, perl = TRUE)
          return(substr(x,  attributes(y)$capture.start, attributes(y)$capture.start + attributes(y)$capture.length - 1))
        }, "")
        temp1[[names(Roolz)[i]]] <- gsub("^ +| +$", "", temp1[[names(Roolz)[i]]])
      }
    } else { # TAIR-specific behaviour
      if (IDs_only) {
        temp1$"Protein ID" <- gsub(" \\| (Symbols: )?.*", "", gsub("^>", "",temp1$Header))
      } else {
        temp1[, c("Protein ID", "Gene", "Common Name", "Chromosome")] <- proteoCraft::Isapply(strsplit(gsub("^>", "",temp1$Header), " \\| (Symbols: )?"),  unlist)
        temp1 <- temp1[, c("Protein ID", "Common Name", "Gene", "Chromosome", "Header")]
        temp1$"Common Name" <- gsub(" / ", ";", gsub("; ", ", ", gsub("  ?", " ", temp1$"Common Name")))
        temp1$Gene <- gsub(", +", ";", temp1$Gene)
        temp1$Organism_Full <- "Arabidopsis thaliana"
        nms.list <- c("Common Name", "Gene", "Protein ID")
      }
    }
    if (("Name" %in% colnames(temp1))&&(trimName)) {
      temp1$Name <- gsub("_[^_]*$", "", temp1$Name)
    }
    if (!IDs_only) {
      if (mode %in% c("REFSEQCDS", "NCBI", "UNIPROTKB")) {
        temp1$"No Isoforms" <- vapply(temp1$"Protein ID", function(x) { unlist(strsplit(x, "-"))[1] }, "")
      }
      # More Ensembl-specific behaviour
      if (mode == "ENSEMBL") {
        #temp1$"Protein ID" <- temp1$"Full ID" # Should be addressed more generally further down
        temp1$Organism <- gsub("_", " ", unlist(strsplit(file, "\\."))[1])
      }
      # UniProtKB-specific behaviour
      if (mode == "UNIPROTKB") {
        temp1$Organism <- gsub(".*_", "", temp1$"Full ID")
        temp1$Organism_Full <- gsub(" [A-Z]{2}=.+$", "", gsub("^>.+OS=", "", temp1$Header))
        temp1$TaxID <- NA
        w <- grep("OX=", temp1$Header)
        temp1$TaxID[w] <- as.numeric(gsub(" .+", "", vapply(strsplit(temp1$Header[w], "OX="), function(x) { unlist(x)[2] }, "")))
        temp1$"Full Name" <- apply(temp1[, c("Name", "Organism")], 1, FUN = function(x) {paste(x, collapse = "_")})
        a <- strsplit(temp1$Header, " [^ ]+=")
        temp1$"Common Name" <- vapply(1:nrow(temp1), function(x) {
          x1 <- unlist(a[x])[1]
          x2 <- temp1$Organism[x]
          x <- unlist(strsplit(x1, paste0("_", x2, " ?")))
          if (length(x) == 2) { return(x[2]) } else { return("") }
        }, "")
        a <- grep("^[I,i]soform [a-z,A-Z,0-9]+ of ", temp1$"Common Name")
        a1 <- temp1$"Common Name"[a]
        b1 <- gsub("^[I,i]soform [a-z,A-Z,0-9]+ of ", "", a1)
        c1 <- nchar(a1)
        c2 <- nchar(b1)
        b2 <- sapply(1:length(a), function(x) {
          substr(temp1$"Common Name"[a][x], start = 1, stop = c1[x]-c2[x]-4)
        })
        temp1$"Common Name"[a] <- vapply(1:length(a), function(x) { paste(b1[x], ", ", b2[x], sep = "") }, "")
        temp1$Gene <- vapply(strsplit(temp1$Header, "GN="), function(x) {
          x <- unlist(x)
          if (length(x) == 2) {
            x <- x[2]
            return(unlist(strsplit(x, " OS=| PE=| SV="))[1])
          } else { return("") }
        }, "")
        nms.list <- c("Common Name", "Name", "Protein ID", "No Isoforms", "Full ID", "Gene")
      }
      # More RefSeq-specific behaviour
      if (mode %in% c("REFSEQRNA", "REFSEQPROTEIN")) {
        if (mode == "REFSEQRNA") {
          a <- strsplit(temp1$Header, " ")
          temp1[, c("Organism", "Full Name")] <- as.data.frame(t(sapply(a, function(x) {
            x <- unlist(x)
            x1 <- paste(x[2:3], collapse = " ")
            x2 <- paste(x[4:length(x)], collapse = " ")
            return(c(x1, x2))
          })))
        }
        if (mode == "REFSEQPROTEIN") { temp1$"Full Name" <- temp1$Name }
        for (tt in c("transcript variant", "isoform")) { #tt <- "transcript variant"
          a <- grep(tt, temp1$"Full Name", ignore.case = TRUE)
          if (length(a) > 0) {
            u <- unlist(strsplit(tt, " "))
            u[1] <- paste0(toupper(substr(u[1], 1, 1)), substr(u[1], 2, nchar(u[1])))
            u <- paste(u, collapse = ".")
            temp1[[u]] <- ""
            u1 <- paste0(" ?", u, " ?")
            u1 <- paste(c(u1, tolower(u1), toupper(u1)), collapse = "|")
            b <- strsplit(temp1$"Full Name"[a], u1)
            if (tt == "isoform") {
              temp1$"No Isoforms" <- temp1$"Full Name"
              temp1[a, c("No Isoforms", u)] <- proteoCraft::Isapply(b, function(x) {
                x <- x[which(x != "")]
                return(c(x[1], paste(x[2:length(x)], collapse = " ")))
              })
            } else { temp1[a, u] <- vapply(b, function(x) { rev(unlist(x))[1] }, "") }
          } else { if (tt == "isoform") { nms.list <- nms.list[which(nms.list != "No Isoforms")] }}
        }
      }
      if (mode == "REFSEQCDS") {
        tt <- gsub("CCDS:CCDS[0-9]+\\.*[0-9]*,", "", temp1$Header)
        temp1$"Gene ID" <- vapply(tt, function(x) {
          y <- regexpr("^>.*\\[db_xref=GeneID\\:([^\\]]+)\\]", x, perl = TRUE)
          return(substr(x,  attributes(y)$capture.start, attributes(y)$capture.start + attributes(y)$capture.length - 1))
        }, "")
      }
      if (!"Common Name" %in% colnames(temp1)) {
        temp1$"Common Name" <- vapply(strsplit(gsub(" \\[[^\\[]+\\]$", "", temp1$Header), " "), function(x) {
          x <- unlist(x)
          if (length(x) > 1) {
            x <- paste(x[2:length(x)], collapse = " ")
          } else { x <- "" }
          return(x)
        }, "")
        # We want to fill important columns with any value, if none is available
        w <- which(temp1$"Common Name" %in% c("", " ", "NA", NA))
        tmpkol <- nms.list[which(nms.list %in% colnames(temp1))]
        if ((length(w))&&(length(tmpkol))) {
          f0 <- function(x) {
            tt <- which(!x %in% c("", " ", "NA", NA))
            return(x[min(tt)])
          }
          temp1$"Common Name" <- do.call(f0, c(temp1[, c("Common Name", tmpkol)]))
        }
      }
      w <- which(as.character(temp1$"Protein ID") %in% c("", " ", "NA"))
      temp1$"Protein ID"[w] <- temp1$"Full ID"[w]
      hdrs <- c(hdrs, length(btch) + 1)
      temp1$Sequence <- gsub(" ", "", vapply(1:(length(hdrs)-1), function(x) {
        paste(btch[(hdrs[x]+1):(hdrs[x+1]-1)], collapse = "")
      }, ""))
      temp1$Sequence <- gsub("\\*$", "", temp1$Sequence) # Required for some types of databases, such as TAIR
      if (Unique) {
        seq <- unique(temp1$Sequence)
        if (length(seq) < nrow(temp1)) {
          wY <- match(seq, temp1$Sequence)
          wN <- which(!c(1:nrow(temp1)) %in% wY)
          lN <- length(wN)
          cat(paste0(lN, " duplicate sequence", c(" was", "s were")[(lN > 1)+1], " removed.\n"))
          temp1 <- temp1[wY, , drop = FALSE]
        }
      }
      temp1$Sequence <- setNames(temp1$Sequence, temp1$"Protein ID")
    } else { temp1 <- temp1$"Protein ID" }
    rm(list = setdiff(ls(), "temp1"))
    return(temp1)
  }
  if (parallel) {
    environment(F0) <- .GlobalEnv
    parallel::clusterExport(cl, list("Roolz", "IDs_only", "mode", "Unique", "nms.list", "trimName"),
                            envir = environment())
    res <- parallel::parLapply(cl, batChes, F0)
  } else { res <- lapply(batChes, F0) }
  #
  res <- plyr::rbind.fill(res)
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
