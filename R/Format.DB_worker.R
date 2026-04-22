#' .Format.DB_worker
#' 
#' Worker function for Format.DB()
#'
#' @param btch Batch of fasta entries
#' @param Roolz List of regex rules
#' @param mode Mode string
#' @param IDs_only Logical, return only IDs?
#' @param trimName Logical, whether to trim names (anything after the last underscore will be trimmed)

.Format.DB_worker <- function(btch,
                              Roolz,
                              mode,
                              IDs_only,
                              trimName) { #i <- 1L
  #DefArg(.Format.DB_worker);btch <- batChes[[1L]]
  hdrs <- grep("^>", btch)
  temp1 <- data.frame(Header = btch[hdrs])
  # Batches
  #temp1 <- temp1[1L:10000L, , drop = FALSE]
  if (IDs_only) { Roolz <- Roolz["Protein ID"] }
  if (mode != "TAIR") {
    for (i in 1L:length(Roolz)) {
      temp1[[names(Roolz)[i]]] <- vapply(temp1$Header, \(x) { #x <- temp1$Header[1L]
        y <- regexpr(Roolz[[i]], x, perl = TRUE)
        w <- which(attributes(y)$capture.length > 0L)
        if (!length(w)) { return("") }
        return(substr(x,  attributes(y)$capture.start[w], attributes(y)$capture.start[w] + attributes(y)$capture.length[w] - 1L))
      }, "")
      temp1[[names(Roolz)[i]]] <- gsub("^ +| +$", "", temp1[[names(Roolz)[i]]])
    }
  } else { # TAIR-specific behaviour
    if (IDs_only) {
      temp1$"Protein ID" <- gsub(" \\| (Symbols: )?.*", "", gsub("^>", "",temp1$Header))
    } else {
      temp1[, c("Protein ID", "Gene", "Common Name", "Chromosome")] <- Isapply(strsplit(gsub("^>", "",temp1$Header), " \\| (Symbols: )?"),
                                                                               unlist)
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
      temp1$"No Isoforms" <- vapply(temp1$"Protein ID", \(x) { unlist(strsplit(x, "-"))[1L] }, "")
    }
    # UniProtKB-specific behaviour
    if (mode == "UNIPROTKB") {
      temp1$Organism <- temp1$Organism_Full <- ""
      g <- grep("^>.+OS=", temp1$Header)
      temp1$Organism[g] <- temp1$Organism_Full[g] <- gsub(" [A-Z]{2}=.+$", "", gsub("^>.+OS=", "", temp1$Header[g]))
      g <- grep("[^_]+_[^_]+", temp1$"Full ID")
      temp1$Organism[g] <- gsub(".*_", "", temp1$"Full ID"[g])
      unique(temp1$Organism_Full)
      temp1$Organism <- gsub(" +\\(.+", "", temp1$Organism)
      w <- which((nchar(temp1$Organism_Full) == 0L)&(nchar(temp1$Organism) > 0L))
      temp1$Organism_Full[w] <- temp1$Organism[w]
      temp1$TaxID <- NA
      w <- grep("OX=", temp1$Header)
      temp1$TaxID[w] <- as.numeric(gsub(" .+", "", vapply(strsplit(temp1$Header[w], "OX="), \(x) {
        unlist(x)[2L]
      }, "")))
      temp1$"Full Name" <- do.call(paste, c(temp1[, c("Name", "Organism")], sep = "_"))
      a <- strsplit(temp1$Header, " [^ ]+=")
      temp1$"Common Name" <- vapply(1L:nrow(temp1), \(x) {
        x1 <- unlist(a[x])[1L]
        x2 <- temp1$Organism[x]
        x <- unlist(strsplit(x1, paste0("_", x2, " ?")))
        if (length(x) == 2L) { return(x[2L]) } else { return("") }
      }, "")
      a <- grep("^[I,i]soform [a-z,A-Z,0-9]+ of ", temp1$"Common Name")
      if (length(a)) {
        a1 <- temp1$"Common Name"[a]
        b1 <- gsub("^[I,i]soform [a-z,A-Z,0-9]+ of ", "", a1)
        c1 <- nchar(a1)
        c2 <- nchar(b1)
        b2 <- vapply(1L:length(a), \(x) {
          substr(a1[x], 1L, c1[x]-c2[x]-4L)
        }, "")
        temp1$"Common Name"[a] <- vapply(1L:length(a), \(x) { paste(b1[x], ", ", b2[x], sep = "") }, "")
      }
      temp1$Gene <- vapply(strsplit(temp1$Header, "GN="), \(x) {
        x <- unlist(x)
        if (length(x) == 2L) {
          x <- x[2L]
          return(unlist(strsplit(x, " OS=| PE=| SV="))[1L])
        } else { return("") }
      }, "")
      nms.list <- c("Common Name", "Name", "Protein ID", "No Isoforms", "Full ID", "Gene")
    }
    # More RefSeq-specific behaviour
    if (mode %in% c("REFSEQRNA", "REFSEQPROTEIN")) {
      if (mode == "REFSEQRNA") {
        a <- strsplit(temp1$Header, " ")
        temp1[, c("Organism", "Full Name")] <- as.data.frame(t(sapply(a, \(x) {
          x <- unlist(x)
          x1 <- paste(x[2L:3L], collapse = " ")
          x2 <- paste(x[4L:length(x)], collapse = " ")
          return(c(x1, x2))
        })))
      }
      if (mode == "REFSEQPROTEIN") { temp1$"Full Name" <- temp1$Name }
      for (tt in c("transcript variant", "isoform")) { #tt <- "transcript variant"
        a <- grep(tt, temp1$"Full Name", ignore.case = TRUE)
        if (length(a)) {
          u <- unlist(strsplit(tt, " "))
          u[1L] <- paste0(toupper(substr(u[1L], 1L, 1L)), substr(u[1L], 2L, nchar(u[1L])))
          u <- paste(u, collapse = ".")
          temp1[[u]] <- ""
          u1 <- paste0(" ?", u, " ?")
          u1 <- paste(c(u1, tolower(u1), toupper(u1)), collapse = "|")
          b <- strsplit(temp1$"Full Name"[a], u1)
          if (tt == "isoform") {
            temp1$"No Isoforms" <- temp1$"Full Name"
            temp1[a, c("No Isoforms", u)] <- Isapply(b, \(x) {
              x <- x[which(x != "")]
              return(c(x[1L], paste(x[2L:length(x)], collapse = " ")))
            })
          } else { temp1[a, u] <- vapply(b, \(x) { rev(unlist(x))[1L] }, "") }
        } else { if (tt == "isoform") { nms.list <- nms.list[which(nms.list != "No Isoforms")] }}
      }
    }
    if (mode == "REFSEQCDS") {
      tt <- gsub("CCDS:CCDS[0-9]+\\.*[0-9]*,", "", temp1$Header)
      temp1$"Gene ID" <- vapply(tt, \(x) {
        y <- regexpr("^>.*\\[db_xref=GeneID\\:([^\\]]+)\\]", x, perl = TRUE)
        return(substr(x,  attributes(y)$capture.start, attributes(y)$capture.start + attributes(y)$capture.length - 1L))
      }, "")
    }
    if (!"Common Name" %in% colnames(temp1)) {
      temp1$"Common Name" <- vapply(strsplit(gsub(" \\[[^\\[]+\\]$", "", temp1$Header), " "), \(x) {
        x <- unlist(x)
        if (length(x) > 1L) {
          x <- paste(x[2L:length(x)], collapse = " ")
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
    hdrs <- c(hdrs, length(btch) + 1L)
    temp1$Sequence <- gsub(" ", "", vapply(1L:(length(hdrs)-1L), \(x) {
      paste(btch[(hdrs[x]+1L):(hdrs[x+1L]-1L)], collapse = "")
    }, ""))
    temp1$Sequence <- gsub("\\*$", "", temp1$Sequence) # Required for some types of databases, such as TAIR
    temp1$Sequence <- stats::setNames(temp1$Sequence, temp1$"Protein ID")
  } else { temp1 <- temp1$"Protein ID" }
  #rm(list = setdiff(ls(), "temp1"))
  return(temp1)
}
