#' .Format.DB_worker
#' 
#' Worker function for Format.DB()
#'
#' @param btch A batch of fasta entries
#' @param Roolz List of regex rules
#' @param mode Mode string
#' @param IDs_only Logical, return only IDs?
#' @param trimName Logical, whether to trim names (anything after the last underscore will be trimmed)
#' @param file Original input file (or vector of lines); used as a last ditch effort to provide an organism name when all else fails, assuming - perhaps foolishly - that the name is in the file name. As you should expect, this can fail quite spectacularly so should never be trusted!

.Format.DB_worker <- function(btch,
                              Roolz,
                              mode,
                              IDs_only,
                              trimName,
                              file) { #btch <- batChes[[1]]
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
      temp1[, c("Protein ID", "Gene", "Common Name", "Chromosome")] <- Isapply(strsplit(gsub("^>", "",temp1$Header), " \\| (Symbols: )?"),  unlist)
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
            temp1[a, c("No Isoforms", u)] <- Isapply(b, function(x) {
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
    temp1$Sequence <- setNames(temp1$Sequence, temp1$"Protein ID")
  } else { temp1 <- temp1$"Protein ID" }
  #rm(list = setdiff(ls(), "temp1"))
  return(temp1)
}