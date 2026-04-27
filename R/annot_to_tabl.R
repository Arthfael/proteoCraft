#' annot_to_tabl
#'
#' @description
#' A function to convert annotated sequences (MaxQuant format) into a table.
#' Useful for creating a table, of, e.g. sequence with amino-acid modifications, site probabilities, etc...
#' 
#' @param x The annotated sequence character.
#' @param start The tag at the start of each annotation. Default = "("
#' @param end The tag at the end of each annotation. Default = ")"
#' @param Nterm Is there a symbol for the N-terminus? Default = TRUE
#' @param Cterm Is there a symbol for the C-terminus? Default = TRUE
#' @param sep Sometimes (for numeric annotations), there can be several per row, for instance when a peptide has several PSMs. The values will be summarized (mean), so we need to know how they are separated. Default = ";"
#' @param numeric_data Is the data numeric? Default = FALSE
#' 
#' @returns
#' A named list, with for each modified peptide a data.frame with 2 columns:
#'  - "Sequence" = amino acid
#'  - "Annotations" = PTM at this position
#' 
#' @export

annot_to_tabl <- function(x,
                          start = "(",
                          end = ")",
                          Nterm = TRUE,
                          Cterm = TRUE,
                          sep = ";",
                          numeric_data = FALSE) {
  #DefArg(annot_to_tabl)
  mySeq <- x <- as.character(x)
  nc <- nchar(mySeq)
  stopifnot(nc > 0L)
  x1 <- substr(mySeq, 1L, 1L)
  xN <- substr(mySeq, nc, nc)
  if (!sum(!x1 %in% AA)) {
    if (Nterm) {
      warning("No N-term symbol detected but Nterm is TRUE => setting it to FALSE")
      Nterm <- FALSE
    }
  } else {
    if (!Nterm) {
      warning("N-term symbol detected but Nterm is FALSE => setting it to TRUE")
      Nterm <- TRUE
    }
  }
  if (!sum(!xN %in% AA)) {
    if (Cterm) {
      warning("No C-term symbol detected but Cterm is TRUE => setting it to FALSE")
      Cterm <- FALSE
    }
  } else {
    if (!Cterm) {
      warning("C-term symbol detected but Cterm is FALSE => setting it to TRUE")
      Cterm <- TRUE
    }
  }
  # Add N- and C-terminal symbols if missing:
  if (!Nterm) { mySeq <- paste0("_", mySeq) }
  if (!Cterm) { mySeq <- paste0(mySeq, "_") }
  #
  mySeq <- gsub(paste0("^", sep, "|", sep, "$"), "", mySeq)
  gr <- grep(sep, mySeq)
  if (length(gr)) {
    mySeq <- gsub("_+", "_", gsub(sep, paste0("_", sep, "_"), mySeq))
    if (!numeric_data) {
      warning(paste0("Some elements contain multiple values separated by \"",
                     sep,"\", I will need to summarize data so will assume it is numeric!"))
      numeric_data <- TRUE
    }
  }
  if ((nchar(start) > 1L)||(nchar(end) > 1L)) {
    stop("Only a single characteer tag can be used for the \"start\" and \"end\" of annotations!")
  }
  #
  mySeq <- strsplit(mySeq, "")
  res <- lapply(mySeq, \(y) {
    y <- unlist(y)
    if (length(y)) {
      left <- which(y == start); ll <- length(left)
      right <- which(y == end); lr <- length(right)
      if(ll != lr) { stop(paste(c("The data appears to be corrupt! ", y), collapse = "")) }
      temp <- data.frame(Sequence = y, Annotations = "")
      if (ll) {
        range <- sapply(c(1L:ll), \(z) { list(left[z]:right[z]) })
        temp$Annotations[left-1L] <- sapply(range, \(z) { 
          z <- unlist(z)
          z <- z[2L:(length(z)-1L)]
          return(paste(y[z], collapse = ""))
        })
        temp <- temp[-unlist(range),]
      }
      if (numeric_data) { temp$Annotations <- as.numeric(temp$Annotations) }
      w <- which(temp$Sequence == sep); lw <- length(w)+1L
      if (lw > 1L) {
        w <- c(w, nrow(temp)+1L)
        temp1 <- temp[1L:(w[1L]-1L),]; colnames(temp1) <- c("Sequence", "A1")
        for (i in 2L:lw) {
          temp1[[paste0("A", i)]] <- temp$Annotations[(w[i-1L]+1L):(w[i]-1L)]
        }
        temp1$Annotations <- rowMeans(temp1[, paste0("A", as.character(1L:lw))], na.rm = TRUE)
        temp1$Annotations[which(!is.all.good(temp1$Annotations, 2L))] <- NA
        temp <- temp1[, c("Sequence", "Annotations")]
      }
      temp$Names <- ""
      temp$Names[1L] <- "N-terminus"
      temp$Names[nrow(temp)] <- "C-terminus"
      w <- which(temp$Names == "")
      temp$Names[w] <- paste0("AA", as.character(1L:length(w)))
      rownames(temp) <- temp$Names
      temp$Names <- NULL
    } else { temp <- data.frame() }
    return(temp)
  })
  names(res) <- mySeq
  return(res)
}
