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
  #proteoCraft::DefArg(annot_to_tabl)
  mySeq <- x <- as.character(x)
  nc <- nchar(mySeq)
  stopifnot(nc > 0)
  x1 <- substr(mySeq, 1, 1)
  xN <- substr(mySeq, nc, nc)
  if (x1 %in% proteoCraft::AA) {
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
  if (xN %in% proteoCraft::AA) {
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
  if ((nchar(start) > 1)||(nchar(end) > 1)) {
    stop("Only a single characteer tag can be used for the \"start\" and \"end\" of annotations!")
  }
  #
  mySeq <- strsplit(mySeq, "")
  res <- lapply(mySeq, function(y) {
    y <- unlist(y)
    if (length(y)) {
      left <- which(y == start); ll <- length(left)
      right <- which(y == end); lr <- length(right)
      if(ll != lr) { stop(paste(c("The data appears to be corrupt! ", y), collapse = "")) }
      temp <- data.frame(Sequence = y, Annotations = "")
      if (ll) {
        range <- sapply(c(1:ll), function(z) { list(left[z]:right[z]) })
        temp$Annotations[left-1] <- sapply(range, function(z) { 
          z <- unlist(z)
          z <- z[2:(length(z)-1)]
          return(paste(y[z], collapse = ""))
        })
        temp <- temp[-unlist(range),]
      }
      if (numeric_data) { temp$Annotations <- as.numeric(temp$Annotations) }
      w <- which(temp$Sequence == sep); lw <- length(w)+1
      if (lw > 1) {
        w <- c(w, nrow(temp)+1)
        temp1 <- temp[1:(w[1]-1),]; colnames(temp1) <- c("Sequence", "A1")
        for (i in 2:lw) {
          temp1[[paste0("A", i)]] <- temp$Annotations[(w[i-1]+1):(w[i]-1)]
        }
        temp1$Annotations <- rowMeans(temp1[, paste0("A", as.character(1:lw))], na.rm = TRUE)
        temp1$Annotations[which(!proteoCraft::is.all.good(temp1$Annotations, mode = 2))] <- NA
        temp <- temp1[, c("Sequence", "Annotations")]
      }
      temp$Names <- ""
      temp$Names[1] <- "N-terminus"
      temp$Names[nrow(temp)] <- "C-terminus"
      w <- which(temp$Names == "")
      temp$Names[w] <- paste0("AA", as.character(1:length(w)))
      rownames(temp) <- temp$Names
      temp$Names <- NULL
    } else { temp <- data.frame() }
    return(temp)
  })
  names(res) <- mySeq
  return(res)
}
