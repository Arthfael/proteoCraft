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
#' @param sep Sometimes (for numeric annotations), there can be several per row, for instance when a peptide has several evidences. The values will be summarized (mean), so we need to know how they are separated. Default = ";"
#' @param numeric_data Is the data numeric? Default = FALSE
#' 
#' @export

annot_to_tabl <- function(x,
                          start = "(",
                          end = ")",
                          Nterm = TRUE,
                          Cterm = TRUE,
                          sep = ";",
                          numeric_data = FALSE) {
  if (!Nterm) { x <- paste0("_", x) }
  if (!Cterm) { x <- paste0(x, "_") }
  x <- gsub(paste0("^", sep, "|", sep, "$"), "", x)
  x <- gsub(sep, paste0("_", sep, "_"), x)
  gr <- grep(sep, x)
  if ((length(gr) > 0)&&(!numeric_data)) {
    warning(paste("Some elements contain multiple values separated by \"", sep,"\", I will need to summarize data so will assume it is numeric!", sep = ""))
    numeric_data <- TRUE
  }
  if ((nchar(start) > 1)||(nchar(end) > 1)) {
    stop("Only a single characteer tag can be used for the \"start\" and \"end\" of annotations!")
  }
  if (class(x) != "character") {
    warning(paste("Converting argument \"x\" to \"character\" as its class is \"", class(x), "\" instead!", sep = ""))
    x <- as.character(x)
  }
  x <- strsplit(x, "")
  res <- lapply(x, function(y) {
    y <- unlist(y)
    if (length(y) > 0) {
      left <- which(y == start); ll <- length(left)
      right <- which(y == end); lr <- length(right)
      if(ll != lr) { stop(paste(c("The data appears to be corrupt! ", y), collapse = "")) }
      temp <- data.frame(Sequence = y, Annotations = "")
      if (ll > 0) {
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
          temp1[[paste("A", i, sep = "")]] <- temp$Annotations[(w[i-1]+1):(w[i]-1)]
        }
        temp1$Annotations <- rowMeans(temp1[,paste("A", 1:lw, sep = "")], na.rm = TRUE)
        temp1$Annotations[which(!proteoCraft::is.all.good(temp1$Annotations, mode = 2))] <- NA
        temp <- temp1[, c("Sequence", "Annotations")]
      }
      temp$Names <- ""
      temp$Names[1] <- "N-terminus"
      temp$Names[nrow(temp)] <- "C-terminus"
      w <- which(temp$Names == "")
      temp$Names[w] <- paste("AA", c(1:length(w)), sep = "")
      rownames(temp) <- temp$Names; temp$Names <- NULL
    } else { temp <- data.frame() }
    return(temp)
  })
  names(res) <- x
  return(res)
}
