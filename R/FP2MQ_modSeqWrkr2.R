#' .FP2MQ_modSeqWrkr2
#' 
#' Worker function used by FP_to_MQ().
#'
#' @param x Temp data, see FP_to_MQ code...
#' @param openSearch Is this an open search?
#' 
#' @export

.FP2MQ_modSeqWrkr2 <- function(x,
                               openSearch = OpenSearch) {
  #x <- tmp[i,]
  seq <- unlist(x[[1L]])
  l <- length(seq)
  rs <- data.frame(Seq = c("N-term", seq, "C-term"),
                   Mod = "")
  ptms <- x[[2L]]
  if (is.data.frame(ptms)) {
    stopifnot(sum(rs$Seq[ptms$Position + 1L] != ptms$Site) == 0)
    rs$Mod[ptms$Position + 1L] <- ptms$"Full name"
    tmp2 <- c(rs$Mod[1L], rs$Mod[2L])
    tmp2 <- tmp2[which(tmp2 != "")]
    if (length(tmp2)) {
      tmp2 <- stats::aggregate(tmp2,
                               list(tmp2),
                               length)
      tmp2$x <- as.character(tmp2$x)
      tmp2 <- paste(gsub("^1 ", "",
                         do.call(paste, c(tmp2[, c("x", "Group.1")], sep = " "))),
                    collapse = ",")
    } else { tmp2 <- "" }
    rs$Mod[2L] <- tmp2; rs$Mod[1L] <- ""
    rs$Mod[l+1L] <- paste0(rs$Mod[l+1L], rs$Mod[l+2L]); rs$Mod[l+2L] <- ""
  }
  rs$Seq[1L] <- rs$Seq[l+2L] <- "_"
  # NB: although MaxQuant does not label fixed modifications...
  # we will keep those labels because it's useful information and does not hurt
  #
  # Open search:
  # Add a simple catch-all mark for delta mass which will just indicate that one should look into the "Mass error [Da]" column
  if (OpenSearch) {
    obsptms <- unlist(x[[3L]])
    if (length(obsptms)) {
      obsptms <- obsptms[which(obsptms != "")]
      rs$Mod[l+1L] <- gsub("^,", "", paste(c(rs$Mod[l+1L], obsptms), collapse = ","))
    }
  }
  w <- which(rs$Mod != "")
  rs$Mod[w] <- paste0("(", rs$Mod[w], ")")
  rs <- paste(do.call(paste, c(rs, sep = "")), collapse = "")
  return(rs)
}
