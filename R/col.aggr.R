#' col.aggr
#'
#' @description
#' A function to add an aggregated column to an existing data frame.
#' Slated for phasing out and deletion: this is rather old and clumsy.
#' 
#' @param Out The destination data frame.
#' @param In The origin data frame.
#' @param Out.names A vector of n column names. The first is the name of the column in the destination df to match to the column in the origin df; the others are the names of the new columns created.
#' @param In.names A vector of n column names. The first is the name of the column in the origin df to match to the column in the destination df and to aggregate by; the others are the names of the columns to aggregate.
#' @param FUN The function to aggregate by.
#' @param Remove If set to TRUE (default), removes columns in the destination file which may have identical names. If not, a number is appended to the column names such that the new column names are not already present in the destination file.
#' @param Coerce Optional parameter: if set to TRUE, will attempt to output columns in the same format as the input. Default = FALSE.
#' 
#' @examples
#' pep <- col.aggr(Out = pep, In = ev, Out.names = c("Modified.sequence", "PEP"), In.names = c("Modified.sequence", "PEP"), FUN = function(x) {
#'   x <- proteoCraft::is.all.good(as.numeric(x))
#'   if (length(x) > 0) {
#'     return(min(x))
#'   } else {
#'     return(NA)
#'   }
#' })
#' # Result: pep now has a "PEP" column with PEP scores. Takes seconds or minutes where sapply takes hours.
#' 
#' @export

col.aggr <- function(Out,
                     In,
                     Out.names,
                     In.names,
                     FUN,
                     Remove = TRUE,
                     Coerce = FALSE) {
  stopifnot(class(Out) == "data.frame",
            class(In) == "data.frame",
            class(Out.names) == "character",
            class(In.names) == "character",
            length(Out.names) >= 2,
            length(Out.names) == length(In.names),
            class(Remove) == "logical",
            Out.names[1] %in% colnames(Out),
            In.names %in% colnames(In))
  In.IDs <- In.names[1]
  In.names <- In.names[2:length(In.names)]
  Out.IDs <- Out.names[1]
  Out.names <- Out.names[2:length(Out.names)]
  if (Coerce) { klass <- sapply(In.names, function(x) {class(In[,x])}) }
  if (Remove) {
    a <- colnames(Out)[which(!colnames(Out) %in% Out.names)]
    Out <- Out[, a, drop = FALSE]
  }
  X <- aggregate(In[, In.names], list(In[[In.IDs]]), FUN)
  if (!Remove) {
    Out.names2 <- Out.names
    w <- which(Out.names2 %in% colnames(Out))
    k <- 0
    while(length(w) > 0) {
      k <- k+1
      Out.names2[w] <- paste0(Out.names[w], "_", k)
    }
    Out.names <- Out.names2
  }
  colnames(X) <- c(In.IDs, Out.names)
  w <- which(Out[[Out.IDs]] %in% X[[In.IDs]])
  Out[w, Out.names] <- X[match(Out[w, Out.IDs], X[[In.IDs]]), Out.names]
  if (Coerce) {
    for (i in 1:length(Out.names)) {
      if (klass[i] != class(Out[, Out.names[i]])) {
        a <- get(paste("as.", klass[i], sep = ""))
        Out[[Out.names[i]]] <- a(Out[[Out.names[i]]])
      }
    }
  }
  return(Out)
}
