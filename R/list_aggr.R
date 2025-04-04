#' list_aggr
#'
#' @description 
#' A function to create melt lists.
#' 
#' @param origin The origin data frame.
#' @param dest The destination data frame. 
#' @param list.col The (character) column of IDs to use for the melting process. See the "insep" argument.
#' @param other.col The other column(s) of interest to be melted.
#' @param dest.match.col The name of the column in the destination data frame whose IDs will be matched to those listed in list.col.
#' @param insep The separator in list.col; default = ";". Change to the correct separator (e.g. ",") if needed.
#' @param outsep The separator to use in the output.
#' 
#' @examples
#' test <- data.frame("values" = c(1:3))
#' test$ids <- list("a;b", "c", "d;e")
#' test2 <- data.frame(ids = c("a", "b", "c", "d", "e", "f"))
#' list_aggr(test, test2, "ids", "values", "ids", insep = ";", outsep = ";")
#' 
#' @export
#' 

list_aggr <- function(origin, dest, list.col, other.col, dest.match.col, insep = ";", outsep = ";") {
  stopifnot(class(origin) %in% c("data.frame", "matrix"),
            class(dest) %in% c("data.frame", "matrix"),
            sum(c(list.col, other.col) %in% colnames(origin)) == length(c(list.col, other.col)),
            dest.match.col %in% colnames(dest))
  origin <- as.data.frame(origin)
  dest <- as.data.frame(dest)
  temp <- origin[[list.col]]
  temp <- strsplit(as.character(temp), insep)
  if (length(other.col) > 1) {
    namez <- apply(origin[,other.col], 1, function(x) {paste(x, collapse = "_-;-_")})
  } else { namez <- origin[[other.col]] }
  names(temp) <- namez
  temp <- melt(temp)
  temp2 <- aggregate(temp$L1, list(temp$value), function(x) {
    paste(unique(sort(as.numeric(x))), collapse = outsep)
  })
  m <- match(dest[[dest.match.col]], temp2$Group.1)
  x <- c(temp2$x, NA)
  m[which(is.na(m))] <- nrow(temp2)+1
  res <- x[m]
  if (length(other.col) > 1) {
    res <- as.data.frame(t(sapply(strsplit(res, "_-;-_"), function(x) {unlist(x)})))
  }
  return(res)
}
