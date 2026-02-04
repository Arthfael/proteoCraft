# Boolean functions to check parameter values
validLogicPar %<o% function(x,
                            length_1 = TRUE,
                            noNAs = TRUE) {
  stopifnot(is.character(x))
  if (!exists(x, .GlobalEnv)) { return(FALSE) }
  val <- get(x, .GlobalEnv)
  if (!"logical" %in% class(val)) { return(FALSE) }
  l <- length(val)
  if ((length_1)&&(l > 1)) { return(FALSE) }
  if (noNAs) {
    tst <- sum(is.na(val))
    if (tst) { return(FALSE) }
  }
  return(TRUE)
}
validIntegPar %<o% function(x,
                            min = 1,
                            length_1 = TRUE,
                            noNAs = TRUE) {
  stopifnot(is.character(x))
  if (!exists(x, .GlobalEnv)) { return(FALSE) }
  val <- get(x, .GlobalEnv)
  if (sum(c("integer", "numeric") %in% class(val)) == 0) { return(FALSE) }
  l <- length(val)
  if ((length_1)&&(l > 1)) { return(FALSE) }
  if (noNAs) {
    tst <- sum(is.na(val))
    if (tst) { return(FALSE) }
  }
  if (max(val) > 10^18) {
    stop("Nope, not having it. Go read on integers and floating point arithmetic then come back here!")
  }
  return(sum((is.na(val))|(val >= min)) == l)
}
validCharPar %<o% function(x,
                           filt,
                           length_1 = TRUE,
                           noNAs = TRUE) {
  stopifnot(is.character(x))
  x <- gsub("^ +| +$", "", x)
  if (!exists(x, .GlobalEnv)) { return(FALSE) }
  val <- get(x, .GlobalEnv)
  if (!"character" %in% class(val)) { return(FALSE) }
  l <- length(val)
  if ((length_1)&&(l > 1)) { return(FALSE) }
  if (noNAs) {
    tst <- sum(is.na(val))
    if (tst) { return(FALSE) }
  }
  return(sum(!val %in% filt) == 0)
}
