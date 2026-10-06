#' clustCheck
#' 
#' @description
#' Claude.AI-designed function to check our parallel cluster.
#' 
#' @param cl A cluster to check.
#' @param n Expected number of logical cores in the cluster.
#' @param timeout Default timeout = 5(s).
#' 
#' @returns
#' TRUE if the cluster is working, otherwise FALSE.
#' 
#' @examples
#' tst <- clustCheck(parClust)
#' 
#' @export

clustCheck <- function(cl,
                       n = N.clust,
                       timeout = 5L) {
  if ((!inherits(cl, "cluster")) || (length(cl) != n)) { return(FALSE) }
  # 2. sockets still open?
  open_ok <- vapply(cl, \(node) {
    tryCatch(isTRUE(isOpen(node$con)), error = \(e) FALSE)
  }, TRUE)
  if (!all(open_ok)) { return(FALSE) }
  # 3. round trip with a nonce, under a timeout
  token <- paste(sample(c(letters, 0L:9L), 16L, replace = TRUE), collapse = "")
  res <- tryCatch(R.utils::withTimeout(parallel::clusterCall(cl,
                                                             \(tok) { list(tok = tok, pid = Sys.getpid()) },
                                                             token),
                                       timeout = timeout, onTimeout = "error"),
                  error = \(e) NULL)
  if ((!is.list(res)) || (length(res) != n)) { return(FALSE) }
  ok_struct <- vapply(res, \(r) {
    is.list(r) &&
      identical(names(r), c("tok", "pid")) &&
      identical(r$tok, token) &&
      is.numeric(r$pid)
  }, TRUE)
  if (!all(ok_struct)) { return(FALSE) }
  pids <- vapply(res, `[[`, 1, "pid")
  if (anyDuplicated(pids)) { return(FALSE) }
  return(TRUE)
}
