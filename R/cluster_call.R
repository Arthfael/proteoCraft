#' cluster_call
#' 
#' Persistent Parallel Apply on PSOCK Clusters
#'
#' A safe and fast wrapper around the parallel apply family that avoids
#' exporting large closures. Objects are assigned once per call to each worker,
#' preventing scoping bugs and dramatically improving performance when running
#' parallel code inside functions.
#'
#' The package namespace of the caller is automatically loaded on workers so
#' internal helpers and S3 methods are available.
#'
#' @param cl A cluster created with [parallel::makeCluster()].
#' @param X Input object.
#' @param FUN Function to evaluate.
#' @param ... Named objects to assign on workers and extra arguments passed to FUN.
#' @param type One of `"lapply"`, `"sapply"`, or `"apply"`.
#' @param MARGIN For `type="apply"` only. Passed to [base::apply()].
#' @param simplify Logical; simplify result for `"sapply"`.
#' @param USE.NAMES Logical; keep names for `"sapply"`.
#' @param LB Logical; use load-balanced scheduling (`clusterApplyLB`).
#' @param export Optional character vector of objects to export from caller env.
#' @param packages Optional additional packages to load on workers.
#'
#' @return Result identical to the corresponding parallel apply function.
#'
#' @details
#' Workers are persistent R sessions. This function:
#' \enumerate{
#'   \item loads required packages,
#'   \item assigns provided objects once,
#'   \item evaluates a lightweight function in parallel.
#' }
#'
#' This avoids repeated serialization of parent environments, a common cause of
#' extremely slow parallel execution inside functions.
#'
#' @examples
#' \dontrun{
#' cl <- parallel::makeCluster(4)
#'
#' A <- rnorm(1000)
#' f <- function(i) A[i]^2
#'
#' cluster_call(cl, 1:1000, f, A=A)
#' cluster_call(cl, 1:1000, f, A=A, type="sapply")
#'
#' M <- matrix(rnorm(1000), ncol=10)
#' cluster_call(cl, M, rowMeans, type="apply", MARGIN=1)
#'
#' parallel::stopCluster(cl)
#' }
#'
#' @export

cluster_call <- function(cl,
                         X,
                         FUN,
                         ...,
                         type = c("lapply","sapply","apply"),
                         MARGIN = NULL,
                         simplify = TRUE,
                         USE.NAMES = TRUE,
                         LB = FALSE,
                         export = NULL,
                         packages = NULL) {
  #
  type <- match.arg(type)
  dots <- list(...)
  #
  ## force promises
  if (length(dots)) {
    for (nm in names(dots)) { dots[[nm]] <- force(dots[[nm]]) }
  }
  ## load packages on workers
  pkgs <- unique(c(packages, utils::packageName()))
  if (length(pkgs)) {
    parallel::clusterExport(cl, "pkgs", envir = environment())
    parallel::clusterEvalQ(cl, {
      for (p in pkgs) library(p, character.only = TRUE)
      NULL
    })
  }
  
  ## optional exports
  if (!is.null(export)) {
    parallel::clusterExport(cl, export, envir = parent.frame())
  }
  ## assign objects once per call
  if (length(dots)) {
    parallel::clusterCall(cl, function(vals) {
      list2env(vals, envir = .GlobalEnv)
      NULL
    }, dots)
  }
  ## dispatch
  if (type == "lapply") {
    if (LB) {
      return(parallel::clusterApplyLB(cl, X, FUN))
    } else {
      return(parallel::parLapply(cl, X, FUN))
    }
  }
  if (type == "sapply") {
    res <- if (LB) {
      parallel::clusterApplyLB(cl, X, FUN)
    } else {
      parallel::parLapply(cl, X, FUN)
    }
    return(base::simplify2array(res, higher = simplify))
  }
  if (type == "apply") {
    if (is.null(MARGIN)) { stop("MARGIN must be provided when type='apply'") }
    idx <- seq_len(dim(X)[MARGIN])
    wrapper <- function(i) {
      slice <- if (MARGIN == 1) { X[i, , drop = FALSE] } else { X[, i, drop = FALSE] }
      FUN(slice)
    }
    res <- if (LB) {
      parallel::clusterApplyLB(cl, idx, wrapper)
    } else {
      parallel::parLapply(cl, idx, wrapper)
    }
    return(base::simplify2array(res))
  }
}