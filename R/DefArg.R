#' DefArg
#'
#' @description 
#' Call this function for any function you want to troubleshoot.
#' This will create default values of the function's parameters in the main environment to allow easy debugging.
#'
#' I used to have an in-house solution for this, but it was slightly buggy.
#' This is based on more standard R code and is more reliable. 
#' 
#' @param FUN name of the function; can be entered as the function object or as its name (character). 
#' @param remove_NULLs Should arguments with default value "NULL" be excluded?\cr
#' TRUE by default, because the risk is that some valuable variable in the parent environment can be "erased" as a result.
#' 
#' @examples
#' temp <- DefArg(F_test)
#' 
#' @export

DefArg <- function(FUN,
                   remove_NULLs = TRUE) {
  #proteoCraft::DefArg(proteoCraft::DefArg)
  #proteoCraft::DefArg(proteoCraft::annot_to_tabl)
  #proteoCraft::DefArg(proteoCraft::PG_assemble)
  #FUN <- proteoCraft::annot_to_tabl
  #FUN <- proteoCraft::PG_assemble
  #FUN <- proteoCraft::Volcano.plot
  if (class(FUN) == "character") { FUN <- eval(parse(text = FUN)) }
  args <- formals(FUN)
  if (remove_NULLs) {
    args <- args[which(!vapply(args, is.null, TRUE))]
  }
  args <- args[which(!vapply(args, function(x) {
    ("name" %in% class(x))&&(length(x) == 1)&&(as.character(x) == "")
  }, TRUE))]
  args <- lapply(args, eval)
  lapply(names(args), function(x) { try(assign(x, args[[x]], envir = .GlobalEnv), silent = TRUE) })
}
