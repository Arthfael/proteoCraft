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
#' @param silent Logical, default = TRUE. If FALSE, the values of the assigned arguments is printed.\cr
#' TRUE by default, because the risk is that some valuable variable in the parent environment can be "erased" as a result.
#' 
#' @details
#' Call this function as the first step of testing/troubleshooting the code of any function.
#' The default values of any argument of the function of interest will be assigned to the global environment, making debugging/testing easier.
#' 
#' @returns
#' By default the function returns nothing. If silent = FALSE, it returns a named list where each item corresponds to the default value of an argument.
#' 
#' @examples
#' temp <- DefArg(F_test)
#' 
#' @export

DefArg <- function(FUN,
                   remove_NULLs = TRUE,
                   silent = TRUE) {
  #proteoCrafts::DefArg(proteoCraft::DefArg)
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
  invisible({
    silent <- as.logical(silent)[1]
    if (is.na(silent)) { silent <- TRUE }
  })
  if (silent) {
    invisible(setNames(lapply(names(args), function(x) {
      try(assign(x, args[[x]], envir = .GlobalEnv), silent = TRUE)
    }), names(args)))
    return()
  }
  setNames(lapply(names(args), function(x) {
    try(assign(x, args[[x]], envir = .GlobalEnv), silent = TRUE)
  }), names(args))
}
