#' Param_filter
#'
#' @description 
#' A simple function to remove some covariates from semicolon-separated covariate aggregate parameters.
#'
#' @param param The parameter to filter.
#' @param filter Which covariate(s) to filter.
#' 
#' @examples
#' Param$Ratios.Groups.Ref.Aggregate.Level
#' # [1] "Exp;Con;Rep"
#' Param_filter(Param$Ratios.Groups.Ref.Aggregate.Level, "Rep")
#' # [1] "Exp;Con"
#' 
#' @export

Param_filter <- function(param, filter) {
  param <- unlist(strsplit(param, ";"))
  param <- param[which(param != filter)]
  return(paste(param, collapse = ";"))
}
