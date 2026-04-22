#' parse.Param.aggreg
#'
#' @description 
#' Internal. Parse an aggregates-type parameter value.
#' Use for assignments to temporary variables, inside functions, or where the value isn't stored in a parameters object.
#' This one is slightly more flexible than the newer parse.Param.aggreg.2(), which takes a parameter name, not its value.
#' However, parse.Param.aggreg.2 is more deeply integrated in the workflows.
#' 
#' @param parameter The parameter to parse.
#' @param aggregates The annotations (basic aggregate building blocks) as a named vector. Default: Aggregates
#' @param map The reference aggregates map. All existing aggregates in this map must exist. Default: Aggregate.map
#' @param aggr.list The list of aggregates. Default: Aggregate.list
#' @param checkExp Logical, default = FALSE. Should we check the Experiment object, dropping it if it has only 1 level?
#' @param experiment Experiment factor object. Allows for dropping Exp if its length is one. Default = "Exp"
#' 
#' @details
#' This function parses a factors aggregate-type parameter in the parameters object and returns a list.
#' 
#' @returns
#' A named list with the name of the factors aggregate, its possible values (levels), the name of the individual aggregated factors, and the name of the corresponding column in the experiment map.
#' 
#' @examples
#' Ratios.Plot.split <- parse.Param.aggreg(parameter = Param$Ratios.Plot.split)
#' 
#' @export

parse.Param.aggreg <- function(parameter,
                               aggregates = Aggregates,
                               map = Aggregate.map,
                               aggr.list = Aggregate.list,
                               checkExp = FALSE,
                               experiment = Exp) {
  if (parameter == "") {
    return()
  }
  if (is.na(parameter)){
    warning("Invalid parameter value!")
    return()
  }
  paramValue <- unlist(strsplit(parameter, ";"))
  if ((checkExp)&&(length(experiment) == 1L)) { paramValue <- setdiff(myPar, "Exp") }
  m <- match(paramValue, names(aggregates))
  if ((!length(m))||(sum(is.na(m)))) { stop("Invalid parameter value!") }
  aggr <- aggregates[m]
  a <- which((lengths(map$Characteristics) == length(aggr))&(vapply(map$Characteristics, function(x) {
    sum(!x %in% aggr)
  }, 1L) == 0L))
  aggr <- map$Aggregate.Name[a]
  val <- aggr.list[[aggr]]
  nms2 <- nms <- unlist(map$Characteristics[which(map$Aggregate.Name == aggr)])
  kol <- if (length(nms) == 1L) { nms } else { aggr }
  if (length(experiment) == 1L) { nms2 <- setdiff(nms2, "Experiment") }
  res <- list(aggregate = aggr,
              values = val,
              names = nms,
              column = kol,
              limmaCol = paste0(paste(nms2, collapse = "_"), "_._"))
  return(res)
}
