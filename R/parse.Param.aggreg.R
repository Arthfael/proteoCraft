#' parse.Param.aggreg
#'
#' @description 
#' Internal. Parse an aggregates-type parameter.
#' Only use for assignments to temporary variables, otherwise use the newer parse.Param.aggreg.2 which is more deeply integrated in the workflows.
#' 
#' @param parameter The parameter to parse.
#' @param aggregates The annotations (basic aggregate building blocks) as a named vector. Default: Aggregates
#' @param map The reference aggregates map. All existing aggregates in this map must exist. Default: Aggregate.map
#' @param aggr.list The list of aggregates. Default: Aggregate.list
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
                               aggr.list = Aggregate.list) {
  aggr <- aggregates[which(names(aggregates) %in% unlist(strsplit(parameter, ";")))]
  a <- which((sapply(map$Characteristics, length) == length(aggr))&(sapply(map$Characteristics, function(x) {
    sum(!x %in% aggr)
  }) == 0))
  aggr <- map$Aggregate.Name[a]
  val <- aggr.list[[aggr]]
  nms <- unlist(map$Characteristics[which(map$Aggregate.Name == aggr)])
  if (length(nms) == 1) { kol <- nms } else { kol = aggr}
  res <- list(aggregate = aggr, values = val, names = nms, column = kol)
  return(res)
}
