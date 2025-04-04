#' parse.Param.aggreg.2
#'
#' @description 
#' A new function to parse an aggregates-type parameter. Extends the functionalities of the former version.
#' The old function is still valid for parsing to temporary variables.
#' 
#' @param param.nm The name of the parameter to parse.
#' @param parameters The parameters file. Default = Param
#' @param aggregates The annotations (basic aggregate building blocks) as a named vector. Default: Aggregates
#' @param map The reference aggregates map. All existing aggregates in this map must exist. Default: Aggregate.map
#' @param aggr.list The list of aggregates. Default: Aggregate.list
#' @param param.list.nm The name of the vector of parsed parameters. Default = "Param.aggreg"
#' @param parsed.param.nm The name of the parsed parameter variable to create. If not specified, defaults to param.nm
#' @param obj.list.nm The name of the list of objects to never purge from the global environment. Default = ".obj"
#' 
#' @examples
#' parse.Param.aggreg.2("Ratios.Plot.split")
#' 
#' @export

parse.Param.aggreg.2 <- function(param.nm, parameters = Param, aggregates = Aggregates, map = Aggregate.map,
                                 aggr.list = Aggregate.list, param.list.nm = "Param.aggreg", parsed.param.nm,
                                 obj.list.nm = ".obj") {
  #proteoCraft::DefArg(proteoCraft::parse.Param.aggreg.2)
  if (missing(parsed.param.nm)) { parsed.param.nm <- param.nm }
  aggr <- aggregates[which(names(aggregates) %in% unlist(strsplit(parameters[[param.nm]], ";")))]
  a <- which((sapply(map$Characteristics, length) == length(aggr))&(sapply(map$Characteristics, function(x) {sum(!x %in% aggr)}) == 0))
  aggr <- map$Aggregate.Name[a]
  val <- aggr.list[[aggr]]
  nms <- unlist(map$Characteristics[which(map$Aggregate.Name == aggr)])
  if (length(nms) == 1) { kol <- nms } else { kol = aggr}
  res <- list(aggregate = aggr, values = val, names = nms, column = kol)
  assign(parsed.param.nm, res, parent.frame())
  assign(param.list.nm, unique(c(get(param.list.nm), parsed.param.nm)), parent.frame())
  assign(obj.list.nm, unique(c(get(obj.list.nm), parsed.param.nm)), parent.frame())
}
