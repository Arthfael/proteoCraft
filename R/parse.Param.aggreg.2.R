#' parse.Param.aggreg.2
#'
#' @description 
#' Internal. Parse an aggregates-type parameter. Extends the functionalities of the former version.
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
#' @details
#' This function parses a factors aggregate-type parameter in the parameters object and assigns (without returning anything) the corresponding factors aggregate list to the global environment, also adding them to the list of persistent objects.
#' 
#' @returns
#' Does not return anything. For a version returning something, use parse.Param.aggreg()
#' 
#' @examples
#'  parse.Param.aggreg.2("Ratios.Plot.split")
#' # This is equivalent to:
#'  parse.Param.aggreg.2(Param$Ratios.Plot.split)
#' #... but with more background maintenance.
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
  assign(parsed.param.nm, res, #parent.frame()
         .GlobalEnv)
  assign(param.list.nm, unique(c(get(param.list.nm), parsed.param.nm)), #parent.frame()
         .GlobalEnv)
  assign(obj.list.nm, unique(c(get(obj.list.nm), parsed.param.nm)), #parent.frame()
         .GlobalEnv)
  return()
}
