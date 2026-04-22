#' parse.Param.aggreg.2
#'
#' @description 
#' Internal. Parse an aggregates-type parameter. Extends the functionalities of the former version, but is also in a way more restrictive:
#'  - Unlike the former, it takes a parameter's name, not value, thus argument parameters is essential.
#'  - It assigns its output to the global environment, which the other doesn't.
#'  - The resulting object is remanent.
#' The old function is still valid. Which one to use depends on exact context (input, aims).
#' 
#' @param param.nm The name of the parameter to parse.
#' @param dest The name of the vector of parsed parameters. Default = "Param.aggreg"
#' @param filt If provided, calls Param_filter() to remove these factors before parsing.
#' @param parameters The parameters file. Default = Param
#' @param aggregates The annotations (basic aggregate building blocks) as a named vector. Default: Aggregates
#' @param map The reference aggregates map. All existing aggregates in this map must exist. Default: Aggregate.map
#' @param aggr.list The list of aggregates. Default: Aggregate.list
#' @param parsed.param.nm The name of the parsed parameter variable to create. If not specified, defaults to param.nm
#' @param obj.list.nm The name of the list of objects to never purge from the global environment. Default = ".obj"
#' @param checkExp Logical, default = FALSE. Should we check the Experiment object, dropping it if it has only 1 level?
#' @param experiment Experiment factor object. Allows for dropping Exp if its length is one. Default = "Exp"
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

parse.Param.aggreg.2 <- function(param.nm,
                                 dest,
                                 filt = NULL,
                                 parameters = Param,
                                 aggregates = Aggregates,
                                 map = Aggregate.map,
                                 aggr.list = Aggregate.list,
                                 param.list.nm = "Param.aggreg",
                                 obj.list.nm = ".obj",
                                 checkExp = FALSE,
                                 experiment = Exp) {
  #DefArg(parse.Param.aggreg.2)
  if (!nchar(param.nm)) { return() }
  if (!param.nm %in% names(parameters)) {
    stop("Invalid parameter name!")
  }
  if (missing(dest)) { dest <- param.nm }
  paramValue <- parameters[[param.nm]]
  if (paramValue == "") {
    return()
  }
  if (is.na(paramValue)){
    warning("Invalid parameter value!")
    return()
  }
  if (!is.null(filt)) {
    paramValue <- Param_filter(paramValue, filt)
  }
  paramValue <- unlist(strsplit(paramValue, ";"))
  if ((checkExp)&&(length(experiment) == 1L)) { paramValue <- setdiff(paramValue, "Exp") }
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
  assign(dest, res, .GlobalEnv)
  assign(param.list.nm, union(get(param.list.nm, .GlobalEnv), dest), .GlobalEnv)
  assign(obj.list.nm, union(dest, get(obj.list.nm, .GlobalEnv)), .GlobalEnv)
  return()
}
