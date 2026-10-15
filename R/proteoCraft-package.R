#' proteoCraft: Proteomics Data Analysis
#'
#' Scripts and functions for processing, normalizing, quantifying and statistically analyzing proteomics data, starting from PSMs level.
#'
#'
#' @keywords internal
#' @importFrom grDevices dev.list dev.off graphics.off rainbow
#' @importFrom graphics layout
#' @importFrom methods as getFunction
#' @importFrom stats aggregate as.dendrogram as.formula dist hclust lm.fit mad median p.adjust prcomp pt qnorm qt sd setNames time weighted.mean
#' @importFrom utils capture.output data installed.packages object.size packageVersion read.csv read.delim sessionInfo write.csv write.table
"_PACKAGE"

utils::globalVariables(".obj")
