#' %<c%
#'
#' @description 
#' These two custom assignment operators are meant to be used within the data analysis scripts.
#' Both scripts regularly purge their environment of temporary files, to avoid unnecessary inflation of data.
#' Both operators are meant to create variables which are proof against this purge.
#' 
#' Each operator:
#'  - assigns the value on the right to the variable on the left, creating the variable if necessary,...
#'  - ... and in addition adds the variable's name to a hidden variable, .obj, which lists all variables which should escape environment purges.
#'  - If .obj does not exist, it is created.
#'  
#'  Both operators differ slightly in usage:
#'  - "\%<o\%" takes the following syntax:
#'          variable \%<o\% value
#'  - "\%<c\%" takes the following syntax:
#'          "variable" \%<o\% value
#' 
#' Careful!!!
#' Those operators have the following limitations:
#'  - repeated assignment within the same line do not work!!!
#'  - if the value on the right is a logical test, it should be surrounded by parentheses!!!
#' 
#' @param x Variable to create, see description for details.
#' @param y Value to assign.
#' 
#' @examples
#' exists(".obj")
#' exists("a")
#' exists("b")
#' a %<o% "test 1"
#' "b" %<c% "test 2"
#' exists("a")
#' exists("b")
#' a
#' b
#' rm(list = ls()[which(!ls() %in% .obj)])
#' a
#' b
#' .obj
#'  
#' @export

"%<c%" <- function(x, y) { # When the variable is created by name stored as the value of x
  nm <- x
  assign(nm, y, envir = .GlobalEnv)
  if (!exists(".obj")) { .obj <<- c(".obj", "%<o%", "%<c%") }
  .obj <<- unique(c(.obj, nm))
}
