#' AA_bias
#' 
#' @description
#' A little function to see if there is a bias against certain amino acids in the dataset, versus parent search database.
#' 
#' @param Ev The evidence data frame ; must contain a "Sequence" column.
#' @param DB The search database, processed as a data.frame by Format.DB; must contain a "Sequence" column.
#' 
#' @examples
#' AA.bias <- AA_bias(Ev = ev, DB = db)
#' 
#' @import data.table
#' @export

AA_bias <- function(Ev,
                    DB) {
  E <- unlist(strsplit(toupper(Ev$Sequence), ""))
  AA <- data.table::data.table(E = E, dumdum = as.integer(rep(1, length(E))))
  AA <- AA[, list(`Occurrences - dataset` = sum(dumdum)), keyby = list(AA = E)]
  AA <- as.data.frame(AA)
  AA$"% - dataset" <- 100*AA$`Occurrences - dataset`/sum(AA$`Occurrences - dataset`)
  D <- unlist(strsplit(toupper(DB$Sequence), ""))
  AA_db <- data.table::data.table(D = D, dumdum = as.integer(rep(1, length(D))))
  AA_db <- AA_db[, list(x = sum(dumdum)), keyby = list(Group.1 = D)]
  AA_db <- as.data.frame(AA_db)
  AA$"Occurrences - database" <- AA_db$x[match(AA$AA, AA_db$Group.1)]
  AA$"% - database" <- 100*AA$`Occurrences - database`/sum(AA$`Occurrences - database`)
  AA$Ratio <- signif(AA$"% - dataset"/AA$"% - database", 3)
  AA$"% - dataset" <- signif(AA$"% - dataset", 3)
  AA$"% - database" <- signif(AA$"% - database", 3)
  return(AA)
}
