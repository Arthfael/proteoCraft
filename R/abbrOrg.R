#' abbrOrg
#'
#' @description
#' A utility function to abbreviate strings of Organism names.
#'  e.g. c("Homo sapiens;Mus musculus (C57Bl6)", "Escherichia coli") -> c("H. sapiens;M. musculus", "E. coli")
#' 
#' @param x String to abbreviate.
#' @param sep Name separator, default = ";"
#' 
#' @returns
#' String of abbreviated names.
#' 
#' @examples
#' x <- c("Homo sapiens;Mus musculus (C57Bl6)", "Escherichia coli")
#' abbrOrg(x)
#' 
#' @export

abbrOrg <- function(x,
                    sep = ";") { #x <- tstOrg2
  x <- strsplit(x, sep)
  return(vapply(x, \(y) { #y <- x[[1]]
    y <- strsplit(gsub(" *\\([^\\)]*\\)", "", y), "")
    y <- vapply(y, \(z) { #z <- y[[1]]
      g1 <- grep("[A-Z]", z)
      g2 <- grep("[a-z]", z)
      w1 <- g2[which((g2-1L) %in% g1)]
      z[w1] <- "."
      gsub("\\.[^ ]+", ".", paste(z, collapse = ""))
    }, "")
    paste(y, collapse = sep)
  }, ""))
}
