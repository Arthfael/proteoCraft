#' GO_annotate2
#'
#' @description 
#' A function to annotate a vector of proteins with GO terms. Can also do KEGG and Reactome pathways, despite its more restrictive name.
#' This was meant to be an improved version of the previous one, which did not work well at all...
#' but I think does not work now, at least I do not think it is being used.
#' Candidate for deletion.
#' 
#' Returns the vector of protein IDs as data.frame with added GO columns.
#' 
#' @param proteins The character vector of UniProt protein IDs to annotate.
#' @param species The species of the proteins.
#' @param chunks Optional, but set by default to 6. Size of the chunks to attempt to download. Emergency argument, normally data should be querried in one batch but this currently fails. May be set to missing by default in future.
#' @param tries When working by chunks, how many times should I attempt each chunk before skipping?
#' @param kegg Default = FALSE. Should KEGG annotations be downloaded too?
#' @param reactome Default = FALSE. Should Reactome annotations be downloaded too?
#' 
#' @examples
#' test <- GO_annotate2(db$Protein.ID, "Homo sapiens")
#' 
#' @export

GO_annotate2 <- function(proteins, species, chunks = 6, tries = 5, kegg = FALSE, reactome = FALSE) {
  sp_info <- availableUniprotSpecies(pattern = paste("^", species, "$", sep = ""), n = Inf)
  if (nrow(sp_info) != 1) {
    if (nrow(sp_info) > 1) {
      stop(c("Ambiguous species name provided! Potential matches: ", paste(sp_info$`Species name`, collapse = ", "), "!"))
    } else {
      stop("The requested species was not found on UniProt! Did you provide a name in the \"Genus species\" form?")
    }
  } else {
    tx <- as.numeric(sp_info$`taxon ID`)
    up <- UniProt.ws::UniProt.ws(taxId = tx)
    no.iso <- gsub("-[0-9]+$", "", proteins)
    prot <- unique(no.iso)
    if (missing("chunks")) {
      Res <- UniProt.ws::select(up, keys = prot, columns = c("GO","GO-ID"), keytype = "UNIPROTKB")
    } else {
      stopifnot(is.numeric(chunks), chunks == round(chunks))
      s <- c(0:ceiling(length(prot)/chunks))*chunks+1
      s[length(s)] <- length(prot)+1
      Res <- data.frame("UNIPROTKB" = prot, "GO" = rep(NA, length(prot)), I_hate_this = rep(NA, length(prot)))
      # Seriously, the way data frames names handle special characters can be frustrating...
      colnames(Res)[3] <- "GO-ID"
      Col <- c("GO","GO-ID")
      if (kegg) {
        Res$KEGG <- rep(NA, length(prot))
        Col <- c(Col, "KEGG")
      }
      if (reactome) {
        Res$Reactome <- rep(NA, length(prot))
        Col <- c(Col, "REACTOME")
      }
      for (i in 1:(length(s)-1)) {
        #i <- 1
        range <- s[i]:(s[i+1]-1); l <- length(range)
        res <- try(UniProt.ws::select(up, keys = prot[range], columns = Col, keytype = "UNIPROTKB"),
                   silent = TRUE)
        if (class(res) == "try-error") {
          res <- data.frame("UNIPROTKB" = prot[range], "GO" = rep(NA, l),
                            No_really_I_hate_this = rep(NA, l))
          # Same issue as above...
          colnames(res)[3] <- "GO-ID"
          if (kegg) { res$KEGG <- rep(NA, l) }
          if (reactome) { res$Reactome <- rep(NA, l) }
          for (j in 1:l) {
            k <- 0; ok <- FALSE
            while ((!ok)&&(k < tries)) {
              r <- try(UniProt.ws::select(up, keys = prot[range[j]], columns = Col, keytype = "UNIPROTKB"),
                       silent = TRUE)
              if (class(r) != "try-error") { ok <- TRUE; res[j,] <- r }
              k <- k+1
            }
          }
        }
        Res[range,] <- res
      }
    }
    Res <- Res[match(no.iso, Res$UNIPROTKB),]
    Res$UNIPROTKB <- proteins
    return(Res)
  }
}
