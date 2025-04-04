#' GO_annotate_DB
#'
#' @description 
#' A function to annotate a vector of proteins with GO terms.
#' 
#' @param proteins The character vector of UniProt protein IDs to annotate.
#' @param species The species of the proteins.
#' @param type The type of organism.
#' 
#' @examples
#' annot <- GO_annotate_DB(db$Protein.ID, "Homo sapiens")
#' 
#' @export

GO_annotate_DB <- function (proteins,
                            species,
                            type = "Standard") {
  library(GO.db)
  proteins <- data.frame(`Protein ID` = proteins, check.names = FALSE)
  bm.host <- c("www.ensembl.org",
               "plants.ensembl.org",
               "protists.ensembl.org")[match(type, c("Standard",
                                                     "Plant",
                                                     "Protist"))]
  bm.mart <- c("ENSEMBL_MART_ENSEMBL", "plants_mart", "protist_mart")[which(c("Standard",  "Plant", "Protist") == type)]
  a <- biomaRt::listDatasets(useMart(bm.mart, host = bm.host))
  b <- tolower(unlist(strsplit(species, " ")))
  b <- paste0(substr(b[1], 1, 1), b[2])
  bm.dataset <- grep(b, a[, 1], value = TRUE)
  a <- biomaRt::listAttributes(biomaRt::useDataset(bm.dataset,
                                                   biomaRt::useMart(bm.mart,
                                                                    host = bm.host)))
  mart <- biomaRt::useDataset(bm.dataset,
                              biomaRt::useMartuseMart(bm.mart,
                                                      host = bm.host))
  q <- c(unique(a$name[grep("UNIPROT", toupper(a$name))]), 
         unique(a$name[grep("GO_", toupper(a$name))]))
  GO <- setNames(c("BP", "CC", "MF"),
                 c("biological_process", "cellular_component",  "molecular_function"))
  res <- list()
  if (type == "Plant") {
    temp <- biomaRt::getBM(q,
                           "uniprot_swissprot_accession",
                           proteins$"Protein ID",
                           mart)
    if (nrow(temp) > 0) {
      for (i in 1:length(GO)) {
        temp2 <- temp[match(names(GO)[i], temp$go_namespace_1003),]
        temp2 <- temp2[which(temp2$uniprot_swissprot_accession %in%  proteins$"Protein ID"), ]
        proteins <- proteoCraft::col.aggr(proteins,
                                    temp2,
                                    c("Protein ID", paste0("go_", GO[i], c("", "_name"))), 
                                    c("uniprot_swissprot_accession",
                                      "go_accession",
                                      "go_name_1006"),
                                    function(x) { paste(sort(unique(x)), collapse = ";") },
                                    Coerce = TRUE)
        res[[paste0("GO_", GO[i])]] <- temp2
      }
      res$Annotated_Proteins <- proteins
    }
    else {
      res <- "No annotations retrieved!"
      warning(res)
    }
  } else {
    temp <- biomaRt::getBM(q,
                           "uniprotsptrembl",
                           proteins$"Protein ID",
                           mart)
    if (nrow(temp) > 0) {
      xx <- AnnotationDbi::as.list(GO.db::GOTERM)
      xx <- xx[which(names(xx) %in% temp$go_id)]
      temp2 <- as.data.frame(t(sapply(temp$go_id, function(x) {
        x <- which(names(xx) == x)
        if (length(x) == 0) { return(c(NA, NA)) } else {
          return(c(Term(xx[[x]]), Ontology(xx[[x]])))
        }
      })))
      temp[, c("Term", "Ontology")] <- temp2
      temp1 <- temp[, c("uniprotsptrembl", "go_id", "Term",  "Ontology")]
      colnames(temp1)[match("uniprotsptrembl", colnames(temp1))] <- "Protein ID"
      for (i in c("uniprot_gn", "uniprotswissprot")) {
        if (i %in% colnames(temp)) {
          temp2 <- temp[, c(i, "go_id", "Term", "Ontology")]
          colnames(temp2)[match(i, colnames(temp2))] <- "Protein ID"
          temp1 <- rbind(temp1, temp2)
        }
      }
      temp <- temp1
      for (go in GO) {
        temp2 <- temp[which(temp$Ontology == go), ]
        temp2 <- temp2[which(temp2$"Protein ID" %in% proteins$"Protein ID"), ]
        proteins <- proteoCraft::col.aggr(proteins,
                                    temp2,
                                    c("Protein.ID", paste0("go_", go, c("", "_name"))),
                                    c("Protein ID", "go_id", "Term"),
                                    function(x) { paste(sort(unique(x)), collapse = ";") },
                                    Coerce = TRUE)
        res[[paste0("GO_", go)]] <- temp2
      }
      res$Annotated_Proteins <- proteins
    }
    else {
      res <- "No annotations retrieved!"
      warning(res)
    }
  }
  return(res)
}

