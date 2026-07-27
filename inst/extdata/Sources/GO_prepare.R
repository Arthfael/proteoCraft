# Prepare GO term maps
if (Annotate) {
  Ontologies %<o% setNames(c("BP", "CC", "MF"), c("Biological Process", "Cell Compartment", "Molecular Function"))
  packs <- c("annotate", "GO.db", "AnnotationDbi")
  for (pck in packs) { if (!require(pck, character.only = TRUE)) { pak::pak(pck) } }
  for (pck in packs) { library(pck, character.only = TRUE) }
  if (file.exists("GO_terms.RData")) { loadFun("GO_terms.RData") }
  if (!exists("GO_terms")) {
    # Terms from AnnotationDbi
    GO_terms <- AnnotationDbi::select(GO.db,
                                keys = keys(GO.db),
                                columns = c("TERM", "ONTOLOGY"),
                                keytype = "GOID")
    GO_terms$TERM <- paste0(do.call(paste, c(GO_terms[, c("TERM", "GOID")], sep = " [")), "]")
    GO_terms <- GO_terms[, c("GOID", "TERM", "ONTOLOGY")]
    colnames(GO_terms) <- c("ID", "Term", "Ontology")
    # Also get terms from our own annotations
    go <- data.frame(ID = gsub(".+ \\[|\\]", "", unique(unlist(strsplit(db$GO, ";")))),
                           Term = unique(unlist(strsplit(db$GO, ";"))))
    go <- go[which(!is.na(GO_terms$ID)),]
    go <- go[grep("^GO:[0-9]{7}$", go$ID),]
    go <- go[which(!go$ID %in% GO_terms$ID),]
    if (nrow(go)) {
      go$Ontology <- NA_character_
      GO_terms <- rbind(go, GO_terms)
    }
    GO_terms$Offspring <- list(NA_character_)
    for (ont in Ontologies) { #ont <- Ontologies[1L]
      wo <- which(filterGOByOntology(GO_terms$ID, ont))
      if (length(wo)) {
        #cat(paste0("Getting offspring for ", length(wo), " ", ont, " terms...\n"))
        GO_terms$Ontology[wo] <- ont
        Offspr <- get(paste0("GO", ont, "OFFSPRING"))
        # Finally I managed to rewrite the tedious topGO code much faster,
        # and without even using parallel!!!
        Offspr <- toTable(Offspr)
        colnames(Offspr) <- c("Offspring", "Parent")
        Offspr <- as.data.table(Offspr)
        Offspr <- Offspr[, .(Offspring = list(Offspring)), by = .(Parent = Parent)]
        Offspr <- as.data.frame(Offspr)
        GO_terms$Offspring[wo] <- Offspr$Offspring[match(GO_terms$ID[wo], Offspr$Parent)]
      }
    }
    GO_terms$Offspring <- parLapply(parClust, GO_terms$Offspring, \(x) {
      x[which(!is.na(x))]
    })
    GO_terms <- GO_terms[which(!is.na(GO_terms$Ontology)),]
    GO_terms <- GO_terms[order(GO_terms$Ontology),]
    saveFun(GO_terms, file = "GO_terms.RData")
  }
  GO_terms %<o% GO_terms
  #loadFun("GO_terms.RData")
  if (file.exists("GO_mappings.RData")) { loadFun("GO_mappings.RData") }
  if (!exists("GO_mappings")) {
    packs <- c("GO.db", "topGO")
    for (pck in packs) { if (!require(pck, character.only = TRUE)) { pak::pak(pck) } }
    for (pck in packs) { library(pck, character.only = TRUE) }
    GO_mappings <- GO_map(db, cl = parClust)$Mappings
    saveFun(GO_mappings, file = "GO_mappings.RData")
  }
  GO_mappings %<o% GO_mappings
}
