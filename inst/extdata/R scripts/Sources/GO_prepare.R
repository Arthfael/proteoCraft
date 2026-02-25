# Prepare GO term maps
if (Annotate) {
  Ontologies %<o% setNames(c("BP", "CC", "MF"), c("Biological Process", "Cell Compartment", "Molecular Function"))
  packs <- c("annotate", "GO.db")
  for (pack in packs) { biocInstall(pack) }
  for (pack in packs) { library(pack, character.only = TRUE) }
  if (file.exists("GO_terms.RData")) { loadFun("GO_terms.RData") }
  if (!exists("GO_terms")) {
    GO_terms <- data.frame(ID = gsub(".+ \\[|\\]", "", unique(unlist(strsplit(db$GO, ";")))),
                           Term = unique(unlist(strsplit(db$GO, ";"))))
    GO_terms <- GO_terms[which(!is.na(GO_terms$ID)),]
    GO_terms <- GO_terms[grep("^GO:[0-9]{7}$", GO_terms$ID),]
    GO_terms$Ontology <- NA
    GO_terms$Offspring <- list(NA)
    for (ont in Ontologies) { #ont <- Ontologies[1]
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
        Offspr <- Offspr[, list(Offspring = list(Offspring)), by = list(Parent = Parent)]
        Offspr <- as.data.frame(Offspr)
        GO_terms$Offspring[wo] <- Offspr$Offspring[match(GO_terms$ID[wo], Offspr$Parent)]
      }
    }
    GO_terms$Offspring <- parLapply(parClust, GO_terms$Offspring, function(x) {
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
    for (pack in packs) { biocInstall(pack) }
    for (pack in packs) { library(pack, character.only = TRUE) }
    GO_mappings <- GO_map(db, cl = parClust)$Mappings
    saveFun(GO_mappings, file = "GO_mappings.RData")
  }
  GO_mappings %<o% GO_mappings
}
