intPrtFst %<o% paste0(wd, "/Proteins of interest.fasta")
if (!exists("prot.list")) { prot.list %<o% c() }
if (!exists("prot.list_pep")) { prot.list_pep %<o% c() }
if ((exists("db"))&&("Protein of interest" %in% colnames(db))) {
  prot.list <- unique(c(prot.list, db$`Protein ID`[which(db$`Protein of interest`)]))
}
if ((exists("Param"))&&("Prot.list" %in% colnames(Param))) {
  prot.list <- unique(c(prot.list, unlist(strsplit(Param$Prot.list, ";"))))
}
if ((exists("Param"))&&("Prot.list_pep" %in% colnames(Param))) {
  prot.list_pep <- unlist(strsplit(Param$Prot.list_pep, ";"))
  prot.list <- unique(c(prot.list, prot.list_pep))
}
if ((exists("Param"))&&("Norma.Prot.Ratio.to.proteins" %in% colnames(Param))) {
  prot.list <- unique(c(prot.list, unlist(strsplit(Param$Norma.Prot.Ratio.to.proteins, ";"))))
}
if (exists("TargetProteins")) { prot.list <- unique(c(prot.list, TargetProteins)) }
prot.list %<o% prot.list
# Update backup fasta
if (length(prot.list)) { writeFasta(db[match(prot.list, db$`Protein ID`),], intPrtFst) }
