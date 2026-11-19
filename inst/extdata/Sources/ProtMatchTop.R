if ((!exists("Update_Prot_matches")) || (!validLogicPar("Update_Prot_matches"))) { Update_Prot_matches <- TRUE }
if (Update_Prot_matches) {
  cat("Re-mapping peptide sequences to database proteins...\n")
  I_eq_L %<o% TRUE
  p2pFl <- paste0(wd, "/evmatch.RDS")
  if (file.exists(p2pFl)) { loadFun(p2pFl) } else {
    Seq <- unique(ev$Sequence)
    DB <- db
    Src <- paste0(libPath, "/extdata/Sources/ProtMatch.R")
    #rstudioapi::documentOpen(Src)
    source(Src)
    saveFun(evmatch, p2pFl)
  }
  ev$Search_engine_proteins <- ev$Proteins
  ev$Proteins <- evmatch$Proteins[match(ev$Sequence, evmatch$Sequence)]
  cat("Done!\n\n")
}
