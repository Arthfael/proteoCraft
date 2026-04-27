#### Code chunk - peptide tables for visualizing the coverage of proteins of interest in 3D using SCV
if ((!is.null(prot.list))&&(length(prot.list))) {
  # From https://stackoverflow.com/questions/52911812/check-if-url-exists-in-r
  valid_url <- function(url_in, t = 2L){
    con <- url(url_in)
    check <- suppressWarnings(try(open.connection(con, open = "rt", timeout = t), silent = TRUE)[1L])
    suppressWarnings(try(close.connection(con), silent = TRUE))
    ifelse(is.null(check), TRUE, FALSE)
  }
  # Mods and their mass shifts
  SCV_PTMs <- TRUE
  if ((!"Mass shift" %in% colnames(Modifs))||(sum(is.na(Modifs$`Mass shift`)))) {
    if ("UniMod" %in% colnames(Modifs)) {
      if (!require("PTMods", quietly = TRUE)) { pak::pak("rformassspectrometry/PTMods") }
      require(PTMods)
      data(modifications, package = "PTMods")
      UniMod <- modifications
      Modifs$"Mass shift" <- UniMod$MonoMass[match(as.integer(gsub("^UniMod:", "", Modifs$UniMod)), UniMod$UnimodId)]
    } else {
      kols <- c("Mass delta", "Delta mass")
      w <- which(vapply(kols, \(k) { (k %in% colnames(Modifs))&&(is.numeric(Modifs[[k]])) }, TRUE))
      if (length(w)) {
        Modifs$"Mass shift" <- Modifs[[kols[w]]]
      } else {
        warning("Could not map PTMs to mass shifts, these will be ignored from the SCV visualisations.")
        SCV_PTMs <- FALSE
      }
    }
  }
  g1 <- grsep2(prot.list, pep$Proteins)
  prVect <- pep$Proteins[g1]
  modSq <- pep$"Modified sequence"[g1]
  if (scrptType == "withReps") {
    intVect <- setNames(lapply(RSA$value, \(x) {
      pep[g1, paste0(pep.ref["Original"], x)]
    }), RSA$value)
  }
  if (scrptType == "noReps") {
    intVect <- setNames(lapply(Exp, \(x) {
      pep[g1, paste0(int.cols["Original"], " - ", x)]
    }), Exp)
  }
  PDB_in_DB <- ("PDB" %in% colnames(db))
  g2 <- grsep2(prot.list, db$`Protein ID`)
  dbPDB <- db$PDB[g2]
  dbPID <- db$"Protein ID"[g2]
  source(parSrc)
  clusterExport(parClust,
                list("grsep", "grsep2", "listMelt", "cov3D", "annot_to_tabl",
                     "prVect", "intVect", "wd", "PDB_in_DB", "SCV_PTMs", "dbPDB", "dbPID", "Modifs", "valid_url", "modSq",
                     "AA_table", "AA"),
                envir = environment())
  tst3D <- parSapply(parClust, prot.list, \(plp) { #plp <- prot.list[1L]
    pdbFls <- c()
    grs <- grsep2(plp, prVect)
    if (!length(grs)) { return(FALSE) }
    dir <- paste0(wd, "/Coverage/", plp)
    if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
    # For each protein we want to download:
    # - all models for all fragments.
    # - latest version only!
    # PDB: we get PDB IDs from parsing the txt file
    if (PDB_in_DB) {
      tmp <- unlist(strsplit(dbPDB[match(plp, dbPID)], ";"))
      tmp <- tmp[which(tmp != "")]
      if (length(tmp)) {
        for (i in tmp) {
          fl <- paste0(dir, "/PDB-", i, ".pdb")
          url <- paste0("https://files.rcsb.org/download/", i, ".pdb")
          tst <- try(utils::download.file(url, fl), silent = TRUE)
          if (!inherits(tst, " try-error")) {
            pdbFls <- union(pdbFls, fl)
          }
        }
      }
    }
    # Then AlphaFold
    mdlNms <- unlist(lapply(1L:3L, \(x) { # (We will test up to 3 fragments)
      # Currently (04/2026), alphaFold reaches up to version 6
      paste0("AF-", plp, "-F", as.character(x), "-model_v", as.character(1L:6L), ".pdb")
    }))
    alphURLsTbl <- data.frame(Name = mdlNms,
                              url = paste0("https://alphafold.ebi.ac.uk/files/", mdlNms),
                              Frag = unlist(lapply(1L:3L, \(x) { rep(x, 6L) })),
                              vers = 1L:6L)
    alphURLsTbl$Valid <- vapply(alphURLsTbl$url, valid_url, TRUE)
    if (sum(alphURLsTbl$Valid)) {
      alphURLsTbl <- alphURLsTbl[which(alphURLsTbl$Valid),]
      wh <- aggregate(1L:nrow(alphURLsTbl), list(alphURLsTbl$Frag), \(x) {
        v <- alphURLsTbl$vers[x]
        x[which(v == max(v))]
      })$x
      alphURLs <- setNames(alphURLsTbl$url[wh], paste0(dir, "/", alphURLsTbl$Name[wh]))
      lapply(names(alphURLs), \(x) {
        utils::download.file(alphURLs[x], x)
      })
      pdbFls <- union(pdbFls, names(alphURLs))
    }
    if (!length(pdbFls)) { return(FALSE) }
    # We have found at least one model which can be used to visualize coverage
    # Let's write peptidoforms
    seq1 <- gsub("_", "", modSq[grs])
    if (SCV_PTMs) {
      seq2 <- gsub("\\)", "]_", gsub("\\(", "_[", seq1))
      seq2 <- strsplit(seq2, "_")
      seq2 <- vapply(seq2, \(x) { #x <- seq2[1L]
        x <- unlist(x)
        w <- grep("\\[.+\\]", x)
        x[w] <- paste0("[", round(Modifs$`Mass shifts`[match(x[w], paste0("[", Modifs$Mark, "]"))], 0L), "]")
        return(paste(x, collapse = ""))
      }, "")
    } else { seq2 <- gsub("[^A-Z]", "", seq2) }
    write(seq2, paste0(dir, "/SCV - observed peptides.txt"))
    lapply(pdbFls, \(fl) { #fl <- pdbFls[1L]
      lapply(names(intVect), \(x) { #x <- names(intVect)[1L]
        nm <- gsub(".*/|\\.pdb$", "", fl)
        pth <- gsub("\\.pdb$", ".html", fl)
        cov3D(fl, seq1, path = pth, ttl = nm, intensities = intVect[[x]][grs], display = FALSE)
      })
    })
    return(TRUE)
  })
  # If this worked, we write a guide in the Coverage folder
  if (sum(tst3D)) {
    Guide <- c("Visualing protein coverage in 3D using SCV",
               "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~",
               "",
               "If, for protein accession of interest \"PROTEIN\", a 3D model is available, then the following files are created in its subfolder:",
               " - \"peptidoforms.txt\": sequences of identified peptides",
               " - For each fragment number i for which a structure is available (ranging from 1 to n), \"AF-PROTEIN-F1-model_v#.pdb\", where \"#\" is the latest valid version.",
               "",
               "To visualise the peptides onto the folded protein structure, navigate to https://scv.lab.gy and:",
               " - paste the peptides into the \"PSM/peptide list\" field",
               " - load the pdb ",
               "")
    write(Guide, paste0(wd, "/Coverage/SCV - how to visualise protein coverage in 3D.txt"))
  }
}
