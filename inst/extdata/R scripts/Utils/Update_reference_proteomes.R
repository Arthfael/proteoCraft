
###############################################################################################################
#                                                                                                             #
# Download fresh versions of all reference proteomes for the organisms we have in our Fasta_databases folder. #
#                                                                                                             #
###############################################################################################################

# You could add new species simply by adding a new empty folder named after them and running this script
RPath <- as.data.frame(library()$results)
RPath <- normalizePath(RPath$LibPath[match("proteoCraft", RPath$Package)], winslash = "/")
libPath <- paste0(RPath, "/proteoCraft")

today <- gsub("-", "", as.character(Sys.Date()))

# Create parallel processing cluster
require(parallel)
N.clust <- detectCores()-1
a <- 1
tst <- try(clusterExport(parClust, "a", envir = environment()), silent = TRUE)
if ("try-error" %in% class(tst)) {
  try(stopCluster(parClust), silent = TRUE)
  parClust <- makeCluster(N.clust, type = "SOCK")
}

dbDir <- "D:/Fasta_databases"
stopifnot(dir.exists(dbDir))
orgDirs <- grep("/", list.dirs(dbDir, full.names = FALSE), invert = TRUE, value = TRUE)
orgDirs <- orgDirs[which(nchar(orgDirs) > 0)]
orgDirs <- grep("[Cc]ontaminants", orgDirs, invert = TRUE, value = TRUE)
Orgs <- orgDirs
Orgs <- gsub("_aka_.*", "", Orgs)
Orgs <- gsub("\\([^\\)]+\\)*", "", Orgs)
orgDirs <- setNames(orgDirs, Orgs)
if (length(Orgs) != length(unique(Orgs))) {
  warning("It looks like you have multiple folders for the same organism!")
}
if (interactive()) {
  Orgs <- svDialogs::dlg_list(Orgs, Orgs, TRUE, "Select species to update")$res
}
pack <- "myTAI"
if (!require(pack, character.only = TRUE)) {
  try(pak::pkg_install("drostlab/myTAI", ask = FALSE, upgrade = TRUE, dependencies = TRUE), silent = TRUE) 
}
if (!require(pack, character.only = TRUE)) {
  Src2 <- paste0(libPath, "/extdata/R scripts/Sources/taxonomy.R")
  source(Src2, local = FALSE)
}
# Map taxonomy
cat("Mapping species names to taxIDs...\n")
taxTst <- try(setNames(lapply(Orgs, function(x) {
  suppressMessages(taxonomy(organism = x, db = "ncbi", output = "classification"))
}), Orgs), silent = TRUE)
if ("try-error" %in% class(taxTst)) {
  kount <- 0
  while ((kount < 20)&&("try-error" %in% class(taxTst))) {
    taxTst <- try(setNames(lapply(Orgs, function(x) {
      suppressMessages(taxonomy(organism = x, db = "ncbi", output = "classification"))
    }), Orgs), silent = TRUE)
  } 
  kount <- kount + 1
}
stopifnot(!"try-error" %in% class(taxTst))
tst <- sapply(taxTst, ncol)
Orgs <- names(taxTst)[which(tst == 3)] # Failures are data.frames with 1 column only!
taxTst <- taxTst[Orgs]
stopifnot(length(Orgs) > 0)
domainS <- setNames(vapply(taxTst, function(x) { x$name[match("domain", x$rank)] }, ""), Orgs)
taxIDs <- setNames(gsub(" ", "_", vapply(taxTst, function(x) { rev(x$id)[1] }, "")), Orgs)
taxLins <- setNames(vapply(taxTst, function(x) { paste(x$name, collapse = ", ") }, ""), Orgs)
#
# Get the info again from UniProt:
# We need some additional id which I could not get from taxonomy()
taxURLs <- paste0("https://rest.uniprot.org/proteomes/stream?&query=reference:true+taxonomy_id:", taxIDs,
                  "&fields=upid,lineage,organism_id&format=tsv")
taxURLs <- URLencode(taxURLs)
up_org_IDs <- setNames(lapply(1:length(Orgs), function(i) { #i <- 1 #i <- 12 #i <- 17
  #print(Orgs[i])
  url <- taxURLs[i]
  a <- try({
    res <- httr::GET(url)
    logFl <- paste0("tmp", i, ".txt")
    if (file.exists(logFl)) { unlink(logFl) }
    tst <- try(writeLines(httr::content(res, encoding = "UTF-8"), logFl), silent = TRUE)
    if (file.exists(logFl)) {
      tmp <- read.table(logFl, TRUE, sep = "\t", check.names = FALSE)
      if (nrow(tmp)) {
        w <- which(tmp$`Taxonomic lineage` == taxLins[i])
        if (length(w)) { w <- w[1] } else { w <- 1 }      
        x <- tmp[w, c("Proteome Id", "Organism Id")]
        x$Organism <- Orgs[i]
      }
      unlink(logFl)
      return(x[which(!is.na(x$"Proteome Id")),])
    } else {
      return()
    }
  }, silent = TRUE)
  if ("try-error" %in% class(a)) { return() }
}), Orgs)
up_org_IDs <- up_org_IDs[which(vapply(up_org_IDs, function(x) {
  ("data.frame" %in% class(x))&&(nrow(x) > 0)
}, TRUE))]
up_org_IDs <- do.call(plyr::rbind.fill, up_org_IDs)
nopeOrgs <- Orgs[which(!Orgs %in% up_org_IDs$Organism)]
Orgs <- up_org_IDs$Organism
#
# Fasta files
Urls_Orgs <- paste0("https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/",
                    domainS[Orgs], "/", up_org_IDs$`Proteome Id`, "/", up_org_IDs$`Proteome Id`, "_",
                    as.character(up_org_IDs$`Organism Id`), ".fasta.gz")
nopeOrgs_taxIDs <- vapply(nopeOrgs, function(org) { rev(taxTst[[org]]$id)[1] }, "")
Urls_nopeOrgs <- paste0("https://rest.uniprot.org/uniprotkb/stream?compressed=false&download=true&format=fasta&includeIsoform=true&query=%28%28taxonomy_id%3A",
                        nopeOrgs_taxIDs, "%29%29")

allOrgs <- c(Orgs, nopeOrgs)
Urls <- setNames(URLencode(c(Urls_Orgs, Urls_nopeOrgs)), allOrgs)
dstFls <- setNames(paste0(dbDir, "/", orgDirs[allOrgs], "/", allOrgs, "_UP_", today, "_Iso.fasta.gz"),
                   allOrgs)
stopifnot(sum(!dir.exists(dirname(dstFls))) == 0)
# Annotation txt files
txt_Urls_Orgs <- gsub("\\.fasta.gz$", ".dat.gz",  Urls_Orgs)
txt_Urls_nopeOrgs <- paste0("https://rest.uniprot.org/uniprotkb/stream?download=true&format=txt&query=%28%28taxonomy_id%3A",
                            nopeOrgs_taxIDs, "%29%29")
txt_Urls <- setNames(URLencode(c(txt_Urls_Orgs, txt_Urls_nopeOrgs)), allOrgs)
txt_dstFls <- c(gsub("\\.fasta.gz$", ".txt.gz",  dstFls[Orgs]),
                gsub("\\.fasta.gz$", ".txt",  dstFls[nopeOrgs]))
#
#unlink(dstFls)
wY <- which(file.exists(dstFls))
lY <- length(wY)
ovrwrt <- FALSE
if (lY) {
  if (interactive()) {
    msg <- paste0(c("Some", "All")[(lY == length(dstFls))+1], " of the fasta files you want to download already exist.\nOverwrite?")
    ovrwrt <- c(TRUE, FALSE)[match(svDialogs::dlg_message(msg, "yesno")$res, c("yes", "no"))]
  } else { ovrwrt <- TRUE }
  if (ovrwrt) { unlink(dstFls[wY]) }
}
wN <- which(!file.exists(dstFls))
if (length(wN)) { dlTst <- curl::multi_download(Urls[wN], dstFls[wN], multi_timeout = Inf) }
wN <- which(!file.exists(dstFls))
lN <- length(wN)
if (lN) { warning(paste0(lN, " downloads failed, trying alternate approach...")) }
kount <- 0
clusterExport(parClust, list("Urls", "dstFls"), envir = environment())
msgs <- c("but I'm a hopeful little script and I'm not giving up...",
          "but I know I will make it",
          "it will all end well, right? Right...?",
          "the grind is real... one last try...",
          "and I'm too tired of this **** to go on trying, so download your stoopid files yourself!\n")
while ((lN)&&(kount < 5)) {
  dlTst2 <- parLapply(parClust, wN, function(i) { #i <- wN[1]
    rs <- try(curl::curl_download(Urls[i], dstFls[i], mode = "w"), silent = TRUE)
    if ("try-error" %in% class(rs)) { rs <- list(Outcome = FALSE) } else {
      rs <- list(Outcome = TRUE, Res = rs)
    }
    return(rs)
  })
  wN <- which(!file.exists(dstFls))
  lN <- length(wN)
  kount <- kount + 1
  if (lN) { warning(paste0(lN, " downloads failed again, ", msgs[kount], "\n")) }
}
if (length(wN)) {
  warning(paste(sapply(wN, function(x) { paste0("  -> ", allOrgs[x], " = ", Urls[x], "\n") }), collapse = ""))
}
wY <- which(file.exists(dstFls))
if (length(wY)) {
  cat("Post-processing files (removing duplicates, adding contaminants)\n")
  clusterExprot(parClust, list("dstFls", "allOrgs", "dbDir"), envir = environment)
  parLapply(parClust, wY, function(i) { #i <- wY[1]
    origFl <- gsub("\\.fasta.gz$", ".fasta", dstFls[i])
    if (file.exists(dstFls[i])) {
      if (file.exists(origFl)) { unlink(origFl) }
      R.utils::gunzip(dstFls[i])
      if (file.exists(dstFls[i])) { unlink(dstFls[i]) }
    }
    nuFl <- gsub("\\.fasta$", "_noDupl_cont.fasta.gz", origFl)
    if (!file.exists(nuFl)) {
      cat(" - ", allOrgs[i], "\n")
      # It would be more efficient to read the source only once...
      # but somehow it gets corrupted everytime I loop
      Src <- readLines(paste0(dbDir, "/Remove_duplicate_protein_sequences_and_fragments.R"))
      Src <- grep("openwd\\(", Src, invert = TRUE, value = TRUE)
      g <- grep("^dbFl +<- +", Src)
      stopifnot(length(g) == 1)
      Src[g] <- paste0("dbFl <- \"", origFl, "\"")
      tmpFl <- tempfile(fileext = ".txt")
      write(Src, tmpFl)
      #system(paste0("open \"", tmpFl, "\""))
      #cat(paste(Src, collapse = "\n"))
      #xprs <- parse(text = Src)
      #eval(xprs, envir = .GlobalEnv)
      source(tmpFl)
      unlink(tmpFl)
    }
  })
}

# Also get the txt files!
wY <- which(file.exists(txt_dstFls))
lY <- length(wY)
if (ovrwrt) { unlink(txt_dstFls[wY]) }
wN <- which(!file.exists(txt_dstFls))
lN <- length(wN)
if (lN) { dlTst <- curl::multi_download(txt_Urls[wN], txt_dstFls[wN], multi_timeout = Inf) }
wN <- which(!file.exists(txt_dstFls))
lN <- length(wN)
if (lN) {
  warning(paste0(lN, " downloads failed, trying alternate approach..."))
  kount <- 0
  clusterExport(parClust, list("txt_Urls", "txt_dstFls"), envir = environment())
  msgs <- c("but I'm a hopeful little script and I'm not giving up...",
            "but I know I will make it",
            "it will all end well, right? Right...?",
            "the grind is real... one last try...",
            "and I'm too tired of this **** to go on trying, so download your stoopid files yourself!\n")
}
while ((lN)&&(kount < 5)) {
  dlTst2 <- parLapply(parClust, wN, function(i) { #i <- wN[1]
    rs <- try(curl::curl_download(txt_Urls[i], txt_dstFls[i], mode = "w"), silent = TRUE)
    if ("try-error" %in% class(rs)) { rs <- list(Outcome = FALSE) } else {
      rs <- list(Outcome = TRUE, Res = rs)
    }
    return(rs)
  })
  wN <- which(!file.exists(txt_dstFls))
  lN <- length(wN)
  kount <- kount + 1
  if (lN) { warning(paste0(lN, " downloads failed again, ", msgs[kount], "\n")) }
}
if (length(wN)) {
  warning(paste(sapply(wN, function(x) { paste0("  -> ", allOrgs[x], " = ", txt_Urls[x], "\n") }), collapse = ""))
}
wY <- which(file.exists(txt_dstFls))
clusterExprot(parClust, list("txt_dstFls"), envir = environment)
parLapply(parClust, wY, function(i) { #i <- wY[1]
  origFl <- gsub("\\.txt.gz$", ".txt", txt_dstFls[i])
  if (file.exists(txt_dstFls[i])) {
    if (file.exists(origFl)) { unlink(origFl) }
    R.utils::gunzip(txt_dstFls[i])
    if (file.exists(txt_dstFls[i])) { unlink(txt_dstFls[i]) }
  }
})
cat("Done!")
