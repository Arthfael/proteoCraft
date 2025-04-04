
###############################################################################################################
#                                                                                                             #
# Download fresh versions of all reference proteomes for the organisms we have in our Fasta_databases folder. #
#                                                                                                             #
###############################################################################################################

# You could add new species simply by adding a new empty folder named after them and running this script
require(parallel)
today <- gsub("-", "", as.character(Sys.Date()))
# Create parallel processing cluster
N.clust <- detectCores()-1
a <- 1
tst <- try(clusterExport(parClust, "a", envir = environment()), silent = TRUE)
if ("try-error" %in% class(tst)) {
  try(stopCluster(parClust), silent = TRUE)
  parClust <- makeCluster(N.clust, type = "SOCK")
}
dbDir <- "D:/Fasta_databases"
stopifnot(dir.exists(dbDir))
Orgs <- grep("/", list.dirs(dbDir, full.names = FALSE), invert = TRUE, value = TRUE)
Orgs <- Orgs[which(nchar(Orgs) > 0)]
Orgs <- grep("[Cc]ontaminants", Orgs, invert = TRUE, value = TRUE)
if (interactive()) {
  Orgs <- svDialogs::dlg_list(Orgs, Orgs, TRUE, "Select species to update")$res
}
pack <- "myTAI"
if (!require(pack, character.only = TRUE)) {
  pak::pkg_install(pack, ask = FALSE, upgrade = TRUE, dependencies = TRUE)
}
if (!require(pack, character.only = TRUE)) {
  source(paste0(libPath, "/extdata/R scripts/Sources/taxonomy.R"))
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
#
# Fasta files
baseUrl <- c("https://rest.uniprot.org/uniprotkb/stream?compressed=false&download=true&format=fasta&includeIsoform=true&query=%28%28taxonomy_id%3A", "%29%29")
Urls <- paste0(baseUrl[1], sapply(Orgs, function(org) { rev(taxTst[[org]]$id)[1] }), baseUrl[2])
dstFls <- paste0(dbDir, "/", Orgs, "/", Orgs, "_UP_", today, "_Iso.fasta")
# Annotation txt files
txt_baseUrl <- c("https://rest.uniprot.org/uniprotkb/stream?download=true&format=txt&query=%28%28taxonomy_id%3A", "%29%29")
txt_Urls <- paste0(txt_baseUrl[1], sapply(Orgs, function(org) { rev(taxTst[[org]]$id)[1] }), txt_baseUrl[2])
txt_dstFls <- paste0(dbDir, "/", Orgs, "/", Orgs, "_UP_", today, "_Iso.txt")
#
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
lN <- length(wN)
if (lN) { dlTst <- curl::multi_download(Urls[wN], dstFls[wN], multi_timeout = Inf) }
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
  warning(paste(sapply(wN, function(x) { paste0("  -> ", Orgs[x], " = ", Urls[x], "\n") }), collapse = ""))
}
wY <- which(file.exists(dstFls))
if (length(wY)) {
  cat("Post-processing files (removing duplicates, adding contaminants)\n")
  for (i in wY) {
    nuFl <- gsub("\\.fasta$", "_noDupl_cont.fasta", dstFls[i])
    if (!file.exists(nuFl)) {
      cat(" - ", Orgs[i], "\n")
      # It would be more efficient to read the source only once...
      # but somehow it gets corrupted everytime I loop
      Src <- readLines(paste0(dbDir, "/Remove_duplicate_protein_sequences_and_fragments.R"))
      Src <- grep("openwd\\(", Src, invert = TRUE, value = TRUE)
      g <- grep("^dbFl +<- +", Src)
      stopifnot(length(g) == 1)
      Src[g] <- paste0("dbFl <- \"", dstFls[i], "\"")
      tmpFl <- tempfile(fileext = ".txt")
      write(Src, tmpFl)
      #system(paste0("open \"", tmpFl, "\""))
      #cat(paste(Src, collapse = "\n"))
      #xprs <- parse(text = Src)
      #eval(xprs, envir = .GlobalEnv)
      source(tmpFl)
      unlink(tmpFl)
    }
  }
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
if (lN) { warning(paste0(lN, " downloads failed, trying alternate approach...")) }
kount <- 0
clusterExport(parClust, list("txt_Urls", "txt_dstFls"), envir = environment())
msgs <- c("but I'm a hopeful little script and I'm not giving up...",
          "but I know I will make it",
          "it will all end well, right? Right...?",
          "the grind is real... one last try...",
          "and I'm too tired of this **** to go on trying, so download your stoopid files yourself!\n")
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
  warning(paste(sapply(wN, function(x) { paste0("  -> ", Orgs[x], " = ", txt_Urls[x], "\n") }), collapse = ""))
}
cat("Done!")
