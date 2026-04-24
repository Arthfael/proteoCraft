#### Load and process search database(s) (or reload a pre-processed version if available)
SSH_on %<o% FALSE
fastas <- unique(fastas)
#fastas <- fastas$Full
#
# Function to only display the useful part of paths
pathAbbr <- \(paths) {
  #pathAbbr(c("C:/file.txt", "D:/file.txt"))
  #pathAbbr(c("C:/file.txt", "C:/file.gff"))
  #pathAbbr(c("C:/file.txt", "C:/file2.txt"))
  #pathAbbr(c("C:/file.txt", "C:/User/file.txt"))
  l1 <- length(paths)
  res <- gsub(".*/", "", paths)
  l2 <- length(unique(res))
  if ((l1 == 1L)||(l1 == l2)) { return(res) }
  tst <- strsplit(gsub("/[^/]+$", "", paths), "/")
  l <- lengths(tst)
  m <- min(l)
  tst2 <- vapply(1L:m, \(i) {
    length(unique(vapply(tst, \(x) { x[[i]] }, ""))) == 1L
  }, TRUE)
  wN <- suppressWarnings(min(which(!tst2)))
  wY <- which(tst2)
  if (length(wY)) { wY <- max(wY[which(wY < wN)]) }
  if (length(wY)) {
    rmv <- paste(tst[[1L]][wY], collapse = "/")
    nc <- nchar(rmv) + 1L
    res <- paste0("...", substring(paths, nc))
  } else {
    res <- paths
  }
  return(res)
}
fastasTbl %<o% data.frame(Full = fastas,
                          Name = basename(fastas),
                          Dir = dirname(fastas),
                          Contaminant = FALSE)
if (exists("isContaminant")) { 
  fastasTbl$Contaminant <- as.logical(isContaminant[fastasTbl$Full])
}
fastasTbl$Contaminant[which(is.na(fastasTbl$Contaminant))] <- FALSE
fastasTbl$Loc <- c("Local", "Cluster")[grepl("^/+nfs/", fastasTbl$Dir)+1L]
fastasTbl$Exists <- file.exists(fastasTbl$Full)
fastasTbl$ExistsHere <- file.exists(paste0(wd, "/", fastasTbl$Name))
fastasTbl$TXT <- gsub("\\.fa((s(ta(\\.fas)?)?)|(a))?$", ".txt", fastasTbl$Full) # Should catch .fasta, .fa, .faa, .fas, .fasta.fas
fastasTbl$TXTHere <- paste0(wd, "/", gsub("\\.fa((s(ta(\\.fas)?)?)|(a))?$", ".txt", fastasTbl$Name))
fastasTbl$TXTExists <- grepl("\\.txt$", fastasTbl$TXT)&(file.exists(fastasTbl$TXT))
fastasTbl$TXTExistsHere <- grepl("\\.txt$", fastasTbl$TXTHere)&(file.exists(fastasTbl$TXTHere))
w1 <- which((fastasTbl$ExistsHere)&(fastasTbl$Loc == "Cluster")) # Conflict between here and cluster
if (length(w1)) {
  require(tools)
  if (!SSH_on) {
    require(ssh)
    msg <- "SSH session required\nEnter host name!"
    sshost <- dlg_input(msg)$res
    sshsess <- ssh_connect(sshost)
    kount <- 1L
    while ((kount < 5L)&&(!inherits(sshsess, "ssh_session"))) {
      sshsess <- ssh_connect(sshost)
      kount <- kount + 1L
    }
    #print(sshsess)
    SSH_on <- inherits(sshsess, "ssh_session")
  }
  if (SSH_on) {
    tst1 <- md5sum(fastasTbl$Name[w1])
    tst2 <- vapply(fastasTbl$Full[w1], \(x) { #x <- fastasTbl$Full[w1[1L]]
      gsub(" .*", "", capture.output(ssh_exec_wait(sshsess, paste0("md5sum \"", x, "\"")))[1L])
    }, "")
    w1y <- w1[which(tst1 == tst2)]
    w1n <- w1[which(tst1 != tst2)]
    fastasTbl$Exists[w1n] <- FALSE
    fastasTbl$ExistsHere[w1y] <- fastasTbl$Exists[w1y] <- TRUE
    fastasTbl$Full[w1y] <- paste0(wd, "/", fastasTbl$Name[w1y])
    fastasTbl$Dir[w1y] <- wd
    fastasTbl$Loc[w1y] <- "Local"
  }
}
w2 <- which((!fastasTbl$Exists)&(fastasTbl$Loc == "Local")) # Local but not found...
if (length(w2)) {
  w2o <- w2[which(file.exists(paste0(wd, "/", fastasTbl$Name[w2])))] # ... and are currently in the data processing folder...
  if (length(w2o)) {
    tst <- (length(w2o) == 1)+1
    msg <- paste0(c("Several", "One")[tst], " Fasta file", c("s", "")[tst], " used to search the data w", c("ere", "as")[tst],
                  " not found at the original location, but ", c("are", "is")[tst], " present in the temporary data processing folder:\n",
                  paste(paste0(" - \"", fastasTbl$Full[w2o], "\""), collapse = "\n"))
    cat(msg)
    fastasTbl$Dir[w2o] <- wd
    fastasTbl$Full[w2o] <- paste0(wd, "/", fastasTbl$Name[w2o])
    fastasTbl$Exists[w2o] <- fastasTbl$ExistsHere[w2o] <- TRUE
    w2 <- which((!fastasTbl$Exists)&(fastasTbl$Loc == "Local"))
  }
}
if (length(w2)) {
  w2i <- w2[which(file.exists(paste0(indir, "/", fastasTbl$Name[w2])))] # ... or the input folder...
  if (length(w2i)) {
    tst <- (length(w2i) == 1L)+1L
    msg <- paste0(c("Several", "One")[tst], " Fasta file", c("s", "")[tst], " used to search the data w", c("ere", "as")[tst],
                  " not found at the original location, but ", c("are", "is")[tst], " present in the input folder:\n",
                  paste(paste0(" - \"", fastasTbl$Full[w2i], "\""), collapse = "\n"))
    warning(msg)
    fastasTbl$Dir[w2i] <- indir
    fastasTbl$Full[w2i] <- paste0(indir, "/", fastasTbl$Name[w2i])
    fastasTbl$Exists[w2i] <- fastasTbl$ExistsHere[w2i] <- TRUE
    w2 <- which((!fastasTbl$Exists)&(fastasTbl$Loc == "Local"))
  }
}
if (length(w2)) {
  # ... the ones which were not found should be located one by one by the user!!!
  dfltDr <- wd
  for (i in w2) { #i <- 1L
    tmp <- selectFile(paste0("Select missing fasta ", fastasTbl$Name[i]), path = wd)
    tmp <- gsub("^~", normalizePath(Sys.getenv("HOME"), winslash = "/"), tmp)
    fastasTbl$Dir[i] <- dfltDr <- dirname(tmp)
    fastasTbl$Full[i] <- tmp
    fastasTbl$Name[i] <- basename(tmp)
    fastasTbl$Exists[i] <- fastasTbl$ExistsHere[i] <- TRUE
    fastasTbl$Loc[i] <- "Local"
  }
}
w3 <- which((!fastasTbl$ExistsHere)&(fastasTbl$Exists)) # These exist in the current directory but not at the original location stipulated in their path
if (length(w3)) {
  fs::file_copy(fastasTbl$Full[w3], wd)
  fastasTbl$ExistsHere[w3] <- TRUE
  nwFast3 <- paste0(wd, "/", fastasTbl$Name[w3])
  w3m <- which(fastas_map$Actual %in% fastasTbl$Full[w3])
  m3 <- match(fastas_map$Actual[w3m], fastasTbl$Full[w3])
  fastas_map$Actual[w3m] <- nwFast3[m3]
  fastasTbl$Full[w3] <- nwFast3
  fastasTbl$Dir[w3] <- wd
  w3b <- which((!fastasTbl$TXTExistsHere[w3])&(fastasTbl$TXTExists[w3]))
  if (length(w3b)) {
    file_copy(fastasTbl$TXT[w3][w3b], wd)
    fastasTbl$TXTExistsHere[w3][w3b] <- TRUE
  }
}
w4 <- which((!fastasTbl$ExistsHere)&(fastasTbl$Loc == "Cluster")) # These do not exist here but are on the cluster
if (length(w4)) {
  require(tools)
  if (!SSH_on) {
    require(ssh)
    msg <- "SSH session required\nEnter host name!"
    sshost <- dlg_input(msg)$res
    sshsess <- ssh_connect(sshost)
    kount <- 1L
    while ((kount < 5L)&&(!inherits(sshsess, "ssh_session"))) {
      sshsess <- ssh_connect(sshost)
      kount <- kount + 1L
    }
    #print(sshsess)
    SSH_on <- inherits(sshsess, "ssh_session")
  }
  if (SSH_on) {
    tst <- try(scp_download(sshsess, fastasTbl$Full[w4], wd), silent = TRUE)
    w4y <- w4[which(file.exists(paste0(wd, "/", fastasTbl$Name[w4])))]
    w4n <- w4[which(!file.exists(paste0(wd, "/", fastasTbl$Name[w4])))]
    if (length(w4y)) {
      fastasTbl$ExistsHere[w4y] <- fastasTbl$Exists[w4y] <- TRUE
      fastasTbl$Full[w4y] <- paste0(wd, "/", fastasTbl$Name[w4y])
      fastasTbl$Dir[w4y] <- wd
      fastasTbl$Loc[w4y] <- "Local"
    }
    if (length(w4n)) {
      dfltDr <- wd
      for (i in w4n) { #i <- 1L
        tmp <- selectFile(paste0("Select missing fasta ", fastasTbl$Name[i]), path = wd)
        tmp <- gsub("^~", normalizePath(Sys.getenv("HOME"), winslash = "/"), tmp)
        fastasTbl$Dir[i] <- dfltDr <- dirname(tmp)
        fastasTbl$Full[i] <- tmp
        fastasTbl$Name[i] <- basename(tmp)
        fastasTbl$Exists[i] <- fastasTbl$ExistsHere[i] <- TRUE
        fastasTbl$Loc[i] <- "Local"
      }
    }
  }
}
#
fastasTbl <- fastasTbl[match(unique(fastasTbl$Full), fastasTbl$Full),]
fastasTbl$Full_list <- as.list(fastasTbl$Full)
# Sometimes, there may be several paths for a same fasta
tst <- aggregate(fastasTbl$Full, list(fastasTbl$Name), list)
tst$L <- lengths(tst$x)
if (max(tst$L) > 1L) {
  msg <- paste0("Possible duplicate fasta databases detected:", 
                paste(paste0("\n\n", unlist(tst$x[which(tst$L > 1L)])), collapse = ""),
                "\n\n\nAre they really duplicates?\n")
  simplFasta <- c(TRUE, FALSE)[match(dlg_message(msg, "yesno", rstudio = TRUE)$res, c("yes", "no"))]
  if (simplFasta) {
    fastasTbl2 <- tst[, c("Group.1", "x")]
    colnames(fastasTbl2) <- c("Name", "Full_list")
    fastasTbl2$Full <- fastasTbl$Full[match(fastasTbl2$Name, fastasTbl$Name)]
    fastasTbl2$Dir <- pathAbbr(fastasTbl2$Full)
    kol <- c("Loc", "Exists", "ExistsHere", "TXT", "TXTHere", "TXTExists", "TXTExistsHere", "Contaminant")
    tmp <- Isapply(fastasTbl2$Full, \(x) {
      w <- which(fastasTbl$Full == x)
      if (length(w) > 1L) {
        w <- which((fastasTbl$Full == x)&(fastasTbl$ExistsHere))
      }
      if (length(w) > 1L) { w <- w[1L] }
      return(fastasTbl[w, kol])
    })
    fastasTbl2[, kol] <- tmp
    for (k in c("Exists", "ExistsHere", "TXTExists", "TXTExistsHere", "Contaminant")) {
      fastasTbl2[[k]] <- as.logical(fastasTbl2[[k]])
    }
    fastasTbl <- fastasTbl2; rm(fastasTbl2)
  }
}
fasta_types %<o% c("UniprotKB", "ENSEMBL", "REFSEQPROTEIN", "NCBI", "TAIR")
opt <- vapply(fasta_types, \(x) { paste(c(x, rep(" ", 400L-nchar(x))), collapse = "") }, "")
cont_types <- setNames(c("No", "Yes"), c("", "+"))
opt2 <- vapply(cont_types, \(x) { paste(c(x, rep(" ", 400L-nchar(x))), collapse = "") }, "")
fastasTbl$Data <- fastasTbl$Headers <- NA
fastasTbl$Species <- rep("prompt user", nrow(fastasTbl))
fastasTbl$Type <- rep("unknown", nrow(fastasTbl))
fastasTbl$"Regexes test" <- ""
tstL <- (nrow(fastasTbl) > 1L)+1L
whFnd <- which(fastasTbl$Exists)
stopifnot(length(whFnd) > 0L)
fastasTbl$Data[whFnd] <- lapply(whFnd, \(x) {
  rs <- readLines(fastasTbl$Full[x])
  if(length(rs) == 0L) { stop(paste0("Fasta database \"", fastasTbl$Full[x], "\" appears to be empty!")) }
  return(rs)
})
fastasTbl$Headers[whFnd] <- lapply(fastasTbl$Data[whFnd], \(x) { grep("^>", x, value = TRUE) })
fastasTbl$Type[whFnd] <- vapply(fastasTbl$Headers[whFnd], \(x) { #x <- fastasTbl$Headers[whFnd[1L]]
  x <- unlist(x)[1L]
  # Note: the , in {2,} below is a concession to Bella who always, always, uses more than two characters here.
  res <- if (grepl("^>(rev_)?[a-z]{2,}\\|[^\\|]+\\|[^ ]+ ?", x)) { "UniprotKB" } else {
    if (grepl("^>(rev_)?AT1", x)) { "TAIR" } else { c("unknown", "NCBI")[grepl("^>(rev_)?[WYN]P_", x) + 1L] }
  }
  return(res)
}, "")
fastasTbl$Species[whFnd] <- vapply(fastasTbl$Headers[whFnd], \(hdrs) { #hdrs <- fastasTbl$Headers[whFnd[1L]]
  hdrs <- unlist(hdrs)
  orgpat <- c(" +OS=", "\\[[^\\[]+\\]")
  tst <- vapply(orgpat, \(x) { length(grep(x, hdrs)) }, 1L)
  tst <- which(tst == max(tst))[1L]
  if (tst == 1L) { pat <- "^.* +OS=| +[A-Z]{2}=.+$" }
  if (tst == 2L) { pat <- "^[^\\[]+\\[|\\].*$" }
  g <- grep(pat, hdrs)
  hdrs <- gsub(pat, "", hdrs[g])
  hdrs <- hdrs[which(nchar(hdrs) > 0L)]
  hdrs <- grep("^[A-Z][a-z]+( [a-z,A-Z]+( .*)?)", hdrs, value = TRUE)
  if (length(hdrs)) {
    hdrs <- data.table(Species = hdrs, Species2 = hdrs)
    hdrs <- hdrs[, .(Count = length(c(Species))), by = .(Species = Species2)]
    hdrs <- as.data.frame(hdrs)
    hdrs <- hdrs[order(hdrs$Count, decreasing = TRUE),]
    hdrs <- if (hdrs$Count[1L] > sum(hdrs$Count)*0.95) { hdrs$Species[1L] } else { "mixed" }
  } else { hdrs <- "prompt user" }
  return(hdrs)
}, "")
optSrc <- unique(c(fasta_types, "unknown"))
optSpc <- unique(c(fastasTbl$Species, "prompt user"))
if (!"Contaminants regex" %in% colnames(fastasTbl)) { fastasTbl$"Contaminants regex" <- "^CON__" }
if (!"Reverse regex" %in% colnames(fastasTbl)) { fastasTbl$"Reverse regex" <- "^rev_" }
# Below: commented code for automated calculation of number of reverse and contaminant hits for a given database
# Impractically slow for inclusion in the app to test regexes
# May be useful for later though...
# fastasTbl$"Regexes test"[whFnd] <- apply(fastasTbl[whFnd, c("Headers", "Type", "Contaminants regex", "Reverse regex")], 1, \(x) {
#   #x <- fastasTbl[whFnd[1L], c("Headers", "Type", "Contaminants regex", "Reverse regex")]
#   mode <- toupper(x[[2L]])
#   regex1 <- x[[3L]]
#   regex2 <- x[[4L]]
#   x <- unlist(x[[1L]])
#   if (!mode %in% toupper(optSrc)) {
#     warning(paste0("The value provided for \"mode\" (\"", mode, "\") could not be recognized, defaulting to \"custom\"!"))
#     mode <- "CUSTOM"
#   }
#   PAT2 <- "^>([^ ]*)"
#   if (mode %in% c("UNIPROT", "UNIPROTKB", "SWISSPROT", "TREMBL")) { PAT1 <- "^>[a-z]{2}\\|([^\\|]*)" }
#   # Ensembl-/RefSeq-RNA/NCBI specific behaviour
#   if (mode %in% c("ENSEMBL", "REFSEQRNA", "REFSEQPROTEIN", "NCBI")) { PAT1 <- PAT2 }
#   # RefSeq-Protein specific behaviour
#   if (mode == "REFSEQPROTEIN") { PAT1 <- "^>gi\\|[0-9]+\\|[a-z]+\\|([^ ]+(?![ \\|])[^ \\|])" }
#   # RefSeq-CDS specific behaviour
#   if (mode == "REFSEQCDS") { PAT1 <- "^>.*\\[protein=([^\\]]+)\\]" }
#   x1 <- lapply(x, \(y) {
#     z <- regexpr(PAT1, y, perl = TRUE)
#     return(substr(y,  attributes(z)$capture.start, attributes(z)$capture.start + attributes(z)$capture.length - 1L))
#   })
#   x2 <- lapply(x, \(y) {
#     z <- regexpr(PAT2, y, perl = TRUE)
#     return(substr(y,  attributes(z)$capture.start, attributes(z)$capture.start + attributes(z)$capture.length - 1L))
#   })
#   g1 <- grep(regex1, x1)
#   g2 <- grep(regex2, x2)
#   x <- paste0("IDs: ", length(x1), "\nCont: ", length(g1), " (", signif(round(100*length(g1)/length(x1)), 3L), "%)\nRev: ",
#               length(g2), " (",  signif(round(100*length(g2)/length(x1)), 3L), "%)\n")
#   cat(x)
#   return(x)
# })
# Annotation files
locFl <- paste0(homePath, "/Default_locations.xlsx")
dfltLocs <- openxlsx2::read_xlsx(locFl)
annotDir <- fastasDir <- dfltLocs$Path[match("Fasta files", dfltLocs$Folder)]
annotDflt <- paste0(fastasDir, "/.+\\.((txt)|(gtf)|(gff))$")
annotOpt <- setNames(c("txt", "gtf", "gff"),
                     c("UniProtKB .txt", "NCBI .gtf", "NCBI .gff"))
AnnotFls <- if ((!exists("AnnotFls"))||(!is.character(AnnotFls))) { c() } else {
  if (length(AnnotFls)) { AnnotFls[which(file.exists(AnnotFls))] }
}
if ((exists("AnnotFlsTbl"))&&(is.data.frame(AnnotFlsTbl))&&(nrow(AnnotFlsTbl))) {
  AnnotFls <- unique(c(AnnotFls, AnnotFlsTbl$Path))
}
AnnotFls <- AnnotFls[which(!is.na(AnnotFls))]
AnnotFls %<o% AnnotFls
nr0 <- length(AnnotFls)
updt_Type1 <- \(file) {
  vapply(file, \(fl) {
    if (!nchar(fl)) { return(NA) }
    x <- tolower(gsub(".*\\.", "", fl))
    if (!x %in% annotOpt) {
      #warning("Unrecognized file extension!")
      return(NA)
    }
    return(names(annotOpt)[match(x, annotOpt)])
  }, "")
}
annotTbl <- if (nr0) {
  data.frame(Select = "",
             Path = AnnotFls,
             Type = updt_Type1(AnnotFls),
             Remove = "")
} else {
  data.frame(Select = "",
             Path = "",
             Type = names(annotOpt)[1L],
             Remove = "")
}
rownames(annotTbl) <- NULL
nr0 <- nrow(annotTbl)
rng0 <- as.character(1L:nr0)
annotTbl2 <- annotTbl
annotTbl2$Path <- pathAbbr(annotTbl2$Path)
annotTbl2$Select <- vapply(paste0("selectAnnotFl___", rng0), \(id) {
  as.character(shiny::actionButton(id, "Select annotation file"))
}, "")
fSlct1 <- \(i, data = annotTbl2, opt = names(annotOpt)) {
  opt2 <- paste0("<option value=\"", opt, "\"",
                 c("", " selected")[(opt == data$Type[i])+1L],
                 ">", opt, "</option>", collapse = "")
  return(paste0("<select id=\"type___", as.character(i),
                "\" style=\"width:200px;\">", opt2, "</select>"))
}
annotTbl2$Type <- vapply(1L:nr0, fSlct1, "") # not rng0 here!
annotTbl2$Remove <- vapply(paste0("removeAnnotFl___", rng0), \(id) {
  as.character(shiny::actionButton(id, "Remove annotation file"))
}, "")
colnames(annotTbl2)[4L] <- ""
#
nr <- nrow(fastasTbl)
rws <- seq_len(nr)
appNm <- paste0(dtstNm, " - Fasta databases")
ui <- fluidPage(
  useShinyjs(),
  setBackgroundColor( # Doesn't work
    color = c(#"#F8F8FF",
      "#E1E2ED"),
    gradient = "linear",
    direction = "bottom"
  ),
  extendShinyjs(text = jsToggleFS, functions = c("toggleFullScreen")),
  tags$head(tags$style(HTML("table {table-layout: fixed;"))), # So table widths can be properly adjusted!
  titlePanel(tag("u", "Fasta databases"),
             appNm),
  h2(dtstNm),
  br(),
  h5("Indicate the source of fasta databases, which will define how the headers are parsed to extract protein accession and name, gene name, taxonomy etc... as available."),
  br(),
  h4("Regexes should be compatible with R's grep function:"),
  h5(" - Use the \"Contaminants regex\" column when a fasta contains a mixture of non-contaminants and contaminant proteins. Matches to the Protein IDs (accession) column will be marked as contaminants."),
  h5(" - Use the \"Reverse regex\" column when a fasta contains reverse proteins (decoys). Matches to the Full ID column will be removed from the database."),
  br(),
  h5(em("Whilst you have the option to choose from several types of common fasta, we only ever use UniProtKB, so for now use the other types at your own risk:")),
  h5(em("they are likely to break the script at a later stage!")),
  br(),
  h4(strong(em(tag("u", "Fasta databases")))),
  DTOutput("Fastas"),
  br(),
  br(),
  br(),
  h4(strong(em(tag("u", "Functional annotation files")))),
  DTOutput("annotFls"),
  shiny::actionButton("addBtn", "+ add annotation file"),
  br(),
  br(),
  br(),
  actionBttn("saveBtn", "Save", icon = icon("save"), color = "success", style = "pill"),
  br()
)
# Dummy table for display in app
fastasTbl2 <- data.frame("Name" = fastasTbl$Name,
                         "File found?" = c("-", "+")[fastasTbl$Exists+1L],
                         "Species" = fastasTbl$Species,
                         "Type" = fastasTbl$Type,
                         "Contaminants-only" = fastasTbl$Contaminant,
                         check.names = FALSE)
fastasTbl2$"Contaminants-only"[which(is.na(fastasTbl2$"Contaminants-only"))] <- FALSE # Doesn't hurts
fastasTbl2$Type <- shinySelectInput(fastasTbl$Type,
                                    "Type",
                                    optSrc,
                                    paste0(30L*max(c(nchar(optSrc), 2L)), "px"))
fastasTbl2$"Contaminants-only" <- shinyCheckInput(fastasTbl2$"Contaminants-only",
                                                  "ContOnly")
fastasTbl2$Species <- shinySelectInput(fastasTbl$Species,
                                       "Species",
                                       optSpc,
                                       paste0(30L*max(c(nchar(optSpc), 2L)), "px"))
fastasTbl2$"Contaminants regex" <- fastasTbl$"Contaminants regex"
fastasTbl2$"Reverse regex" <- fastasTbl$"Reverse regex"
NmWdth <- paste0(as.character(min(c(80L, max(nchar(fastasTbl2$Name))))*8L), "px")
wTest <- list(list(width = NmWdth, targets = 0L),
              list(width = "100px", targets = 1L:6L))
edith <- list(target = "column",
              enable = list(columns = grep(" regex$", colnames(fastasTbl2))-1))
tmp <- c(0L:(ncol(fastasTbl2)-1L))
tmp <- tmp[which(!tmp %in% edith$enable$columns)]
edith$disable <- list(columns = tmp)
if (exists("fastasTbl3")) { rm(fastasTbl3) }
if (exists("annotTbl3")) { rm(annotTbl3) }
slctXprs <- expression({
  dat <- ANNOTTBL()
  nAnnot <- nrow(dat)
  rg <- 1L:nAnnot
  rg0 <- rg[which(rg != i)]
  fls <- dat$Path[rg0]
  fl <- rstudioapi::selectFile(paste0("Select ", c("", "additional ")[(i>1L)+1L], " functional annotation file"),
                               path = ANNOTDIR())
  if ((length(fl) == 1L)&&(!is.na(fl))&&(file.exists(fl))) {
    if (fl %in% fls) {
      warning("You already selected this file! Ignoring...")
    } else {
      dat2 <- ANNOTTBL2()
      ANNOTDIR(dirname(fl))
      assign("annotDir", dirname(fl), envir = .GlobalEnv)
      dat$Path[i] <- dat2$Path[i] <- fl
      cat("")
      tp <- updt_Type1(fl)
      dat$Type[i] <- tp
      ANNOTTBL(dat)
      assign("annotTbl", dat, envir = .GlobalEnv)
      dat2$Type[i] <- fSlct1(i, dat)
      ANNOTTBL2(dat2)
      dat2$Path <- pathAbbr(dat2$Path)
      assign("annotTbl2", dat2, envir = .GlobalEnv)
      output$annotFls <- updt_AnnotFls()
    }
  }
})
#eval(parse(text = runApp), envir = .GlobalEnv)
typeXprs <- expression({
  dat <- ANNOTTBL()
  dat$Type[i] <- evnt$value
  ANNOTTBL(dat)
  if ("" %in% colnames(dat)) { stop() }
  assign("annotTbl", dat, envir = .GlobalEnv)
  dat2 <- ANNOTTBL2()
  assign("tmp", fSlct1(i, dat), envir = .GlobalEnv)
  dat2$Type[i] <- tmp
  ANNOTTBL2(dat2)
  dat2$Path <- pathAbbr(dat2$Path)
  assign("annotTbl2", dat2, envir = .GlobalEnv)
  output$annotFls <- updt_AnnotFls()
})
#eval(parse(text = runApp), envir = .GlobalEnv)
rmvXprs <- expression({
  dat <- ANNOTTBL()
  nAnnot <- nrow(dat)
  w <- which(1L:nAnnot != i)
  if (length(w)) { # Can't remove all!!!
    dat <- dat[w,]
    ANNOTTBL(dat)
    if ("" %in% colnames(dat)) { stop() }
    assign("annotTbl", dat, envir = .GlobalEnv)
    dat2 <- ANNOTTBL2()[w,]
    nAnnot2 <- nrow(dat)
    chRg2 <- as.character(1L:nAnnot2)
    # Re-generate table with IDs from updated row position
    dat2$Select <- vapply(paste0("selectAnnotFl___", chRg2), \(id) {
      as.character(shiny::actionButton(id, "Select annotation file"))
    }, "")
    dat2$Type <- vapply(1L:nAnnot2, fSlct1, "", dat)
    dat2[[4L]] <- vapply(paste0("removeAnnotFl___", chRg2), \(id) {
      as.character(shiny::actionButton(id, "Remove annotation file"))
    }, "")
    ANNOTTBL2(dat2)
    dat2$Path <- pathAbbr(dat2$Path)
    assign("annotTbl2", dat2, envir = .GlobalEnv)
    output$annotFls <- updt_AnnotFls()
  }
})
server <- \(input, output, session) {
  ANNOTTBL <- shiny::reactiveVal(annotTbl)
  ANNOTTBL2 <- shiny::reactiveVal(annotTbl2)
  ANNOTDIR <- shiny::reactiveVal(annotDir)
  # Outputs
  fastasTbl3 <- fastasTbl
  output$Fastas <- renderDT({ fastasTbl2 },
                            FALSE,
                            escape = FALSE,
                            class = "compact",
                            selection = "none",
                            editable = edith,
                            rownames = FALSE,
                            options = list(dom = 't',
                                           paging = FALSE,
                                           ordering = FALSE,
                                           autowidth = TRUE,
                                           columnDefs = wTest,
                                           scrollX = TRUE),
                            # the callback is essential to capture the inputs in each row
                            callback = JS("table.rows().every(function(i, tab, row) {
        var $this = $(this.node());
        $this.attr('id', this.data()[0]);
        $this.addClass('shiny-input-container');
      });
      Shiny.unbindAll(table.table().node());
      Shiny.bindAll(table.table().node());"))
  # Reactive functions
  updt_AnnotFls <- \(reactive = TRUE) {
    dat2 <- if (reactive) { ANNOTTBL2() } else { annotTbl2 }
    return(DT::renderDT({ dat2 },
                        FALSE,
                        escape = FALSE,
                        class = "compact",
                        selection = "none",
                        rownames = FALSE,
                        editable = FALSE,
                        options = list(dom = "t",
                                       paging = FALSE,
                                       ordering = FALSE,
                                       autowidth = TRUE,
                                       columnDefs = list(list(width = "150px", targets = 0L),
                                                         list(width = "300px", targets = 1L),
                                                         list(width = "200px", targets = 2L),
                                                         list(width = "100px", targets = 3L)),
                                       scrollX = FALSE),
                        # the callback is essential to capture the inputs in each row
                        callback = DT::JS("
// Buttons
table.on('click', 'button', function() {
  var id = this.id;
  Shiny.setInputValue('dt_event', {type: 'button', id: id}, {priority: 'event'});
});
// Dropdowns
table.on('change', 'select', function() {
  var id = this.id;
  var val = this.value;
  Shiny.setInputValue('dt_event', {type: 'select', id: id, value: val}, {priority: 'event'});
});
")))
  }
  output$annotFls <- updt_AnnotFls(FALSE)
  # Observers
  observeEvent(input$Fastas_cell_edit, {
    kol <- colnames(fastasTbl2)[input$Fastas_cell_edit$col+1L]
    if ((length(kol))&&(kol %in% colnames(fastasTbl3))) {
      fastasTbl3[input$Fastas_cell_edit$row, kol] <<- input$Fastas_cell_edit$value
    }
  }, ignoreNULL = FALSE)
  observeEvent(input$dt_event, {
    evnt <- input$dt_event
    id <- evnt$id
    i <- as.integer(gsub(".*___", "", id))
    root <- gsub("___.*", "", id)
    if (evnt$type == "button") {
      # Directory selection
      if (root == "selectAnnotFl") { eval(slctXprs) }
      # Row removal?
      if (root == "removeAnnotFl") { eval(rmvXprs) }
    }
    if (evnt$type == "select") {
      # Software selection
      if(root == "searchSoft") { eval(typeXprs) }
    }
  })
  observeEvent(input$addBtn, {
    dat <- ANNOTTBL()
    dat2 <- ANNOTTBL2()
    fls <- dat$Path
    nri <- length(fls)
    j <- nri+1L
    jChr <- as.character(j)
    datDflt <- data.frame(Select = "",
                          Path = dirname(dat$Path[nri]),
                          Type = dat$Type[nri],
                          Remove = "")
    rownames(datDflt) <- NULL
    dat <- rbind(dat, datDflt)
    ANNOTTBL(dat)
    if ("" %in% colnames(dat)) { stop() }
    assign("annotTbl", dat, envir = .GlobalEnv)
    datDflt$Path <- ""
    datDflt2 <- datDflt
    datDflt2$Select <- as.character(shiny::actionButton(paste0("selectAnnotFl___", jChr),
                                                        "Select annotation file"))
    datDflt2$Type <- fSlct1(i, datDflt)
    datDflt2$Remove <- as.character(shiny::actionButton(paste0("removeAnnotFl___", jChr),
                                                        "Remove annotation file"))
    colnames(datDflt2)[4L] <- ""
    dat2 <- rbind(dat2, datDflt2)
    ANNOTTBL2(dat2)
    dat2$Path <- pathAbbr(dat2$Path)
    assign("annotTbl2", dat2, envir = .GlobalEnv)
    output$annotFls <- updt_AnnotFls()
  })
  observeEvent(input$saveBtn, {
    fastasTbl3$Contaminant <- vapply(rws, \(x) { input[[paste0("ContOnly___", as.character(x))]] }, TRUE)
    fastasTbl3$Type <- vapply(rws, \(x) { input[[paste0("Type___", as.character(x))]] }, "a")
    fastasTbl3$Species <- vapply(rws, \(x) { input[[paste0("Species___", as.character(x))]] }, "a")
    assign("fastasTbl3", fastasTbl3, envir = .GlobalEnv)
    assign("annotTbl3", ANNOTTBL(), envir = .GlobalEnv)
    stopApp()
  })
  #observeEvent(input$cancel, { stopApp() })
  session$onSessionEnded(\() { stopApp() })
}
runKount <- 0L
while ((!runKount)||(!exists("fastasTbl3"))) {
  eval(parse(text = runApp), envir = .GlobalEnv)
  shinyCleanup()
  runKount <- runKount+1L
}
annotTbl3 <- annotTbl3[which(file.exists(annotTbl3$Path)),]
isContaminant <- setNames(fastasTbl3$Contaminant, fastasTbl3$Full)
Annotate %<o% (nrow(annotTbl3) > 0L)
if (Annotate) {
  AnnotFlsTbl %<o% annotTbl3
  AnnotFls %<o% annotTbl3$Path
}
fastasTbl %<o% fastasTbl3
w <- which(fastasTbl$Species == "prompt user")
if (length(w)) {
  for (i in w) {
    fastasTbl$Species[i] <- dlg_input(paste0("What is the main species in Fasta \"",
                                             fastasTbl$Full[i],
                                             "\" (ignore contaminants; write \"mixed\" if not single species)?"))$res
  }
}
Sp <- gsub(" *[\\(|\\[].*", "", fastasTbl$Species)
pack <- "myTAI"
if (!require(pack, character.only = TRUE)) {
  try(pak::pkg_install("drostlab/myTAI", ask = FALSE, upgrade = TRUE, dependencies = TRUE), silent = TRUE) 
}
if (!require(pack, character.only = TRUE)) {
  Src2 <- paste0(libPath, "/extdata/Sources/taxonomy.R")
  source(Src2, local = FALSE)
}
tst %<o% try(setNames(lapply(Sp, \(x) {
  suppressMessages(taxonomy(organism = x, db = "ncbi", output = "classification"))
}), Sp), silent = TRUE)
if (inherits(tst, "try-error")) {
  kount <- 0L
  while ((kount < 20L)&&(inherits(tst, "try-error"))) {
    tst %<o% try(setNames(lapply(Sp, \(x) {
      suppressMessages(taxonomy(organism = x, db = "ncbi", output = "classification"))
    }), Sp), silent = TRUE)
    kount <- kount + 1L
  } 
}
taxTst %<o% !inherits(tst, "try-error")
if (taxTst) {
  Taxonomies <- tst
  w <- which(!is.na(fastasTbl$Species))
  fastasTbl$Kingdom <- fastasTbl$Taxonomy <- NA
  fastasTbl$Taxonomy[w] <- vapply(Taxonomies, \(x) { #x <- Taxonomies[[1L]]
    paste(x$name, collapse = "; ")
  }, "")
  fastasTbl$Kingdom[w] <- vapply(Taxonomies, \(x) { #x <- Taxonomies[[1L]]
    x$name[match("superkingdom", x$rank)]
  }, "")
}
source(parSrc, local = FALSE)
dbs <- setNames(lapply(whFnd, \(i) { #i <- whFnd[1L] #i <- whFnd[2L]
  tmp <- Format.DB(unlist(fastasTbl$Data[[i]]), #file <- unlist(fastasTbl$Data[[i]])
                   in.env = TRUE,
                   mode = fastasTbl$Type[i],
                   parallel = TRUE,
                   cl = parClust)
  tmp$Source <- fastasTbl$Type[i]
  tmp$"Potential contaminant" <- c("", "+")[fastasTbl$Contaminant[i]+1L]
  if (!fastasTbl$Contaminant[i]) {
    g <- grep(fastasTbl$"Contaminants regex"[i], tmp$"Protein ID")
    tmp$"Potential contaminant"[g] <- "+"
    tmp$"Protein ID"[g] <- gsub(fastasTbl$"Contaminants regex"[i], "", tmp$"Protein ID"[g])
  }
  if (nchar(fastasTbl$`Reverse regex`[i])) {
    g <- try(grep(fastasTbl$"Reverse regex"[i], tmp$"Full ID", invert = TRUE), silent = TRUE)
    if (!inherits(g, "try-error")) { tmp <- tmp[g,] }
  }
  tmp$"Protein of interest" <- (fastasTbl$Name[i] == "Proteins of interest.fasta")
  if ((!is.na(fastasTbl$Species[i]))&&(!fastasTbl$Species[i] %in% c("mixed", "prompt user", "unknown"))) {
    tmp$Organism_Full <- fastasTbl$Species[i]
  }
  if (taxTst) {
    tmp$Taxonomy <- fastasTbl$Taxonomy[i]
    tmp$Kingdom <- fastasTbl$Kingdom[i]
  }
  return(tmp)
}), fastasTbl$Full[whFnd])
#
kol <- c("Organism_Full", "Organism")
dbs_Org %<o% setNames(lapply(dbs, \(db) { #db <- dbs[[1L]] #db <- dbs[[2L]] #db <- dbs[[3L]] #db <- dbs[[4L]]
  kol <- kol[which(kol %in% colnames(db))]
  tst <- vapply(kol, \(x) { length(unique(db[which(!as.character(db[[x]]) %in% c("", "NA")), x])) }, 1L)
  kol <- kol[order(tst, decreasing = TRUE)][1L]
  w <- which(db$`Potential contaminant` != "+")
  if (!length(w)) { return() }
  org %<o% aggregate(w, list(db[w, kol]), length)
  colnames(org) <- c("Organism", "Count")
  org <- org[which(org$Count == max(org$Count)[1L]),]
  org$Source <- aggregate(db$Source[which(db[[kol]] %in% org$Organism)], list(db[which(db[[kol]] %in% org$Organism), kol]), \(x) {
    unique(x[which(!is.na(x))])
  })$x
  org$Source[which(is.na(org$Source))] <- ""
  return(org)
}), fastasTbl$Full[whFnd])
tmp <- listMelt(fastasTbl$Full_list, fastasTbl$Full, c("Original", "Final"))
dbs_Txt <- setNames(vapply(fastas_map$Actual, \(x) { #x <- fastas_map$Actual[[1L]]
  m <- unique(c(x, tmp$Final[match(basename(x), basename(tmp$Original))]))
  tbl <- do.call(rbind, c(dbs_Org[m]))
  # tbl <- data.frame(Organism = c("Arabidopsis thaliana", "Arabidopsis thaliana", "Brassica oleracea", "Escherichia coli", "Arabidopsis thaliana", "Homo sapiens"),
  #                   Count = c(10000L, 2000L, 5000L, 300L, 1500L, 25L),
  #                   Source = c("UniprotKB", "UniprotKB", "UniprotKB", "UniprotKB", "TAIR", "UniprotKB"))
  tbl <- aggregate(tbl$Count, list(tbl$Source, tbl$Organism), sum)
  colnames(tbl) <- c("Source", "Organism", "Count")
  tbl <- tbl[order(tbl$Organism),]
  tbl <- tbl[order(tbl$Count, decreasing = TRUE),]
  tbl <- tbl[order(tbl$Source),]
  tbl1 <- aggregate(tbl[, c("Organism", "Count")], list(tbl$Source), list)
  colnames(tbl1)[1L] <- "Source"
  tbl1$charTst <- lapply(tbl1$Organism, \(y) {
    (tolower(substr(y, 1L, 1L)) %in% c("a", "e", "i", "o", "u")) + 1L })
  Txt <- vapply(1L:nrow(tbl1), \(y) { #y <- 1L #y <- 2L
    org <- tbl1$Organism[[y]]
    knt <- tbl1$Count[[y]]
    tst <- tbl1$charTst[[y]]
    l <- length(org)
    if (l > 1L) {
      org <- paste0(paste(org[1L:(l-1L)], collapse = ", "), " and ", org[l], " databases")
      knt <- paste0(paste(knt[1L:(l-1L)], collapse = ", "), " and ", knt[l], " entries, resp.")
    } else {
      org <- paste0("a", c("", "n")[tst], " ", org, " database")
      knt <- paste0(as.character(knt), " entr", c("y", "ies")[(knt > 1L)+1L])
    }
    return(paste0(org, " obtained from ", tbl1$Source[y], " (", knt, ")"))
  }, "")
  l <- length(Txt)
  if (l > 1L) { Txt <- paste0(paste(Txt[1L:(l-1L)], collapse = ", "), " as well as ", Txt[l]) }
  return(paste0(Txt, "."))
}, ""), names(fastas_map$Actual))
# Update MatMet_Search
for (dir in names(MatMet_Search)) { #dir <- names(MatMet_Search)[1L]
  MatMet_Search[[dir]] <- gsub("TEMPLATELIBTEXT", dbs_Txt[dir], MatMet_Search[[dir]])
}
MatMet_Search <- paste(MatMet_Search, collapse = "\n")
#cat(MatMet_Search)
#
Org %<o% do.call(rbind, dbs_Org)
Org <- aggregate(Org$Count, list(Org$Organism), sum)
colnames(Org) <- c("Organism", "Count")
Org <- Org[order(Org$Count, decreasing = TRUE),]
#tstOrgNm %<o% c("", "n")[(tolower(substr(Org$Organism[1L], 1L, 1L)) %in% c("a", "e", "i", "o", "u"))+1]
#
db %<o% plyr::rbind.fill(dbs)
#
# In case we are assembling several databases into one, we may have redundant sequences.
# Remove them, keeping in priority those matching our accessions of interest.
wInt <- which(db$"Protein of interest")
wRst <- which(!db$"Protein of interest")
db <- db[c(wInt, wRst),]
tst <- data.table::data.table(Row = 1L:nrow(db), Seq = db$Sequence)
tst <- tst[, list(Row = min(Row)), by = list(Seq)]
n <- nrow(db)-nrow(tst)
if (n) {
  cat("(Removing", n, "redundant protein sequences from the", c("", "combined")[(length(dbs) > 1L)+1L], "database.)\n")
  db <- db[tst$Row,]
}
#
w2 <- which((!db$`Protein ID` %in% db$`Protein ID`[w])|(db$`Protein of interest`))
db <- db[w2,]
db <- db[order(db$"Protein of interest", decreasing = TRUE),]
for (i in 1L:nrow(fastasTbl)) { if (!file.exists(paste0(wd, "/", fastasTbl$Name[i]))) {
  fs::file_copy(fastasTbl$Full[i], wd)
} }
if (scrptType == "noReps") { AnalysisParam$fastasTbl <- list(fastasTbl$Full) }

# Filter for repeats, mark as contaminants
tmp <- as.data.table(db[, c("Potential contaminant", "Protein ID")])
tst <- tmp[, list(x = c(`Potential contaminant`)), by = list(Group.1 = `Protein ID`)]
tst <- as.data.frame(tst)
u <- tst$Group.1
tst$L <- lengths(tst$x)
w <- which(tst$L > 1L)
if (length(w)) {
  tst <- tst[w,]
  tst$Cont <- vapply(tst$x, \(x) { c("", "+")[("+" %in% x)+1L] }, "")
  w <- which(db$`Protein ID` %in% tst$Group.1)
  db$`Potential contaminant`[w] <- tst$Cont[match(db$`Protein ID`[w], tst$Group.1)]
  db <- db[match(u, db$`Protein ID`),]
}
#
for (i in 1L:nrow(fastasTbl)) { if (!file.exists(paste0(wd, "/", fastasTbl$Name[i]))) {
  fs::file_copy(fastasTbl$Full[i], wd)
} }
# Remove reverse entries
db <- db[grep("^>rev_", db$Header, invert = TRUE),]

# Also load and append contaminants database
setwd(wd)
if (sum(c("MAXQUANT", "DIANN", "FRAGPIPE") %in% SearchSoft)) {
  # NB:
  # DIA-NN does not provide an inbuilt contaminants database.
  # But we are providing one, slightly modified from the CCP's cRAPome.fasta, which should normally be used for DiaNN searches.
  # FragPipe can add the CCP's cRAPome.fasta to the search and we recommend to do so.
  # For MaxQuant, use the contaminants.fasta which is also copied with the package.
  contDBFls <- paste0(libPath, "/extData/", c("CCP_cRAPome.fasta",
                                              "contaminants.fasta")[c(sum(c("DIANN", "FRAGPIPE") %in% SearchSoft),
                                                                      "MAXQUANT" %in% SearchSoft)])
  contDBFls <- contDBFls[which(file.exists(contDBFls))]
  if (length(contDBFls)) {
    contDB <- lapply(contDBFls, \(contDBFl) {
      x <- readLines(contDBFl)
      x <- Format.DB(x, in.env = TRUE, cl = parClust)
      return(x)
    })
    contDB <- plyr::rbind.fill(contDB)
  }
  if ("MAXQUANT" %in% SearchSoft) {
    contCsv <- paste0(wd, "/contaminants.csv")
    if (file.exists(contCsv)) { contDB <- read.csv(contCsv, check.names = FALSE) } else {
      # This is a bastard Fasta...
      dr <- paste0(fastasDir, "/Contaminants")
      if (!dir.exists(dr)) { dir.create(dr, recursive = TRUE) }
      clusterExport(parClust, "dr", envir = environment())
      tst <- parLapply(parClust, contDB$`Full ID`, \(x) { #x <- contDB$`Full ID`[1L]
        if (!grepl(":", x)) {
          tmp <- paste0("https://rest.uniprot.org/uniprotkb/", x, ".fasta")
          dest <- paste0(dr, "/", x, ".fasta")
          if ((!file.exists(dest))||(base::file.size(dest) == 0)) {
            try(utils::download.file(tmp, dest), silent = TRUE)
            
          }
        }
      })
      w <- which(file.exists(paste0(dr, "/", contDB$`Full ID`, ".fasta")))
      tst <- parSapply(parClust, paste0(dr, "/", contDB$`Full ID`[w], ".fasta"), \(x) { #x <- paste0(dr, "/", contDB$`Full ID`[w], ".fasta")
        tst <- try(proteoCraft::Format.DB(x), silent = TRUE)
        tst <- if (inherits(tst, "try-error")) { list(Outcome = FALSE) } else {
          list(Outcome = TRUE,
               tbl = tst)
        }
        return(tst)
      })
      tst <- tst[which(vapply(tst, \(x) { x$Outcome }, TRUE))]
      tst <- lapply(tst, \(x) { x$tbl })
      tst <- plyr::rbind.fill(tst)
      contDB <- plyr::rbind.fill(tst,
                                 contDB[which(!contDB$`Full ID` %in% tst$`Protein ID`),])
      w1 <- which(grepl("^>(sp)|(tr)\\|", contDB$Header))
      org <- unique(contDB$Organism_Full[w1])
      library(UniProt.ws)
      txIDs <- setNames(lapply(org, \(x) { c() }), org)
      w <- which(org %in% names(Taxonomies))
      if (length(w)) { txIDs[w] <- Taxonomies[w] }
      w <- which(!org %in% names(Taxonomies))
      if (length(w)) {
        txIDs[w] <- setNames(parLapply(parClust, org[w], \(x) {
          tst <- try(myTAI::taxonomy(organism = x, db = "ncbi", output = "classification"), silent = TRUE)
          rs <- if (inherits(rs, "try-error")) { list(Outcome = FALSE) } else {
            list(Outcome = TRUE,
                 Res = tst) 
          }
        }), org[w])
        wY <- w[which(vapply(txIDs[w], \(x) { x$Outcome }, TRUE))]
        wN <- w[which(!vapply(txIDs[w], \(x) { x$Outcome }, TRUE))]
        txIDs[wN] <- c()
        txIDs[wY] <- lapply(txIDs[wY], \(x) { x$Res })
        Taxonomies[org[wY]] <- txIDs[org[wY]]
      }
      txIDs <- Taxonomies[which(vapply(Taxonomies, \(x) { is.data.frame(x) }, TRUE))]
      txIDs <- vapply(txIDs, \(x) { as.character(x$id[which(x$rank == "species")]) }, "")
      kount <- 0L
      w2 <- grep("^>(sp)|(tr)\\|", contDB$Header, invert = TRUE)
      while ((length(w2))&&(kount < length(txIDs))) {
        kount <- kount + 1L
        sp <- gsub(" ", "+", txIDs[kount])
        #tmp <- paste0("https://www.uniprot.org/uniprot/?query=organism_name:", sp, "&format=fasta") # Old 
        tmp <- paste0("https://rest.uniprot.org/uniprotkb/search?query=organism_id:", sp, "&format=fasta")
        dest <- paste0("D:/Fasta_databases/", gsub(" ", "_", sp), "_-_uniprot-all-accessions_", gsub("-", "", Sys.Date()), ".fasta")
        tmp2 <- if (!file.exists(dest)) { try(download.file(tmp, dest), silent = TRUE) } else { 0L }
        if ((!inherits(tmp2, "try-error"))&&(!tmp2)) {
          tmp2 <- Format.DB(dest)
          wY <- which(tmp2$Sequence %in% contDB$Sequence[w2])
          if (length(wY)) {
            contDB <- plyr::rbind.fill(contDB[which(!contDB$Sequence %in% tmp2$Sequence),],
                                       tmp2[wY,])
            w2 <- grep("^>(sp)|(tr)\\|", contDB$Header, invert = TRUE)
          }
        }
      }
      contDB <- contDB[grep("^>(sp)|(tr)\\|", contDB$Header),] # We do not want anything unusable
      contDB$Organism_Full <- "Contaminant"
      contDB$Organism <- "Contaminant"
      contDB$"Protein ID" <- paste0("CON__", gsub("^CON__", "", contDB$"Protein ID"))
      contDB$"Potential contaminant" <- "+"
      write.csv(contDB, contCsv, row.names = FALSE)
    }
  }
  #db <- db[which(db$"Potential contaminant" != "+"),]
  w <- which(db$Sequence %in% contDB$Sequence)
  db$"Potential contaminant"[w] <- "+"
  for (kol in c("Organism", "Organism_Full")) { if (!kol %in% colnames(db)) { db[[kol]] <- NA } }
  db[w, c("Organism", "Organism_Full")] <- "Contaminant"
  db$"Protein ID"[w] <- paste0("CON__", gsub("^CON__", "", db$"Protein ID"[w]))
  contDB %<o% contDB[which(!contDB$Sequence %in% db$Sequence),]
  if (nrow(contDB)) {
    contDB$"Potential contaminant" <- "+"
    db <- rbind.fill(db, contDB)
  }
} else { if (SearchSoft == "PROTEOMEDISCOVERER") { stop("This part has not yet been re-written for PD!") } }

# Gene column!
# This column can be missing, yet it is critical for some functions (e.g. GO enrichment)
if (!"Gene" %in% colnames(db)) { db$Gene <- "" }
wY <- which((is.na(db$Gene))|(db$Gene == ""))  # "Yes", as in "yes it's broken, fix it!"
wN <- which((!is.na(db$Gene))&(db$Gene != "")) # "No", as in "nope, everything's fine here, carry on"
lY <- length(wY)
if (lY) {
  tmp <- as.character(1L:lY)
  tst <- max(c(max(nchar(tmp)), 6L))
  gns <- paste0("Gn", vapply(tst-nchar(tmp[1L:lY]), \(x) {
    paste(rep("0", x), collapse = "")
  }, ""), tmp)
  wY2 <- which(gns %in% db$Gene[wN])
  wN2 <- which(!gns %in% db$Gene[wN])
  lY2 <- length(wY2)
  while (lY2) {
    lY2 <- length(wY2)
    tmp2 <- as.character(lY+(1L:lY2))
    tst2 <- max(c(max(nchar(tmp2)), 6L))
    gns2 <- paste0("Gn", vapply(tst2-nchar(tmp2[1L:lY2]), \(x) {
      paste(rep("0", x), collapse = "")
    }, ""), tmp)
    lY <- lY+lY2
    gns <- c(gns[wN2], gns2)
    wY2 <- which(gns %in% db$Gene[wN])
    wN2 <- which(!gns %in% db$Gene[wN])
  }
  db$Gene[wY] <- gns
}

#
w <- which((is.na(db$Organism))|(db$Organism == ""))
if (length(w)) {
  w1 <- w[which(db$`Potential contaminant`[w] == "+")]
  w2 <- w[which(db$`Potential contaminant`[w] != "+")]
  db$Organism[w1] <- "Contaminant"
  db$Organism[w2] <- "Unknown"
}
w <- which((is.na(db$Organism_Full))|(db$Organism_Full == ""))
if (length(w)) { db$Organism_Full[w] <- db$Organism[w] }

# Organism columns
w <- which(c("Organism_Full", "Organism") %in% colnames(db))
tstOrg <- (length(w) > 0L)
if (tstOrg) {
  dbOrgKol <- c("Organism_Full", "Organism")[w[1L]]
  tst <- gsub(" *(\\(|\\[).*", "", db[[dbOrgKol]])
  tst <- aggregate(tst, list(tst), length)
  tst <- tst[order(tst$x, decreasing = TRUE),]
  mainOrg %<o% tst$Group.1[1L]
} else {
  dbOrgKol <- c()
  mainOrg <- "[UNKNOWN_ORGANISM]"
}
dbOrgKol %<o% dbOrgKol
tstOrg %<o% tstOrg
mainOrg %<o% mainOrg

# We will not need the super large "Headers" and "Data" columns anymore,
# so let's remove them to make the object smaller.
fastasTbl$Data <- NULL
fastasTbl$Headers <- NULL
