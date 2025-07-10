#### Code chunk - Load and process search database(s) (or reload a pre-processed version if available)
SSH_on %<o% FALSE
fastas <- unique(fastas)
#fastas <- fastas$Full
#
fastasTbl %<o% data.frame(Full = fastas,
                          Name = basename(fastas),
                          Dir = dirname(fastas),
                          Contaminant = FALSE)
fastasTbl$Loc <- c("Local", "Cluster")[sapply(fastasTbl$Dir, function(x) { grepl("^/+nfs/", x) }) + 1 ]
fastasTbl$Exists <- file.exists(fastasTbl$Full)
fastasTbl$ExistsHere <- file.exists(fastasTbl$Name)
fastasTbl$TXT <- gsub("\\.fa(s(ta(\\.fas)?)?)?$", ".txt", fastasTbl$Full)
fastasTbl$TXTHere <- gsub("\\.fa(s(ta(\\.fas)?)?)?$", ".txt", fastasTbl$Name)
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
    kount <- 1
    while ((kount < 5)&&(!"ssh_session" %in% class(sshsess))) {
      sshsess <- ssh_connect(sshost)
      kount <- kount + 1
    }
    #print(sshsess)
    SSH_on <- "ssh_session" %in% class(sshsess)
  }
  if (SSH_on) {
    tst1 <- md5sum(fastasTbl$Name[w1])
    tst2 <- sapply(fastasTbl$Full[w1], function(x) { #x <- fastasTbl$Full[w1[1]]
      gsub(" .*", "", capture.output(ssh_exec_wait(sshsess, paste0("md5sum \"", x, "\"")))[1])
    })
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
    tst <- (length(w2i) == 1)+1
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
  for (i in w2) { #i <- 1
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
  fastasTbl$Full[w3] <- paste0(wd, "/", fastasTbl$Name[w3])
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
    kount <- 1
    while ((kount < 5)&&(!"ssh_session" %in% class(sshsess))) {
      sshsess <- ssh_connect(sshost)
      kount <- kount + 1
    }
    #print(sshsess)
    SSH_on <- "ssh_session" %in% class(sshsess)
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
      for (i in w4n) { #i <- 1
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
# Sometimes, there may be several paths for a same fasta
tst <- aggregate(fastasTbl$Full, list(fastasTbl$Name), list)
tst$L <- sapply(tst$x, length)
if (max(tst$L) > 1) {
  msg <- paste0("Possible duplicate fasta databases detected:", 
                paste(paste0("\n\n", unlist(tst$x[which(tst$L > 1)])), collapse = ""),
                "\n\n\nAre they really duplicates?\n")
  simplFasta <- c(TRUE, FALSE)[match(dlg_message(msg, "yesno", rstudio = FALSE)$res, c("yes", "no"))]
  if (simplFasta) {
    fastasTbl2 <- data.frame(Full = unique(fastasTbl$Full))
    fastasTbl2$Name <- gsub(".*/", "", fastasTbl2$Full)
    fastasTbl2$Dir <- gsub("/[^/]+$", "", fastasTbl2$Full)
    kol <- c("Loc", "Exists", "ExistsHere", "TXT", "TXTHere", "TXTExists", "TXTExistsHere")
    tmp <- Isapply(fastasTbl2$Full, function(x) {
      w <- which(fastasTbl$Full == x)
      if (length(w) > 1) {
        w <- which((fastasTbl$Full == x)&(fastasTbl$ExistsHere))
      }
      if (length(w) > 1) { w <- w[1] }
      return(fastasTbl[w, kol])
    })
    fastasTbl2[, kol] <- tmp
    for (k in c("Exists", "ExistsHere", "TXTExists", "TXTExistsHere")) {
      fastasTbl2[[k]] <- as.logical(fastasTbl2[[k]])
    }
    fastasTbl <- fastasTbl2; rm(fastasTbl2)
  }
}
fasta_types %<o% c("UniprotKB", "ENSEMBL", "REFSEQPROTEIN", "NCBI", "TAIR")
opt <- sapply(fasta_types, function(x) { paste(c(x, rep(" ", 400-nchar(x))), collapse = "") })
cont_types <- setNames(c("No", "Yes"), c("", "+"))
opt2 <- sapply(cont_types, function(x) { paste(c(x, rep(" ", 400-nchar(x))), collapse = "") })
fastasTbl$Data <- fastasTbl$Headers <- NA
fastasTbl$Species <- rep("prompt user", nrow(fastasTbl))
fastasTbl$Type <- rep("unknown", nrow(fastasTbl))
fastasTbl$"Regexes test" <- ""
tstL <- (nrow(fastasTbl) > 1)+1
whFnd <- which(fastasTbl$Exists)
stopifnot(length(whFnd) > 0)
fastasTbl$Data[whFnd] <- lapply(whFnd, function(x) {
  rs <- readLines(fastasTbl$Full[x])
  if(length(rs) == 0) { stop(paste0("Fasta database \"", fastasTbl$Full[x], "\" appears to be empty!")) }
  return(rs)
})
fastasTbl$Headers[whFnd] <- lapply(fastasTbl$Data[whFnd], function(x) { grep("^>", x, value = TRUE) })
fastasTbl$Type[whFnd] <- sapply(fastasTbl$Headers[whFnd], function(x) { #x <- fastasTbl$Headers[whFnd[1]]
  x <- unlist(x)[1]
  if (grepl("^>(rev_)?[a-z]{2}\\|[^\\|]+\\|[^ ]+ ", x)) { res <- "UniprotKB" } else {
    if (grepl("^>(rev_)?ATl", x)) { res <- "TAIR" } else {
      if (grepl("^>(rev_)?[WYN]P_", x)) { res <- "NCBI" } else { res <- "unknown" }
    }
  }
  return(res)
})
fastasTbl$Species[whFnd] <- sapply(fastasTbl$Headers[whFnd], function(hdrs) { #hdrs <- fastasTbl$Headers[whFnd[1]]
  hdrs <- unlist(hdrs)
  orgpat <- c(" +OS=", "\\[[^\\[]+\\]")
  tst <- sapply(orgpat, function(x) { length(grep(x, hdrs)) })
  tst <- which(tst == max(tst))[1]
  if (tst == 1) { pat <- "^.* +OS=| +[A-Z]{2}=.+$" }
  if (tst == 2) { pat <- "^[^\\[]+\\[|\\].*$" }
  g <- grep(pat, hdrs)
  hdrs <- gsub(pat, "", hdrs[g])
  hdrs <- hdrs[which(nchar(hdrs) > 0)]
  hdrs <- grep("^[A-Z][a-z]+( [a-z,A-Z]+( .*)?)", hdrs, value = TRUE)
  if (length(hdrs)) {
    hdrs <- aggregate(hdrs, list(hdrs), length)
    hdrs <- hdrs[order(hdrs$x, decreasing = TRUE),]
    hdrs <- hdrs$Group.1[1]
  } else { hdrs <- "prompt user" }
  return(hdrs)
})
optSrc <- unique(c(fasta_types, "unknown"))
optSpc <- unique(c(fastasTbl$Species, "prompt user"))
if (!"Contaminants regex" %in% colnames(fastasTbl)) { fastasTbl$"Contaminants regex" <- "^CON__" }
if (!"Reverse regex" %in% colnames(fastasTbl)) { fastasTbl$"Reverse regex" <- "^rev_" }
# Below: commented code for automated calculation of number of reverse and contaminant hits for a given database
# Impractically slow for inclusion in the app to test regexes
# May be useful for later though...
# fastasTbl$"Regexes test"[whFnd] <- apply(fastasTbl[whFnd, c("Headers", "Type", "Contaminants regex", "Reverse regex")], 1, function(x) {
#   #x <- fastasTbl[whFnd[1], c("Headers", "Type", "Contaminants regex", "Reverse regex")]
#   mode <- toupper(x[[2]])
#   regex1 <- x[[3]]
#   regex2 <- x[[4]]
#   x <- unlist(x[[1]])
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
#   x1 <- lapply(x, function(y) {
#     z <- regexpr(PAT1, y, perl = TRUE)
#     return(substr(y,  attributes(z)$capture.start, attributes(z)$capture.start + attributes(z)$capture.length - 1))
#   })
#   x2 <- lapply(x, function(y) {
#     z <- regexpr(PAT2, y, perl = TRUE)
#     return(substr(y,  attributes(z)$capture.start, attributes(z)$capture.start + attributes(z)$capture.length - 1))
#   })
#   g1 <- grep(regex1, x1)
#   g2 <- grep(regex2, x2)
#   x <- paste0("IDs: ", length(x1), "\nCont: ", length(g1), " (", signif(round(100*length(g1)/length(x1)), 3), "%)\nRev: ",
#               length(g2), " (",  signif(round(100*length(g2)/length(x1)), 3), "%)\n")
#   cat(x)
#   return(x)
# })
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
  h5(" - Use the \"Reverse regex\" column when a fasta contains reverse proteins (decoys). Matches to the Full ID column will be removed from the database."), br(),
  actionBttn("saveBtn", "Save", icon = icon("save"), color = "success", style = "pill"),
  DTOutput("Fastas")
)
# Dummy table for display in app
fastasTbl2 <- data.frame("Name" = fastasTbl$Name,
                         "File found?" = c("-", "+")[fastasTbl$Exists+1],
                         "Species" = fastasTbl$Species,
                         "Type" = fastasTbl$Type,
                         "Contaminants-only" = fastasTbl$Contaminant,
                         check.names = FALSE)
fastasTbl2$"Contaminants-only"[which(is.na(fastasTbl2$"Contaminants-only"))] <- FALSE # Doesn't hurts
fastasTbl2$Type <- shinySelectInput(fastasTbl$Type,
                                    "Type",
                                    optSrc,
                                    paste0(30*max(c(nchar(optSrc), 2)), "px"))
fastasTbl2$"Contaminants-only" <- shinyCheckInput(fastasTbl2$"Contaminants-only",
                                                  "ContOnly")
fastasTbl2$Species <- shinySelectInput(fastasTbl$Species,
                                       "Species",
                                       optSpc,
                                       paste0(30*max(c(nchar(optSpc), 2)), "px"))
fastasTbl2$"Contaminants regex" <- fastasTbl$"Contaminants regex"
fastasTbl2$"Reverse regex" <- fastasTbl$"Reverse regex"
NmWdth <- paste0(as.character(min(c(80, max(nchar(fastasTbl2$Name))))*8), "px")
wTest <- list(list(width = NmWdth, targets = 0),
              list(width = "100px", targets = 1:6))
edith <- list(target = "column",
              enable = list(columns = grep(" regex$", colnames(fastasTbl2))-1))
tmp <- c(0:(ncol(fastasTbl2)-1))
tmp <- tmp[which(!tmp %in% edith$enable$columns)]
edith$disable <- list(columns = tmp)
if (exists("fastasTbl3")) { rm(fastasTbl3) }
server <- function(input, output, session) {
  # Output
  fastasTbl3 <- fastasTbl
  output$Fastas <- renderDT({ fastasTbl2 },
                            FALSE,
                            escape = FALSE,
                            selection = "none",
                            editable = edith,
                            rownames = FALSE,
                            options = list(dom = 't',
                                           paging = FALSE,
                                           ordering = FALSE,
                                           autowidth = TRUE,
                                           columnDefs = wTest,
                                           scrollX = TRUE),
                            ,
                            # the callback is essential to capture the inputs in each row
                            callback = JS("table.rows().every(function(i, tab, row) {
        var $this = $(this.node());
        $this.attr('id', this.data()[0]);
        $this.addClass('shiny-input-container');
      });
      Shiny.unbindAll(table.table().node());
      Shiny.bindAll(table.table().node());"))
  observeEvent(input$Fastas_cell_edit, {
    kol <- colnames(fastasTbl2)[input$Fastas_cell_edit$col+1]
    if ((length(kol))&&(kol %in% colnames(fastasTbl3))) {
      fastasTbl3[input$Fastas_cell_edit$row, kol] <<- input$Fastas_cell_edit$value
    }
  }, ignoreNULL = FALSE)
  observeEvent(input$saveBtn, {
    fastasTbl3$Contaminant <- vapply(rws, function(x) { input[[paste0("ContOnly___", as.character(x))]] }, TRUE)
    fastasTbl3$Type <- vapply(rws, function(x) { input[[paste0("Type___", as.character(x))]] }, "a")
    fastasTbl3$Species <- vapply(rws, function(x) { input[[paste0("Species___", as.character(x))]] }, "a")
    assign("fastasTbl3", fastasTbl3, envir = .GlobalEnv)
    stopApp()
  })
  #observeEvent(input$cancel, { stopApp() })
  session$onSessionEnded(function() { stopApp() })
}
runKount <- 0
while ((!runKount)||(!exists("fastasTbl3"))) {
  eval(parse(text = runApp), envir = .GlobalEnv)
  runKount <- runKount+1
}
fastasTbl %<o% fastasTbl3
#
w <- which(fastasTbl$Species == "prompt user")
if (length(w)) {
  for (i in w) {
    fastasTbl$Species[i] <- dlg_input(paste0("What is the main species in Fasta \"", fastasTbl$Full[i], "\" (ignore contaminants; write \"mixed\" if not single species)?"))$res
  }
}
Sp <- gsub(" *[\\(|\\[].*", "", fastasTbl$Species)
pack <- "myTAI"
if (!require(pack, character.only = TRUE)) {
  pak::pkg_install(pack, ask = FALSE, upgrade = TRUE, dependencies = TRUE)
}
if (!require(pack, character.only = TRUE)) {
  Src2 <- paste0(libPath, "/extdata/R scripts/Sources/taxonomy.R")
  source(Src2, local = FALSE)
}
tst %<o% try(setNames(lapply(Sp, function(x) {
  suppressMessages(taxonomy(organism = x, db = "ncbi", output = "classification"))
}), Sp), silent = TRUE)
if ("try-error" %in% class(tst)) {
  kount <- 0
  while ((kount < 20)&&("try-error" %in% class(tst))) {
    tst %<o% try(setNames(lapply(Sp, function(x) {
      suppressMessages(taxonomy(organism = x, db = "ncbi", output = "classification"))
    }), Sp), silent = TRUE)
  } 
  kount <- kount + 1
}
taxTst %<o% !"try-error" %in% class(tst)
if (taxTst) {
  Taxonomies <- tst
  w <- which(!is.na(fastasTbl$Species))
  fastasTbl$Kingdom <- fastasTbl$Taxonomy <- NA
  fastasTbl$Taxonomy[w] <- sapply(Taxonomies, function(x) { #x <- Taxonomies[[1]]
    paste(x$name, collapse = "; ")
  })
  fastasTbl$Kingdom[w] <- sapply(Taxonomies, function(x) { #x <- Taxonomies[[1]]
    x$name[match("superkingdom", x$rank)]
  })
}
source(parSrc, local = FALSE)
dbs <- lapply(whFnd, function(i) { #i <- whFnd[1] #i <- whFnd[2]
  tmp <- Format.DB(unlist(fastasTbl$Data[[i]]), in.env = TRUE, mode = fastasTbl$Type[i], parallel = TRUE, cl = parClust)
  tmp$Source <- fastasTbl$Type[i]
  tmp$"Potential contaminant" <- c("", "+")[fastasTbl$Contaminant[i]+1]
  if (!fastasTbl$Contaminant[i]) {
    g <- grep(fastasTbl$"Contaminants regex"[i], tmp$"Protein ID")
    tmp$"Potential contaminant"[g] <- "+"
    tmp$"Protein ID"[g] <- gsub(fastasTbl$"Contaminants regex"[i], "", tmp$"Protein ID"[g])
  }
  if (nchar(fastasTbl$`Reverse regex`[i])) {
    g <- try(grep(fastasTbl$"Reverse regex"[i], tmp$"Full ID", invert = TRUE), silent = TRUE)
    if (!"try-error" %in% class(g)) { tmp <- tmp[g,] }
  }
  tmp$"Protein of interest" <- (fastasTbl$Name[i] == "Proteins of interest.fasta")
  if (taxTst) {
    tmp$Taxonomy <- fastasTbl$Taxonomy[i]
    tmp$Kingdom <- fastasTbl$Kingdom[i]
  }
  return(tmp)
})
db %<o% plyr::rbind.fill(dbs)
w <- which(db$"Protein of interest")
w2 <- which((!db$`Protein ID` %in% db$`Protein ID`[w])|(db$`Protein of interest`))
db <- db[w2,]
db <- db[order(db$"Protein of interest", decreasing = TRUE),]
for (i in 1:nrow(fastasTbl)) { if (!file.exists(paste0(wd, "/", fastasTbl$Name[i]))) {
  fs::file_copy(fastasTbl$Full[i], wd)
} }
if (scrptType == "noReps") { AnalysisParam$fastasTbl <- list(fastasTbl$Full) }

# Filter for repeats, mark as contaminants
tmp <- as.data.table(db[, c("Potential contaminant", "Protein ID")])
tst <- tmp[, list(x = c(`Potential contaminant`)), by = list(Group.1 = `Protein ID`)]
tst <- as.data.frame(tst)
u <- tst$Group.1
tst$L <- vapply(tst$x, length, 1)
w <- which(tst$L > 1)
if (length(w)) {
  tst <- tst[w,]
  tst$Cont <- vapply(tst$x, function(x) { c("", "+")[("+" %in% x)+1] }, "")
  w <- which(db$`Protein ID` %in% tst$Group.1)
  db$`Potential contaminant`[w] <- tst$Cont[match(db$`Protein ID`[w], tst$Group.1)]
  db <- db[match(u, db$`Protein ID`),]
}
#
for (i in 1:nrow(fastasTbl)) { if (!file.exists(paste0(wd, "/", fastasTbl$Name[i]))) { fs::file_copy(fastasTbl$Full[i], wd) } }
if (SearchSoft == "FRAGPIPE") { # Remove reverse entries
  db <- db[grep("^>rev_", db$Header, invert = TRUE),]
}

# Also load and append contaminants database
setwd(wd)
if (SearchSoft %in% c("MAXQUANT", "DIANN", "FRAGPIPE")) {
  # NB:
  # DIA-NN does not provide an inbuilt contaminants database.
  # But we are providing one, slightly modified from the CCP's cRAPome.fasta, which should normally be used for DiaNN searches.
  # FragPipe can add the CCP's cRAPome.fasta to the search and we recommend to do so.
  # For MaxQuant, use the contaminants.fasta which is also copied with the package.
  contDBFl <- paste0(libPath, "/extData/", c("CCP_cRAPome.fasta",
                                             "contaminants.fasta")[(SearchSoft == "MAXQUANT")+1])
  contDB <- readLines(contDBFl)
  contDB <- Format.DB(contDB, in.env = TRUE)
  if (SearchSoft == "MAXQUANT") {
    contCsv <- paste0(wd, "/contaminants.csv")
    if (file.exists(contCsv)) { contDB <- read.csv(contCsv, check.names = FALSE) } else {
      # This is a bastard Fasta...
      fl <- paste0(homePath, "/Default_locations.xlsx")
      dflts <- openxlsx2::read_xlsx(fl)
      dr <- paste0(dflts$Path[which(dflts$Folder == "Fasta files")], "/Contaminants")
      if (!dir.exists(dr)) { dir.create(dr, recursive = TRUE) }
      clusterExport(parClust, "dr", envir = environment())
      tst <- parLapply(parClust, contDB$`Full ID`, function(x) { #x <- contDB$`Full ID`[1]
        if (!grepl(":", x)) {
          tmp <- paste0("https://rest.uniprot.org/uniprotkb/", x, ".fasta")
          dest <- paste0(dr, "/", x, ".fasta")
          if ((!file.exists(dest))||(file.size(dest) == 0)) {
            try(utils::download.file(tmp, dest), silent = TRUE)
            
          }
        }
      })
      w <- which(file.exists(paste0(dr, "/", contDB$`Full ID`, ".fasta")))
      tst <- parSapply(parClust, paste0(dr, "/", contDB$`Full ID`[w], ".fasta"), function(x) { #x <- paste0(dr, "/", contDB$`Full ID`[w], ".fasta")
        tst <- try(proteoCraft::Format.DB(x), silent = TRUE)
        if ("try-error" %in% class(tst)) { tst <- list(Outcome = FALSE) } else {
          tst <- list(Outcome = TRUE,
                      tbl = tst)
        }
        return(tst)
      })
      tst <- tst[which(sapply(tst, function(x) { x$Outcome }))]
      tst <- lapply(tst, function(x) { x$tbl })
      tst <- plyr::rbind.fill(tst)
      contDB <- plyr::rbind.fill(tst,
                                 contDB[which(!contDB$`Full ID` %in% tst$`Protein ID`),])
      w1 <- which(grepl("^>(sp)|(tr)\\|", contDB$Header))
      org <- unique(contDB$Organism_Full[w1])
      library(UniProt.ws)
      txIDs <- setNames(lapply(org, function(x) { c() }), org)
      w <- which(org %in% names(Taxonomies))
      if (length(w)) { txIDs[w] <- Taxonomies[w] }
      w <- which(!org %in% names(Taxonomies))
      if (length(w)) {
        txIDs[w] <- setNames(parLapply(parClust, org[w], function(x) {
          tst <- try(myTAI::taxonomy(organism = x, db = "ncbi", output = "classification"), silent = TRUE)
          if ("try-error" %in% class(tst)) { rs <- list(Outcome = FALSE) } else {
            rs <- list(Outcome = TRUE,
                       Res = tst) 
          }
        }), org[w])
        wY <- w[which(sapply(txIDs[w], function(x) { x$Outcome }))]
        wN <- w[which(!sapply(txIDs[w], function(x) { x$Outcome }))]
        txIDs[wN] <- c()
        txIDs[wY] <- lapply(txIDs[wY], function(x) { x$Res })
        Taxonomies[org[wY]] <- txIDs[org[wY]]
      }
      txIDs <- Taxonomies[which(sapply(Taxonomies, function(x) { "data.frame" %in% class(x) }))]
      txIDs <- sapply(txIDs, function(x) { x$id[which(x$rank == "species")] })
      kount <- 0
      w2 <- grep("^>(sp)|(tr)\\|", contDB$Header, invert = TRUE)
      while ((length(w2))&&(kount < length(txIDs))) {
        kount <- kount + 1
        sp <- gsub(" ", "+", txIDs[kount])
        #tmp <- paste0("https://www.uniprot.org/uniprot/?query=organism_name:", sp, "&format=fasta") # Old 
        tmp <- paste0("https://rest.uniprot.org/uniprotkb/search?query=organism_id:", sp, "&format=fasta")
        dest <- paste0("D:/Fasta_databases/", gsub(" ", "_", sp), "_-_uniprot-all-accessions_", gsub("-", "", Sys.Date()), ".fasta")
        if (!file.exists(dest)) { tmp2 <- try(download.file(tmp, dest), silent = TRUE) } else { tmp2 <- 0 }
        if ((!"try-error" %in% class(tmp2))&&(tmp2 == 0)) {
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
  db <- rbind.fill(db, contDB)
} else { if (SearchSoft == "PROTEOMEDISCOVERER") { stop("This part has not yet been re-written for PD!") } }

# Gene column!
# This column can be missing, yet it is critical for some functions (e.g. GO enrichment)
if (!"Gene" %in% colnames(db)) { db$Gene <- "" }
wY <- which((is.na(db$Gene))|(db$Gene == ""))
wN <- which((!is.na(db$Gene))&(db$Gene != ""))
k <- length(wY)
if (length(wY)) {
  tmp <- as.character(1:k)
  tst <- max(c(max(nchar(tmp)), 6))
  gns <- paste0("Gn", sapply(1:k, function(i) {
    paste(rep("0", tst-nchar(tmp[i])), collapse = "")
  }), tmp)
  wY2 <- which(gns %in% db$Gene[wN])
  wN2 <- which(!gns %in% db$Gene[wN])
  while (length(wY2)) {
    l <- length(wY2)
    tmp2 <- as.character(k+(1:l))
    tst2 <- max(c(max(nchar(tmp2)), 6))
    gns2 <- paste0("Gn", sapply(1:l, function(i) {
      paste(rep("0", tst2-nchar(tmp2[i])), collapse = "")
    }), tmp2)
    k <- k+l
    gns <- c(gns[wN2], gns2)
    wY2 <- which(gns %in% db$Gene[wN])
    wN2 <- which(!gns %in% db$Gene[wN])
  }
  db$Gene[wY] <- gns
}

# Organism columns
w <- which(c("Organism_Full", "Organism") %in% colnames(db))
tstorg %<o% (length(w) > 0)
if (tstorg) {
  dbOrgKol %<o% c("Organism_Full", "Organism")[w[1]]
  tst <- gsub(" *(\\(|\\[).*", "", db[[dbOrgKol]])
  tst <- aggregate(tst, list(tst), length)
  tst <- tst[order(tst$x, decreasing = TRUE),]
  mainOrg %<o% tst$Group.1[1]
}

# We will not need the super large "Headers" and "Data" columns anymore,
# so let's remove them to make the object smaller.
fastasTbl$Data <- NULL
fastasTbl$Headers <- NULL
