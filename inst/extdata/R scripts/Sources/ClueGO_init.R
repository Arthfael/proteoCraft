# Initialize ClueGO
clueGOahead %<o% FALSE
if (CytoScape) {
  txt_2_df %<o% \(table.text) {
    table <- NULL
    rows <- unlist(strsplit(table.text, "\n"))
    header <- t(unlist(strsplit(rows[1L], "\t")))
    if (length(rows) > 1L) {
      for (i in 2:length(rows)) {
        if (is.null(table)) {
          table <- t(unlist(strsplit(rows[i], "\t")))
        } else {
          table <- rbind(table, t(unlist(strsplit(rows[i], "\t"))))
        }
      }
      table <- as.data.frame(table, col.names = header)
    } else {
      warning("Empty table!")
      table <- as.data.frame(matrix(nrow = 1L, ncol = length(header)))
    }
    names(table) <- header
    return(table)
  }
  clueGO_Home %<o% gsub("/Documents$", "", normalizePath(Sys.getenv("HOME"), winslash = "/"))
  clueGO_Home2 %<o% paste0(clueGO_Home, "/ClueGOConfiguration")
  clueGO_Vers %<o% grep("v", list.dirs(clueGO_Home2, recursive = FALSE, full.names = FALSE), value = TRUE)
  if (length(clueGO_Vers) > 1L) {
    tst <- as.data.frame(t(as.data.frame(strsplit(gsub("^[Vv]", "", clueGO_Vers), "\\."))))
    tst$Vers <- clueGO_Vers
    clueGO_Vers <- tst$Vers[order(tst$V1, tst$V2, tst$V3)[1L]]
  }
  clueGO_Home3 %<o% paste(clueGO_Home2, clueGO_Vers, sep = "/")
  myPort %<o% 1234L
  myHost %<o% "localhost"
  cytoscp_URL %<o% paste0("http://", myHost, ":", toString(myPort), "/v1")
  clueGO_URL %<o% paste(cytoscp_URL, "apps", "cluego", "cluego-manager", sep = "/")
  #### Start up ClueGO
  response <- httr::POST(url = paste(cytoscp_URL, "apps", "cluego", "start-up-cluego", sep = "/"), encode = "json")
  if (httr::http_status(response)$category != "Success") {
    Sys.sleep(2L)
    httr::http_status(response)$category != "Success"
  }
  msg <- ""
  if (httr::http_status(response)$category != "Success") {
    msg <- "It appears clueGO isn't installed in Cytoscape! Install it and try again!"
  }
  #
  #### Select the ClueGO Organism to analyze
  rqst <- paste0(clueGO_URL, "/organisms/get-all-organism-info")
  response <- httr::GET(rqst, encode = "json")
  if (httr::http_status(response)$category != "Success") {
    stop(msg)
  } else { msg <- "" }
  clueGOrgs <- httr::content(response, encode = "json")
  myOrgNm %<o% Org$Organism[1L]
  myOrgNm <- unique(c(myOrgNm, sapply(strsplit(myOrgNm, " "), \(x) { #x <- strsplit(myOrgNm, " ")
    # I have to do this because ClueGO is inconsistent with capitals in species names!
    x <- unlist(x)
    l <- length(x)
    if (length(x) >= 2L) {
      prt2 <- tolower(x[2L])
      nc <- nchar(prt2)
      prt2 <- c(prt2, paste0(toupper(substr(prt2, 1L, 1L)), substr(prt2, 2L, nc)))
    }
    rs <- paste0(x[1L], " ", prt2)
    if (l > 2L) { rs <- paste(c(rs, x[3L:l]), collapse = " ") }
    return(rs)
  })))
  orgNm <- myOrgNm[which(myOrgNm %in% names(clueGOrgs))]
  clueGOahead <- (length(orgNm) > 0L)
  if (!clueGOahead) {
    msg <- "I could not identify the organism automatically, select organism from this list"
    dlNew <- " -> download new organism..."
    optNms <- c(names(clueGOrgs), dlNew)
    opt <- setNames(sapply(optNms, \(x) { paste(c(x, rep(" ", 250L-nchar(x))), collapse = "") }), optNms)
    tmp <- dlg_list(opt, opt[dlNew], title = msg)$res
    orgNm <- optNms[match(tmp, opt)]
    if ((length(orgNm) == 1L)&&(orgNm == dlNew)) {
      msg <- paste0("Download the new organism in ClueGO from within the CytoScape GUI\n(click on the leftmost blue paw icon)\nthen click \"ok\" to continue.")
      dlg_message(msg, "ok")
      rqst <- paste0(clueGO_URL, "/organisms/get-all-organism-info")
      response <- httr::GET(rqst, encode = "json")
      clueGOrgs <- httr::content(response, encode = "json")
      orgNm <- myOrgNm[which(myOrgNm %in% names(clueGOrgs))]
      if (!length(orgNm)) {
        msg <- "I could not identify the organism automatically, select organism from this list"
        optNms <- names(clueGOrgs)
        opt <- setNames(sapply(optNms, \(x) { paste(c(x, rep(" ", 250L-nchar(x))), collapse = "") }), optNms)
        tmp <- dlg_list(opt, opt[dlNew], title = msg)$res
        orgNm <- optNms[match(tmp, opt)]
      }
    }
    clueGOahead <- (length(orgNm) == 1L)
  }
  if (clueGOahead) {
    orgNm <- orgNm[1L]
    #
    # Number of sets
    #N_sets %<o% 2
    N_sets %<o% 1L
    myClust %<o% 1L # Regulated only, Background is now provided as custom .txt file
    nodeShape %<o% "Ellipse" # ("Ellipse", "Diamond", "Hexagon", "Octagon", "Parallelogram", "Rectangle", "Round Rectangle", "Triangle", "V")
    clustCols %<o% c("#00FF00", "#FF0000") # The color in hex, e.g. #F3A455
    clueGO_types <- c("Enrichment (Right-sided hypergeometric test)",
                      "Enrichment/Depletion (Two-sided hypergeometric test)")
    noRestr %<o% "false" # # so min #/% genes are not ignored!
    visStyles %<o% setNames(c("ShowClusterDifference", "ShowGroupDifference", "ShowSignificanceDifference"),
                            c("clustDiff", "groupDiff", "signifDiff"))
    myOnt %<o% RJSONIO::toJSON(c("3;Ellipse", "8;Triangle", "9;Rectangle")) # (run "3.1 Get all available Ontologies" to get all options)
    minLvl %<o% 5L
    maxLvl %<o% 6L
    allLvl %<o% "false" # so the 2 parameters above are not ignored!
    clueGO_opt %<o% "Continue analysis" # ("Continue analysis", "Skip the grouping", "Cancel and refine selection")  -> Analysis option in case there are more than 1000 terms found!
    #
    rqst <- paste(clueGO_URL, "organisms", "set-organism", URLencode(orgNm), sep = "/")
    response <- httr::PUT(rqst, encode = "json")
    if (httr::http_status(response)$category != "Success") { Sys.sleep(2L) }
    rqst <- paste(clueGO_URL, "cluster", "max-input-panel", N_sets, sep = "/")
    response <- httr::PUT(rqst, encode = "json")
    if (httr::http_status(response)$category != "Success") { Sys.sleep(2L) }
    cat(paste0("Analysis Parameters for ", dtstNm,
               "\nUser Home Folder: ", clueGO_Home,
               "\nCytoscape Base URL: ", cytoscp_URL,
               "\nClueGO Base URL: ", clueGO_URL,
               "\nOrganism: ", orgNm, "\n"))
    #
  }
}
