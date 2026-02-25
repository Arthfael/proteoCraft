########################################
# Source for running a ClueGO analysis #
########################################
#
# Should be run immediately after GO_enrich.R, as it utilizes some of the temporary variables created by the former!
#
if (clueGOahead) {
  cat("\n   ClueGO analysis\n   ---------------\n")
  nuWidth <- 1800
  custRefFl <- paste0(subfolder, "/ref.txt")
  m <- match(clueGO_type, clueGO_types)
  Kappa <- c(0.5, 0.8)[m]
  minGenesPerTerm %<o% c(3, 10)[m]
  minPercGenesMapped %<o% c(4, 10)[m]
  cat(paste0("   | Parameters:",
             "\n   |> min. genes per term = ", minGenesPerTerm,
             "\n   |> min. % of genes mapped = ", minPercGenesMapped,
             "\n   |> min. GO level = ", minLvl,
             "\n   |> max. GO level = ", maxLvl, 
             "\n   |> Kappa level = ", Kappa, " ]\n"))
  for (nm in names(filters)) { #nm <- names(filters)[1]
    nm1 <- gsub("___", "_", nm)
    nm2 <- proteoCraft::cleanNms(nm)
    cat(paste0("      - Filter = ", nm2, "\n"))
    clueTST <- try({
      # Request again the correct organism, just to be sure!
      rqst <- paste(clueGO_URL, "organisms", "set-organism", URLencode(orgNm), sep = "/")
      response <- httr::PUT(rqst, encode = "json")
      if (httr::http_status(response)$category != "Success") { Sys.sleep(2) }
      #
      for (ii in 1:2) { #ii <- 1 #ii <- 2
        # (ii == 1 is for the filter, ii == 2 is the reference)
        #
        # Useful guide to the Cytoscape RESTful API:
        # http://127.0.0.1:1234/v1/swaggerUI/swagger-ui/index.html?url=http://127.0.0.1:1234/v1/swagger.json#/Apps5832ClueGO
        # (Includes a section on ClueGO)
        #
        # Create cluster
        if (ii == 1) {
          rqst <- paste(clueGO_URL, "cluster", "set-analysis-properties", ii, URLencode(nodeShape),
                        URLencode(clustCols[ii], reserved = TRUE), 
                        minGenesPerTerm, minPercGenesMapped, noRestr, sep = "/")
          response <- httr::PUT(rqst, encode = "json")
          if (httr::http_status(response)$category != "Success") { Sys.sleep(2) }
        }
        #
        # Note to self:
        #  - use httr::http_status(response) to inspect the outcome, print(response) fails (there does not seem to be a print method for this class of object)
        #  - using print by error will return an error even though the call actually worked!
        #
        # Set gene lists
        if (ii == 1) {
          myFlt <- myFlt1 <- filters[[nm]]
          tmpDat <- Prot[myFlt,]
        }
        if (ii == 2) {
          myFlt <- myFlt2 <- ref.filters[[nm]]
          tmpDat <- parentData[myFlt,]
        }
        if ("Potential contaminant" %in% colnames(tmpDat)) {
          w <- which((is.na(tmpDat$"Potential contaminant"))|(tmpDat$"Potential contaminant" != "+"))
          tmpDat <- tmpDat[w,]
        }
        if (ii == 1) {
          IDs_lst <- unique(unlist(strsplit(tmpDat[myFlt, ID_col], ";"))) # inherited from GO_enrich.R, which should always be run before!
        }
        if (ii == 2) {
          IDs_lst <- unique(unlist(strsplit(tmpDat[, parentCol], ";"))) # inherited from GO_enrich.R, which should always be run before!
        }
        #
        IDs_lst <- grep("^((cRAP)|(CON__))", IDs_lst, value = TRUE, invert = TRUE)
        IDs_lst <- IDs_lst[which(IDs_lst != "NA")]
        #writeClipboard(IDs_lst)
        if (ii == 1) {
          IDs_lst1 <- IDs_lst
          IDs_lst <- RJSONIO::toJSON(IDs_lst)
          rqst <- paste(clueGO_URL, "cluster", "upload-ids-list", URLencode(as.character(ii)), sep = "/")
          response <- httr::PUT(rqst, body = IDs_lst, encode = "json", httr::content_type_json())
          if (httr::http_status(response)$category != "Success") { Sys.sleep(2) }
        }
        if (ii == 2) {
          IDs_lst2 <- IDs_lst
          stopifnot(sum(!IDs_lst1 %in% IDs_lst2) == 0)
          IDs_lst <- RJSONIO::toJSON(IDs_lst)
          # We are instead saving locally and passing it later as reference list
          #rqst <- paste(clueGO_URL, "cluster", "upload-ids-list", URLencode(as.character(ii)), sep = "/")
          #response <- httr::PUT(rqst, body = IDs_lst, encode = "json", httr::content_type_json())
          #if (httr::http_status(response)$category != "Success") { Sys.sleep(2) }
          write(IDs_lst, custRefFl)
        }
      }
      rqst <- paste(clueGO_URL, "ontologies", "set-ontologies", sep = "/")
      response <- httr::PUT(rqst, body = myOnt, encode = "json", httr::content_type_json())
      if (httr::http_status(response)$category != "Success") { Sys.sleep(2) }
      rqst <- paste(clueGO_URL, "ontologies", "set-min-max-levels", minLvl, maxLvl, allLvl, sep = "/")
      response <- httr::PUT(rqst, encode = "json")
      if (httr::http_status(response)$category != "Success") { Sys.sleep(2) }
      rqst <- paste0(clueGO_URL, "/ontologies/get-kappa-score-level/", Kappa)
      response <- httr::PUT(rqst, encode = "json")
      if (httr::http_status(response)$category != "Success") { Sys.sleep(2) }
      #rqst <- paste0(clueGO_URL, "/stats/", URLencode(clueGO_type), "/Benjamini-Hochberg/false/false/false")
      rqst <- paste0(clueGO_URL, "/stats/", URLencode(clueGO_type), "/Benjamini-Hochberg/false/false/true/",
                     URLencode(custRefFl, TRUE))
      response <- httr::PUT(rqst, encode = "json")
      if (httr::http_status(response)$category != "Success") { Sys.sleep(2) }
      #
      # Run the analysis and save log file
      rqst <- paste(clueGO_URL, URLencode(nm1), URLencode(clueGO_opt), sep = "/")
      response <- httr::GET(rqst, httr::timeout(3600)) # 1h seems absolutely ginormous!
      if (httr::http_status(response)$category != "Success") { Sys.sleep(2) }
      logFl <- paste0(clueGO_outDir, "/", nm2, "-log.txt")
      msg <- response$message
      tst <- try(writeLines(httr::content(response, encoding = "UTF-8"), logFl), silent = TRUE)
      if (!"try-error" %in% class(tst)) {
        #print(httr::content(response, encode = "text"))
        #
        # Get network id (SUID) (CyRest function from Cytoscape)
        rqst <- paste(cytoscp_URL, "networks/currentNetwork", sep = "/")
        response <- httr::GET(rqst)
        if (httr::http_status(response)$category != "Success") { Sys.sleep(2) }
        currNtwrk_SUID <- httr::content(response, encode = "json")$data$networkSUID
        #print(currNtwrk_SUID)
        #
        # Fit the first available Network View  to the current window.
        rqst <- paste(cytoscp_URL, "apply/fit", currNtwrk_SUID, sep = "/")
        response <- httr::GET(rqst)
        if (httr::http_status(response)$category != "Success") { Sys.sleep(2) }
        #
        # Save network
        #
        cxFl <- paste0(clueGO_outDir, "/", nm2, ".cx")
        for (j in 2:3) {
          rqst <- paste(clueGO_URL, "cluster", "select-visual-style", visStyles[j], sep = "/")
          response <- httr::PUT(rqst, encode = "json")
          if (httr::http_status(response)$category != "Success") { Sys.sleep(2) }
          #
          # Save network
          if (j == 2) {
            exportNetwork(cxFl, "CX", currNtwrk_SUID, cytoscp_URL, TRUE)
          }
          #
          # Get network graphics (CyRest function from Cytoscape)
          Sys.sleep(2) # Precaution
          imgType <- "svg" # png, pdf
          rqst <- paste(cytoscp_URL, "networks", currNtwrk_SUID, "views", paste0("first.", imgType), sep = "/")
          response <- httr::GET(rqst)
          if (httr::http_status(response)$category != "Success") { Sys.sleep(2) }
          imgFl <- gsub("\\.cx$", paste0("_", names(visStyles)[j], ".", imgType), cxFl)
          writeBin(httr::content(response, encode = "raw"), imgFl)
          # Adjust svg's default display size - (not sure this works)
          a <- readLines(imgFl)
          l <- length(a)
          g1 <- grep("^ *<svg *", a)
          g2 <- grep(">", a[g1:l])[1]
          rg <- g1:(g2+g1-1)
          gw <- rg[grep("^ *width=\"", a[rg])[1]]
          gh <- rg[grep("^ *height=\"", a[rg])[1]]
          # NB: do not edit the dimensions in viewBox as they define how the internal dimensions will be scaled to fit the global width and height!
          #gv <- rg[grep("^ *viewBox=\"", a[rg])[1]]
          wTxt <- a[gw]
          hTxt <- a[gh]
          #vTxt <- a[gv]
          w <- as.numeric(gsub("^ *width=\"|px\".*", "", wTxt))
          h <- as.numeric(gsub("^ *height=\"|px\".*", "", hTxt))
          #v <- as.numeric(unlist(strsplit(gsub("^ *viewBox=\"|\".*", "", vTxt), " ")))
          corr <- nuWidth/w
          wTxt <- unlist(strsplit(wTxt, "width=\"[0-9]+(\\.[0-9]+)?px"))
          hTxt <- unlist(strsplit(hTxt, "height=\"[0-9]+(\\.[0-9]+)?px"))
          #vTxt <- unlist(strsplit(hTxt, "viewBox=\"([0-9]+(\\.[0-9]+)? +){3}([0-9]+(\\.[0-9]+)?)"))
          wTxt <- paste0(wTxt[1], "width=\"", ceiling(w*corr), "px", wTxt[2])
          hTxt <- paste0(hTxt[1], "height=\"", ceiling(h*corr), "px", hTxt[2])
          #vTxt <- paste0(vTxt[1], "viewBox=\"", paste(c(0, 0, ceiling(w*corr), ceiling(h*corr)), collapse = " "), vTxt[2])
          a[gw] <- wTxt
          a[gh] <- hTxt
          #a[gv] <- vTxt
          writeLines(a, imgFl)
          #system(paste0("open \"", imgFl, "\""))
        }
      }
    }, silent = TRUE)
    if ("try-error" %in% class(clueTST)) {
      warning(paste0("ClueGO analysis failed for ", nm2))
    }
  }
  if (file.exists(custRefFl)) { unlink(custRefFl) }
  cat("     Done!\n")
}
