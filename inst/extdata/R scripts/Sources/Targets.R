if ("Target" %in% colnames(Exp.map)) {
  TargetProteins %<o% unique(unlist(strsplit(Exp.map$Target[which(!is.na(Exp.map$Target))], ";")))
  TargetProteins <- TargetProteins[which(!TargetProteins %in% c("NA", "", "Control"))]
  w <- which(!TargetProteins %in% db$`Protein ID`)
  if (length(w)) {
    targProt <- setNames(rep(NA, length(w)), TargetProteins[w])
    protHeads3 <- c(protHeads2[prot.list[which(prot.list %in% names(protHeads2))]], "Other")
    names(protHeads3) <- NULL
    appNm <- paste0(dtstNm, " - Baits")
    dfltProt <- setNames(lapply(w, function(x) {
      x <- targProt[x]
      if (is.na(x)) { x <- NULL } else { x <- protHeads3[x] }
      return(x)
    }), TargetProteins[w])
    dfltProt2 <- setNames(lapply(w, function(x) {
      x <- targProt[x]
      if (is.na(x)) { x <- NULL }
      return(x)
    }), TargetProteins[w])
    tmpList <- list(lapply(TargetProteins[w], function(id) {
      list(list(br()),
           fluidRow(column(5, pickerInput(id,
                                          id,
                                          protHeads3,
                                          dfltProt[[id]],
                                          inline = TRUE,
                                          width = "300px",
                                          options = pickerOptions(title = "Search me",
                                                                  `live-search` = TRUE,
                                                                  actionsBox = TRUE,
                                                                  deselectAllText = "Clear search",
                                                                  showTick = TRUE))),
                    column(1, " or... "),
                    column(6, textInput(paste0(id, "_free"),
                                        "",
                                        dfltProt2[[id]],
                                        placeholder = "... enter UniProt accession"))),
           list(br()))
    }))
    ui <- fluidPage(
      useShinyjs(),
      extendShinyjs(text = jsToggleFS, functions = c("toggleFullScreen")),
      titlePanel(tag("u", "Map pull-down baits to protein Accessions"),
                 appNm),
      br(),
      withSpinner(uiOutput("Baits")),
      actionButton("saveBtn", "Save"),
      br(),
      br()
    )
    server <- function(input, output, session) {
      TARGPROT <- reactiveVal(targProt)
      output$Baits <- renderUI({ tmpList })
      sapply(TargetProteins[w], function(id) {
        id2 <- paste0(id, "_free")
        observeEvent(input[[id]],
                     {
                       tmp <- TARGPROT()
                       tmp[id] <- input[[id]]
                       TARGPROT(tmp)
                     })
        observeEvent(input[[id2]],
                     {
                       tmp <- TARGPROT()
                       tmp[id] <- input[[id2]]
                       TARGPROT(tmp)
                     })
      })
      observeEvent(input$saveBtn, {
        sapply(TargetProteins[w], function(id) {
          targProt <<- TARGPROT()
        })
        stopApp()
      })
      #observeEvent(input$cancel, { stopApp() })
      session$onSessionEnded(function() { stopApp() })
    }
    # Modify App so that any remaining NAs turn off saving/closing the App!
    eval(parse(text = runApp), envir = .GlobalEnv)
    w1 <- which((nchar(targProt) > 0)&(!targProt %in% protHeads3))
    w2 <- which(targProt %in% protHeads3)
    if (length(w1)) {
      tmp <- lapply(paste0("https://rest.uniprot.org/uniprotkb/", targProt[w1], ".fasta"), function(x) {
        x <- try(proteoCraft::Format.DB(readLines(file(x)), in.env = TRUE), silent = TRUE)
        res <- list(Outcome = (!"try-error" %in% class(x)))
        if (res$Outcome) { res$Output <- x }
        return(res)
      })
      tmp <- tmp[which(sapply(tmp, function(x) { x$Outcome }))]
      if (length(tmp)) {
        tmp <- plyr::rbind.fill(tmp)
        db <- plyr::rbind.fill(db, tmp)
      }
    }
    if (length(w2)) {
      tmp <- targProt[w2]
      tmp <- tmp[which(tmp != "Other")]
      m <- match(tmp, protHeads2)
      w3 <- which(!is.na(m))
      targProt[w2[w3]] <- names(protHeads2)[m[w3]]
    }
    w <- which(targProt == "")
    targProt[w] <- names(targProt)[w]
    Exp.map$"Bait (aka Target) name" <- Exp.map$Target # Alias as backup
    Exp.map$Target <- setNames(targProt[Exp.map$Target], NULL)
  }
  Tar <- unique(Exp.map$Target)
  FactorsLevels$Target <- Tar
  tmp <- list(Factors = Factors, Levels = FactorsLevels)
  save(tmp, file = "Factors.RData")
  # Modify app so that any remaining NAs turn off saving!
}
