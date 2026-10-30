if ("Target" %in% colnames(Exp.map)) {
  TargetProteins %<o% unique(unlist(strsplit(Exp.map$Target[!is.na(Exp.map$Target)], ";")))
  TargetProteins <- base::setdiff(TargetProteins, c("NA", "", "Control"))
  w <- which(!TargetProteins %in% db$`Protein ID`)
  if (length(w)) {
    targProt <- setNames(rep(NA, length(w)), TargetProteins[w])
    protHeads3 <- c(protHeads2[base::intersect(prot.list, names(protHeads2))], "Other")
    names(protHeads3) <- NULL
    appNm <- paste0(dtstNm, " - Baits")
    dfltProt <- setNames(lapply(w, \(x) {
      x <- targProt[x]
      if (is.na(x)) { x <- NULL } else { x <- protHeads3[x] }
      return(x)
    }), TargetProteins[w])
    dfltProt2 <- setNames(lapply(w, \(x) {
      x <- targProt[x]
      if (is.na(x)) { x <- NULL }
      return(x)
    }), TargetProteins[w])
    tmpList <- list(lapply(TargetProteins[w], \(id) {
      list(list(shiny::br()),
           shiny::fluidRow(shiny::column(5L,
                                         shinyWidgets::pickerInput(id,
                                                                   id,
                                                                   protHeads3,
                                                                   dfltProt[[id]],
                                                                   inline = TRUE,
                                                                   width = "300px",
                                                                   options = shinyWidgets::pickerOptions(title = "Search me",
                                                                                                         `live-search` = TRUE,
                                                                                                         actionsBox = TRUE,
                                                                                                         deselectAllText = "Clear search",
                                                                                                         showTick = TRUE))),
                           shiny::column(1L,
                                         " or... "),
                           shiny::column(6L,
                                         shiny::textInput(paste0(id, "_free"),
                                                          "",
                                                          dfltProt2[[id]],
                                                          placeholder = "... enter UniProt accession"))),
           list(shiny::br()))
    }))
    ui <- shiny::fluidPage(
      shinyjs::useShinyjs(),
      shinyjs::extendShinyjs(text = jsFun$toggle_full_screen,
                             functions = c("toggleFullScreen")),
      shiny::titlePanel(shiny::tag("u", "Map pull-down baits to protein Accessions"),
                        appNm),
      shiny::br(),
      shiny::h5("Map Target names to accessions in the search database."),
      shiny::h5(shiny::em("(Ignore cases, such as Control/Ctrl/IgGminus/etc... where this is not feasible.)")),
      shinycssloaders::withSpinner(shiny::uiOutput("Baits")),
      shiny::actionButton("saveBtn", "Save"),
      shiny::br(),
      shiny::br()
    )
    server <- \(input, output, session) {
      TARGPROT <- shiny::reactiveVal(targProt)
      output$Baits <- shiny::renderUI({ tmpList })
      sapply(TargetProteins[w], \(id) {
        id2 <- paste0(id, "_free")
        shiny::observeEvent(input[[id]],
                            {
                              tmp <- TARGPROT()
                              tmp[id] <- input[[id]]
                              TARGPROT(tmp)
                              if (sum(is.na(tmp))) {
                                shinyjs::disable("saveBtn")
                              } else {
                                shinyjs::enable("saveBtn")
                              }
                            })
        shiny::observeEvent(input[[id2]],
                            {
                              tmp <- TARGPROT()
                              tmp[id] <- input[[id2]]
                              TARGPROT(tmp)
                              if (sum(is.na(tmp))) {
                                shinyjs::disable("saveBtn")
                              } else {
                                shinyjs::enable("saveBtn")
                              }
                            })
      })
      shiny::observeEvent(input$saveBtn, {
        sapply(TargetProteins[w], \(id) {
          assign("targProt", TARGPROT(), envir = .GlobalEnv)
        })
        shiny::stopApp()
      })
      #observeEvent(input$cancel, { shiny::stopApp() })
      session$onSessionEnded(\() { shiny::stopApp() })
    }
    # Modify App so that any remaining NAs turn off saving/closing the App!
    g <- shiny:::.globals
    g$appState <- NULL
    eval(parse(text = run_App), envir = .GlobalEnv)
    shinyCleanup()
    w1 <- which((nchar(targProt) > 0L) & (!targProt %in% protHeads3))
    w2 <- which(targProt %in% protHeads3)
    if (length(w1)) {
      tmp <- lapply(paste0("https://rest.uniprot.org/uniprotkb/", targProt[w1], ".fasta"), \(x) {
        x <- try(Format.DB(readLines(file(x)), in.env = TRUE), silent = TRUE)
        res <- list(Outcome = !inherits(x, "try-error"))
        if (res$Outcome) { res$Output <- x }
        return(res)
      })
      tmp <- tmp[vapply(tmp, \(x) { x$Outcome }, TRUE)]
      if (length(tmp)) {
        tmp <- plyr::rbind.fill(tmp)
        db <- plyr::rbind.fill(db, tmp)
      }
    }
    if (length(w2)) {
      tmp <- targProt[w2]
      tmp <- base::setdiff(tmp, "Other")
      m <- match(tmp, protHeads2)
      w3 <- which(!is.na(m))
      targProt[w2[w3]] <- names(protHeads2)[m[w3]]
    }
    w <- which(targProt == "")
    targProt[w] <- names(targProt)[w]
    Exp.map$"Bait (aka Target) name" <- Exp.map$Target # Alias as backup
    w <- which(Exp.map$Target %in% names(targProt))
    if (length(w)) {
      Exp.map$Target[w] <- setNames(targProt[Exp.map$Target[w]], NULL)
    }
  }
  Tar <- unique(Exp.map$Target)
  FactorsLevels$Target <- Tar
  tmp <- list(Factors = Factors, Levels = FactorsLevels)
  save(tmp, file = "Factors.RData")
}
