# Define Contrasts of interest

# expMap: limma-compatible version of Exp.map
hasBlocks %<o% FALSE
origCoeff <- union(RSA$names, RRG$names)
if ((!is.na(Param$Blocking.factors)) && (Param$Blocking.factors != "")) {
  parse.Param.aggreg.2("Blocking.factors")
  origCoeff <- union(origCoeff, Blocking.factors$names)
  hasBlocks <- TRUE
}
if ((!is.na(Param$Batch.effect)) && (Param$Batch.effect != "")) {
  parse.Param.aggreg.2("Batch.effect")
  origCoeff <- union(origCoeff, Batch.effect$names)
}
Coefficients %<o% origCoeff
expMap %<o% Exp.map[, origCoeff]
rownames(expMap) <- Exp.map$Ref.Sample.Aggregate
#
#
# Replace hyphens by dots and prepend root to numeric factor levels to avoid issues with evaluating contrasts
# - > This chunk below must remain:
# Yes, in the Experimental_Factors_editor.R source, we already dealt with hyphens...
# except for the Target factor (where values are protein accessions and may include hyphens, e.g. for UniProtKB isoforms).
# Also, numeric factors such as replicates cannot be used as column names in the design matrix... and for good reason.
# Other exceptions may come up...
# This is also an additional reason to keep Exp.map and expMap as two separate representations of the same experimental structure.
# This also creates a factor for each coefficient.
for (Coeff in origCoeff) { #Coeff <- "Replicate" #Coeff <- origCoeff[2L]
  assign(Coeff, factor(expMap[[Coeff]], FactorsLevels[[Coeff]]))
  nuCoeff <- paste0(Coeff, "_._")
  stopifnot(!nuCoeff %in% colnames(Exp.map)) # Not expMap in case we are re-running a small chunk
  expMap[[nuCoeff]] <- expMap[[Coeff]]
  l <- length(grep("-", Exp.map[[Coeff]])) # Not expMap in case we are re-running a small chunk
  if (l) {
    expMap[[nuCoeff]] <- gsub("-", ".", expMap[[Coeff]])
    assign(nuCoeff, factor(expMap[[nuCoeff]], gsub("-", ".", FactorsLevels[[Coeff]])))
  }
  w <- which(suppressWarnings(as.character(as.numeric(as.character(Exp.map[[Coeff]]))) == as.character(Exp.map[[Coeff]])))
  if (length(w)) {
    expMap[w, nuCoeff] <- paste0(Coeff, "_", as.character(expMap[w, Coeff]))
    assign(nuCoeff, factor(expMap[[nuCoeff]], paste0(Coeff, "_", as.character(FactorsLevels[[Coeff]]))))
  }
  Coefficients[match(Coeff, origCoeff)] <- nuCoeff
  assign(nuCoeff, expMap[[nuCoeff]])
  nuCoeff %<o% nuCoeff
  expMap[[Coeff]] <- NULL
}
#
for (aggrNm in c("RSA", "VPAL", "RG", "RRG", "Blocking.factors", "Batch.effect")) { # Update in case we had to do some fixing
  #aggrNm <- "RSA"
  if (exists(aggrNm)) {
    aggr <- get(aggrNm)
    tmp <- aggr$names
    tmp <- Coefficients[match(tmp, origCoeff)]
    nm <- paste(tmp, collapse = "_")
    if (!grepl("_\\._$", nm)) { nm <- paste0(nm, "_._") }
    aggr$limmaCol <- nm
    expMap[[nm]] <- do.call(paste, c(expMap[, tmp, drop = FALSE], sep = "___"))
    expMap[[nm]] <- factor(expMap[[nm]], levels = unique(expMap[[nm]]))
    nm %<c% as.factor(expMap[[nm]])
    if (length(Exp) == 1L) { tmp <- setdiff(tmp, "Experiment_._") }
    if (length(tmp)) {
      nm <- paste(tmp, collapse = "_")
      if (!grepl("_\\._$", nm)) { nm <- paste0(nm, "_._") }
      aggr$limmaCol <- nm
      if (!nm %in% colnames(expMap)) {
        expMap[[nm]] <- do.call(paste, c(expMap[, tmp, drop = FALSE], sep = "___"))
        expMap[[nm]] <- factor(expMap[[nm]], levels = unique(expMap[[nm]]))
        nm %<c% factor(expMap[[nm]])
      }
    }
    assign(aggrNm, aggr)
  }
}
#
#
#
expMap <- expMap[order(#expMap[[RRG$limmaCol]], # Do not use RRG here
  expMap[[RG$limmaCol]], # Cf. below: safe because each RG contains at least one ref samples group
  # (Otherwise every gets confusing and  my head starts hurting...)
  expMap[[VPAL$limmaCol]],
  expMap$Replicate_._),]

# Make contrasts
ratGrps <- unique(expMap[[RG$limmaCol]])
contrBlocks <- setNames(lapply(ratGrps, \(grp) {
  x <- as.character(unique(expMap[which(expMap[[RG$limmaCol]] == grp), VPAL$limmaCol]))
  setNames(x, cleanNms(x, rep = " "))
}), ratGrps)
contrBlocks2 <- stack(contrBlocks)
colnames(contrBlocks2) <- c("Samples group", "Comparison group")
contrBlocks2$Name <- rownames(contrBlocks2)
rownames(contrBlocks2) <- NULL
tmp1 <- cleanNms(VPAL$values, rep = " ")
tmp2 <- cleanNms(VPAL$values, rep = "_")
contrBlocks2$Name_ <- tmp2[match(contrBlocks2$Name, tmp1)]
dfltContr_Opt <- setNames(lapply(ratGrps, \(grp) { #grp <- ratGrps[1L]
  x <- as.data.frame(gtools::permutations(length(contrBlocks[[grp]]), 2L, names(contrBlocks[[grp]])))
  colnames(x) <- c("A", "B")
  do.call(paste, c(x[, c("A", "B")], sep = " - "))
}), ratGrps)
alsoDouble <- (length(ratGrps) > 1L) || (length(dfltContr_Opt[[1L]]) > 3L) # We need 4 different things for an interaction double contrast
dfltContr_Opt <- stack(dfltContr_Opt)
colnames(dfltContr_Opt) <- c("Contrast", "Comparison group")
dfltContr_Opt[, c("A", "B")] <- do.call(rbind, strsplit(dfltContr_Opt$Contrast, " - "))
#
fullContrFun <- \(prim, sec) {
  if (missing(sec) || is.na(sec) || (length(sec) != 1L) || (sec == "")) { return(prim) }
  return(paste0("(", prim, ") - (", sec, ")"))
}
contrastsFl %<o% paste0(wd, "/Contrasts.rds")
if (file.exists(contrastsFl)) {
  myContrasts <- readr::read_rds(contrastsFl)
  g <- grep("\\) - \\(", myContrasts$Contrast)
  if (length(g) && (!alsoDouble)) {
    warning("Invalid contrasts reloaded!")
    myContrasts <- myContrasts[-g,]
  }
  myContrasts$Primary <- myContrasts$Contrast
  if (alsoDouble) {
    myContrasts$Secondary <- ""
    g <- grep("\\) - \\(", myContrasts$Contrast)
    if (length(g)) {
      myContrasts$Primary[g] <- sub("^\\(", "", sub("\\) - \\(.*", "", myContrasts$Contrast[g]))
      myContrasts$Secondary[g] <- sub("\\)$", "", sub(".*\\) - \\(", "", myContrasts$Contrast[g]))
    }
  }
} else {
  myContrasts <- data.frame("Primary" = character(),
                            "Secondary" = character(),
                            "Contrast" = character(),
                            "Up-only" = logical(),
                            check.names = FALSE)
  if (!alsoDouble) { myContrasts$Secondary <- NULL }
}
myContrasts %<o% myContrasts
myContrasts2 <- myContrasts[, c("Contrast", "Up-only")]
nr2 <- nrow(myContrasts2)
myContrasts2$Remove <- if (nr2) {
  vapply(1L:nr2, \(i) {
    iChr <- as.character(i)
    as.character(actionButton(paste0("rmvBtn_", iChr), "remove contrast"))
  }, "")
} else { character() }
AB <- dfltContr_Opt[match(dfltContr_Opt$Contrast[1L], dfltContr_Opt$Contrast), c("A","B")]
tmp <- c("", dfltContr_Opt$Contrast[(!dfltContr_Opt$A %in% AB)&(!dfltContr_Opt$B %in% AB)])
makeContr <- data.frame("Contrast" = as.character(selectInput("Primary", "", dfltContr_Opt$Contrast)),
                        "(opt. secondary contrast)" = as.character(selectInput("Secondary", "", tmp)),
                        "Up-regulated only?" = as.character(checkboxInput("upOnly", "", FALSE)),
                        "add contrast" = as.character(actionButton("addContr", "add contrast")),
                        check.names = FALSE)
colDefs1 <- list(list(width = "250px", targets = 0L),
                 list(width = "100px", targets = 1L+alsoDouble),
                 list(width = "50px", targets = 2L+alsoDouble))
colDefs2 <- list(list(width = "500px", targets = 0L),
                 list(width = "100px", targets = 1L+alsoDouble))
if (alsoDouble) {
  colDefs1 <- append(colDefs1,
                     list(list(width = "250px", targets = 1L)))
  colDefs2 <- append(colDefs2,
                     list(list(width = "100px", targets = 1L)))
} else {
  makeContr$"(opt. secondary contrast)" <- NULL
}
appNm <- "Define contrasts"
make_ui0 <- \() {
  shinyUI(
    fluidPage(
      shinyjs::useShinyjs(),
      shinyWidgets::setBackgroundColor( # Doesn't work
        color = c(#"#F8F8FF",
          "#EBEFF7"),
        gradient = "linear",
        direction = "bottom"
      ),
      shinyjs::extendShinyjs(text = jsToggleFS, functions = c("toggleFullScreen")),
      tags$head(tags$script(HTML("
Shiny.addCustomMessageHandler('updateSecondary', function(data) {
  var sel = document.getElementById('Secondary');
  if (!sel) return;

  sel.innerHTML = '';

  data.choices.forEach(function(choice) {
    var opt = document.createElement('option');
    opt.value = choice;
    opt.text = choice;
    sel.appendChild(opt);
  });
});
"))),
      tags$head(tags$style(HTML("table {table-layout: fixed;"))), # So table widths can be properly adjusted!
      titlePanel(tag("u", "Define contrasts of interest"),
                        appNm),
      mainPanel(
        em("Select each \"A - B\" contrast of interest, then click \"Add\" to add it to the analysis."),
        br(),
        em("(NB: \"A - B\" means \"A versus B\", i.e. FC = A/B and logFC = A-B)"),
        br(),
        br(),
        em("Tick \"one-sided\" for \"up-regulated\" (for \"down-regulated\", just select the reverse contrast first)."),
        br(),
        if (alsoDouble) {
          em(HTML("<br>You may also add double contrasts (interaction contrasts) of the form \"(A - B) - (C - D)\" (where all 4 are distinct)."))
        },
        if (alsoDouble) {
          em(HTML("<br><br>For now we are restricting users to \"sensible contrasts\" - this may change!"))
        },
        if (alsoDouble) {
          br()
        },
        em(HTML("<br><br>Once you are finished, click&nbsp;")),
        shinyWidgets::actionBttn("saveBtn", "Save", icon = icon("save"), color = "success", style = "pill"),
        em(HTML("&nbsp;to continue.")),
        br(),
        em("Saving is only possible once at least one contrast has been entered"),
        em("(otherwise why even run this workflow?)"),
        br(),
        br(),
        DT::DTOutput("makeContrasts"),
        uiOutput("msg"),
        br(),
        br(),
        br(),
        em("-> your current contrasts:"),
        DT::DTOutput("myContrasts"),
        br(),
        br()
      )))
}
server0 <- shinyServer(\(input, output, session) {
  CONTRASTSTBL <- reactiveVal(myContrasts2)
  MSG <- reactiveVal("")
  formVALS <- reactiveValues(Primary = dfltContr_Opt$Contrast[1L],
                             Secondary = "",
                             upOnly = FALSE)
  output$makeContrasts <- DT::renderDT({ makeContr },
                                       FALSE,
                                       escape = FALSE,
                                       class = "compact",
                                       selection = "none",
                                       rownames = FALSE,
                                       editable = FALSE,
                                       options = list(dom = "t",
                                                      paging = FALSE,
                                                      ordering = FALSE,
                                                      autoWidth = TRUE,
                                                      columnDefs = colDefs1,
                                                      scrollX = FALSE),
                                       callback = DT::JS("
// Buttons
table.on('click', 'button', function() {
  var id = this.id;
  Shiny.setInputValue('dt1_event', {type: 'button', id: id}, {priority: 'event'});
});
// Checkboxes
table.on('change', 'input[type=\"checkbox\"]', function() {
  var id = this.id;
  var val = this.checked;
  Shiny.setInputValue('dt1_event', {type: 'checkbox', id: id, value: val}, {priority: 'event'});
});
// Dropdowns
table.on('change', 'select', function() {
  var id = this.id;
  var val = this.value;
  Shiny.setInputValue('dt1_event', {type: 'select', id: id, value: val}, {priority: 'event'});
});
"))
  # Reactive functions to update UI\
  updt_ContrTbl <- \(reactive = TRUE) {
    dat <- if (reactive) { CONTRASTSTBL() } else { myContrasts2 }
    return(DT::renderDT({ dat },
                        FALSE,
                        escape = FALSE,
                        class = "compact",
                        selection = "none",
                        rownames = FALSE,
                        editable = FALSE,
                        options = list(dom = "t",
                                       paging = FALSE,
                                       ordering = FALSE,
                                       autoWidth = TRUE,
                                       columnDefs = colDefs2,
                                       scrollX = FALSE),
                        callback = DT::JS("
// Buttons
table.on('click', 'button', function() {
  var id = this.id;
  Shiny.setInputValue('dt2_event', {type: 'button', id: id}, {priority: 'event'});
});
")))
  }
  updt_Msg <- \(reactive = TRUE) {
    msg <- if (reactive) { MSG() } else { "" }
    renderUI(h5(strong(em(msg,
                          style = "color:red",
                          .noWS = "outside"))))
  }
  #
  output$myContrasts <- updt_ContrTbl(FALSE)
  output$msg <- updt_Msg(FALSE)
  #
  # Event observers
  observeEvent(input$dt1_event, {
    info <- input$dt1_event
    if (info$type == "select") {
      formVALS[[info$id]] <- info$value
      if (info$id == "Primary") {
        AB <- dfltContr_Opt[match(info$value, dfltContr_Opt$Contrast), c("A","B")]
        tmp <- c("", dfltContr_Opt$Contrast[(!dfltContr_Opt$A %in% AB)&(!dfltContr_Opt$B %in% AB)])
        session$sendCustomMessage("updateSecondary",
                                  list(choices = tmp))
        formVALS$Secondary <- ""
      }
    }
    if (info$type == "checkbox") {
      formVALS[[info$id]] <- as.logical(info$value)
    }
  })
  observeEvent(input$Primary, {
    #AB <- dfltContr_Opt[match(dfltContr_Opt$Contrast[1L], dfltContr_Opt$Contrast), c("A", "B")]
    AB <- dfltContr_Opt[match(input$Primary, dfltContr_Opt$Contrast), c("A", "B")]
    tmp <- c("", dfltContr_Opt$Contrast[which((!dfltContr_Opt$A %in% AB)&(!dfltContr_Opt$B %in% AB))])
    updateSelectInput(inputId = "Secondary",
                      choices = tmp,
                      selected = "")
  })
  observeEvent(input$dt1_event, {
    info <- input$dt1_event
    if ((info$type == "button") && (info$id == "addContr")) {
      dat <- CONTRASTSTBL()
      nr <- nrow(dat)
      contr <- if (alsoDouble) { fullContrFun(formVALS$Primary, formVALS$Secondary) } else { formVALS$Primary }
      m <- match(contr, dat$Contrast)
      shinyjs::enable("saveBtn")
      if (!is.na(m)) {
        if (formVALS$upOnly != dat$"Up-only"[m]) {
          MSG("")
          output$msg <- updt_Msg()
          dat$"Up-only"[m] <- formVALS$upOnly
          CONTRASTSTBL(dat)
          output$myContrasts <- updt_ContrTbl()
        } else {
          MSG("Contrast already added!")
          output$msg <- updt_Msg()
        }
      } else {
        wTst <- c()
        if (nr) {
          indiv <- unlist(strsplit(gsub("^\\(|\\)$", "", gsub("\\) - \\(", " - ", contr)), " - "))
          indiv2 <- lapply(dat$Contrast, \(x) {
            unlist(strsplit(gsub("^\\(|\\)$", "", gsub("\\) - \\(", " - ", x)), " - "))
          })
          wTst <- which(vapply(indiv2, \(x) {
            (length(x) == length(indiv))&(sum(!indiv %in% x) == 0L)
          }, TRUE))
        }
        if (length(wTst)) {
          MSG("Warning: similar contrast already added (only the order of terms changes)! Do you really want this one included too?")
        } else {
          MSG("")
        }
        output$msg <- updt_Msg()
        i <- nr + 1L
        iChr <- as.character(i)
        tmpDF <- data.frame("Contrast" = contr,
                            "Up-only" = formVALS$upOnly,
                            "Remove" = as.character(actionButton(paste0("rmvBtn_", iChr), "remove contrast")),
                            check.names = FALSE)
        dat <- rbind(dat, tmpDF)
        CONTRASTSTBL(dat)
        output$myContrasts <- updt_ContrTbl()
      }
    }
  })
  observeEvent(input$dt2_event, {
    info <- input$dt2_event
    if ((info$type == "button") && (grepl("^rmvBtn_", info$id))) {
      i <- as.integer(sub("rmvBtn_", "", info$id))
      dat <- CONTRASTSTBL()
      nr <- nrow(dat)
      if (nr) {
        if (i <= nr) {
          dat <- dat[-i, , drop = FALSE]
          if (nrow(dat)) {
            shinyjs::enable("saveBtn")
            dat$Remove <- vapply(seq_len(nrow(dat)), \(j) {
              as.character(actionButton(paste0("rmvBtn_", j), "remove contrast"))
            }, "")
          } else {
            shinyjs::disable("saveBtn")
          }
        } else {
          shinyjs::enable("saveBtn")
        }
        CONTRASTSTBL(dat)
        output$myContrasts <- updt_ContrTbl()
      }
    }
  })
  observeEvent(input$saveBtn, {
    myContrasts2 <- CONTRASTSTBL()[, c("Contrast", "Up-only")]
    assign("myContrasts2", myContrasts2, envir = .GlobalEnv)
    readr::write_rds(myContrasts2, contrastsFl)
    assign("appRunTest", TRUE, envir = .GlobalEnv)
    stopApp()
  })
  observeEvent(input$cancel, { stopApp() })
  session$onSessionEnded(\() { stopApp() })
})
if (exists("appRunTest")) { rm(appRunTest) }
appTxt0 <- sub("myApp", "myApp0", sub("\\(ui", "(ui0", sub(", server", ", server0", runApp)))
runKount <- 0L
while ((!runKount) || (!exists("appRunTest"))) {
  ui0 <- make_ui0() # Update ui with current values
  eval(parse(text = appTxt0), envir = .GlobalEnv)
  shinyCleanup()
  if (file.exists(contrastsFl)) { myContrasts <- readr::read_rds(contrastsFl) }
  runKount <- runKount+1L
}
# Some post-processing
# We used spaces in our app for visibility, but limma doesn't like spaces, so...
myContrasts$Primary <- myContrasts$Contrast
if (alsoDouble) {
  myContrasts$Secondary <- ""
  g <- grep("\\) - \\(", myContrasts$Contrast)
  if (length(g)) {
    myContrasts$Primary[g] <- sub("^\\(", "", sub("\\) - \\(.*", "", myContrasts$Contrast[g]))
    myContrasts$Secondary[g] <- sub("\\)$", "", sub(".*\\) - \\(", "", myContrasts$Contrast[g]))
  }
}
tmp1 <- strsplit(myContrasts$Primary, " - ")
tmp1 <- lapply(tmp1, \(x) {
  contrBlocks2$Name_[match(x, contrBlocks2$Name)]
})
myContrasts$Contrast <- myContrasts$Primary <- vapply(tmp1, paste, "", collapse = " - ")
myContrasts$A <- contrBlocks2$`Samples group`[match(sub(" - .*", "", myContrasts$Primary), contrBlocks2$Name_)]
myContrasts$B <- contrBlocks2$`Samples group`[match(sub(".* - ", "", myContrasts$Primary), contrBlocks2$Name_)]
if (alsoDouble) {
  myContrasts$C <- myContrasts$D <- ""
  tmp2 <- strsplit(myContrasts$Secondary[g], " - ")
  tmp2 <- lapply(tmp2, \(x) {
    contrBlocks2$Name_[match(x, contrBlocks2$Name)]
  })
  myContrasts$Secondary[g] <- vapply(tmp2, paste, "", collapse = " - ")
  myContrasts$Contrast[g] <- paste0("(", myContrasts$Contrast[g], ") - (", myContrasts$Secondary[g], ")")
  myContrasts$C[g] <- contrBlocks2$`Samples group`[match(sub(" - .*", "", myContrasts$Secondary[g]), contrBlocks2$Name_)]
  myContrasts$D[g] <- contrBlocks2$`Samples group`[match(sub(".* - ", "", myContrasts$Secondary[g]), contrBlocks2$Name_)]
}
kol <- c("A", "B", "C", "D")
kol <- intersect(kol, colnames(myContrasts))
for (i in kol) {
  myContrasts[[paste0(i, "_samples")]] <- lapply(myContrasts[[i]], \(x) {
    rownames(expMap)[which(expMap[[VPAL$limmaCol]] == x)]
  })
}
if (!"Secondary" %in% colnames(myContrasts)) {
  myContrasts$Secondary <- ""
}
myContrasts$isDouble <- myContrasts$Secondary != "" # For convenience

# Make limma type designMatr and contrMatr
#
# Design matrix
# -------------
# We need to update limmaForm because we use different names for factors in expMap!
limmaForm <- paste0("~ 0 + ", VPAL$limmaCol)
tmpForm <- VPAL$limmaCol
if ((!is.na(Param$Batch.effect)) && (Param$Batch.effect != "")) {
  # Let's first make a design matrix without batch effect for use by ComBat
  designMatr_noBatch %<o% model.matrix(as.formula(limmaForm)#,
                                       #data = expMap
                                       )
  # Test before we edit column names
  tst <- lapply(tmpForm, \(x) { grep(topattern(x), colnames(designMatr_noBatch)) })
  l <- length(tmpForm)
  if (l > 1L) {
    for (i in 2L:l) {
      stopifnot(sum(tst[[i]] %in% unlist(tst[[1L:(i-1L)]])) == 0L)
    }
  }
  # Edit - using dimnames here to avoid stripping other attributes
  for (i in 1L:l) {
    dimnames(designMatr_noBatch)[[2L]][tst[[i]]] <- sub(topattern(tmpForm[i]), "", dimnames(designMatr_noBatch)[[2L]][tst[[i]]])
  }
  dimnames(designMatr_noBatch)[[2L]] <- gsub("___", "_", dimnames(designMatr_noBatch)[[2L]])
  dimnames(designMatr_noBatch)[[1L]] <- gsub("___", "_", as.character(expMap[[RSA$limmaCol]]))
  #
  # Now we add the batch effect
  limmaForm <- paste0(limmaForm, " + ", Batch.effect$limmaCol)
  tmpForm <- c(tmpForm, Batch.effect$limmaCol)
}
designMatr %<o% model.matrix(as.formula(limmaForm)#,
                             #data = expMap
                             )
# Test before we edit column names
tst <- lapply(tmpForm, \(x) { grep(topattern(x), colnames(designMatr)) })
l <- length(tmpForm)
if (l > 1L) {
  for (i in 2L:l) {
    stopifnot(sum(tst[[i]] %in% unlist(tst[[1L:(i-1L)]])) == 0L)
  }
}
# Edit - using dimnames here to avoid stripping other attributes
for (i in 1L:l) {
  dimnames(designMatr)[[2L]][tst[[i]]] <- sub(topattern(tmpForm[i]), "", dimnames(designMatr)[[2L]][tst[[i]]])
}
dimnames(designMatr)[[2L]] <- gsub("___", "_", dimnames(designMatr)[[2L]])
dimnames(designMatr)[[1L]] <- gsub("___", "_", as.character(expMap[[RSA$limmaCol]]))
expMap$designMatr_rowNames <- dimnames(designMatr)[[1L]]
#
# Contrasts matrix
# ----------------
contrMatr %<o% makeContrasts(contrasts = myContrasts$Contrast,
                             levels = designMatr)

# Reference column
# ----------------
#
# Although we have gotten rid of the Reference-centered approach to stats (to shift to contrasts),
# there are still cases where we need one:
#  - for SAINTexpress (pull-downs)... but this is addressed there, and better done there than here
#  - for Amica
# What we will do then is detect References from the contrasts table.
# If we have several references per comparison group, we will:
#  - aim for one with NA/control Target (bait) protein
#  - failing that, pick the first
# If this doesn't work, consider asking here/later for user input.
if ((!"Reference" %in% colnames(Exp.map)) || (!is.logical(Exp.map$Reference)) || (length(unique(Exp.map$Reference[which(!is.na(Exp.map$Reference))])) != 2L)) {
  allRefs <- setdiff(union(myContrasts$B, myContrasts$D), "")
  Rfs <- rownames(expMap)[which(expMap[[VPAL$limmaCol]] %in% allRefs)]
  kol <- c(VPAL$column, RG$column)
  pullDwnTst <- WorkFlow %in% c("PULLDOWN", "BIOID")
  if (pullDwnTst) {
    kol <- c(kol, "Target")
  }
  em <- Exp.map[match(Rfs, Exp.map[[RSA$column]]), kol]
  nc <- ncol(em)
  nr <- nrow(em)
  colnames(em) <- c("samplesGroup", "compGroup", "Bait")[1L:nc]
  kol2 <- c("samplesGroup", "Bait")[1L:(pullDwnTst+1L)]
  tst <- aggregate(1L:nr, list(em$compGroup), \(x) {
    m <- match(unique(em$samplesGroup[x]), em$samplesGroup)
    return(set_colnames(em[m, kol2, drop = FALSE], kol2))
  })
  tst$samplesGroup <- lapply(1L:nrow(tst), \(i) { tst$x[i]$samplesGroup })
  if (pullDwnTst) {
    tst$Bait <- lapply(1L:nrow(tst), \(i) { tst$x[i]$Bait })
  }
  tst$x <- NULL
  # - Step 1: give priority to NA baits
  if (pullDwnTst) {
    tst$L <- lengths(tst$samplesGroup)
    w <- which(tst$L > 1L)
    if (length(w)) {
      tst$samplesGroup[w] <- lapply(w, \(x) {
        w <- which(is.na(tst$Bait[[x]]))
        x <- tst$samplesGroup[[x]]
        if (length(w)) { x <- x[w] }
        return(x)
      })
    }
  }
  # - Step 2: if that fails, take all
  #tst$samplesGroup <- vapply(tst$samplesGroup, \(x) { x[[1L]] }, "")
  Exp.map$Reference <- FALSE
  for (i in 1L:nrow(tst)) { #i <- 1L
    w1 <- which((Exp.map[[RG$column]] == tst$Group.1[i]))
    w2 <- which((Exp.map[[RG$column]] == tst$Group.1[i])
                &(Exp.map[[VPAL$column]] %in% tst$samplesGroup[[i]]))
    if (length(w2) < length(w1)) {
      rfVal <- opt <- setNames(vapply(cleanNms(VPAL$values), \(x) {
        paste(c(x, rep(" ", 250L-nchar(x))), collapse = "")
      }, ""), VPAL$values)
      while (!length(setdiff(opt, rfVal))) {
        rfVal <- dlg_list(opt, opt[1L], TRUE,
                          title = paste0(tst$Group.1[i], " -> select one or more reference levels"))$res
      }
      Exp.map$Reference <- Exp.map[[VPAL$column]] %in% names(opt)[match(rfVal, opt)]
    } else {
      Exp.map$Reference[w2] <- TRUE
    }
  }
}
