# Define Contrasts of interest

# expMap: limma-compatible version of apap
hasBlocks %<o% FALSE
origCoeff <- union(RSA$names, RRG$names)
if ((!is.na(Param$Blocking.factors))&&(Param$Blocking.factors != "")) {
  parse.Param.aggreg.2("Blocking.factors")
  origCoeff <- union(origCoeff, Blocking.factors$names)
  hasBlocks <- TRUE
}
if ((!is.na(Param$Batch.effect))&&(Param$Batch.effect != "")) {
  parse.Param.aggreg.2("Batch.effect")
  origCoeff <- union(origCoeff, Batch.effect$names)
}
Coefficients %<o% origCoeff
expMap %<o% Exp.map[, origCoeff]
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
for (Coeff in origCoeff) { #Coeff <- "Replicate"
  assign(Coeff, factor(expMap[[Coeff]], FactorsLevels[[Coeff]]))
  l <- length(grep("-", Exp.map[[Coeff]])) # Not expMap in case we are re-running a small chunk
  if (l) {
    nuCoeff <- paste0(Coeff, "_._")
    stopifnot(!nuCoeff %in% colnames(Exp.map)) # Not expMap in case we are re-running a small chunk
    Coefficients[match(Coeff, origCoeff)] <- nuCoeff
    expMap[[nuCoeff]] <- gsub("-", ".", expMap[[Coeff]])
    expMap[[Coeff]] <- NULL
    assign(nuCoeff, factor(expMap[[nuCoeff]], gsub("-", ".", FactorsLevels[[Coeff]])))
  }
  w <- which(suppressWarnings(as.character(as.numeric(as.character(Exp.map[[Coeff]]))) == as.character(Exp.map[[Coeff]])))
  if (length(w)) {
    nuCoeff <- paste0(Coeff, "_._")
    stopifnot(!nuCoeff %in% colnames(Exp.map)) # Not expMap in case we are re-running a small chunk
    Coefficients[match(Coeff, origCoeff)] <- nuCoeff
    expMap[[nuCoeff]] <- expMap[[Coeff]]
    expMap[w, nuCoeff] <- paste0(Coeff, "_", as.character(expMap[w, Coeff]))
    expMap[[Coeff]] <- NULL
    assign(nuCoeff, factor(expMap[[nuCoeff]], paste0(Coeff, "_", as.character(FactorsLevels[[Coeff]]))))
  }
}
#
for (aggrNm in c("RSA", "VPAL", "RG", "RRG", "Blocking.factors", "Batch.effect")) { # Update in case we had to do some fixing
  if (exists(aggrNm)) {
    aggr <- get(aggrNm)
    tmp <- aggr$names
    tmp <- Coefficients[match(tmp, origCoeff)]
    nm <- paste0(paste(tmp, collapse = "_"), "_._")
    aggr$limmaCol <- nm
    expMap[[nm]] <- do.call(paste, c(expMap[, tmp, drop = FALSE], sep = "___"))
    expMap[[nm]] <- factor(expMap[[nm]], levels = unique(expMap[[nm]]))
    nm %<c% as.factor(expMap[[nm]])
    if (length(Exp) == 1) {
      tmp <- setdiff(tmp, "Experiment")
      nm <- paste0(paste(tmp, collapse = "_"), "_._")
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
contrBlocks <- setNames(lapply(ratGrps, function(grp) {
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
dfltContr_Opt <- setNames(lapply(ratGrps, function(grp) { #grp <- ratGrps[1]
  x <- as.data.frame(gtools::#combinations
                       permutations(2, length(contrBlocks[[grp]]), names(contrBlocks[[grp]])))
  colnames(x) <- c("A", "B")
  do.call(paste, c(x[, c("A", "B")], sep = " - "))
}), ratGrps)
alsoDouble <- length(ratGrps > 1)||(length(dfltContr_Opt[[1]]) > 3) # We need 4 different things for an interaction double contrast
dfltContr_Opt <- stack(dfltContr_Opt)
colnames(dfltContr_Opt) <- c("Contrast", "Comparison group")
dfltContr_Opt[, c("A", "B")] <- do.call(rbind, strsplit(dfltContr_Opt$Contrast, " - "))
#
fullContrFun <- function(prim, sec) {
  if ((missing(sec))||(is.na(sec))||(length(sec) != 1)||(sec == "")) { return(prim) }
  return(paste0("(", prim, ") - (", sec, ")"))
}
contrastsFl %<o% paste0(wd, "/Contrasts.rds")
if (file.exists(contrastsFl)) {
  myContrasts <- readr::read_rds(contrastsFl)
  g <- grep("\\) - \\(", myContrasts$Contrast)
  if ((length(g))&&(!alsoDouble)) {
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
  vapply(1:nr2, function(i) {
    iChr <- as.character(i)
    as.character(actionButton(paste0("rmvBtn_", iChr), "remove contrast"))
  }, "")
} else { character() }
AB <- dfltContr_Opt[match(dfltContr_Opt$Contrast[1], dfltContr_Opt$Contrast), c("A","B")]
tmp <- c("", dfltContr_Opt$Contrast[(!dfltContr_Opt$A %in% AB)&(!dfltContr_Opt$B %in% AB)])
makeContr <- data.frame("Contrast" = as.character(selectInput("Primary", "", dfltContr_Opt$Contrast)),
                        "(opt. secondary contrast)" = as.character(selectInput("Secondary", "", tmp)),
                        "Up-regulated only?" = as.character(checkboxInput("upOnly", "", FALSE)),
                        "add contrast" = as.character(actionButton("addContr", "add contrast")),
                        check.names = FALSE)
if (!alsoDouble) { makeContr$"(opt. secondary contrast)" <- NULL }
colDefs1 <- list(list(width = "250px", targets = 0),
                 if (alsoDouble) { list(width = "250px", targets = 1) },
                 list(width = "100px", targets = 2),
                 list(width = "50px", targets = 3))
colDefs2 <- list(list(width = "500px", targets = 0),
                if (alsoDouble) { list(width = "100px", targets = 1) },
                list(width = "100px", targets = 2))
appNm <- "Define contrasts"
make_ui0 <- function() {
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
        em("Tick \"one-sided\" for \"up-regulated\" (for \"down-regulated\", just select the reverse contrast first)."),
        br(),
        if (alsoDouble) {
          em("You may also add double contrasts (interaction contrasts) of the form \"(A - B) - (C - D)\" (where all 4 are distinct).")
        },
        if (alsoDouble) {
          em("(for now we are restricting users to \"sensible contrasts\" - this may change)")
        },
        if (alsoDouble) {
          br()
        },
        em(HTML("Once you are finished, click&nbsp")),
        shinyWidgets::actionBttn("saveBtn", "Save", icon = icon("save"), color = "success", style = "pill"),
        em(HTML("&nbspto continue.")),
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
server0 <- shinyServer(function(input, output, session) {
  CONTRASTSTBL <- reactiveVal(myContrasts2)
  MSG <- reactiveVal("")
  formVALS <- reactiveValues(Primary = dfltContr_Opt$Contrast[1],
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
                                                      autowidth = TRUE,
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
  updt_ContrTbl <- function(reactive = TRUE) {
    if (reactive) { dat <- CONTRASTSTBL() } else { dat <- myContrasts2 }
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
                                       autowidth = TRUE,
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
  updt_Msg <- function(reactive = TRUE) {
    if (reactive) { msg <- MSG() } else { msg <- "" }
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
    #AB <- dfltContr_Opt[match(dfltContr_Opt$Contrast[1], dfltContr_Opt$Contrast), c("A", "B")]
    AB <- dfltContr_Opt[match(input$Primary, dfltContr_Opt$Contrast), c("A", "B")]
    tmp <- c("", dfltContr_Opt$Contrast[which((!dfltContr_Opt$A %in% AB)&(!dfltContr_Opt$B %in% AB))])
    updateSelectInput(inputId = "Secondary",
                      choices = tmp,
                      selected = "")
  })
  observeEvent(input$dt1_event, {
    info <- input$dt1_event
    if ((info$type == "button")&&(info$id == "addContr")) {
      dat <- CONTRASTSTBL()
      nr <- nrow(dat)
      contr <- if (alsoDouble) { fullContrFun(formVALS$Primary, formVALS$Secondary) } else { formVALS$Primary }
      m <- match(contr, dat$Contrast)
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
          indiv2 <- lapply(dat$Contrast, function(x) {
            unlist(strsplit(gsub("^\\(|\\)$", "", gsub("\\) - \\(", " - ", x)), " - "))
          })
          wTst <- which(vapply(indiv2, function(x) {
            (length(x) == length(indiv))&(sum(!indiv %in% x) == 0)
          }, TRUE))
        }
        if (length(wTst)) {
          MSG("Warning: similar contrast already added (only the order of terms changes)! Do you really want this one included too?")
        } else {
          MSG("")
        }
        output$msg <- updt_Msg()
        i <- nr + 1
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
    if ((info$type == "button")&&(grepl("^rmvBtn_", info$id))) {
      i <- as.integer(sub("rmvBtn_", "", info$id))
      dat <- CONTRASTSTBL()
      nr <- nrow(dat)
      if (nr &&(i <= nr)) {
        dat <- dat[-i, , drop = FALSE]
        if (nrow(dat)) {
          shinyjs::enable("saveBtn")
          dat$Remove <- vapply(seq_len(nrow(dat)), function(j) {
            as.character(actionButton(paste0("rmvBtn_", j), "remove contrast"))
          }, "")
        } else {
          shinyjs::disable("saveBtn")
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
  session$onSessionEnded(function() { stopApp() })
})
if (exists("appRunTest")) { rm(appRunTest) }
appTxt0 <- sub("myApp", "myApp0", sub("\\(ui", "(ui0", sub(", server", ", server0", runApp)))
runKount <- 0
while ((!runKount)||(!exists("appRunTest"))) {
  ui0 <- make_ui0() # Update ui with current values
  eval(parse(text = appTxt0), envir = .GlobalEnv)
  myContrasts <- readr::read_rds(contrastsFl)
  shinyCleanup()
  runKount <- runKount+1
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
tmp1 <- lapply(tmp1, function(x) {
  contrBlocks2$Name_[match(x, contrBlocks2$Name)]
})
myContrasts$Contrast <- myContrasts$Primary <- vapply(tmp1, paste, "", collapse = " - ")
myContrasts$A <- contrBlocks2$`Samples group`[match(sub(" - .*", "", myContrasts$Primary), contrBlocks2$Name_)]
myContrasts$B <- contrBlocks2$`Samples group`[match(sub(".* - ", "", myContrasts$Primary), contrBlocks2$Name_)]
if (alsoDouble) {
  myContrasts$C <- myContrasts$D <- ""
  tmp2 <- strsplit(myContrasts$Secondary[g], " - ")
  tmp2 <- lapply(tmp2, function(x) {
    contrBlocks2$Name_[match(x, contrBlocks2$Name)]
  })
  myContrasts$Secondary[g] <- vapply(tmp2, paste, "", collapse = " - ")
  myContrasts$Contrast[g] <- paste0("(", myContrasts$Contrast[g], ") - (", myContrasts$Secondary[g], ")")
  myContrasts$C[g] <- contrBlocks2$`Samples group`[match(sub(" - .*", "", myContrasts$Secondary[g]), contrBlocks2$Name_)]
  myContrasts$D[g] <- contrBlocks2$`Samples group`[match(sub(".* - ", "", myContrasts$Secondary[g]), contrBlocks2$Name_)]
}
kol <- VPAL$column
if (length(Exp) == 1) { kol <- sub("Exp", "", kol) }
if (nchar(kol) == 3) { kol <- Factors[kol] }
kol2 <- c("A", "B", "C", "D")
kol2 <- intersect(kol2, colnames(myContrasts))
for (i in kol2) {
  myContrasts[[paste0(i, "_samples")]] <- lapply(myContrasts[[i]], function(x) {
    Exp.map$Ref.Sample.Aggregate[which(Exp.map[[kol]] == x)]
  })
}


# Make limma type designMatr and contrMatr
#
# Design matrix
tmpForm <- unlist(strsplit(limmaForm, "\\+"))
tmpForm <- tmpForm[2:length(tmpForm)]
tmpForm <- paste0(tmpForm, "_._")
tmpForm2 <- paste0("~0+", paste(tmpForm, collapse = "+"))
designMatr %<o% model.matrix(as.formula(tmpForm2))
# Test before we edit column names
tst <- lapply(tmpForm, function(x) {
  grep(topattern(x), colnames(designMatr))
})
l <- length(tmpForm)
if (l > 1) {
  for (i in 2:l) {
    stopifnot(sum(tst[[i]] %in% unlist(tst[[1:(i-1)]])) == 0)
  }
}
# Edit
for (i in 1:l) {
  colnames(designMatr)[tst[[i]]] <- sub(topattern(tmpForm[i]), "", colnames(designMatr)[tst[[i]]])
}
colnames(designMatr) <- gsub("___", "_", colnames(designMatr))
rownames(designMatr) <- gsub("___", "_", as.character(expMap[[RSA$limmaCol]]))
#
# Contrasts matrix
contrMatr %<o% makeContrasts(contrasts = myContrasts$Contrast, levels = designMatr)
