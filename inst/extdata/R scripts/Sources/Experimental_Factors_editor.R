# Create Experimental Factors
if (exists("tstXpFct")) { rm(tstXpFct) }
minFact <- c("Experiment", "Replicate")
minFactDesc <- setNames(c("grouping of sample groups to compare", "maximum number of replicates per sample group"),
                        minFact)
if (file.exists("Factors.RData")) {
  load("Factors.RData")
  Factors <- setNames(tmp$Factors, substr(tmp$Factors, 1, 3))
  Factors <- Factors[which(!is.na(Factors))]
  if (length(Factors)) { FactorsLevels %<o% tmp$Levels[Factors] } else { rm(Factors) }
}
if (WorkFlow == "PULLDOWN") {
  minFact <- c(minFact, "Target")
  minFactDesc["Target"] <- "ID (e.g. UniProtKB accession) in the search database(s) of the bait protein"
}
if (WorkFlow == "TIMECOURSE") {
  minFact <- c(minFact, "Time.point")
  minFactDesc["Time.point"] <- "number without unit"
}
if (WorkFlow == "LOCALISATION") {
  minFact <- c(minFact, "Compartment")
  minFactDesc["Compartment"] <- "compartment/subcellular fractionation"
}
if (LabelType == "Isobaric") {
  minFact <- c(minFact, "Isobaric.set")
  minFactDesc["Isobaric.set"] <- "grouping of individual isobarically-labelled samples pooled together (integer)"
}
if (!exists("Factors")) {
  Factors <- minFact
  FactorsLevels %<o% setNames(lapply(Factors, function(x) { c("") }), Factors)
} else { Factors %<o% unique(c(Factors, minFact)) }
if (WorkFlow == "BIOID") {
  Factors <- unique(c(Factors, "Target")) # In this case target is not obligatory
  minFactDesc["Target"] <- "ID (e.g. UniProtKB accession) in the search database(s) of the bait protein"
}
if (!exists("FactorsLevels")) {
  FactorsLevels %<o% setNames(lapply(Factors, function(x) { c("") }), Factors)
  if (LabelType == "Isobaric") {
    FactorsLevels["Isobaric.set"] <- sort(unique(FracMap$Isobaric.set))
  }
}
w <- which(!Factors %in% names(FactorsLevels))
if (length(w)) {
  FactorsLevels[Factors[w]] <- c()
}
Factors %<o% Factors[which(!is.na(Factors))]
if (LabelType == "Isobaric") {
  for (Fct in Factors[which(!Factors %in% c("Experiment", "Replicate", "Isobaric.set"))]) {
    FactorsLevels[[Fct]] <- unique(c(FactorsLevels[[Fct]], "Mixed_IRS"))
  }
}
FactorsLevels %<o% FactorsLevels[Factors]
#rm(Factors, FactorsLevels)
# Do not use my usual meta-coding approach!!! It is flexible... but not reactive, and almost as bad as shiny to debug!
appNm <- paste0(dtstNm, " - Exp. Factors")
if (exists("runTst")) { rm(runTst) }
ui <- fluidPage(
  useShinyjs(),
  setBackgroundColor( # Doesn't work
     color = c(#"#F8F8FF",
       "#F2F0FA"),
     gradient = "linear",
     direction = "bottom"
  ),
  extendShinyjs(text = jsToggleFS, functions = c("toggleFullScreen")),
  titlePanel(tag("u", "Experimental Factors editor"),
             appNm),
  h2(dtstNm), 
  h3("Enter the names and levels of all the Factors describing the samples."),
  br(),
  h4("The following Factors are included by defaults and cannot be removed:"),
  uiOutput("minFact"),
  br(),
  sidebarLayout(
    sidebarPanel(
      textInput("nuFact", "Experimental Factor name(s)", ""),
      actionBttn("addFactor", "Create Factor(s)", color = "primary", size = "xs", style = "pill"),
      actionBttn("rmvFactor", "Remove Factor(s)", color = "primary", size = "xs", style = "pill"),
      h5(em("Factor names must be:")),
      h5(em(" - at least 3 characters long")),
      h5(em(" - start with a capital letter,")),
      h5(em(" - followed by only lower case letters, numbers or dots.")),
      br(),
      br(),
    ),
    mainPanel(
      h5(em("Add Factor levels:")),
      h5(em(" - use underscores, no spaces or hyphens within level names; sequences of 3 or more consecutive underscores are forbidden!\"\"")),
      h5(em(" - use spaces to add multiple levels at a time")),
      br(),
      span(uiOutput("Message"), style = "color:red"),
      h3("Factor levels:"),
      uiOutput("Factors"),
      br(),
      br(),
      actionBttn("saveBtn", "Save", icon = icon("save"), color = "success", style = "pill"),
    )
  ),
  br(),
  br()
)
server <- function(input, output, session) {
  # Initialize reactive variables
  FACT <- reactiveVal(Factors)
  FACTLevels <- reactiveVal(FactorsLevels)
  intFact <- c("Replicate", "Isobaric.set")
  dfltInt <- c(2, 1)
  # Create function to update output$Factors for UI
  updtFactUI <- function(reactive = TRUE) {
    if (reactive) {
      FAKT <- FACT()
      FAKTLevels <- FACTLevels()
    } else {
      FAKT <- Factors
      FAKTLevels <- FactorsLevels
    }
    L <- length(FAKT)
    # Update UI
    return(renderUI({
      lst <- vector("list", L)
      for (i in 1:L) {
        Fact <- FAKT[i]
        if (Fact %in% intFact) {
          miN <- dfltInt[match(Fact, intFact)]
          dflt <- max(c(miN, FAKTLevels[[Fact]]))
          txt <- as.character(dflt)
          if (!length(txt)) { txt <- "" }
          lst[[i]] <- list(fluidRow(numericInput(paste0(Fact, "_lev"),
                                                 paste0("N. of ", gsub("\\.", " ", Fact), "s = ", txt), dflt, miN, step = 1)))
        }
        if (Fact == "Time.point") {
          lst[[i]] <- list(fluidRow(textInput(paste0(Fact, "_lev"),
                                              paste0(Fact, ", levels (numeric(s) only!) = ", paste(FAKTLevels[[Fact]], collapse = " / ")), "")))
        }
        if (!Fact %in% c(intFact, "Time.point")) {
          lst[[i]] <- list(fluidRow(textInput(paste0(Fact, "_lev"),
                                              paste0(Fact, ", levels = ", paste(FAKTLevels[[Fact]], collapse = " / ")), "")))
        }
        lst[[i]] <- append(lst[[i]],
                           list(fluidRow(column(3, actionBttn(paste0(Fact, "_levAdd"), "Add level(s)", color = "primary", size = "xs", style = "pill")),
                                         column(3, actionBttn(paste0(Fact, "_levRmv"), "Remove level(s)", color = "primary", size = "xs", style = "pill"))),
                                tags$hr(style = "border-color: grey;"),
                                br()))
      }
      return(list(lst))
    }))
  }
  #
  output$minFact <- renderUI({ HTML(paste0(" - ", minFact, ": ", minFactDesc, collapse = "<br>")) })
  output$Message <- renderUI({ em(" ") })
  # Initialize
  output$Factors <- updtFactUI(reactive = FALSE)
  #
  # Observers for already extent factors
  #  - Add new Factor level
  sapply(Factors, function(Fact) {
    observeEvent(input[[paste0(Fact, "_levAdd")]], {
      vals <- input[[paste0(Fact, "_lev")]]
      vals <- vals[which(!is.na(vals))]
      if (length(vals)) {
        if ((is.character(vals))&&(Fact != "Target")) { vals <- gsub("-", ".", unlist(strsplit(vals, " "))) }
        tmp2 <- FACTLevels()
        if (Fact %in% intFact) {
          tmp2[[Fact]] <- 1:as.integer(max(c(vals, dfltInt[match(Fact, intFact)])))
        } else {
          if (Fact == "Time.point") {
            vals <- suppressWarnings(as.numeric(vals))
            vals <- vals[which(!is.na(vals))]
          }
          tmp <- unique(c(FACTLevels()[[Fact]], vals))
          tmp <- tmp[which(tmp != "")]
          tmp2[[Fact]] <- tmp
        }
        FACTLevels(tmp2)
        output$Factors <- updtFactUI()
      }
    })
  })
  #  - Remove Factor level
  sapply(Factors, function(Fact) {
    observeEvent(input[[paste0(Fact, "_levRmv")]], {
      vals <- input[[paste0(Fact, "_lev")]]
      vals <- vals[which(!is.na(vals))]
      if (length(vals)) {
        if ((is.character(vals))&&(Fact != "Target")) { vals <- gsub("-", ".", unlist(strsplit(vals, " "))) }
        tmp2 <- FACTLevels()
        if (Fact %in% intFact) {
          tmp2[[Fact]] <- 1:(max(c(as.integer(vals)-1, dfltInt[match(Fact, intFact)])))
        } else {
          if (Fact == "Time.point") {
            vals <- suppressWarnings(as.numeric(vals))
            vals <- vals[which(!is.na(vals))]
          }
          tmp2[[Fact]] <- tmp2[[Fact]][which(!tmp2[[Fact]] %in% vals)]
        }
        FACTLevels(tmp2)
        output$Factors <- updtFactUI()
      }
    })
  })
  #
  # Create new Factor(s)
  observeEvent(input$addFactor, {
    msg <- " "
    Facts <- gsub("[^A-Z,a-z,0-9]", "\\.", unlist(strsplit(input$nuFact, " +")))
    sapply(Facts, function(Fact) {
      if ((nchar(Fact) < 3)||(grepl("^[0-9]", Fact))||(substr(Fact, 1, 3) %in% substr(FACT(), 1, 3))) {
        msg <- "Invalid Factor name! Must be at least 3 characters long and start with a capital letter! The first 3 characters must be unique to this Factor!"
      } else {
        Fact <- paste0(toupper(substr(Fact, 1, 1)), tolower(substr(Fact, 2, nchar(Fact))))
        if (!Fact %in% FACT()) {
          FACT(c(FACT(), Fact))
          tmp <- FACTLevels()
          tmp[Fact] <- list(NULL)
          FACTLevels(tmp)
          # Also, ABSOLUTELY crucial: create new observers for level addition/removal!!!
          observeEvent(input[[paste0(Fact, "_levAdd")]], {
            vals <- input[[paste0(Fact, "_lev")]]
            vals <- vals[which(!is.na(vals))]
            if (length(vals)) {
              if ((is.character(vals))&&(Fact != "Target")) { vals <- gsub("-", ".", unlist(strsplit(vals, " "))) }
              tmp2 <- FACTLevels()
              if (Fact %in% intFact) {
                tmp2[[Fact]] <- 1:as.integer(max(c(vals, dfltInt[match(Fact, intFact)])))
              } else {
                if (Fact == "Time.point") {
                  vals <- suppressWarnings(as.numeric(vals))
                  vals <- vals[which(!is.na(vals))]
                }
                tmp <- unique(c(FACTLevels()[[Fact]], vals))
                tmp <- tmp[which(tmp != "")]
                tmp2[[Fact]] <- tmp
              }
              FACTLevels(tmp2)
              output$Factors <- updtFactUI()
            }
          })
          observeEvent(input[[paste0(Fact, "_levRmv")]], {
            vals <- input[[paste0(Fact, "_lev")]]
            vals <- vals[which(!is.na(vals))]
            if (length(vals)) {
              if ((is.character(vals))&&(Fact != "Target")) { vals <- gsub("-", ".", unlist(strsplit(vals, " "))) }
              tmp2 <- FACTLevels()
              if (Fact %in% intFact) {
                tmp2[[Fact]] <- 1:(max(c(as.integer(vals)-1, dfltInt[match(Fact, intFact)])))
              } else {
                if (Fact == "Time.point") {
                  vals <- suppressWarnings(as.numeric(vals))
                  vals <- vals[which(!is.na(vals))]
                }
                tmp <- unique(c(FACTLevels()[[Fact]], vals))
                tmp2[[Fact]] <- tmp2[[Fact]][which(!tmp2[[Fact]] %in% vals)]
              }
              FACTLevels(tmp2)
              output$Factors <- updtFactUI()
            }
          })
        }
      }
    })
    output$Message <- renderUI({ em(msg) })
    output$Factors <- updtFactUI()
  })
  # Remove extant Factor
  observeEvent(input$rmvFactor, {
    msg <- " "
    tmp <- gsub("[^A-Z,a-z,0-9]", "\\.", input$nuFact)
    if ((nchar(tmp) < 3)||(grepl("^[0-9]", tmp))) {
      msg <- "Invalid Factor name! Must be at least 3 characters long and start with a capital letter!"
    } else {
      tmp <- paste0(toupper(substr(tmp, 1, 1)), tolower(substr(tmp, 2, nchar(tmp))))
      if (!tmp %in% FACT()) {
        msg <- "Cannot remove a non-existent Factor!"
      } else {
        if (tmp %in% minFact) {
          msg <- paste0("Factor ", tmp, " is included by default in this workflow and cannot be removed!")
        } else {
          tmp2 <- FACT()
          tmp2 <- tmp2[which(tmp2 != tmp)]
          FACT(tmp2)
          FACTLevels(FACTLevels()[FACT()])
        }
      }
    }
    output$Message <- renderUI({ em(msg) })
    output$Factors <- updtFactUI()
  })
  #
  observeEvent(input$saveBtn, {
    assign("Factors", FACT(), envir = .GlobalEnv)
    assign("FactorsLevels", FACTLevels(), envir = .GlobalEnv)
    assign("runTst", TRUE, envir = .GlobalEnv)
    stopApp()
  })
  #observeEvent(input$cancel, { stopApp() })
  session$onSessionEnded(function() { stopApp() })
}
runKount <- 0
while ((!runKount)||(!exists("runTst"))) {
  eval(parse(text = runApp), envir = .GlobalEnv)
  runKount <- runKount+1
}
#
if ("Target" %in% Factors) {
  FactorsLevels$Target <- FactorsLevels$Target[which(!is.na(FactorsLevels$Target))]
  if (length(FactorsLevels$Target) == 1) {
    FactorsLevels$Target <- c(FactorsLevels$Target, "Control")
    if (length(FactorsLevels$Target) == 1) {
      FactorsLevels$Target <- c(FactorsLevels$Target, "Ctrl")
    }
  }
}
FactorsLevels <- setNames(lapply(Factors, function(fct) {
  x <- FactorsLevels[[fct]]
  x[which(!is.na(x))]
}), Factors)
Factors <- Factors[which(sapply(FactorsLevels[Factors], length) > 0)]
names(Factors) <- substr(Factors, 1, 3)
Factors <- Factors[c("Exp", names(Factors)[which(!names(Factors) %in% c("Exp", "Rep"))], "Rep")]
FactorsLevels <- FactorsLevels[Factors]
tmp <- list(Factors = Factors, Levels = FactorsLevels)
save(tmp, file = "Factors.RData")

tstXpFct <- TRUE
