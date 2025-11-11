# Statistical tests
#if (!Param$Param_suppress_UI) {
#dlg_message("Remember: replace this app with a \"Statistics: contrasts editor\" app!", "ok")
#
expMap <- Exp.map # temporary aliases
expMap$Row <- 1:nrow(Exp.map)
vpal <- unlist(strsplit(Param$Volcano.plots.Aggregate.Level, ";"))
#
# We should select a reference level per ratio group, not ratio reference group!
#rrg <- unlist(strsplit(Param$Ratios.Ref.Groups, ";"))
rg <- unlist(strsplit(Param$Ratios.Groups, ";"))
if (length(unique(expMap$Experiment)) == 1) {
  if (length(vpal) > 1)  { vpal <- vpal[which(vpal != "Exp")] }
  #if (length(rrg) > 1)  { rrg <- rrg[which(rrg != "Exp")] }
  if (length(rg) > 1)  { rg <- rg[which(rg != "Exp")] }
}
kolVPAL <- Factors[vpal]
#kolRRG <- Factors[rrg]
kolRG <- Factors[rg]
factLevComb1 <- do.call(paste, c(expMap[, kolVPAL, drop = FALSE], sep = " ")) # Sample group levels
factLevComb2 <- do.call(paste, c(expMap[, kolRG, drop = FALSE], sep = " ")) # Ratio group levels
#factLevComb3 <- do.call(paste, c(expMap[, kolRRG, drop = FALSE], sep = " ")) # Ratio reference group levels - currently not used
refTst <- ("Reference" %in% colnames(expMap))&&("logical" %in% class(expMap$Reference))&&
  (sum(expMap$Reference))&&(sum(!expMap$Reference))
rfLev <- data.frame("Group" = unique(factLevComb2))
rfLev$"All levels" <- lapply(rfLev$Group, function(x) {
  unique(factLevComb1[which(factLevComb2 == x)])
})
# Remove groups with one sample group - no comparison possible!
tstL <- vapply(rfLev$"All levels", length, 1)
w1 <- which(tstL == 1)
l1 <- length(w1)
if (l1) {
  tmp <- rfLev$Group[w1]
  t1 <- l1 > 1
  if (t1) { tmp <- paste0(paste(rfLev$Group[w1[1:(l1-1)]], collapse = ", "), " and ", rfLev$Group[w1[l1]]) }
  msg <- paste0("Ratio group", c("", "s")[(l1 > 1)+1], " ", tmp, " contain", c("", "s")[t1+1], " only one samples group, skipping!")
  warning(msg)
  lev2Rmv <- unique(unlist(rfLev$`All levels`[w1]))
  w <- which(!factLevComb1 %in% lev2Rmv)
  expMap <- expMap[w,] # We want to keep those in the original, Exp.map, e.g. for mixed channels for TMT IRS normalization!!!
  factLevComb1 <- factLevComb1[w]
  factLevComb2 <- factLevComb2[w] # Ratio group levels
  #factLevComb3 <- factLevComb3[w] # Ratio reference group levels - currently not used
  rfLev <- rfLev[-w1,]
}
rfLev$"Reference level" <- vapply(rfLev$"All levels", function(x) { #x <- unlist(rfLev$"All levels"[2])
  if (refTst) {
    w <- which(expMap$Reference[match(x, factLevComb1)])
    if (length(w)) { w <- w[1] } else { w <- 1 }
    rs <- x[w]
  } else {
    rs <- factLevComb1[w[1]]
  }
  return(rs)
}, "")
makeSec <- FALSE
# if (Param$F.test) {
#   Factors3 <- Factors2[which(vapply(Factors2, function(Fact) {
#     length(FactorsLevels[[Fact]]) > 1
#   }, TRUE))]
#   if (WorkFlow == "TIMECOURSE") { Factors3 <- unique(c(Factors3, "Time.point")) }
#   makeSec <- length(Factors3) > 1
# }
wdth <- paste0(30*max(c(nchar(unlist(rfLev$`All levels`)), 2)), "px")
appNm <- paste0(dtstNm, " - Stat-tests")
ui2 <- fluidPage(
  useShinyjs(),
  setBackgroundColor( # Doesn't work
    color = c(#"#F8F8FF",
      "#ECE6F7"),
    gradient = "linear",
    direction = "bottom"
  ),
  extendShinyjs(text = jsToggleFS, functions = c("toggleFullScreen")),
  titlePanel(tag("u", "Statistical test - define comparisons"), 
             appNm),
  h2(dtstNm), 
  br(),
  actionBttn("saveBtn", "Save", icon = icon("save"), color = "success", style = "pill"),
  br(),
  fluidPage(
    h3(" -> T-tests: define reference factors combination per comparisons group."),
    h5(em("Select a single reference level sample group per comparisons group.")),
    h5(em("Samples from this samples group will serve as the denominator of ratios and the controls for all t-tests performed within the comparisons group.")),
    br(),
    #withSpinner(
    DTOutput("refLevels")#)
    ,
    br(),
    br()
  ),
  # if (Param$F.test) {
  #   sidebarLayout(
  #     sidebarPanel(
  #       h3(" -> F-tests: design analyses"),
  #       h5(em("For F-tests, you can define multiple analyses, each with simple and/or double contrasts to different references.")),
  #       h5(em("For each analysis, select reference levels for each Factor relevant to the desired contrasts (comparisons) then add Analysis.")),
  #       h5(em("N.B.: Only combinations of levels corresponding to at least 1 valid sample group are allowed.")),
  #       br(),
  #       actionButton("addAnalysis", "Add new analysis"),
  #       withSpinner(uiOutput("NewAnalysis")),
  #       span(uiOutput("Msg"), style = "color:red"),
  #     ),
  #     mainPanel(
  #       h3("F-test analyses"),
  #       br(),
  #       uiOutput("Analyses")
  #     )
  #   )
  # },
  br()
)
#
NuAn_tmplt <- data.frame(Analysis = "", Primary = "", Primary_Ref = "")
if (makeSec) { NuAn_tmplt$Secondary <- NuAn_tmplt$Secondary_Ref <- "" }
rfLev2 <- rfLev[, c("Group", "Reference level")]
rfLev2$"Reference level" <- sapply(1:nrow(rfLev2), function(i) {
  as.character(selectInput(paste0("Ref_", as.character(i)),
                           "",
                           rfLev$`All levels`[[i]],
                           rfLev$`Reference level`[i],
                           selectize = FALSE,
                           width = wdth))
})
if (exists("appRunTest")) { rm(appRunTest) }
server2 <- function(input, output, session) {
  appRunTest <- TRUE
  output$refLevels <- renderDT({ rfLev2 },
                               FALSE,
                               escape = FALSE,
                               class = "compact",
                               selection = "none",
                               editable = TRUE,
                               rownames = FALSE,
                               options = list(dom = 't',
                                              paging = FALSE,
                                              ordering = FALSE
                               ),
                               callback = JS("table.rows().every(function(i, tab, row) {
        var $this = $(this.node());
        $this.attr('id', this.data()[0]);
        $this.addClass('shiny-input-container');
      });
      Shiny.unbindAll(table.table().node());
      Shiny.bindAll(table.table().node());"))
  # sapply(1:nrow(rfLev2), function(x) {
  #   id <- paste0("Ref_", as.character(x))
  #   observeEvent(input[[id]],
  #                {
  #                  rfLev$`Reference level`[x] <- input[[id]]
  #                  assign("rfLev", rfLev, envir = .GlobalEnv)
  #                })
  # })
  # if (Param$F.test) { # Initialize variables to create in main environment
  #   #
  #   # Parse F-test Parameters
  #   fFact <- fRef <- FALSE
  #   #Fact <- Factors2[which(!Factors2 %in% c("Experiment", "Replicate"))]
  #   #w <- which(vapply(Fact, function(x) { length(FactorsLevels[[x]]) > 1 }, TRUE))
  #   #Fact[w]
  #   #tmp <- setNames(substr(Fact, 1, 3), Fact)
  #   #Param$F.test_factors <- paste0(paste(tmp[1:2], collapse = "___"), "_;_", tmp[3])
  #   #Param$F.test_factors_ref <- paste0(paste(sapply(names(tmp)[1:2], function(x) { FactorsLevels[[x]][1] }), collapse = "___"), "_;_", FactorsLevels[[names(tmp)[3]]][1])
  #   if ("F.test_factors" %in% colnames(Param)) {
  #     F_Fact <- unique(unlist(strsplit(Param$F.test_factors, "_\\|_")))
  #     F_Fact <- F_Fact[which(F_Fact != "")]
  #     if (length(F_Fact)) {
  #       fFact <- TRUE
  #     }
  #   }
  #   if ("F.test_factors_ref" %in% colnames(Param)) {
  #     F_Ref <- unique(unlist(strsplit(Param$F.test_factors_ref, "_\\|_")))
  #     F_Ref <- F_Ref[which(F_Ref != "")]
  #     if (length(F_Ref)) {
  #       fRef <- TRUE
  #     }
  #   }
  #   if (fFact && fRef) {
  #     F_Fact <- strsplit(F_Fact, "_;_")
  #     stopifnot(max(vapply(F_Fact, length, 1)) <= 3) # We can have at most 3 factors aggregates: factors for primary and secondary contrasts, and blocking factors for nested designs.
  #     F_Fact <- lapply(F_Fact, function(x) { x[1:min(c(2, length(x)))] }) # We will add nesting later for simplicity using the UI. 
  #     F_Fact <- lapply(F_Fact, function(x) {
  #       x <- strsplit(x, "___")
  #       x <- x[1:2]
  #       return(x)
  #     })
  #     F_Fact <- setNames(F_Fact, paste0("Analysis_", 1:length(F_Fact)))
  #     F_Ref <- strsplit(F_Ref, "_;_")
  #     stopifnot(max(vapply(F_Ref, length, 1)) <= 2) # We can have at most 2 here: no blocking factors.
  #     F_Ref <- lapply(F_Ref, function(x) { x[1:min(c(2, length(x)))] })
  #     F_Ref <- lapply(F_Ref, function(x) {
  #       x <- strsplit(x, "___")
  #       x <- x[1:2] 
  #       return(x)
  #     })
  #     F_Ref <- setNames(F_Ref, paste0("Analysis_", 1:length(F_Ref)))
  #     tst <- length(F_Fact) == length(F_Ref)
  #     if (tst) {
  #       tst <- sum(vapply(1:length(F_Fact), function(x) { length(F_Fact[[x]]) != length(F_Ref[[x]]) }, TRUE)) == 0
  #       if (tst) {
  #         tmp <- substr(Factors, 1, 3)
  #         tst <- sum(unlist(sapply(1:length(F_Fact), function(x) {
  #           sapply(1:length(F_Ref[[x]]), function(y) {
  #             if (length(F_Fact[[x]][[y]]) > 0) {
  #               vapply(1:length(F_Fact[[x]][[y]]), function(z) {
  #                 sum(!F_Ref[[x]][[y]][[z]] %in% FactorsLevels[[Factors[match(substr(F_Fact[[x]][[y]][[z]],1 , 3), tmp)]]])
  #               }, 1)  
  #             } else { 0 }
  #           })
  #         }))) == 0
  #         if (!tst) { fFact <- fRef <- FALSE }
  #       } else { fFact <- fRef <- FALSE }
  #     } else { fFact <- fRef <- FALSE }
  #   }
  #   if (fFact && fRef) {
  #     nms <- names(F_Fact)
  #     F_Analyses <- data.frame(Analysis = nms, row.names = nms)
  #     F_Analyses$Primary <- lapply(F_Fact, function(x) { x[[1]] })
  #     F_Analyses$Primary_Ref <- lapply(F_Ref, function(x) { x[[1]] })
  #     if (makeSec) {
  #       F_Analyses$Secondary <- lapply(F_Fact, function(x) { x[[2]] })
  #       F_Analyses$Secondary_Ref <- lapply(F_Ref, function(x) { x[[2]] })
  #     }
  #   } else {
  #     F_Analyses <- data.frame(Analysis = "", Primary = "", Primary_Ref = "")
  #     if (makeSec) { F_Analyses$Secondary <- F_Analyses$Secondary_Ref <- "" }
  #     F_Analyses <- F_Analyses[integer(0),]
  #   }
  #   #
  #   F_An <- reactiveVal(F_Analyses)
  #   nr <- reactiveVal(nrow(F_Analyses))
  #   mxNr <- reactiveVal(nrow(F_Analyses))
  #   #
  #   updtAnUI <- function(reactive = TRUE) {
  #     # Update UI
  #     if (reactive) { FA <- F_An() } else { FA <- F_Analyses }
  #     return(renderUI({
  #       if (nrow(FA)) {
  #         txts <- apply(FA, 1, function(r) {
  #           #r <- FA[1,]
  #           pr <- unlist(r[[2]])
  #           prRf <- unlist(r[[3]])
  #           txt <- paste0("-> Primary: ", paste(sapply(1:length(pr), function(x) {
  #             paste0(pr[x], " = ", prRf[x])
  #           }), collapse = ", "))
  #           if (makeSec) {
  #             sc <- unlist(r[[4]])
  #             if (length(sc)) {
  #               scRf <- unlist(r[[5]])
  #               txt <- paste0(txt, "_|_",
  #                             paste0("-> Secondary: ", paste(sapply(1:length(sc), function(x) {
  #                               paste0(sc[x], " = ", scRf[x])
  #                             }), collapse = ", ")))
  #             }
  #           }
  #           return(txt)
  #         })
  #         txts <- strsplit(txts, "_\\|_")
  #         lst <- vector("list", length(txts)*2)
  #         for (r in 1:length(txts)) {
  #           txt <- txts[[r]]
  #           if (length(txt) == 1) {
  #             lst[[r*2-1]] <- list(tags$table(
  #               tags$tr(width = "100%", tags$td(width = "100%", br())),
  #               tags$tr(width = "100%", tags$td(width = "100%", HTML(paste0("Analysis ", r, ":")))),
  #               tags$tr(width = "100%", tags$td(width = "100%", HTML(txt[[1]])))
  #             ))
  #           } 
  #           if (length(txt) == 2) {
  #             lst[[r*2-1]] <- list(tags$table(
  #               tags$tr(width = "100%", tags$td(width = "100%", br())),
  #               tags$tr(width = "100%", tags$td(width = "100%", HTML(paste0("Analysis ", r, ":")))),
  #               tags$tr(width = "100%", tags$td(width = "100%", HTML(txt[[1]]))),
  #               tags$tr(width = "100%", tags$td(width = "100%", HTML(txt[[2]])))
  #             ))
  #           }
  #           lst[[r*2]] <- list(actionButton(paste0("Rmv", r), "Remove"))
  #         }
  #       } else { lst <- list(em("No analysis yet")) }
  #       return(lst)
  #     }))
  #   }
  #   # Create original remove analysis observers
  #   if (nrow(F_Analyses)) {
  #     sapply(1:nrow(F_Analyses), function(r) {
  #       observeEvent(input[[paste0("Rmv", r)]], {
  #         FA <- F_An()
  #         rws <- 1:nrow(FA)
  #         rws <- rws[which(rws != r)]
  #         FA <- FA[rws,]
  #         nr(nrow(FA))
  #         if (nr()) { FA$Analysis <- paste0("Analysis_", 1:nr()) } # Rename analyses
  #         F_An(FA)
  #         output$Analyses <- updtAnUI()
  #       })
  #     })
  #     output$Analyses <- updtAnUI(reactive = FALSE)
  #   }
  #   #
  #   # Box to define a new analysis
  #   dfltsPrim <- reactive(setNames(rep("Not used", length(Factors3)), Factors3))
  #   if (makeSec) { dfltsSec <- reactive(setNames(rep("Not used", length(Factors3)), Factors3)) }
  #   output$NewAnalysis <- renderUI({
  #     tags$table(
  #       tags$tr(width = "100%",
  #               tags$td(width = "20%", div(em(""))),
  #               tags$td(width = "40%", div(em("Primary contrasts"))),
  #               if (makeSec) { tags$td(width = "40%", div(em("Secondary contrasts"))) }
  #       ),
  #       lapply(Factors3, function(Fact){
  #         tags$tr(width = "100%",
  #                 tags$td(width = "20%", div(strong(Fact))),
  #                 tags$td(width = "40%", div(selectInput(paste0("Prim", Fact), "",
  #                                                        c("Not used", FactorsLevels[[Fact]]), dfltsPrim()[[Fact]]))),
  #                 if (makeSec) { tags$td(width = "40%", div(selectInput(paste0("Sec", Fact), "",
  #                                                                       c("Not used", FactorsLevels[[Fact]]),
  #                                                                       dfltsSec()[[Fact]]))) }
  #         )
  #       })
  #     )
  #   })
  #   # Observers for factor reference selection
  #   # These will change the selected value for other factors in the column, if incompatible
  #   EM <- expMap
  #   for (Fact in Factors3) { EM[which(is.na(EM[[Fact]])), Fact] <- "NA" }
  #   output$Msg <- renderUI({ em(" ") })
  #   sapply(Factors3, function(Fact) {
  #     sapply(c("Prim", "Sec")[1:(makeSec+1)], function(i) {
  #       observeEvent(input[[paste0(i, Fact)]], {
  #         m <- EM[which(as.character(EM[[Fact]]) == input[[paste0(i, Fact)]]),]
  #         facTrs <- Factors3[which(Factors3 != Fact)]
  #         k <- 0
  #         for (fcT in facTrs) {
  #           if (!input[[paste0(i, fcT)]] %in% c("Not used", m[[fcT]])) {
  #             k <- k+1
  #             updateTextInput(inputId = paste0(i, fcT), value = "Not used")
  #           }
  #         }
  #         output$Msg <- renderUI({ em(c(" ", "Invalid levels combination!")[(k > 0)+1]) })
  #       })
  #     })
  #   })
  #   # Parse new analysis + update UI
  #   # Template new Analysis row
  #   kol <- c("Primary", "Primary_Ref", "Secondary", "Secondary_Ref")[1:(2^(1+makeSec))]
  #   observeEvent(input$addAnalysis, {
  #     Prim <- setNames(sapply(Factors3, function(x) { input[[paste0("Prim", x)]] }), Factors3)
  #     #Prim <- setNames(rep(NA, length(Factors3)), Factors3) ; Prim[Factors3[1]] <- FactorsLevels[[Factors3[1]]][1]
  #     #
  #     Prim <- Prim[which(Prim != "Not used")]
  #     if (length(Prim)) {
  #       NuAn <- NuAn_tmplt
  #       NuAn$Primary <- list(names(Prim))
  #       NuAn$Primary_Ref <- list(Prim)
  #       FA <- F_An()
  #       #FA <- F_Analyses
  #       if (makeSec) {
  #         Sec <- setNames(sapply(Factors3, function(x) { input[[paste0("Sec", x)]] }), Factors3)
  #         #Sec <- setNames(rep(NA, length(Factors3)), Factors3) ; Sec[Factors3[2]] <- FactorsLevels[[Factors3[2]]][1]
  #         Sec <- Sec[which(Sec != "Not used")]
  #         NuAn$Secondary <- list(names(Sec))
  #         NuAn$Secondary_Ref <- list(Sec)
  #       }
  #       tst1 <- apply(NuAn[, kol], 1, function(x) {
  #         x <- unlist(x)
  #         names(x) <- NULL
  #         return(paste(x, collapse = "@"))
  #       })
  #       tst2 <- apply(FA[, kol], 1, function(x) {
  #         x <- unlist(x)
  #         names(x) <- NULL
  #         return(paste(x, collapse = "@"))
  #       })
  #       if (!tst1 %in% tst2) {
  #         FA <- rbind(FA, NuAn)
  #         nr(nrow(FA))
  #         FA$Analysis <- paste0("Analysis_", 1:nr()) #FA$Analysis <- paste0("Analysis_", 1:nrow(FA))
  #         F_An(FA)
  #         # Create new analysis remove observers
  #         if (nr() > mxNr()) {
  #           sapply((mxNr()+1):nr(), function(r) {
  #             observeEvent(input[[paste0("Rmv", r)]], {
  #               FA <- F_An()
  #               rws <- 1:nrow(FA)
  #               rws <- rws[which(rws != r)]
  #               FA <- FA[rws,]
  #               nr(nrow(FA))
  #               if (nr()) { FA$Analysis <- paste0("Analysis_", 1:nr()) } # Rename analyses
  #               F_An(FA)
  #               output$Analyses <- updtAnUI()
  #             })
  #           })
  #           mxNr(max(c(mxNr()), nr()))
  #         }
  #         # Update UI
  #         output$Analyses <- updtAnUI()
  #       }
  #     }
  #   })
  #   #
  # }
  observeEvent(input$saveBtn, {
    # if (Param$F.test) {
    #   F_Analyses <<- F_An()
    #   #print(F_Analyses)
    #   # Update parameters here!
    #   if (nrow(F_Analyses)) {
    #     PAR <- Param
    #     tmpFct <- data.frame(Prim = vapply(F_Analyses$Primary, function(x) { paste(substr(x, 1, 3), collapse = "___") }, ""))
    #     tmpRf <- data.frame(Prim = vapply(F_Analyses$Primary_Ref, paste, "", collapse = "___"))
    #     if (makeSec) {
    #       tmpFct$Sec <- vapply(F_Analyses$Secondary, function(x) { paste(substr(x, 1, 3), collapse = "___") }, "")
    #       tmpRf$Sec <- vapply(F_Analyses$Secondary_Ref, paste, "", collapse = "___")
    #     } else {
    #       tmpFct$Sec <- ""
    #       tmpRf$Sec <- ""
    #     }
    #     if (Param$Ratios.Groups_Nested) { tmpFct$Block <- "Rep" }
    #     PAR$F.test_factors <- paste(gsub("_;_$", "", do.call(paste, c(tmpFct, sep = "_;_"))), collapse = "_|_")
    #     PAR$F.test_factors_ref <- paste(gsub("_;_$", "", do.call(paste, c(tmpRf, sep = "_;_"))), collapse = "_|_")
    #     Param <<- PAR
    #   } else { Param$F.test <- FALSE }
    # }
    rfLev$`Reference level` <- sapply(1:nrow(rfLev2), function(x) {
      input[[paste0("Ref_", as.character(x))]]
    })
    assign("rfLev", rfLev, envir = .GlobalEnv)
    assign("appRunTest", appRunTest, envir = .GlobalEnv)
    stopApp()
  })
  #observeEvent(input$cancel, { stopApp() })
  session$onSessionEnded(function() {
    assign("appRunTest", appRunTest, envir = .GlobalEnv)
    stopApp()
  })
}
appTxt2 <- gsub("myApp", "myApp2", gsub("\\(ui", "(ui2", gsub(", server", ", server2", runApp)))
runKount <- 0
while ((!runKount)||(!exists("appRunTest"))) {
  eval(parse(text = appTxt2), envir = .GlobalEnv)
  runKount <- runKount+1
}
#
wRf <- lapply(1:nrow(rfLev), function(i) { #i <- 1 #i <- 2
  expMap$Row[which((factLevComb2 == rfLev$Group[i])&(factLevComb1 == rfLev$"Reference level"[i]))]
}) # Row indices valid for Exp.map, not necessarily for expMap!
Exp.map$Reference <- 1:nrow(Exp.map) %in% unlist(wRf)
tmp <- Exp.map
tmp$MQ.Exp <- vapply(tmp$MQ.Exp, paste, "", collapse = ";")
tst <- try(write.csv(tmp, file = ExpMapPath, row.names = FALSE), silent = TRUE)
while (("try-error" %in% class(tst))&&(grepl("cannot open the connection", tst[1]))) {
  dlg_message(paste0("File \"", ExpMapPath, "\" appears to be locked for editing, close the file then click ok..."), "ok")
  tst <- try(write.csv(tmp, file = ExpMapPath, row.names = FALSE), silent = TRUE)
}

# In case we have enriched for a PTM, it helps to check how good the enrichment was:
# This should be before any PSMs are filtered out - we want to look at the data "straight out of the MS"
tstEnrich <- unique(FracMap$`PTM-enriched`)
tstEnrich <- tstEnrich[which((!is.na(tstEnrich))&(tstEnrich != "NA"))]
if (length(tstEnrich)) {
  dir <- paste0(wd, "/Summary plots")
  dirlist <- unique(c(dirlist, dir))
  if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
  for (Mod in tstEnrich) { #Mod <- tstEnrich[1]
    pat <- paste0("\\(", Modifs$Mark[match(Mod, Modifs$`Full name`)], "\\)")
    ev[[Mod]] <- grepl(pat, ev$"Modified sequence")
    for (Type in c("PSMs", "Pep")) { #Type <- "PSMs"
      if (Type == "PSMs") {
        tst <- data.table(Mod = ev[[Mod]],
                          MS_file = ev$"Raw file path",
                          temp = 1)
        tst <- tst[, list(Count = sum(temp)), by = list(MS_file = MS_file, Mod = Mod)]
        tst <- as.data.frame(tst)
        Root <- "PSM"
      }
      if (Type == "Pep") {
        tst <- data.table(Mod = ev[[Mod]],
                          MS_file = ev$"Raw file path",
                          ModSeq = ev$`Modified sequence`)
        tst <- tst[, list(Count = length(unique(ModSeq))), by = list(MS_file = MS_file, Mod = Mod)]
        tst <- as.data.frame(tst)
        Root <- "peptidoform"
      }
      colnames(tst)[2] <- Mod
      m <- match(tst$MS_file, FracMap$"Raw file")
      tst$Sample <- FracMap$MQ.Exp[m]
      tst <- tst[which(!is.na(tst$Sample)),]
      tst <- tst[which(tst$Sample %in% unlist(Exp.map$MQ.Exp)),]
      tst <- tst[which(tst$Sample %in% unlist(Exp.map$MQ.Exp)),]
      #which(vapply(Exp.map$MQ.Exp, function(y) { x %in% unlist(y) }, TRUE))
      tst2 <- reshape2::melt(tst)
      colnames(tst2) <- gsub("\\(|\\)|\\[|\\]", "", gsub(" ", "_", colnames(tst2)))
      frml <- as.formula(paste0("MS_file ~ `", gsub("\\(|\\)|\\[|\\]", "", gsub(" ", "_", Mod)), "`"))
      tst2 <- cast(tst2, frml, fun.aggregate = sum)
      kN <- paste0("Non-", Mod, "-modified")
      kY <- paste0(Mod, "-modified")
      colnames(tst2)[which(colnames(tst2) == "FALSE")] <- kN
      colnames(tst2)[which(colnames(tst2) == "TRUE")] <- kY
      tst2[[paste0(Mod, " [%]")]] <- signif(100*tst2[[kY]]/(tst2[[kY]]+tst2[[kN]]), 3)
      tst2$Sample <- tst$Sample[match(tst2$MS_file, tst$MS_file)]
      tst2 <- tst2[, c("Sample", "MS_file", kN, kY, paste0(Mod, " [%]"))]
      write.csv(tst2, paste0(dir, "/", Mod, "-", Root, "s per MS file.csv"), row.names = FALSE)
      rw <- unique(tst$MS_file)
      for (i in c(TRUE, FALSE)) {
        w <- which(tst[[Mod]] == i)
        w2 <- which(!rw %in% tst$MS_file[w])
        if (length(w2)) {
          tmp <- tst[which((tst$MS_file %in% rw[w2])&(tst[[Mod]] == !i)),]
          tmp[[Mod]] <- i
          tmp$Count <- 0
          tst <- rbind(tst, tmp)
        }
      }
      tst[[Mod]] <- factor(c("-", "+")[tst[[Mod]]+1], levels = c("-", "+"))
      ttl <- paste0(Mod, "-", Root, "s per MS file")
      plot <- ggplot(tst) + geom_bar(stat = "identity", position = "dodge",
                                     aes(x = MS_file, y = Count, fill = .data[[Mod]])) +
        theme_bw() + scale_fill_viridis(discrete = TRUE) + ggtitle(ttl) +
        theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1, size = 5),
              plot.margin = unit(c(0, 0, 0, 3), "in"))
      #poplot(plot, 12, 20)
      ggsave(paste0(dir, "/", ttl, ".jpg"), plot, dpi = 300)
      ggsave(paste0(dir, "/", ttl, ".pdf"), plot, dpi = 300)
      ReportCalls <- AddPlot2Report()
    }
  }
}

# Write PTMs table
temp <- Modifs
w <- which(vapply(colnames(Modifs), function(x) { "list" %in% class(Modifs[[x]]) }, TRUE))
for (i in w) { temp[[i]] <- vapply(temp[[i]], paste, "", collapse = ", ") }
dir <- paste0(wd, "/Workflow control")
dirlist <- unique(c(dirlist, dir))
if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
write.csv(temp, paste0(dir, "/Modifications.csv"), row.names = FALSE)
ReportCalls$Calls <- append(ReportCalls$Calls, "body_add_fpar(Report, fpar(ftext(\"PTMs table:\", prop = WrdFrmt$Body_text_ital), fp_p = WrdFrmt$just))")
ReportCalls$Calls <- append(ReportCalls$Calls, "body_add_table(Report, ReportCalls$Objects$AABiases)")
ReportCalls$Objects$AABiases <- temp
ReportCalls <- AddSpace2Report()

#
if (LabelType == "LFQ") {
  aggrCol <- "RSA"
  tstMQXp2 <- setNames(vapply(MQ.Exp, function(x) {
    unique(Exp.map$Ref.Sample.Aggregate[tstMQXp[[x]]])
  }, ""), MQ.Exp) # Here for LFQ we should not have multiples, this should fail if it is the case!!!
  ev[[aggrCol]] <- tstMQXp2[ev$MQ.Exp]
  tmp <- ev[, c("Proteins", "Intensity", aggrCol)]
}
if (LabelType == "Isobaric") {
  # Isobaric case - subtle difference
  aggrCol <- "Iso"
  tstMQXp2 <- setNames(vapply(MQ.Exp, function(x) {
    unique(Exp.map$Isobaric.set[tstMQXp[[x]]])
  }, 1), MQ.Exp)
  tmp <- ev[, c("Proteins", "Intensity")]
  tmp[[aggrCol]] <- tstMQXp2[ev$MQ.Exp]
}
# Plot of contamination levels per sample
dir <- paste0(wd, "/Summary plots")
dirlist <- unique(c(dirlist, dir))
if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
tmp2 <- listMelt(strsplit(tmp$Proteins, ";"), 1:nrow(tmp), c("Protein", "Row"))
m <- match(tmp2$Protein, db$`Protein ID`)
w <- which(!is.na(m))
tmp2 <- tmp2[w,]; m <- m[w]
# For our purpose here we must match contaminant proteins.
tmp2$Cont <- db$`Potential contaminant`[m]
if (tstOrg) {
  tmp2$Organism <- db[m, dbOrgKol]
  tmp2 <- as.data.table(tmp2)
  f0 <- function(x) { c("Target", "Contaminant")[("Contaminant" %in% x)+1] }
  tmp2 <- tmp2[, list(x = f0(Organism)), by = list(Group.1 = Row)]
  tmp2 <- as.data.frame(tmp2)
} else {
  tmp2 <- as.data.table(tmp2)
  f0 <- function(x) { c("Target", "Contaminant")[("+" %in% x)+1] }
  tmp2 <- tmp2[, list(x = f0(Cont)), by = list(Group.1 = Row)]
  tmp2 <- as.data.frame(tmp2)
}
tmp$Organism <- tmp2$x[match(1:nrow(tmp), tmp2$Group.1)]
tmp <- tmp[which(!is.na(tmp$Organism)),]
tmp$Organism <- factor(tmp$Organism, levels = c("Contaminant", "Target"))
tmp$Intensity <- as.numeric(tmp$Intensity)
tmp <- aggregate(tmp$Intensity, list(tmp[[aggrCol]], tmp$Organism), sum, na.rm = TRUE)
if (LabelType == "LFQ") {
  k <- "Sample"
  colnames(tmp) <- c(k, "Organism", "Total intensity")
  tmp[[k]] <- factor(cleanNms(tmp[[k]]), levels = cleanNms(unique(Exp.map$Ref.Sample.Aggregate)))
}
if (LabelType == "Isobaric") {
  k <- "Isobaric.set"
  colnames(tmp) <- c(k, "Organism", "Total intensity")
  tmp[[k]] <- factor(cleanNms(tmp[[k]]), levels = Iso)
}
ttl <- "Contributions to TIC"
plot <- ggplot(tmp) +
  geom_bar(stat = "identity", aes(x = .data[[k]], y = `Total intensity`, fill = Organism)) +
  theme_bw() + scale_fill_viridis(discrete = TRUE, begin = 0.8, end = 0.2) +
  ggtitle(ttl, subtitle = "Summed TIC for each class of identified peptides") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(plot) # This type of QC plot does not need to pop up, the side panel is fine
ggsave(paste0(wd, "/Summary plots/", ttl, ".jpeg"), plot, dpi = 150, width = 10, height = 10, units = "in")
ggsave(paste0(wd, "/Summary plots/", ttl, ".pdf"), plot, dpi = 150, width = 10, height = 10, units = "in")

# Time points
if (exists("Tim")) {
  if (("Time.Points" %in% colnames(Param))&&(Param$Time.Points != "")) {
    tp %<o% as.character(sort(as.numeric(unlist(strsplit(Param$Time.Points, ";")))))
    if (("Time.Point.Names" %in% colnames(Param))&&(Param$Time.Point.Names != "")) {
      names(tp) <- gsub(" ", ".", unlist(strsplit(Param$Time.Point.Names, ";")))
    } else {
      names(tp) <- tp
    }
    if (sum(tp != Tim) > 0) {
      stop("Review your time points!")
    } else {
      Tim <- tp
    }
  } else {
    Tim <- as.character(sort(as.numeric(Tim)))
  }
}

rm(list = ls()[which(!ls() %in% .obj)])
Script <- readLines(ScriptPath)
saveImgFun(BckUpFl)
#loadFun(BckUpFl)
source(parSrc, local = FALSE)
