### A sourced script to edit the Experiment map (table defining experimental structure)
#
library(shiny)
library(shinyjs)
library(DT)
#
expKl <- c("MQ.Exp", "Parent sample")
expKl <- expKl[which(expKl %in% colnames(FracMap))[1]]

Factors2 %<o% Factors[which(!Factors %in% c("Experiment",
                                            "Replicate",
                                            "Isobaric.set",
                                            "Time.point"))]
#
labelMode <- match(LabelType, c("LFQ", "Isobaric"))
if ((exists("Exp.map"))&&(nrow(Exp.map))) {
  tst <- (sum(!c("MQ.Exp", "Sample name") %in% colnames(Exp.map)) == 0)&&
    (sum(!FracMap[[expKl]] %in% Exp.map$MQ.Exp) == 0)
  if (tst) { ExpMap <- Exp.map }
}
if ((exists("ExpMap"))&&(nrow(ExpMap))) {
  tst <- (sum(!c("MQ.Exp", "Sample name") %in% colnames(ExpMap)) == 0)&&
    (sum(!FracMap[[expKl]] %in% ExpMap$MQ.Exp) == 0)
  if (!tst) { rm(ExpMap) }
}
if ((!exists("ExpMap"))||(!nrow(ExpMap))) {
  if (LabelType == "LFQ") {
    ExpMap <- data.frame("Experiment" = "?",
                         "Replicate" = "?",
                         "MQ.Exp" = MQ.Exp,
                         "Reference" = FALSE,
                         "Sample name" = MQ.Exp,
                         "Use" = TRUE,
                         check.names = FALSE)
  }
  if (LabelType == "Isobaric") {
    ExpMap <- data.frame("Experiment" = "?",
                         "Replicate" = "?",
                         "MQ.Exp" = as.character(unlist(sapply(MQ.Exp, function(x) { rep(x, length(get(IsobarLab))) }))),
                         "Reference" = FALSE,
                         "Sample name" = as.character(unlist(sapply(MQ.Exp, function(x) { paste0(x, "_", get(IsobarLab)) }))),
                         "Isobaric label" = rep(get(IsobarLab), length(MQ.Exp)),
                         "Isobaric label details" = rep(IsobarLabDet, length(MQ.Exp)),
                         "Isobaric.set" = "?",
                         "Use" = TRUE,
                         check.names = FALSE)
  }
}
if (!"list" %in% class(ExpMap$MQ.Exp)) {
  ExpMap$MQ.Exp <- strsplit(ExpMap$MQ.Exp, ";")
}
if ((LocAnalysis)&&(!"Proportion" %in% colnames(ExpMap))) { ExpMap$Proportion <- 1 }
for (Fact in Factors2) { if (!Fact %in% colnames(ExpMap)) { ExpMap[[Fact]] <- "?" } }
tmp <- aggregate(FracMap$Fraction, list(FracMap[[expKl]]), function(x) { paste(sort(unique(x)), collapse = ";") })
tmp2 <- listMelt(ExpMap$MQ.Exp, 1:nrow(ExpMap))
tmp2$Fractions <- tmp$x[match(tmp2$value, tmp$Group.1)]
tmp2 <- aggregate(tmp2$Fractions, list(tmp2$L1), function(x) { paste(as.character(sort(unique(as.integer(x)))), collapse = ";") })
ExpMap$Fractions <- tmp2$x[match(1:nrow(ExpMap), tmp2$Group.1)]
if (LabelType == "LFQ") {
  ExpMap$Use <- as.logical(vapply(ExpMap$MQ.Exp, function(x) { max(FracMap$Use[match(x, FracMap[[expKl]])]) }, 1))
}
tmpTbl <- ExpMap
tst <- lapply(colnames(tmpTbl), function(x) { typeof(tmpTbl[[x]]) })
w <- which(tst == "list")
if (length(w)) { for (i in w) { tmpTbl[[i]] <- vapply(tmpTbl[[i]], paste, "", collapse = ";") }}
tst <- try(write.csv(tmpTbl, file = ExpMapPath, row.names = FALSE), silent = TRUE)
while (("try-error" %in% class(tst))&&(grepl("cannot open the connection", tst[1]))) {
  dlg_message(paste0("File \"", ExpMapPath, "\" appears to be locked for editing, close the file then click ok..."), "ok")
  tst <- try(write.csv(tmpTbl, file = ExpMapPath, row.names = FALSE), silent = TRUE)
}
ExpMap <- ExpMap[which(vapply(ExpMap$MQ.Exp, function(x) { sum(x %in% FracMap[[expKl]]) > 0 }, TRUE)),]
#
# Edit map
ExpData <- read.csv(ExpMapPath, check.names = FALSE)
ExpData <- ExpData[which(vapply(ExpData$MQ.Exp, function(x) { sum(x %in% FracMap[[expKl]]) }, 1) > 0),]
for (Fact in Factors[which(!Factors %in% colnames(ExpData))]) { ExpData[[Fact]] <- "?" }
if (LabelType == "LFQ") {
  ExpData$Use <- as.logical(vapply(ExpData$MQ.Exp, function(x) { max(FracMap$Use[match(x, FracMap[[expKl]])]) }, 1))
}
if (LocAnalysis) {
  if (!"Proportion" %in% colnames(ExpData)) {
    ExpData$Proportion <- 1
  } else {
    ExpData$Proportion <- suppressWarnings(as.numeric(ExpData$Proportion))
    w <- which(is.na(ExpData$Proportion))
    if (length(w)) {
      warning("Invalid values in \"Proportion\" column!")
      ExpData$Proportion[w] <- 1
    }
  }
}
tst <- vapply(FactorsLevels, length, 1)
Fact1 <- Factors[which(tst == 1)]
Fact2 <- Factors[which(tst > 1)]
Others <- c("MQ.Exp", "Sample name")
Others <- Others[which(Others %in% colnames(ExpData))]
nr <- nrow(ExpData)
rws <- seq_len(nr)
OtherIDs <- setNames(lapply(Others, function(x) { paste0(x, "___", rws)} ), Others)
Fact2IDs <- setNames(lapply(Fact2, function(x) { paste0(x, "___", rws)} ), Fact2)
AllIDs <- append(OtherIDs, Fact2IDs)
AllIDs$Use <- paste0("Use___", rws)
ALLIDS <- setNames(unlist(AllIDs), NULL)
for (fct in Fact1) { #fct <- Fact1[1]
  # (Since the table can be manually edited too, we want to make sure to correct any typos in our single level factor columns)
  ExpData[[fct]] <- rep(FactorsLevels[[fct]], nr)
}
L <- length(FactorsLevels[["Replicate"]])
dflt_Rpl <- 1:(nr+L) %% L
dflt_Rpl[which(dflt_Rpl == 0)] <- L
dflt_Rpl <- FactorsLevels[["Replicate"]][dflt_Rpl] # Easy template
k1 <- c("MQ.Exp", "Experiment", "Replicate")
k2 <- colnames(ExpData)
k2 <- k2[which(!k2 %in% c(k1, "Sample name", "Use"))]
ExpData <- ExpData[which(vapply(ExpData$MQ.Exp, function(x) { sum(x %in% MQ.Exp) }, 1) > 0), ]
# Test order
tst <- suppressWarnings(as.integer(substr(ExpData$MQ.Exp, 1, 1)))
if (!sum(is.na(tst))) {
  nc1 <- nchar(ExpData$MQ.Exp)
  nc2 <- nchar(gsub("^[0-9]+", "", ExpData$MQ.Exp))
  ord <- as.integer(substr(ExpData$MQ.Exp, 1, nc1-nc2))
  if (length(unique(ord)) == length(ord)) {
    ExpData <- ExpData[order(ord, decreasing = FALSE),]
  }
}
#
if (LabelType == "Isobaric") {
  ExpData$"Parent sample" <- do.call(paste, c(ExpData[, c("MQ.Exp", "Isobaric label details")], sep = "_"))
} else {
  ExpData[["Parent sample"]] <- ExpData$MQ.Exp
}
# Original table column widths
wTest0 <- setNames(vapply(colnames(ExpData), function(k) { #k <- colnames(ExpData)[1]
  tmp <- ExpData[[k]]
  if ("logical" %in% class(tmp)) { tmp <- as.integer(tmp) }
  tmp <- as.character(tmp)
  tst <- k %in% Fact2
  x <- max(nchar(c(k, tmp)) + 3 + 3*tst, na.rm = TRUE)
  if (tst) {
    x <- max(c(x, nchar(FactorsLevels[[k]]) + 6), na.rm = TRUE)
  }
  x <- x*10
  if (is.na(x)) { x <- 15 } else { x <- max(c(ceiling(x/10)*10, 30)) }
  return(x)
}, 1), colnames(ExpData))
#
# Dummy table for app
tst <- vapply(FactorsLevels, length, 1)
Fact1 <- Factors[which(tst == 1)]
Fact2 <- Factors[which(tst > 1)]
#Fact2 <- "Developmental.stage"
myKol <- c(Others, "Sample name")
if (LabelType == "Isobaric") {
  myKol <- c(myKol, "Isobaric label", "Isobaric label details")
}
myKol <- unique(c(myKol, Factors, "Sample name"))
myKol2 <- c("Parent sample", myKol[which(myKol != "Parent sample")])
ExpData2 <- ExpData[, myKol2]
kol <- c()
for (fct in Fact2) { #fct <- Fact2[1]
  IDs <- Fact2IDs[[fct]]
  lvls <- FactorsLevels[[fct]]
  lvls2 <- lvls[which(!is.na(lvls))]
  ExpData2[[fct]] <- shinySelectInput(ExpData[[fct]],
                                      fct,
                                      lvls,
                                      paste0(wTest0[fct], "px"),
                                      TRUE)
  fdNm <- paste0(fct, "___FD")
  wTest0[fdNm] <- 15
  ExpData2[[fdNm]] <- shinyFDInput(fct, nr, TRUE, paste0(wTest0[fdNm], "px"))
  kol <- c(kol, fct, fdNm)
  if (fct == "Replicate") {
    incrNm <- paste0(fct, "___INCR")
    wTest0[incrNm] <- 15
    ExpData2[[incrNm]] <- shinyPlusFD(fct, nr, width = paste0(wTest0[incrNm], "px"))
    kol <- c(kol, incrNm)
  }
}
if (LocAnalysis) {
  ExpData2$Proportion <- shinyPropInput(ExpData$Proportion)
  fdNm <- "Proportion___FD"
  wTest0[fdNm] <- 15
  ExpData2[[fdNm]] <- shinyFDInput("Proportion", nr, TRUE, paste0(wTest0[fdNm], "px"))
  kol <- c(kol, "Proportion", fdNm)
  ALLIDS <- c(ALLIDS, paste0("Proportion___", rws))
}
idsL <- length(ALLIDS)
smplWdth <- paste0(as.character(max(nchar(ExpData2$"Parent sample"))*10), "px")
kol2 <- colnames(ExpData2)[which(!colnames(ExpData2) %in% c(kol, "Use", "Sample name"))]
ExpData2$Use <- shinyCheckInput(ExpData$Use,
                                "Use")
ExpData2$Use___FD <- shinyFDInput("Use", nr, TRUE)
ExpData2 <- ExpData2[, c(kol2, kol, "Use", "Use___FD", "Sample name")]
# Estimate table column widths
wTest1 <- vapply(colnames(ExpData2), function(k) { #k <- colnames(ExpData2)[1]
  if (k == "Parent sample") { k <- "MQ.Exp" }
  if (k %in% names(wTest0)) { x <- wTest0[k] } else { x <- 30 }
  return(x)
}, 1)
wTest2 <- sum(wTest1) + 15 + ncol(ExpData2)*5
wTest1 <- paste0(as.character(wTest1), "px")
wTest1 <- aggregate((1:length(wTest1))-1, list(wTest1), c)
wTest1 <- apply(wTest1, 1, function(x) {
  x2 <- as.integer(x[[2]])
  list(width = x[[1]],
       targets = x2,
       names = colnames(ExpData2)[x2+1])
})
#
g <- grep("___((FD)|(INCR))$", colnames(ExpData2))
colnames(ExpData2)[g] <- ""
#
facLevels2 <- lapply(FactorsLevels, function(x) {
  unique(c(x, NA))
})
facLevels2$Use <- c(TRUE, FALSE)
#
edith <- list(target = "column",
              disable = list(columns = c(0, match(Fact1, colnames(ExpData2))-1)))
tmp <- c(0:(ncol(ExpData2)-1))
tmp <- tmp[which(!tmp %in% edith$disable$columns)]
edith$enable <- list(columns = tmp)
isoMsg <- ""
if (LabelType == "Isobaric") {
  isoMsg <- "This is an isobarically-labelled experiment: assign value \"Mixed_IRS\", if and as relevant, where the value is available, otherwise use NA. For replicates assign an arbitrary replicate number."
}
#
ExpData2$MQ.Exp <- vapply(ExpData2$MQ.Exp, paste, "", collapse = ";")
#
appNm <- paste0(dtstNm, " - Exp. map")
ui <- fluidPage(
  useShinyjs(),
  setBackgroundColor( # Doesn't work
    color = c(#"#F8F8FF",
      "#EBEFF7"),
    gradient = "linear",
    direction = "bottom"
  ),
  extendShinyjs(text = jsToggleFS, functions = c("toggleFullScreen")),
  # Dummy, hidden div to load arrow-down icon
  # tags$div(
  #   class = "container",
  #   style = "display: none;",
  #   tags$div(
  #     style = "margin-top: 50xp;",
  #     actionButton(
  #       "add_thing",
  #       label = "do it",
  #       class = "btn-success",
  #       icon = icon("arrow-down")
  #     )
  #   )
  # ),
  #
  tags$head(tags$style(HTML("table {table-layout: fixed;"))), # So table widths can be properly adjusted!
  titlePanel(tag("u", ExpMapNm),
             appNm),
  h2(dtstNm), 
  h3("Define the relationship between biological samples and experimental factors."),
  br(),
  h4(em(tags$div("\"Parent sample\" is equivalent with \"Experiments\" in the sense of MaxQuant's \"Raw data\" or FragPipe's \"Workflow\" tabs.", tags$br(),
                 "Note that this is different from the meaning of \"Experiment\" we use here (that of a group of samples to process together and compare to each other)."))),
  h4(strong(em(isoMsg, style = "color:red", .noWS = "outside"))),
  br(),
  actionBttn("saveBtn", "Save", icon = icon("save"), color = "success", style = "pill"),
  withSpinner(DT::DTOutput("ExpTbl", width = wTest2))
)
if (exists("ExpData3")) { rm(ExpData3) }
server <- function(input, output, session) {
  # Initialize output table
  ExpData3 <- ExpData # Output table
  #
  # Render dummy table
  output$ExpTbl <- DT::renderDT({ ExpData2 },
                                FALSE,
                                escape = FALSE,
                                class = "compact",
                                selection = "none",
                                rownames = FALSE,
                                editable = edith,
                                options = list(dom = "t",
                                               paging = FALSE,
                                               ordering = FALSE,
                                               autowidth = TRUE,
                                               columnDefs = wTest1,
                                               scrollX = FALSE),
                                # the callback is essential to capture the inputs in each row
                                callback = JS("table.rows().every(function(i, tab, row) {
  var $this = $(this.node());
  $this.attr('id', this.data()[0]);
  $this.addClass('shiny-input-container');
});
Shiny.unbindAll(table.table().node());
Shiny.bindAll(table.table().node());"))
  #
  # Fill down
  sapply(1:idsL, function(x) {
    id1 <- ALLIDS[x]
    id2 <- paste0(ALLIDS[x], "___FD")
    tmp <- unlist(strsplit(id1, "___"))
    fct <- tmp[[1]]
    i <- as.integer(tmp[[2]])
    if (i < nr) {
      observeEvent(input[[id2]],
                   {
                     x <- input[[id1]]
                     # if (fct %in% names(facLevels2)) {
                     #   tp <- typeof(facLevels2[[fct]])
                     #   if (typeof(x) != tp) { x <- get(paste0("as.", tp))(x) }
                     # }
                     for (k in (i+1):nr) {
                       idK <- paste0(fct, "___", as.character(k))
                       if (fct == "Use") {
                         updateCheckboxInput(session, idK, NULL, x)
                       } else {
                         if (fct == "Proportion") {
                           updateNumericInput(session, idK, NULL, x, 0, 1, 0.001)
                         } else {
                           updateSelectInput(session, idK, NULL, facLevels2[[fct]], x)
                         }
                       }
                     }
                   })
    }
  })
  # Incremental fill-down for replicates
  sapply(rws, function(i) {
    if (i < nr) {
      id1 <- paste0("Replicate___", i)
      id2 <- paste0("Replicate___", i, "___INCR")
      observeEvent(input[[id2]],
                   {
                     x <- input[[id1]]
                     #tp <- typeof(facLevels2[[fct]])
                     #if (typeof(x) != tp) { x <- get(paste0("as.", tp))(x) }
                     rplRg <- (i+1):nr
                     l <- length(rplRg)
                     m <- match(x, dflt_Rpl)+1
                     rplVal <- dflt_Rpl[m:(m+l-1)]
                     for (k in rplRg) {
                       y <- rplVal[k-i]
                       idK <- paste0("Replicate___", as.character(k))
                       updateSelectInput(session, idK, NULL, facLevels2[["Replicate"]], y)
                     }
                     
                   }
      )
    }
  })
  # Manual cell edit (sample names)
  observeEvent(input$ExpTbl_cell_edit, {
    kl <- colnames(ExpData2)[input$ExpTbl_cell_edit$col+1]
    if (kl %in% colnames(ExpData3)) {
      ExpData3[input$ExpTbl_cell_edit$row,
               kl] <- input$ExpTbl_cell_edit$value
    } else {
      warning(paste0("Could not find column ", kl, " from dummy table in the real table!"))
    }
  })
  # Save
  observeEvent(input$saveBtn, {
    kls <- c(Fact2, "Use")
    if (LocAnalysis) { kls <- c(kls, "Proportion") }
    for (kl in kls) {
      ExpData3[[kl]] <- sapply(rws, function(i) {
        input[[paste0(kl, "___", i)]]
      })
      if (kl == "Use") {
        ExpData3[[kl]] <- as.logical(ExpData3[[kl]])
      } else {
        if (kl == "Proportion") {
          ExpData3[[kl]] <- as.numeric(ExpData3[[kl]])
        } else {
          typ <- typeof(FactorsLevels[[kl]])
          ExpData3[[kl]] <- get(paste0("as.", typ))(ExpData3[[kl]])
          # Consider here detecting if a factor needs conversion to numeric/integer...
        }
      }
    }
    ExpData3$"Sample name" <- ExpData2$"Sample name"
    ExpData3$MQ.Exp <- ExpData2$MQ.Exp
    assign("ExpData3", ExpData3, envir = .GlobalEnv)
    stopApp()
  })
  #observeEvent(input$cancel, { stopApp() })
  session$onSessionEnded(function() { stopApp() })
}
runKount <- 0
while ((!runKount)||(!exists("ExpData3"))) {
  eval(parse(text = runApp), envir = .GlobalEnv)
  runKount <- runKount+1
}
#
Exp.map %<o% ExpData3
k0 <- colnames(ExpMap)
k0 <- k0[which(!k0 %in% colnames(Exp.map))]
if (length(k0)) {
  e1 <- vapply(Exp.map$MQ.Exp, paste, "", collapse = ";")
  e3 <- vapply(ExpData3$MQ.Exp, paste, "", collapse = ";")
  Exp.map[, k0] <- ExpData3[match(e1, e3), k0]
}
Exp.map$Use <- as.logical(Exp.map$Use)
Exp.map$Use[which(is.na(Exp.map$Use))] <- FALSE
#sum(!Exp.map$Use)
#sum(Exp.map$Use)
tmpTbl <- Exp.map
tst <- lapply(colnames(tmpTbl), function(x) { typeof(tmpTbl[[x]]) })
w <- which(tst == "list")
if (length(w)) { for (i in w) { tmpTbl[[i]] <- vapply(tmpTbl[[i]], paste, "", collapse = ";") }}
tst <- try(write.csv(tmpTbl, file = ExpMapPath, row.names = FALSE), silent = TRUE)
while (("try-error" %in% class(tst))&&(grepl("cannot open the connection", tst[1]))) {
  dlg_message(paste0("File \"", ExpMapPath, "\" appears to be locked for editing, close the file then click ok..."), "ok")
  tst <- try(write.csv(tmpTbl, file = ExpMapPath, row.names = FALSE), silent = TRUE)
}
#
#system(paste0("open \"", wd, "/", ExpMapNm, ".csv\""))
#Exp.map2 %<o% read.csv(ExpMapPath, check.names = FALSE)
Exp.map <- Exp.map[which(Exp.map$Use),]
kol <- colnames(Exp.map)
kol1 <- c("MQ.Exp", "Parent sample", "Sample name", "Use")
if (LabelType == "Isobaric") {
  kol1 <- unique(c(kol1, "Isobaric label", "Isobaric label details"))
}
kol2 <- c("Reference", "Fractions", "Proportion", Factors)
kol2 <- kol2[which(kol2 %in% colnames(Exp.map))]
Exp.map2 <- aggregate(Exp.map[, kol1],
                     lapply(kol2, function(k) {
                       x <- Exp.map[[k]]
                       w <- which(is.na(x))
                       if (length(w)) { x[w] <- "NA" }
                       return(x)
                     }), function(x) { paste(unique(x), collapse = ";") })
colnames(Exp.map2) <- c(kol2, kol1)
Exp.map2 <- Exp.map2[, c(kol1, kol2)]
for (k in kol2) {
  x <- Exp.map2[[k]]
  w <- which(x == "NA")
  if (length(w)) { Exp.map2[w, k] <- NA }
}
if (!"list" %in% class(Exp.map2$MQ.Exp)) { Exp.map2$MQ.Exp <- strsplit(Exp.map2$MQ.Exp, ";") }
Exp.map2$tempName <- Exp.map2$`Sample name`
call <- paste0("Exp.map2 <- arrange(Exp.map2, ", paste(c("tempName", Factors), collapse = ", "), ")")
#cat(call)
eval(parse(text = call))
Exp.map2$tempName <- NULL
Exp.map <- Exp.map2
kl <- c("Level", "Count")
tst <- setNames(lapply(Factors, function(Fact) {
  magrittr::set_colnames(aggregate(Exp.map[[Fact]], list(Exp.map[[Fact]]), length), kl)
}), Factors)
tst <- lapply(Factors, function(x) { #x <- 1 #x <- 2
  dat <- tst[[x]]
  rg <- (1:nrow(dat))+1
  NC <- setNames(lapply(kl, function(k) {
    nchar(c(k, dat[[k]]))
  }), kl)
  maxNC <- setNames(vapply(NC, max, 1), kl)
  dat <- as.data.frame(do.call(cbind, lapply(kl, function(k) { #k <- kl[1]
    dt <- as.character(dat[[k]])
    nc <- NC[[k]][rg]
    w <- which(nc < maxNC[k])
    if (length(w)) {
      dt[w] <- vapply(w, function(z) { paste(c(dt[z], rep(" ", maxNC[k] - nc[z])), collapse = "") }, "")
    }
    return(dt)
  })))
  colnames(dat) <- vapply(kl, function(k) {
    if (NC[[k]][1] < maxNC[k]) { k <- paste(c(k, rep(" ", maxNC[k] - NC[[k]][1])), collapse = "") }
    return(k)
  }, "")
  for (i in 1:ncol(dat)) {
    dat[[colnames(dat)[i]]] <- format(dat[[colnames(dat)[i]]], width = ceiling(maxNC[[kl[i]]]*1.5), justify = "left")
  }
  return(dat)
})

msg2 <- msg <- paste(c("Check the number of samples per factor level below. Is everything ok? If not, click \"no\" to get back to editing the table.\n\n   -----\n", unlist(lapply(names(tst), function(nm) {
  c(paste0("-> ", nm),
    paste0("     ",
           c(paste(colnames(tst[[nm]]), collapse = "          "),
                   apply(tst[[nm]], 1, paste, collapse = "          "))),
    c("\n   -----\n"))
})), "\n"), collapse = "\n")
if (nchar(msg) > 1000) {
  cat(msg)
  msg2 <- paste0(substr(msg, 1, 996), "...")
}
tstXpMp <- c(TRUE, FALSE)[match(dlg_message(msg2, "yesno", rstudio = TRUE)$res, c("yes", "no"))]
