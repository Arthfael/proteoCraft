#### Define analysis parameters
Exp %<o% unique(SamplesMap$Experiment)
moreThan1Exp %<o% (length(Exp) > 1)
SpeciesTst %<o% "Unspecified"
if ("Taxonomy" %in% colnames(db)) {
  SpeciesTst <- unique(db$Taxonomy[which(gsub(" *(\\(|\\[).*", "", db[[dbOrgKol]]) == mainOrg)])
  SpeciesTst <- SpeciesTst[which(as.character(SpeciesTst) != "NA")][1]
}
if ("Kingdom" %in% colnames(db)) {
  KingdomTst %<o% aggregate(db$Kingdom, list(db$Kingdom), length)
  KingdomTst <- KingdomTst[order(KingdomTst$x, decreasing = TRUE),]
  KingdomTst <- KingdomTst$Group.1[1]
  isEukaLike <- (KingdomTst %in% c("Eukaryota", "Archaea"))
} else {
  isEukaLike <- TRUE # Default - reasonable guess
}
isEukaLike %<o% isEukaLike
if (file.exists("AnalysisParam.RData")) {
  tmp <- AnalysisParam
  loadFun("AnalysisParam.RData")
  for (nm in names(AnalysisParam)) {
    if (!nm %in% names(tmp)) { tmp[[nm]] <- AnalysisParam[[nm]] }
  }
  AnalysisParam <- tmp
}
ptmDflt1 <- grep("^[Pp]hospho", Modifs$`Full name`, value = TRUE, invert = TRUE)
ptmDflt2 <- paste(grep("^[Pp]hospho", Modifs$`Full name`, value = TRUE), collapse = ";")
Mod4Quant %<o% Modifs$Mark[match(ptmDflt1, Modifs$`Full name`)]
if ("PTMs eligible for quantitation" %in% names(AnalysisParam)) {
  Mod4Quant <- AnalysisParam$"PTMs eligible for quantitation"
}
Mod2Xclud %<o% set_colnames(Modifs[which(!Modifs$Mark %in% Mod4Quant), c("Mark", "AA")],
                            c("Mark", "Where"))
if ("PTMs to exclude from quantitation" %in% names(AnalysisParam)) {
  Mod2Xclud <- AnalysisParam$"PTMs to exclude from quantitation"
}
klustChoices %<o% c("K-means", "hierarchical")
KlustMeth %<o% 2
#
if (Annotate) {
  allGO <- unique(unlist(strsplit(db$GO[which(!is.na(db$GO))], ";")))
  allGO2 <- paste0("GO:", gsub(".* \\[GO:|\\]$", "", allGO))
}
if (!"GO.terms.for.proteins.of.interest" %in% names(AnalysisParam)) { AnalysisParam$GO.terms.for.proteins.of.interest <- FALSE }
if (!"Custom.PGs" %in% names(AnalysisParam)) { AnalysisParam$Custom.PGs <- "" }
if (!"CRAPome_file" %in% names(AnalysisParam)) { AnalysisParam$CRAPome_file <- "" }
# shpDflt <- c("none", "vsn", "loess")
# names(shpDflt) <- gsub("^none$", "FALSE", shpDflt)
#shpDflt <- shpDflt[match(as.character(AnalysisParam$Norma.Pep.Intens.Shape), names(shpDflt))]
#if (is.na(shpDflt)) { shpDflt <- "none" }
pr <- c("Norma.Ev.Intens", "Norma.Pep.Intens", "Adv.Norma.Pep.Intens", "Norma.Pep.Intens.IRS", "Norma.Prot.Ratio", "Adv.Norma.Prot.Intens")
for (p in pr) { if ((!p %in% names(AnalysisParam))||(!is.logical(AnalysisParam[[p]]))||(is.na(AnalysisParam[[p]]))) { AnalysisParam[[p]] <- TRUE } }
p <- "Adv.Norma.Ev.Intens"
if ((!p %in% names(AnalysisParam))||(!is.logical(AnalysisParam[[p]]))||(is.na(AnalysisParam[[p]]))) { AnalysisParam[[p]] <- length(unique(FracMap$Fraction)) > 1 } 
if ((!is.null(AnalysisParam$Norma.Pep.Intens.Shape))&&(!toupper(as.character(AnalysisParam$Norma.Pep.Intens.Shape)) %in% c("FALSE", "VSN", "LOESS"))) { Norma.Pep.Intens.Shape <- FALSE }
pr <- c("Norma.Pep.Ratio", "Adv.Norma.Pep.Ratio", "Norma.Prot.Ratio.to.Biot")
for (p in pr) { if ((!is.logical(AnalysisParam[[p]]))||(is.na(AnalysisParam[[p]]))) { AnalysisParam[[p]] <- FALSE } }
QuantMethods %<o% setNames(c("Prot.Quant", "Prot.Quant + weights", "Prot.Quant.Unique", "Prot.Quant.Unique + weights",
                             "Prot.Quant2 + weights", "Prot.Quant2", "IQ_MaxLFQ", "Top3", "Top1"),
                           c(paste0("Profile_avg.", c("", ", weights = -log10(PEP)/CV", c(", unique peptides in priority", ", weights = -log10(PEP)/CV, unique peptides in priority"))),
                             paste0("Profile_avg.v2", c(", weights = -log10(PEP)/CV", "")), "MaxLFQ (iq)", "Top3", "Top1"))
QMdef <- "Prot.Quant.Unique"
if (("QuantMeth" %in% names(AnalysisParam))&&(AnalysisParam$QuantMeth %in% QuantMethods)) {
  QMdef <- AnalysisParam$QuantMeth
}
AnalysisParam$QuantMeth <- QMdef
#
if (!QMdef %in% QuantMethods[1:6]) { QMdef <- "Prot.Quant.Unique" }
QMdefnm <- names(QuantMethods)[match(QMdef, QuantMethods)]
if ("Proteome ruler calculated" %in% names(AnalysisParam)) { # Correct typo in older versions
  AnalysisParam$"Proteomic ruler calculated" <- AnalysisParam$"Proteome ruler calculated"
  AnalysisParam$"Proteome ruler calculated" <- NULL
}
if (("Proteomic ruler calculated" %in% names(AnalysisParam))&&
    (is.logical(AnalysisParam$"Proteomic ruler calculated"))&&
    (!is.na(AnalysisParam$"Proteomic ruler calculated"))) {
  protrul %<o% AnalysisParam$"Proteomic ruler calculated"
} else {
  protrul %<o% (WorkFlow %in% c("Discovery", "Regulation"))
  # Archaea and Eukaryotes have introns and histones, Bacteria do not
  protrul <- c(protrul, FALSE)[(!isEukaLike)+1]
}
AnalysisParam$"Proteomic ruler calculated" <- protrul
#
if (("ProtRulNuclL" %in% names(AnalysisParam))&&(!is.na(as.integer(AnalysisParam$ProtRulNuclL)))) {
  ProtRulNuclL <- as.integer(AnalysisParam$ProtRulNuclL)
} else {
  ProtRulNuclL <- 196
}
AnalysisParam$ProtRulNuclL <- ProtRulNuclL
#
Update_Prot_matches %<o% TRUE # See https://github.com/vdemichev/DiaNN/discussions/1631
if ("Update_Prot_matches" %in% names(AnalysisParam)) {
  Update_Prot_matches <- as.logical(AnalysisParam$Update_Prot_matches)
  if ((is.na(Update_Prot_matches))||(is.null(Update_Prot_matches))) { ImputeMissData <- TRUE }
}
AnalysisParam$Update_Prot_matches <- Update_Prot_matches
#
ImputeMissData %<o% FALSE
if ("Pep.Impute" %in% names(AnalysisParam)) {
  warning("Parameter \"Pep.Impute\" is currently not used, use \"ImputeMissData\" instead!")
  AnalysisParam$Pep.Impute <- NULL
}
if ("ImputeMissData" %in% names(AnalysisParam)) {
  ImputeMissData <- as.logical(AnalysisParam$ImputeMissData)
  if ((is.na(ImputeMissData))||(is.null(ImputeMissData))||(!moreThan1Exp)) { ImputeMissData <- FALSE }
}
AnalysisParam$ImputeMissData <- ImputeMissData
#
PepFoundInAtLeast %<o% 1
if ("PepFoundInAtLeast" %in% names(AnalysisParam)) {
  PepFoundInAtLeast <- suppressWarnings(as.integer(AnalysisParam$PepFoundInAtLeast))
  if ((is.na(PepFoundInAtLeast))||(PepFoundInAtLeast < 1)||(!moreThan1Exp)||(PepFoundInAtLeast > length(Exp))) {
    PepFoundInAtLeast <- 1
  }
}
AnalysisParam$PepFoundInAtLeast <- PepFoundInAtLeast
#
dbOrd <- 1:nrow(db)
protDflt <- NULL
tmp <- c()
if ("Protein of interest" %in% names(db)) {
  tmp <- unique(c(tmp, db$`Protein ID`[which(db$"Protein of interest")]))
}
if ("Prot.list" %in% names(AnalysisParam)) {
  tmp <- unique(c(tmp, unlist(strsplit(AnalysisParam$Prot.list, ";"))))
}
if ("Prot.list_pep" %in% names(AnalysisParam)) {
  tmp <- unique(c(tmp, unlist(strsplit(AnalysisParam$Prot.list_pep, ";"))))
}
if (length(tmp)) {
  m <- match(tmp, db$`Protein ID`)
  m <- m[which(!is.na(m))]
  dbOrd <- c(m, which(!db$`Protein ID` %in% tmp))
}
protHeads <- gsub("^>", "", db$Header[dbOrd])
if (length(tmp)) { protDflt <- protHeads[1:length(m)] }
#
tstAdvOpt <- try(sum(file.exists(AnalysisParam$Custom.PGs, AnalysisParam$CRAPome_file)) > 0)
if ("try-error" %in% class(tstAdvOpt)) { tstAdvOpt <- FALSE }
if (MakeRatios) {
  RatiosThresh %<o% 1
  if ("Ratios analysis - threshold" %in% names(AnalysisParam)) {
    RatiosThresh <- AnalysisParam$"Ratios analysis - threshold"
  } else {
    AnalysisParam$"Ratios analysis - threshold" <-  RatiosThresh
  }
  if ((is.na(RatiosThresh))||(!is.numeric(RatiosThresh))||(!RatiosThresh < 0)) {
    RatiosThresh <- 1
  }
  RatiosThresh_2sided %<o% c(TRUE, FALSE)[(WorkFlow == "Pull-down")+1]
  if ("Ratios analysis - threshold is two-sided" %in% names(AnalysisParam)) {
    RatiosThresh_2sided <- AnalysisParam$"Ratios analysis - threshold is two-sided"
  } else {
    AnalysisParam$"Ratios analysis - threshold is two-sided" <- RatiosThresh_2sided
  }
  RatiosThresh_2sided <- as.logical(RatiosThresh_2sided)
  if ((is.na(RatiosThresh_2sided))||(is.null(RatiosThresh_2sided))) {
    RatiosThresh_2sided <- TRUE
  }
}
NPep %<o% 2
if ("N. of peptidoforms for quantitation" %in% names(AnalysisParam)) {
  NPep <- AnalysisParam$"N. of peptidoforms for quantitation"
} else { AnalysisParam$"N. of peptidoforms for quantitation" <- NPep }
NPep <- as.integer(NPep)
if ((is.na(NPep))||(NPep < 1)) { NPep <- 2 }
#
removeMBR %<o% FALSE
if ("Proteins list: remove match-between-runs" %in% names(AnalysisParam)) {
  removeMBR <- as.logical(AnalysisParam$"Proteins list: remove match-between-runs")
  if ((is.na(removeMBR))||(is.null(removeMBR))||(!moreThan1Exp)) { removeMBR <- FALSE }
}
AnalysisParam$"Proteins list: remove match-between-runs" <- removeMBR
#
Pep4QuantOpt %<o% c("Unique peptide IDs", "Razor peptide IDs", "Peptide IDs")
Pep4Quant %<o% Pep4QuantOpt[c(1, 2)[isEukaLike+1]]
if (("Peptide classes eligible for quantitation" %in% names(AnalysisParam))&&
    (length(AnalysisParam$"Peptide classes eligible for quantitation"))&&
    (AnalysisParam$"Peptide classes eligible for quantitation" %in% Pep4QuantOpt)) {
  Pep4Quant <- AnalysisParam$"Peptide classes eligible for quantitation"
}
AnalysisParam$"Peptide classes eligible for quantitation" <- Pep4Quant
#
NormalizePG %<o% moreThan1Exp
if ("NormalizePG" %in% colnames(AnalysisParam)) {
  NormalizePG <- as.logical(AnalysisParam$NormalizePG)
  if ((is.na(NormalizePG))||(is.null(NormalizePG))||(!moreThan1Exp)) { NormalizePG <- moreThan1Exp }
}
AnalysisParam$NormalizePG <- NormalizePG
#
# GO enrichment analysis
# - vs total proteome
globalGO_dflt <- (Annotate&(WorkFlow %in% c("Discovery", "Regulation", "Pull-down")))
if ((exists("globalGO"))&&(length(globalGO) == 1)&&(is.logical(globalGO))&&(!is.na(globalGO))) {
  globalGO_dflt <- globalGO
} else { globalGO <- globalGO_dflt }
if ("globalGO" %in% colnames(AnalysisParam)) {
  globalGO_dflt <- globalGO <- as.logical(AnalysisParam$globalGO)
}
if (is.na(globalGO)) { globalGO <- FALSE }
if (is.na(globalGO_dflt)) { globalGO_dflt <- FALSE }
globalGO %<o% globalGO
# - vs reference sample
enrichGO %<o% (globalGO&MakeRatios)
if (Annotate) {
  if (!exists("GO_filter")) { GO_filter %<o% c() }
  if ("GO terms of interest" %in% names(AnalysisParam)) {
    GO_filter <- AnalysisParam$"GO terms of interest"
  } else { AnalysisParam$"GO terms of interest" <- GO_filter }
  GO_filter1 <- allGO[match(GO_filter, allGO2)]
}
#
Mod2Write %<o% c()
if ("Mod2Write" %in% names(AnalysisParam)) {
  Mod2Write <- AnalysisParam$Mod2Write
} else { AnalysisParam$Mod2Write <- Mod2Write }
if ((length(Mod2Write))&&(!Mod2Write %in% Modifs$Mark)) { Mod2Write <- c() }
#
# Venn diagrams
Venn_Obs %<o% moreThan1Exp
if ("Venn diagrams: observed" %in% names(AnalysisParam)) {
  Venn_Obs <- AnalysisParam$"Venn diagrams: observed"
}
Venn_Obs <- moreThan1Exp # No matter the parameters, if there is only one sample, we DO NOT draw Venn diagrams, period
Venn_Obs <- as.logical(Venn_Obs)
if ((is.na(Venn_Obs))||(is.null(Venn_Obs))) { Venn_Obs <- FALSE }
AnalysisParam$"Venn diagrams: observed" <- Venn_Obs
Venn_Obs %<o% Venn_Obs
#
# Cytoscape
CytoScExe %<o% c()
tmp <- grep("cytoscape", list.dirs("C:/PROGRA~1", recursive = FALSE), value = TRUE, ignore.case = TRUE)
CytoScape %<o% (length(tmp) > 0)
if (length(tmp)) {
  CytoScExe <- unlist(lapply(tmp, function(x) { grep("/Cytoscape\\.exe$", list.files(x, recursive = TRUE, full.names = TRUE), value = TRUE) }))
  if (length(CytoScExe) > 1) {
    tst <- sapply(CytoScExe, function(x) { file.info(x)$mtime })
    CytoScExe <- CytoScExe[order(tst, decreasing = TRUE)]
  }
} else {
  msg <- "Could not locate Cytoscape executable!"
  ReportCalls <- AddMsg2Report()
  CytoScape <- FALSE
}
CytoScExe <- CytoScExe[1]
#
# Some default parameters which nicely follow the same structure
myPar <- c("Pepper", "GSEA", "ProfPlots", "RankAbundPlots")
for (parI in myPar) {
  parNm <- paste0("run", parI)
  # Lowest level default: defined by context
  par_dflt <- par_dflt2 <- c(FALSE,
                             (WorkFlow %in% c("Discovery", "Regulation", "Pull-down"))&Annotate,
                             moreThan1Exp,
                             TRUE)[match(parI, myPar)]
  # Middle level default: defined by existing value
  parOK <- (exists(parNm))
  if (parOK) {
    tmpPar <- get(parNm)
    parOK <- (length(tmpPar) == 1)&&(is.logical(tmpPar))&&(!is.na(tmpPar))
  }
  if (parOK) { par_dflt <- tmpPar }
  # Top level default: defined by AnalysisParam
  if (parNm %in% colnames(AnalysisParam)) {
    par_dflt <- as.logical(AnalysisParam[[parNm]])
  }
  # Backup value
  if (is.na(par_dflt)) { par_dflt <- par_dflt2 }
  #
  if (!parOK) { assign(parNm, par_dflt) }
  assign(paste0(parNm, "_dflt"), par_dflt)
  AnalysisParam[[parNm]] <- par_dflt
  .obj <- unique(c(parNm, .obj))
}
#
appNm <- paste0(dtstNm, " - Parameters")
ui <- shiny::fluidPage(
  shinyjs::useShinyjs(),
  shinyWidgets::setBackgroundColor( # Doesn't work
    color = c(#"#F8F8FF",
      "#E6F7F4"),
    gradient = "linear",
    direction = "bottom"
  ),
  shiny::tags$head(
    shiny::tags$style(
      '.inner.open {
            overflow-y: hidden !important;
        }'
    )
  ),
  shinyjs::extendShinyjs(text = jsToggleFS, functions = c("toggleFullScreen")),
  shiny::titlePanel(shiny::tag("u", "Parameters"),
                    appNm),
  shiny::h2(dtstNm), 
  shiny::br(),
  shiny::tags$hr(style = "border-color: black;"),
  shiny::h4("Proteins of interest"),
  shinyWidgets::pickerInput("IntProt", NULL, protHeads, protDflt, TRUE,
                            shinyWidgets::pickerOptions(title = "Search me",
                                                        `live-search` = TRUE,
                                                        actionsBox = TRUE,
                                                        deselectAllText = "Clear search")),
  shiny::br(),
  shiny::tags$hr(style = "border-color: black;"),
  shiny::h4("Data processing"),
  shiny::fluidRow(
    if (moreThan1Exp) {
      shiny::column(2,
                    shiny::checkboxInput("Impute", "Impute missing peptides-level values?",
                                         Impute, "100%"))
    },
    shiny::column(2,
                  shiny::checkboxInput("Update_Prot_matches",
                                       paste0("Update ", names(SearchSoft),
                                              "'s original protein-to-peptides assignments?"),
                                       Update_Prot_matches, "100%"),
                  shinyBS::bsTooltip("Update_Prot_matches",
                                     "Checking assignments may result in removal of some identifications. It is nonetheless recommended because we have observed occasional inconsistent peptides-to-protein assignments with some search software.",
                                     placement = "right", trigger = "hover",
                                     options = list(container = "body"))#,
                  #shinycssloaders::withSpinner(shiny::uiOutput("ReloadMatches"))
    ),
    if (moreThan1Exp) {
      shiny::column(2, shiny::checkboxInput("prtNorm", "Normalize data",
                                            NormalizePG, "100%")) 
    },
    shiny::column(2, shiny::checkboxInput("runPepper", "Run Pepper ML-based PSMs intensity correction?",
                                          runPepper_dflt, "100%"))
  ),
  shiny::tags$hr(style = "border-color: black;"),
  shiny::br(),
  # Quantitation
  ## Choice of algorithm + Proteomics ruler
  shiny::h4("Protein Groups quantitation"),
  shiny::fluidRow(
    #shiny::column(2,
    shiny::column(2,
                  shiny::selectInput("QuantMeth", "Protein Groups-level quantitation algorithm:",
                                     names(QuantMethods)[1:6], QMdefnm, width = "100%")),
    shiny::column(2,
                  shiny::selectInput("Pep4Quant", "Peptides eligible for quantitation (where available, unique peptides will be prioritized)",
                                     Pep4QuantOpt, Pep4Quant, width = "100%")),
    shiny::column(2,
                  shiny::numericInput("NPep", "Min. number of peptidoforms for discovery/quantitation",
                                      NPep, 1, width = "100%")),
    if (moreThan1Exp) {
      shiny::column(2,
                    shiny::numericInput("PepFoundInAtLeast",
                                        "Use only peptidoforms found in at least how many samples?",
                                        PepFoundInAtLeast, 1, length(Exp), 1, "100%"))
    }
  ),
  shiny::br(),
  shiny::fluidRow(
    shiny::column(2,
                  shiny::checkboxInput("ProtRul",
                                       "Apply Proteomic Ruler to estimate copy numbers per cell? (uses signal from all histones as reference; assumes inter-nucleosomal space = 196 bp, do not use if this assumption does not hold!)",
                                       protrul, "100%"),
                  shiny::numericInput("ProtRulNuclL", "Use inter-nucleosome length = ? (kb)",
                                      ProtRulNuclL, 1, Inf, 1, "100%")),
    if (moreThan1Exp) {
      shiny::column(2,
                    shiny::checkboxInput("removeMBR",
                                         "Exclude match-between-runs peptides from coverage analysis?",
                                         removeMBR, "100%")) 
    }
  ),
  shinycssloaders::withSpinner(uiOutput("Ratios")),
  # Note to self: I am for now excluding some methods, because I need to add code to calculate some columns for those, namely ratios.
  # This should be remedied asap, especially since there include such community favourites as IQ (= MaxLFQ) and Top3!!!
  shiny::br(),
  shiny::tags$hr(style = "border-color: black;"),
  shiny::fluidRow(
      shiny::column(2,
                    if (moreThan1Exp) {
                      shiny::checkboxInput("runProfPlots", "Draw protein profile plots?",
                                           runProfPlots_dflt, "100%")
                    },
                    shiny::checkboxInput("runRankAbundPlots", "Draw protein ranked abundance plots?",
                                                     runRankAbundPlots_dflt, "100%")),
      if (moreThan1Exp) {
        shiny::column(2,
                      shiny::radioButtons("Clustering", "Clustering method",
                                          klustChoices, klustChoices[1], TRUE, "100%")) 
      },
      if (Annotate) {
        shiny::column(2, shiny::checkboxInput("runGSEA", "Run Gene Set Enrichment Analysis (GSEA)?",
                                              runGSEA_dflt, "100%"))
      }
  ),
  if (moreThan1Exp) {
    shiny::checkboxInput("Venn_Obs", "Draw Venn diagrams?", Venn_Obs, "100%")
  },
  shiny::br(),
  shiny::tags$hr(style = "border-color: black;"),
  shinycssloaders::withSpinner(uiOutput("GO")),
  shiny::h4("Post-translational modifications (PTMs)"),
  shiny::fluidRow(shiny::column(2,
                                shinyWidgets::pickerInput("PTMsQuant", "Select PTM(s) eligible for use for Protein Groups quantitation:",
                                                          Modifs$`Full name`, ptmDflt1, TRUE,
                                                          shinyWidgets::pickerOptions(title = "Search me",
                                                                                      `live-search` = TRUE,
                                                                                      actionsBox = TRUE,
                                                                                      deselectAllText = "Clear search"))),
                  shiny::column(2, shinyWidgets::pickerInput("Mod2Write", "Select PTM(s) for which to write a specific tab in the report",
                                                             Modifs$`Full name`, ptmDflt2, TRUE,
                                                             shinyWidgets::pickerOptions(title = "Search me",
                                                                                         `live-search` = TRUE,
                                                                                         actionsBox = TRUE,
                                                                                         deselectAllText = "Clear search"))),
                  #column(2, checkboxInput("PTMsReNorm", "Re-normalize modified peptides ratios to those of parent Protein Group(s)?", TRUE, "100%"))
  ),
  shiny::br(),
  shiny::tags$hr(style = "border-color: black;"),
  shiny::checkboxInput("AdvOptOn", "Advanced options", tstAdvOpt),
  shinycssloaders::withSpinner(uiOutput("AdvOpt")),
  shiny::br(),
  shinyWidgets::actionBttn("saveBtn", "Save", icon = icon("save"), color = "success", style = "pill"),
  shiny::br(),
  shiny::br()
)
server <- function(input, output, session) {
  # Initialize variables to create in main environment
  m4Quant <- shiny::reactiveVal(Mod4Quant)
  m2Xclud <- shiny::reactiveVal(Mod2Xclud)
  ADVOPT <- shiny::reactiveVal(tstAdvOpt)
  PARAM <- shiny::reactiveVal(AnalysisParam)
  #
  # Ratios
  output$Ratios <- shiny::renderUI({
    if (MakeRatios) {
      lst <- list(list(h4("Ratios analysis")),
                  list(shiny::fluidRow(shiny::column(2,
                                                     shiny::numericInput("RatiosThresh", "Fold change threshold (log2)",
                                                                         RatiosThresh, 0, width = "100%")),
                                       shiny::column(2,
                                                     shiny::checkboxInput("RatiosThresh_2sided",
                                                                          "Are you interested in down-regulated protein groups as well?",
                                                                          RatiosThresh_2sided, "100%")))))
    } else { lst <- list() }
    return(lst)
  })
  #
  # GO
  output$GO <- shiny::renderUI({
    lst <- list(list(br()))
    if (Annotate) {
      lst <- list(
        list(h4("GO terms enrichment"),
             shiny::fluidRow(shiny::column(1,
                                           shiny::checkboxInput("GOenrich", "GO enrichment", globalGO_dflt, "100%")),
                             shiny::column(2,
                                           shinyWidgets::pickerInput("GO.tabs", "GO terms of interest",
                                                                     allGO, GO_filter1, TRUE,
                                                                     shinyWidgets::pickerOptions(title = "Search me",
                                                                                                 `live-search` = TRUE,
                                                                                                 actionsBox = TRUE,
                                                                                                 deselectAllText = "Clear search"))),
                             shiny::column(1,
                                           shiny::checkboxInput("GO2Int", "Use GO terms to define list of proteins of interest?",
                                                                AnalysisParam$GO.terms.for.proteins.of.interest , "100%"))))
      )
    }
    return(lst)
  })
  #
  # output$ReloadMatches <- shiny::renderUI({
  #   if ("evmatch.RData" %in% list.files(wd)) {
  #     msg <- "Peptide-to-protein matches backup detected in folder: do you want to reload it?\n"
  #     lst <- list(list(list(br()),
  #                      shiny::tags$table(
  #                        shiny::tags$tr(width = "100%", shiny::tags$td(width = "55%", shiny::checkboxInput("Reuse_Prot_matches", msg, TRUE)))
  #                      )))
  #   } else {
  #     shiny::em(" ")
  #   }
  # })
  #
  #
  # Event observers
  # Optional input files
  updtOptOn <- function(reactive = TRUE) {
    if (reactive) { tst <- ADVOPT() } else { tst <- tstAdvOpt }
    if (tst) {
      lst <- vector("list", 1)
      lst[[1]] <- list(shiny::fluidRow(shiny::column(2,
                                                     shiny::em("->"), 
                                                     shinyFiles::shinyFilesButton("CustPG",
                                                                                  shiny::em("Custom Protein Groups"), "", FALSE),
                                                     shiny::br(),
                                                     shiny::em("Allows \"cheating\" with the naive Protein Groups assembly algorithm."), shiny::br(),
                                                     shiny::em("Useful when e.g. samples express from a custom construct to which matching peptides should be assigned in priority."), shiny::br(),
                                                     shiny::em("Should be a table with two columns: \"Leading protein IDs\" (\";\"-separated) and \"Priority\" (integer)."),
                                                     shiny::br(), shiny::br(),
                                                     shiny::em("Current selection = "),
                                                     shiny::span(AnalysisParam$Custom.PGs, style = "color:blue", .noWS = "outside"),
                                                     shiny::br()
      ),
      shiny::column(2,
                    shiny::em("->"),
                    shinyFiles::shinyFilesButton("CRAPome", shiny::em("CRAPome filter"), "", FALSE), shiny::br(),
                    shiny::em("CRAPome-like filter: 1 column table of protein accessions to mark as contaminants."), shiny::br(),
                    shiny::em("Column name = \"Protein ID\" or \"Protein IDs\", use \";\" if including more than one ID per row)."),
                    shiny::br(), shiny::br(),
                    shiny::em("Current selection = "),
                    shiny::span(AnalysisParam$CRAPome_file, style = "color:blue", .noWS = "outside"),
                    shiny::br()
      )
      ))
    } else { lst <- list(list(em(""))) }
    shiny::renderUI(lst)
  }
  output$AdvOpt <- updtOptOn(FALSE)
  shiny::observeEvent(input$AdvOptOn, {
    ADVOPT(input$AdvOptOn)
    output$AdvOpt <- updtOptOn()
  })
  shiny::observe({ shinyFiles::shinyFileChoose(input, "CustPG", roots = shinyFiles::getVolumes(), filetypes = "csv")
    {
      tmp <- input$CustPG
      if ((!is.null(tmp))&&(is.list(tmp))) {
        tmp <- shinyFiles::parseFilePaths(shinyFiles::getVolumes(), tmp)$datapath
        Par <- PARAM()
        Par$Custom.PGs <- normalizePath(tmp, winslash = "/")
        PARAM(Par)
      }
  }
  })
  shiny::observe({ shinyFiles::shinyFileChoose(input, "CRAPome", roots = shinyFiles::getVolumes(), filetypes = "csv" )
    {
      tmp <- input$CRAPome
      if ((!is.null(tmp))&&(is.list(tmp))) {
        tmp <- shinyFiles::parseFilePaths(getVolumes(), tmp)$datapath
        Par <- PARAM()
        Par$CRAPome_file <- normalizePath(tmp, winslash = "/")
        PARAM(Par)
      }
  }
  })
  #
  # Pep4Quant
  shiny::observeEvent(input$Pep4Quant, {
    Pep4Quant <<- input$Pep4Quant
    Par <- PARAM()
    Par$"Peptide classes eligible for quantitation" <- Pep4Quant
    PARAM(Par)
  })
  #
  shiny::observeEvent(input$RatiosThresh, {
    RatiosThresh <<- input$RatiosThresh
    Par <- PARAM()
    Par$"Ratios analysis - threshold" <- RatiosThresh
    PARAM(Par)
  })
  shiny::observeEvent(input$RatiosThresh_2sided, {
    RatiosThresh_2sided <<- input$RatiosThresh_2sided
    Par <- PARAM()
    Par$"Ratios analysis - threshold is two-sided" <- RatiosThresh_2sided
    PARAM(Par)
  })
  # Impute?
  shiny::observeEvent(input$Impute, {
    ImputeMissData <<- input$Impute
    Par <- PARAM()
    Par$ImputeMissData <- ImputeMissData
    PARAM(Par)
  })
  # Update PSM-to-Protein matches?
  shiny::observeEvent(input$Update_Prot_matches, {
    Update_Prot_matches <<- input$Update_Prot_matches
    Par <- PARAM()
    Par$Update_Prot_matches <- Update_Prot_matches
    PARAM(Par)
    # if (input$Update_Prot_matches) { shinyjs::enable("Reuse_Prot_matches") }
    # if (!input$Update_Prot_matches) { shinyjs::disable("Reuse_Prot_matches") }
  })
  # observeEvent(input[["Reuse_Prot_matches"]], {
  #   Reuse_Prot_matches <<- input$Reuse_Prot_matches
  #   Par <- PARAM()
  #   Par$Reuse_Prot_matches <- Reuse_Prot_matches
  #   PARAM(Par)
  # })
  # Clustering method
  shiny::observeEvent(input$Clustering, {
    assign("KlustMeth", match(input$Clustering, klustChoices), envir = .GlobalEnv)
  })
  # Are analyses Two-sided?
  shiny::observeEvent(input$TwoSided, {
    Par <- PARAM()
    Par$Two.sided <- input$TwoSided == "Both directions"
    PARAM(Par)
  })
  # Normalizations
  shiny::observeEvent(input$prtNorm, {
    assign("NormalizePG", as.logical(input$prtNorm), envir = .GlobalEnv)
    Par <- PARAM()
    Par$NormalizePG <- NormalizePG
    PARAM(Par)
  })
  # Pepper
  shiny::observeEvent(input$runPepper, {
    assign("runPepper", as.logical(input$runPepper), envir = .GlobalEnv)
    Par <- PARAM()
    Par$runPepper <- runPepper
    PARAM(Par)
  })
  # GSEA
  shiny::observeEvent(input$runGSEA, {
    assign("runGSEA", as.logical(input$runGSEA), envir = .GlobalEnv)
    Par <- PARAM()
    Par$runGSEA <- runGSEA
    PARAM(Par)
  })
  # Profile plots
  shiny::observeEvent(input$runProfPlots, {
    assign("runProfPlots", as.logical(input$runProfPlots), envir = .GlobalEnv)
    Par <- PARAM()
    Par$runProfPlots <- runProfPlots
    PARAM(Par)
  })
  # Ranked abundance plots
  shiny::observeEvent(input$runRankAbundPlots, {
    assign("runRankAbundPlots", as.logical(input$runRankAbundPlots), envir = .GlobalEnv)
    Par <- PARAM()
    Par$runRankAbundPlots <- runRankAbundPlots
    PARAM(Par)
  })
  # Quantitation
  # observeEvent(input$QuantMeth, {
  #   Par <- PARAM()
  #   Par$QuantMeth <- QuantMethods[match(input$QuantMeth, names(QuantMethods))]
  #   PARAM(Par)
  # })
  shiny::observeEvent(input$ProtRul, {
    assign("protrul", input$ProtRul, envir = .GlobalEnv)
    Par <- PARAM()
    Par$"Proteomic ruler calculated" <- protrul
    PARAM(Par)
  })
  shiny::observeEvent(input$ProtRulNuclL, {
    assign("ProtRulNuclL", as.integer(input$ProtRulNuclL), envir = .GlobalEnv)
    Par <- PARAM()
    Par$ProtRulNuclL <- ProtRulNuclL
    PARAM(Par)
  })
  shiny::observeEvent(input$removeMBR, {
    assign("removeMBR", input$ProtRul, envir = .GlobalEnv)
    Par <- PARAM()
    Par$"Proteins list: remove match-between-runs" <- removeMBR
    PARAM(Par)
  })
  shiny::observeEvent(input$NPep, {
    assign("NPep", input$NPep, envir = .GlobalEnv)
    Par <- PARAM()
    Par$"N. of peptidoforms for quantitation" <- NPep
    PARAM(Par)
  })
  shiny::observeEvent(input$PepFoundInAtLeast, {
    assign("PepFoundInAtLeast", as.integer(input$PepFoundInAtLeast), envir = .GlobalEnv)
    Par <- PARAM()
    Par$PepFoundInAtLeast <- PepFoundInAtLeast
    PARAM(Par)
  })
  # Proteins of interest
  shiny::observeEvent(input$IntProt, {
    assign("prot.list", db$`Protein ID`[dbOrd][match(input$IntProt, protHeads)], envir = .GlobalEnv)
    Par <- PARAM()
    Par$Prot.list <- prot.list
    PARAM(Par)
  }, ignoreNULL = FALSE)
  # GO enrichment
  shiny::observeEvent(input$GOenrich, {
    assign("globalGO", as.logical(input$GOenrich), envir = .GlobalEnv)
    assign("enrichGO", globalGO&MakeRatios, envir = .GlobalEnv)
    Par <- PARAM()
    Par$"GO terms enrichment analysis" <- globalGO
    PARAM(Par)
  })
  # GO terms of interest
  if (Annotate) {
    shiny::observeEvent(input$GO.tabs, {
      assign("GO_filter", allGO2[match(input$GO.tabs, allGO)], envir = .GlobalEnv)
      if (length(GO_filter1)) {
        mY <- match(GO_filter1, allGO)
        mN <- which(!allGO %in% GO_filter1)
        assign("allGO", allGO[c(mY, mN)], envir = .GlobalEnv)
        assign("allGO2", allGO2[c(mY, mN)], envir = .GlobalEnv)
      }
      Par <- PARAM()
      Par$"GO terms of interest" <- GO_filter
      PARAM(Par)
    }, ignoreNULL = FALSE)
    shiny::observeEvent(input$GO2Int, {
      Par <- PARAM()
      Par$GO.terms.for.proteins.of.interest <- input$GO2Int
      PARAM(Par)
    })
  }
  #
  shiny::observeEvent(input$Venn_Obs, {
    assign("Venn_Obs", input$Venn_Obs, envir = .GlobalEnv)
    Par <- PARAM()
    Par$"Venn diagrams: observed" <- Venn_Obs
    PARAM(Par)
  })
  # PTMs to use for PG Quant
  shiny::observeEvent(input$PTMsQuant, {
    m4Quant(Modifs$Mark[match(unlist(input$PTMsQuant), Modifs$`Full name`)])
    m2Xclud(set_colnames(Modifs[which(!Modifs$Mark %in% Mod4Quant), c("Mark", "AA")],
                         c("Mark", "Where")))
  }, ignoreNULL = FALSE)
  # PTMs to write a tab for
  shiny::observeEvent(input$Mod2Write, {
    assign("Mod2Write", Modifs$Mark[match(input$Mod2Write, Modifs$`Full name`)], envir = .GlobalEnv)
    Par <- PARAM()
    Par$Mod2Write <- Mod2Write
    PARAM(Par)
  }, ignoreNULL = FALSE)
  # Re-normalize PTM peptides
  # observeEvent(input$PTMsReNorm, {
  #   Par <- PARAM()
  #   Par$PTM.analysis_Norm <- input$PTMsReNorm
  #   PARAM(Par)
  # })
  #
  # Save
  shiny::observeEvent(input$saveBtn, {
    Par <- PARAM()
    assign("AnalysisParam", Par, envir = .GlobalEnv)
    assign("Mod4Quant", m4Quant(), envir = .GlobalEnv)
    assign("Mod2Xclud", m2Xclud(), envir = .GlobalEnv)
    stopApp()
  })
  #observeEvent(input$cancel, { stopApp() })
  session$onSessionEnded(function() { shiny::stopApp() })
}
eval(parse(text = runApp), envir = .GlobalEnv)
#
# Post-processing
prot.list %<o% AnalysisParam$Prot.list
if (AnalysisParam$GO.terms.for.proteins.of.interest) {
  tmpGO <- AnalysisParam$"GO terms of interest"
  GO_prot.list <- list()
  GO_prot.list$Offspring <- lapply(tmpGO, function(x) {
    ont <- Ontology(x)
    x <- c(x, get(paste0("GO", ont, "OFFSPRING"))[[x]])
    x <- x[which(!is.na(x))]
    return(x)
  })
  tmpGO2 <- listMelt(strsplit(db$`GO-ID`, ";"), db$`Protein ID`)
  GO_prot.list$Proteins <- lapply(GO_prot.list$Offspring, function(x) {
    unique(tmpGO2$L1[which(tmpGO2$value %in% unlist(x))])
  })
  prot.list <- unique(c(prot.list, unlist(GO_prot.list$Proteins)))
}
prot.list.Cond %<o% (length(prot.list) > 0)
if (prot.list.Cond) {
  m <- match(prot.list, db$`Protein ID`)
  x <- aggregate(db$`Protein ID`[m], list(db$`Common Name`[m]), c)
  IDs.list %<o% setNames(as.list(x$x), x$Group.1)
  prot.names %<o% names(IDs.list)
  db$"Potential contaminant"[which(db$`Protein ID` %in% prot.list)] <- ""
}
#
custPGs_file %<o% AnalysisParam$Custom.PGs
custPGsTst <- (!is.na(custPGs_file))&(file.exists(custPGs_file))
if (custPGsTst) {
  if (dirname(custPGs_file) != wd) {
    file.copy(custPGs_file, wd)
    custPGs_file <- paste0(wd, "/", basename(custPGs_file))
  }
  custPGs %<o% read.delim(custPGs_file, check.names = FALSE)
  if (colnames(custPGs)[1] != "Leading protein IDs") { custPGs <- read.delim(custPGs_file, check.names = FALSE, sep = ",") }
  if (colnames(custPGs)[1] != "Leading protein IDs") { custPGs <- read.csv(custPGs_file, check.names = FALSE) }
  if (colnames(custPGs)[1] != "Leading protein IDs") { custPGs <- read.csv(custPGs_file, check.names = FALSE, sep = "\t") }
  if (colnames(custPGs)[1] != "Leading protein IDs") {
    warning("I could not make sense of that file! Skipping.")
    custPGs <- NA
  }
} else { custPGs %<o% NA }
#
if (prot.list.Cond) {
  temp <- db[which(db$`Protein ID` %in% prot.list),]
  writeFasta(temp, paste0(wd, "/Proteins of interest.fasta"))
  AnalysisParam$"Proteins list: proteins" <- IDs.list
  AnalysisParam$"Proteins list: names" <- prot.names
}

# GO terms of interest (for profile plots, LFQ plots and similar, GO-specific tabs...)
GO_filt %<o% FALSE
if (exists("GO_filter")) { GO_filt <- length(GO_filter) > 0 }
if (GO_filt) {
  library(GO.db)
  AllTerms %<o% unique(unlist(strsplit(db$`GO-ID`, ";")))
  AllTermNames %<o% unique(unlist(strsplit(db$GO, ";")))
}

# Venn diagrams
AnalysisParam$"Venn diagrams: observed" <- Venn_Obs
Venn_Ratios %<o% (Venn_Obs & MakeRatios)
AnalysisParam$"Venn diagrams: regulated (up)" <- Venn_Ratios
AnalysisParam$"Venn diagrams: regulated (down)" <- Venn_Ratios & RatiosThresh_2sided
#
AnalysisParam$"PTMs eligible for quantitation" <- Mod4Quant
AnalysisParam$"PTMs to exclude from quantitation" <- Mod2Xclud

# Experiments
SamplesMapPath %<o% paste0(wd, "/SamplesMap.csv")
tst <- (("Experiment" %in% colnames(ev))&&(sum(is.na(ev$Experiment)) == 0)) 
if (tst) {
  if ("Reference" %in% colnames(SamplesMap)) {
    SamplesMap <- SamplesMap[which(!is.na(SamplesMap$Reference)),]
  }
  Exp %<o% SamplesMap$Experiment
  w1 <- which(ev$Experiment %in% Exp)
  w2 <- which(!ev$Experiment %in% Exp)
  if (length(w2)) {
    warning(paste0("Removing ", length(w2), " evidences from undefined experiments (check \"Reference\" column of Samples map)..."))
    ev <- ev[w1,]
  }
} else {
  Exp %<o% "Exp1"
  ev$Experiment <- Exp
  SamplesMap <- data.frame("Experiment" = Exp,
                           "Ratios group" = 1,
                           "Reference" = TRUE,
                           "Negative Filter" = FALSE,
                           "Use" = TRUE,
                           check.names = FALSE)
  tst <- try(write.csv(SamplesMap, file = SamplesMapPath, row.names = FALSE), silent = TRUE)
  while ("try-error" %in% class(tst)) {
    dlg_message(paste0("File \"", SamplesMapPath, "\" appears to be locked for editing, close the file then click ok..."), "ok")
    tst <- try(write.csv(SamplesMap, file = SamplesMapPath, row.names = FALSE), silent = TRUE)
  }
}
if (identical(c("MQ.Exp", "Experiment") %in% colnames(SamplesMap), c(TRUE, FALSE))) {
  # Here we use Experiment as synonym for MQ.Exp
  colnames(SamplesMap)[which(colnames(SamplesMap) == "MQ.Exp")] <- "Experiment"
}
