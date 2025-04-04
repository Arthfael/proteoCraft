#' MatMet_WetLab
#'
#' @description
#' A function to handle writing the wet-lab part of Materials & Methods, interactively asking relevant questions as it moves down the decision tree.
#' Very much a Rube Goldberg engine, which is fine. Just remember that it creates a WIP which should be checked/corrected before publication!
#' 
#' !!! When editing function, always save with UTF-8 encoding!!!

#' @param File2Reload Optional file to reload the text from so you do not have to rewrite the whole thing. Default = "Materials and methods_WIP.docx".
#' @param exp.map The experiments map.
#' @param exp.map.col Name of the column in exp.map which contains biological sample names. Default = "MQ.Exp".
#' @param Label Were the samples labelled isobarically ("TMT", "TMTPRO" or "ITRAQ"), with "SILAC", or label-free ("LFQ")? This is only relevant for isobaric labels (otherwise ignored) because a) this affects the choice of iST kit, and b) the labelling procedure will then be reported in the samples processing.
#' @param Columns Path to a table of LC columns to pick from.
#'
#' @examples
#' if (LabelType == "LFQ") { MatMet_Text1 <- MatMet_WetLab() }
#' if (LabelType == "Isobaric") { MatMet_Text1 <- MatMet_WetLab(Label = IsobarLab) }
#' 
#' @export

MatMet_WetLab <- function(File2Reload = "Materials and methods_WIP.docx",
                          exp.map = Exp.map,
                          exp.map.col = "MQ.Exp",
                          Label = "LFQ",
                          Columns = "default") {
  #rm(list = ls()[which(!ls() %in% c("wd", "Exp.map"))])
  #proteoCraft::DefArg(proteoCraft::"MatMet_WetLab")
  if (!missing("Label")) { Label <- toupper(gsub("-", "", Label)) }
  # LC columns
  path <- paste0(normalizePath(Sys.getenv("HOME"), winslash = "/"), "/R")
  colTst <- (Columns != "default")&(file.exists(Columns))
  if (!colTst) { Columns <- paste0(path, "/proteoCraft/LC_columns.xlsx") }
  if (colTst) {
    allKolumns <- openxlsx::read.xlsx(Columns, check.names = FALSE)
  } else {
    Kolumns <- c("EasySpray C18 column (2 µm particle size, 75 µm ID x 50 cm length, ThermoFisher Scientific P/N ES903)",
                 "200 cm C18 µPAC column (micro-Pillar Array Column, PharmaFluidics P/N 5525031518210B)",
                 "50 cm C18 µPAC GEN2 column (2nd generation micro-Pillar Array Column, PharmaFluidics P/N COL-nano050G2B)")
    preKolumns <- c("Acclaim PepMap C18 pre-column (5 µm particle size, 0.3 mm ID x 5 mm length, ThermoFisher Scientific P/N 160454)",
                    "C18 µPAC trapping-column (PharmaFluidics P/N 55250200018001)")
    allKolumns <- data.frame(Name = c(Kolumns, preKolumns),
                             Material = "C18",
                             "Length.(cm)" = c(50, 200, 50, 0.5, NA),
                             "ID.(µm)" = c(75, NA, NA, 300, NA),
                             "Particles.size.(µm)" = c(2, NA, NA, 5, NA),
                             "Vendor" = c("ThermoFisher Scientific", "PharmaFluidics", "PharmaFluidics", "ThermoFisher Scientific", "PharmaFluidics"),
                             "P/N" = c("ES903", "5525031518210B", "COL-nano050G2B", "160454", "55250200018001"),
                             Function = c("Analytical", "Analytical", "Analytical", "Trap",  "Trap"),
                             check.names = FALSE)
    Columns <- paste0(path, "/proteoCraft/LC_columns.xlsx")
  }
  if (!"Description" %in% colnames(allKolumns)) { allKolumns$Description <- "" }
  kolDescr <- function(Colonnes) {
    apply(Colonnes[, c("Name", "Material", "Length.(cm)", "ID.(µm)", "Particles.size.(µm)", "Vendor", "P/N"), drop = FALSE], 1, function(x) {
      #x <- Colonnes[1, c("Name", "Material", "Length.(cm)", "ID.(µm)", "Particles.size.(µm)", "Vendor", "P/N")]
      x <- gsub("^ +| +$", "", as.character(unlist(x)))
      dimz <- c(x[3], x[4])
      w <- which(!is.na(dimz))
      dimzTst <- length(w) > 0
      if (dimzTst) { dimz <- paste(paste0(gsub("\\.0+$", "", dimz), c(" cm", " µm ID"))[w], collapse = " * ") } else { dimz <- NA }
      pn <- c(x[6], x[7])
      w <- which(!is.na(pn))
      pnTst <- length(w) > 0
      if (pnTst) { pn <- paste(paste0(c("", "P/N "), pn)[w], collapse = " ") } else { pn <- NA }
      pnTst <- !is.na(x[7])
      if (pnTst) { x[7] <- paste0("P/N ", x[7]) }
      prtcls <- c(x[5], x[2])
      w <- which(!is.na(prtcls))
      prtclsTst <- length(w) > 0
      if (prtclsTst) { prtcls <- paste(paste0(prtcls, c(" µm", "-coated particles"))[w], collapse = " ") } else { prtcls <- NA }
      res <- paste0(x[1], " column (", paste(c(prtcls, dimz, pn)[which(c(prtclsTst, dimzTst, pnTst))], collapse = ", "), ")")
      res <- gsub("column column", "column", res, ignore.case = TRUE)
      return(res)
    })
  }
  #w <- which((nchar(allKolumns$Description) == 0)|(is.na(allKolumns$Description)))
  w <- 1:nrow(allKolumns)
  allKolumns$Description[w] <- kolDescr(allKolumns[w,])
  Kolumns <- unique(c(allKolumns$Description[which(allKolumns$Function == "Fractionation")], "Add new..."))
  AddKolKount  <- 0
  #
  N <- nrow(exp.map)
  if ("Use" %in% colnames(exp.map)) { N <- sum(as.logical(exp.map$Use)) }
  reuseMatMeth <- FALSE
  if ((!is.null(File2Reload))&&(file.exists(File2Reload))) {
    if (dirname(File2Reload) == ".") { File2Reload <- paste0(getwd(), "/", File2Reload) }
    reuseMatMeth <- TRUE
    MatMetTxt <- officer::read_docx(File2Reload)
    MatMetTxt <- officer::docx_summary(MatMetTxt)
    MatMetTxt <- MatMetTxt$text
    MatMetTxt <- MatMetTxt[match("Samples preparation", MatMetTxt) + 1]
    reuseMatMeth <- MatMetTxt != "TEMPLATE"
  }
  if (!reuseMatMeth) {
    SPMethods <- c(paste0("iST", rep(c("", "-NHS"), 3), " kit", c(rep("", 2),
                                                                  paste0(" adapted for ", c(rep("on-paramagnetic beads digest", 2),
                                                                                            rep("tissue samples", 2))))),
                   "FASP", paste0(c("on-beads", paste0("in-", c("gel", "solution"))), " digest"), "none of those")
    SPMeth <- svDialogs::dlg_list(SPMethods,
                                  title = "Select the sample processing method used:")$res
    if (SPMeth != "none of those") {
      SP3 <- FALSE
      PreOmics <- grepl("iST", SPMeth)
      onBeads <- grepl("( |-)beads?( |\\.)?", SPMeth) # The SP3 kit is not compatible with other on-beads digest
      if (grepl("SP3", SPMeth)) { SP3 <- TRUE } else {
        if ((PreOmics)&&(!onBeads)) {
          msg <- "Did you use the SP3 add-on kit?"
          SP3 <- c(TRUE, FALSE)[match(svDialogs::dlg_message(msg, "yesno")$res,
                                      c("yes", "no"))]
        }
      }
      if (SP3) {
        SP3Beads <- svDialogs::dlg_input("Which amount of beads (mg) was used for SP3?",
                                         50)$res
      }
      if (!PreOmics) {
        Reduct <- svDialogs::dlg_list(c("TCEP", "DTT", "β-Mercaptoethanol", "Other"), "TCEP",
                                      title = "How were the samples reduced?")$res
        if (Reduct == "Other") { Reduct <- svDialogs::dlg_input("Enter the name of the reducing agent used:")$res }
        ReductC <- svDialogs::dlg_input("Enter the concentration of the reducing agent used:",
                                        "25 mM")$res
        ReductT <- svDialogs::dlg_input("Enter the time (min) for which the samples were reduced:",
                             30)$res
        ReductTmp <- svDialogs::dlg_input("Enter the temperature (°C) at which the samples were reduced:",
                                          95)$res
        Alk <- svDialogs::dlg_list(c(alkTbl$Full_name, "Other"),
                                   alkTbl$Full_name[1], title = "Which alkylating agent was used?")$res
        if (Alk == "Other") {
          Alk <- svDialogs::dlg_input("Enter the name of the alkylating agent used:")$res
        }
        AlkC <- svDialogs::dlg_input("Enter the concentration of the alkylating agent used:",
                                     "50 mM")$res
        AlkT <- svDialogs::dlg_input("Enter the time (min) for which the samples were alkylated:",
                                     30)$res
      }
      if (SPMeth == "in-gel digest") {
        DestainMeth <- ""
        msg <- "How was the gel stained?"
        StainMethods <- c("InstantBlue", "Other MS-compatible Coomassie blue variation", "MS-incompatible Coomassie blue staining", "Silver nitrate", "SYPRO-Ruby", "Other", "None")
        StainMeth <- svDialogs::dlg_list(StainMethods, StainMethods[1],
                                         title = msg)$res
        if (StainMeth == "Other MS-compatible Coomassie blue variation") {
          StainMeth <- svDialogs::dlg_input("Enter the name of the Coomassie variation used:",
                                            " MS-compatible Coomassie blue")$res
        }
        if (StainMeth == "Other") {
          StainMeth <- svDialogs::dlg_input("Enter the name of the staining method used:")$res
        }
        Destain <- FALSE
        if (StainMeth != "None") {
          Destain <- svDialogs::dlg_list(c("yes, our standard protocol",
                                           "yes, manufacturer's protocol",
                                           "no"),
                              title = "Were the samples destained?")$res
          DestainMethDflt <- "destained by washing x3 in destaining buffer (900 µL 50% acetonitrile = ACN, 50 mM triethylammonium bicarbonate = TEAB) for 15 min with shaking"
          if (Destain == "yes, our standard protocol") {
            DestainMeth <- DestainMethDflt
          }
          if (Destain == "yes, manufacturer's protocol") {
            DestainMeth <- svDialogs::dlg_input("Describe how the samples were destained:",
                                                DestainMethDflt)
          }
          DestainMeth <- paste0(", ", DestainMeth)
        }
      }
      if (grepl("iST", SPMeth)) { Dig <- "Trypsin" } else {
        Dig <- svDialogs::dlg_input("Which enzyme(s) was used for digestion?\nIf multiple, write details as if writing a methods text (e.g. \"Trypsin+Lys-C, Lys-C alone or Chymotrypsin\")",
                                    "Trypsin")$res
      }
      DigT <- svDialogs::dlg_input(paste0("Enter the time (h) for which the samples were ",
                                          Dig, "-digested (or \"ON\" for overnight)"),
                                   c("ON", 3)[PreOmics+1])$res
      if (!is.na(suppressWarnings(as.character(as.numeric(DigT))))) { DigT <- paste0(DigT, " h") }
      if (toupper(DigT) == "ON") { DigT <- "overnight" }
      #LablMethods <- c("Label-free", "SILAC")
      #if (SPMeth != "in-gel digest") {
      #  LablMethods <- c(LablMethods, paste0("TMT", c(paste0("-", c(6, 10, 11), "plex"), paste0("pro-", c(16, 18), "plex"))), "iTRAQ")
      #}
      #
      if (Label %in% c("TMT", "TMTPRO", "ITRAQ")) {
        if (Label == "TMT") { LablMethods <- paste0("TMT", c(paste0("-", c(6, 10, 11), "plex"), paste0("pro-", c(16, 18), "plex"))) }
        if (Label == "TMTPRO") { LablMethods <- paste0("TMTpro-", c(16, 18), "plex") }
        if (Label == "ITRAQ") { LablMethods <- paste0("iTRAQ-", c(4, 8), "plex") }
        LablMeth <- svDialogs::dlg_list(LablMethods,
                                        title = "Select the precise labelling method used:")$res
        LablMeth <- paste0(LablMeth, " (", c(rep("ThermoFisher Scientific", 2), "Sciex")[match(Label, c("TMT", "TMTPRO", "ITRAQ"))], ")")
      }
      test <- (Label %in% c("TMT", "TMTPRO", "ITRAQ"))&(grepl("iST kit", SPMeth))
      while (test) {
        msg <- paste0("The classic iST kit is not compatible with ", Label, " labelling, did you mean to select an iST-NHS variant instead?")
        tst1 <- c(TRUE, FALSE)[match(svDialogs::dlg_message(msg, "yesno")$res,
                                     c("yes", "no"))]
        if (tst1) {
          samplprep <- grep("^iST-NHS ", SPMethods, value = TRUE)
          SPMeth <- svDialogs::dlg_list(samplprep,
                                        title = "Select the sample processing method used:")$res
        } else { stop() }
        test <- (Label %in% c("TMT", "TMTPRO", "ITRAQ"))&(grepl("iST kit", SPMeth))
      }
      EnrichMethods <- c("Not enriched", "Phospho-enrichment (Ti-IMAC)",
                         paste0(c("acetylated", "ubiquitin (GlyGly) remnant"), " peptides enrichment"))
      EnrichMeth <- svDialogs::dlg_list(EnrichMethods, EnrichMethods[1],
                                        title = "Select the peptides-enrichment method(s) used, if and as applicable:")$res
      Enrich <- EnrichMeth != "Not enriched"
      EnrichSepPak <- FALSE
      if (Enrich) {
        EnrichSepPak <- c(TRUE, FALSE)[match(svDialogs::dlg_message("Was enrichment followed by Sep-Pak clean-up?",
                                                                    "yesno")$res, c("yes", "no"))]
      }
      if (SPMeth != "in-gel digest") {
        FracMethods <- c("Not fractionated",
                         "iST-Fractionation Add-on (PreOmics)",
                         "High pH Reversed Phase chromatography",
                         "High pH Reversed Phase fractionation kit (Pierce)")
        FracMeth <- svDialogs::dlg_list(FracMethods, FracMethods[1],
                                        title = "Select the fractionation method used, as applicable:")$res
        if (FracMeth %in% FracMethods[c(2, 4)]) {
          FracMeth <- paste0("the ", FracMeth)
        }
        FracMethods[c(2, 4)] <- paste0("the ", FracMethods[c(2, 4)])
        Frac <- FracMeth != FracMethods[1]
        if (Frac) {
          if (Enrich) {
            FracSources <- c("flow-through from enrichment", "whole sample")
            FracSource <- svDialogs::dlg_list(FracSources, FracSources[1],
                                              title = "What type of input material was used for fractionation?")$res
            if (FracSource == "whole sample") {
              EnrichPerc <- svDialogs::dlg_input("Enter the percentage of the sample used for enrichment", "50%")$res
            }
          } else { FracSource <- "sample" }
          msg <- paste0("How many fractions were run per sample?")
          def <- c(3, 12, 8)[match(FracMeth, FracMethods[2:4])] + 1*((FracMeth != FracMethods[3])&(Label %in% c("TMT", "TMTPRO")))
          NFrac <- suppressWarnings(as.integer(svDialogs::dlg_input(msg, def)$res))
          while ((is.na(NFrac))|(NFrac <= 0)) {
            NFrac <- suppressWarnings(as.integer(svDialogs::dlg_input(msg, def)$res))
          }
          if (FracMeth == FracMethods[3]) {
            FracCol <- svDialogs::dlg_list(Kolumns, Kolumns[1],
                                           title = paste0("Which column was used for offline fractionation?"))$res
            if (FracCol == "Add new...") {
              kol <- svDialogs::dlg_input("Fill in new column's characteristics (\"|\"-separated):",
                                          "Name | Material | Length (cm) | ID (µm) | Particles size (µm) | Vendor | P/N")$res
              kol <- unlist(strsplit(kol, " *\\| *"))
              l <- length(kol)
              if (l < 7) { kol <- c(kol, rep(NA, 7-l)) }
              kol <- kol[1:7]
              tmp <- as.data.frame(t(as.data.frame(kol)))
              colnames(tmp) <- c("Name", "Material", "Length.(cm)", "ID.(µm)", "Particles.size.(µm)", "Vendor", "P/N")
              tmp$Function <- "Fractionation"
              tmp$Description <- kolDescr(tmp)
              k <- colnames(allKolumns)
              k <- k[which(!k %in% colnames(tmp))]
              tmp[,k] <- NA
              tmp[1, which(tmp[1,] == "NA")] <- NA
              allKolumns <- rbind(allKolumns, tmp)
              rownames(allKolumns) <- NULL
              tst <- apply(allKolumns, 1, paste, collapse = "|")
              tst <- aggregate(1:length(tst), list(tst), min)
              allKolumns <- allKolumns[sort(tst$x),]
              require(openxlsx)
              Head <- createStyle(textDecoration = c("bold", "underline"))
              wb <- createWorkbook()
              addWorksheet(wb, "Sheet1")
              writeData(wb, "Sheet1", allKolumns, 1, 1)
              addStyle(wb, "Sheet1", Head, 1, c(1:ncol(allKolumns)))
              tst <- sapply(1:ncol(allKolumns), function(x) {
                x1 <- ceiling(nchar(colnames(allKolumns)[x])*1.2)
                x2 <- ceiling(nchar(allKolumns[[x]])*1.2)
                return(max(c(x1, min(c(60, x2), na.rm = TRUE))))
              })
              setColWidths(wb, "Sheet1", 1:ncol(allKolumns), tst)
              #saveWorkbook(wb, gsub("\\.xlsx$", "2.xlsx", Columns, ""), TRUE)
              try(saveWorkbook(wb, Columns, TRUE), silent = TRUE) # We do not want the function to fail if this fails, it isn't worth the trouble
              #proteoCraft::openwd(dirname(Columns))
            }
            FracCol <- paste0("a", c("", "n")[(tolower(substr(FracCol, 1, 1)) %in% c("a", "e", "i", "o", "u"))+1], " ", FracCol)
            FracMPA <- svDialogs::dlg_input("Enter the composition of solvent A?",
                                            "de-ionized water + 10 mM NH4OH")$res
            FracMPB <- svDialogs::dlg_input("Enter the composition of solvent B?",
                                            "90% LC-grade Acetonitrile + 10 mM NH4OH")$res
            FracGrad <- svDialogs::dlg_input("Describe the gradient's stages:",
                                             "Flow: 0.15 mL/min; 0-4 min: 1% B, 115 min: 25%, 140 min: 40%, 148 min: 75%, maintained for 12 min, followed by 45 min equilibration at 1% B")$res
            NOrigFrac <- 24
            dflTxt <- paste0(NOrigFrac, " fractions of ", round((164-4)*0.15/NOrigFrac, 3),
                             " mL were collected from 4 to 164 min.")
            tst <- NOrigFrac/NFrac
            if (round(tst) == tst) {
              if (tst >= 2) {
                if (tst %% 2 == 0) {
                  Klps <- function(range) {
                    Mn <- min(range)
                    Mx <- max(range)
                    if (Mn == Mx) { return(as.character(Mn)) } else {
                      if (Mn < Mx-1) { coll <- ":" } else { coll <- "-" }
                      return(paste(as.character(c(Mn, Mx)), collapse = coll))
                    }
                  }
                  sq <- 1:(tst/2)
                  repl <- paste0(" and combined ", tst/2, "-by-", tst/2, " at mid-gradient (", Klps(sq), " with ", Klps(sq+NOrigFrac/2), "___",
                                 Klps(sq+tst/2), " with ", Klps(sq+NOrigFrac/2+tst/2), ").")
                  if (max(sq+NOrigFrac/2+tst/2) == NOrigFrac) {
                    repl <- gsub("___", " and ", repl)
                  } else {
                    repl <- gsub("\\)\\.$", paste0("... ", Klps(NOrigFrac/2-rev(sq)+1),  " with ", Klps(NOrigFrac-rev(sq)+1),  ")."), gsub("___", ", ", repl))
                  }
                }
              } else { repl <- paste0(" by ", tst, " into ", NFrac, ".") }
              dflTxt <- gsub("\\.", repl, dflTxt)
            }
            FracCollect <- svDialogs::dlg_input("Describe the fraction collection scheme",
                                                dflTxt)$res
          }
        }
      }
      SepPak <- FALSE
      if (!grepl("iST", SPMeth)) {
        SepPak <- c(TRUE, FALSE)[match(svDialogs::dlg_message("Were the samples Sep-Pak cleaned-up prior to analysis?",
                                                              "yesno")$res, c("yes", "no"))]
      }
      if (PreOmics) {
        tst <- grepl("kit adapted.+$", SPMeth)+1
        kit <- gsub(" kit .+$", " kit", SPMeth)
        protocol <- gsub("iST(-NHS)? kit( adapted)?", "", SPMeth)
        if (SP3) {
          insrt <- paste0(" first cleaned up by SP3 using a commercial kit (PreOmics GmbH, ",
                          SP3Beads, " mg of beads per sample), then")
        } else { insrt <- "" }
        Txt <- paste0(c("The", "All")[(N>1)+1], " sample",  c(" was", "s were")[(N>1)+1], insrt, " processed using the ",
                      kit, " (PreOmics GmbH) ", c("according to the manufacturer's instructions",
                                                  "using the manufacturer's modified protocol")[tst], protocol, ".")
        Txt <- paste0(Txt, " Tryptic digestion was stopped after ", DigT, c("", " incubation")[(DigT == "overnight")+1], ".")
      }
      if (SPMeth == "in-gel digest") {
        tst <- (DestainMeth == "")+1
        Txt <- paste0(c("The", "All")[(N>1)+1], " sample",  c(" was", "s were")[(N>1)+1], " processed by ", SPMeth,
                      ". Briefly, ", c("the ", "")[(N>1)+1], "gel band",  c(" was", "s were")[(N>1)+1], " cut into roughly 1 mm side cubes, washed x1 in dH₂O", DestainMeth,
                      ", de-hydrated in 900 µL 100% ", c("", "acetonitrile (")[tst], "ACN", c("", ")")[tst], " with shaking for 10 min then removal of all liquid,",
                      " then reduced by incubation in ", ReductC, " ", Reduct, " in 100 mM ", c("", "triethylammonium bicarbonate (")[tst], "TEAB", c("", ")")[tst], " at ",
                      ReductTmp, "°C for ", ReductT, " min, de-hydrated again, then alkylated with ", AlkC, " ", Alk, " in 100 mM TEAB for ", AlkT, " min in the dark. The sample",
                      c(" was", "s were")[(N>1)+1], " dehydrated again, then ", Dig, "-digested at 37°C ", c("for ", "")[(DigT == "overnight")+1], DigT, ".",
                      " Digestion was stopped by addition of 100 µL 1% trifluoroacetic acid and incubation at 37°C for 1h with shaking.",
                      " Peptides were extracted in three steps, by incubating the gel bands successively in 100 µL 25% ACN, 100 µL 100% ACN, and 500 µL 100% ACN, each time for 1h at 37°C with shaking, collecting and pooling the supernatant after each step. ",
                      c("The p", "P")[(N>1)+1], "eptides supernatant",  c(" was", "s were")[(N>1)+1])
        if (Enrich) {
          Txt <- paste0(Txt, " subjected to ", EnrichMeth, " then flow-through and enriched fraction were ",
                        c("", "desalted using Sep-Pak tC18 plates (Waters) and ")[EnrichSepPak+1])
        }
        Txt <- paste0(Txt, " vacuum dried then re-dissolved with 10 min sonication in 1% trifluoroacetic acid.")
      }
      if (SPMeth == "FASP") {
        Txt <- paste0(c("The", "All")[(N>1)+1], " sample",c(" was", "s were")[(N>1)+1], " processed by Filter-Aided Sample Preparation (", SPMeth, "). ",
                      c("The ", "Each")[(N>1)+1], " sample was re-dissolved in 100 µL lysis buffer (4% SDS, 150 mM NaCl, 100 mM Tris HCl pH 7.4, 25 mM TCEP, 10% glycerol)",
                      " and denatured and reduced by incubation at ", ReductTmp, "°C for ", ReductT, " min,",
                      " then allowed to cool down, sonicated with a probe sonicator to shear DNA and cleared by spinning for 10 min at 14,000 g.",
                      " 100 µL of UA buffer (8 M urea in 100 mM triethylammonium bicarbonate = TEAB) was loaded onto a 30,000 MWCO Vivacon 500 Hydrosart spin filter unit to pre-wash it.",
                      " At each step, liquid was pushed through the filter by spinning for 5 min at 14,000 g for every 100 µL of liquid (these times are indicative, and may need to be extended if necessary, especially for concentrated samples).",
                      " The sample supernatant was loaded onto the filter unit, spun, then washed x2 with 300 µL UA buffer and alkylated by incubation for ", AlkT, " min in the dark with ", AlkC, " ", Alk, " in UA buffer.",
                      " After spinning, the filter was washed again in 300 µL UA buffer, then x3 in 100 mM TEAB. The filter unit was transferred to a fresh, protein LoBind collection tube, and 1 µg ", Dig, " in 40 µL in 100 mM TEAB added.",
                      " The assembly was sealed with parafilm and incubated at 37°C ", c("for ", "")[(DigT == "overnight")+1], DigT, ".",
                      " The filter was spun for 10 min, and peptides were successively washed away from the filter by addition of 50 µL 100 mM TEAB, then 500 mM NaCl, each time with 10 min spinning.",
                      " Pooled eluted peptides were acidified by addition of 15 µL 1% trifluoroacetic acid",
                      c("", ", desalted using Sep-Pak tC18 plates (Waters)")[SepPak+1],
                      " and vacuum dried.")
      }
      if (SPMeth == "on-beads digest") {
        if (Dig == "Trypsin") { LysC <- c(TRUE, FALSE)[match(svDialogs::dlg_message("Was LysC pre-digestion performed?", "yesno")$res,
                                                             c("yes", "no"))] }
        if (LysC) { LysCT <- svDialogs::dlg_input("For how long (h)?", 4)$res }
        Txt <- paste0("Paramagnetic beads-bound sample",c(" was", "s were")[(N>1)+1], " digested directly on beads. ",
                      " Beads were re-suspended in 2x beads volume (~40-80 µL) 100 mM triethylammonium bicarbonate (TEAB).")
        if (LysC) { Txt <- gsub("\\.$", paste0(" and pre-digested with LysC (Promega; 400 ng per sample, ", LysCT, " h at 37°C)."), Txt) }
        Txt <- paste0(Txt,
                      " Beads were captured with a magnetic concentrator and the supernatant collected into a fresh protein LoBind tube.",
                      " Proteins on the beads were reduced in 250 µL 100 mM TEAB, ", ReductC, " ", Reduct, " for ", ReductT,  " min at ", ReductTmp, "°C with shaking,",
                      " then alkylated at room temperature in the dark in 250 µL 100 mM TEAB, ", AlkC, " ", Alk, "for ", AlkT, " min.",
                      " After removing the supernatant, samples were digested with 1 µg ", Dig, " at 37°C ", c("for ", "")[(DigT == "overnight")+1], DigT, ".",
                      " Digestion was stopped by addition of 10 µL 10% TFA, the digest was collected",
                      c("", ", pooled with the LysC pre-digest")[LysC+1],
                      c("", ", desalted using Sep-Pak tC18 plates (Waters)")[SepPak+1],
                      " and vacuum dried.")
      }
      # Add isobaric labelling
      if (Label %in% c("TMT", "TMTPRO", "ITRAQ")) {
        N <- length(unique(exp.map[[exp.map.col]]))
        Txt <- paste0(Txt,
                      " Peptides were then labelled with ", LablMeth, " according to the manufacturer's instructions, the reaction was quenched,",
                      " samples were combined according to SUPPLEMENTARY_TABLE_LABELLING_SCHEME, clean-up was completed and the combined sample", c(" was", "s were")[(N>1)+1],
                      " vacuum dried.")
      } else { if (PreOmics) { Txt <- gsub("\\.$", paste0(" and ", c("the ", "")[(N>1)+1], "cleaned-up sample", c(" was", "s were")[(N>1)+1], " vacuum dried."), Txt) } }
      # Add fractionation/enrichment text
      if ((PreOmics)||(SPMeth %in% c("EnrichSepPak", "FASP"))) { # Excluding in-gel digest because we do not fractionate after digest for these (but we may enrich)
        if ((Enrich)&&(Frac)&&(FracSource == "whole sample")) {
          Txt <- paste0(Txt, " ", EnrichPerc, " of ", c("the", "each")[(N>1)+1], " sample was used for ", EnrichMeth, ", while the rest was fractionated into ", NFrac,
                        " fractions using ", FracMeth, ".")
          if (FracMeth == FracMethods[3]) {
            Txt <- paste0(Txt, " Peptides were separated on ", FracCol, ". Mobile phase A: ", FracMPA, "; B: ", FracMPB, "; gradient: ", FracGrad, "; ", FracCollect)
          }
          N2 <- N*(2+NFrac)
        } else {
          if (Enrich) {
            Txt <- paste0(Txt, " The sample was subjected to ", EnrichMeth, ".")
            N2 <- N*2
          }
          if (Frac) {
            if (Enrich) { Txt <- gsub("\\.$", paste0(", then the ", FracSource, " was fractionated into ", NFrac, " fractions using ", FracMeth, "."), Txt) } else {
              Txt <- paste0(Txt, " ", c("The", "Each")[(N>1)+1], " whole ", FracSource, " was fractionated into ", NFrac, " fractions using ", FracMeth, ".")
            }
            if (FracMeth == FracMethods[3]) {
              Txt <- paste0(Txt, " Fractionation was done on ", FracCol, ". Mobile phase A: ", FracMPA, "; B: ", FracMPB, "; gradient: ", FracGrad, "; ", FracCollect)
            }
            N2 <- N*(1+NFrac)
          }
        }
        if (SepPak||EnrichSepPak) {
          tmp <- c("enriched fraction", "flow-through", c("sample", "RP fraction")[Frac+1])[which(c(EnrichSepPak, EnrichSepPak, SepPak))]; l <- length(tmp)
          if (l > 1) { tmp <- paste0(paste(tmp[1:(l-1)], collapse = ", "), " and ", tmp[l]) }
          Txt <- paste0(Txt, " ", c("The", "Each")[(N2>1)+1], " ", tmp, " ", c("was", "were")[(N2>1)+1], " desalted using Sep-Pak tC18 plates (Waters).")
        }
        homePath <- paste0(normalizePath(Sys.getenv("HOME"), winslash = "/"), "/R/proteoCraft")
        fl1 <- paste0(homePath, "/Sample_solvents.txt")
        if (file.exists(fl1)) {
          opt1 <- readLines(fl1)
        } else {
          opt1 <- c("a 9:1 mix of the iST kit's LC-LOAD buffer and 1.5% w/v N-Dodecyl-β-D-maltoside (DDM, 0.15% final)",
                    "a 9:1 mix of the iST kit's LC-LOAD buffer and 0.15% w/v N-Dodecyl-β-D-maltoside (DDM, 0.015% final)",
                    "the iST kit's LC LOAD buffer",
                    "1% TFA (Trifluoroacetic Acid)",
                    "1% FA (Formic Acid)",
                    "5% FA (Formic Acid)")
          write(opt1, fl1)
        }
        opt <- unique(c(opt1, "Other..."))
        m <- max(250, nchar(opt))
        opt <- sapply(opt, function(x) { paste(c(x, rep(" ", m-nchar(x))), collapse = "") })
        loadBuf <- svDialogs::dlg_list(opt, opt[PreOmics+1], "Which LCMS loading buffer were the samples re-dissolved in?")$res
        loadBuf <- gsub(" +$", "", loadBuf)
        if (loadBuf == "Other...") {
          loadBuf <- svDialogs::dlg_input("Enter LCMS loading buffer composition:", "")$res
          write(c(loadBuf, opt1), fl1)
        }
        Txt <- paste0(Txt, " Finally, samples were re-dissolved with 10 min sonication in ", loadBuf, ".")
      }
      MatMetTxt <- Txt  
    } else { MatMetTxt <- "TEMPLATE" }
  }
  return(MatMetTxt)
}
