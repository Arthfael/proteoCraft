# GO_enrich
# ---------
# Source for GO terms enrichment analysis using the topGO package.
# Formerly a function, converted to source to speed it up (because of parallel shenanigans) and to simplify troubleshooting.
#
# Below are the definitions of the former arguments:
# - DB: The formatted protein database, with protein IDs and GO terms. Also serves as background list if Mode = "dataset". Importantly, it MUST contain a "GO" column (containing go term in the form "Term name [GO:ID]", multiple values per row, separated by ";" or "; ").
# - db_ID_col: Name of the column of Protein IDs (one per row) in the database. Default = "Protein ID"
# - Prot: Protein groups file, with at least protein IDs, protein group IDs, and optionally fold change. Importantly, it MUST contain a "GO-ID" column (multiple IDs per row, separated by ";" or "; ").
# - Prot_is_Pep: FALSE by default. If TRUE, this means the "protein groups" table provided is actually a peptides table. The names of some columns generated will be changed accordingly.
# - ID_col: Name of the column listing all Leading Protein IDs in the protein groups file (separated by ";"); default = "Leading protein IDs".
# - ID_col2: Name of the column of group IDs in the protein groups file; default = "id". If missing or not found in the table's columns, will be ignored.
# - db_Gene_col: Name of the column of Gene IDs/names in the database. Default = "Gene". Will be used to create columns of gated Genes per filter.
# - Prot_FC_root: Name of the root of the name of the fold change (or ratio) column in the protein groups file. In "dataset" Mode, will not be extended. In "regulated" Mode, this root will be appended a suffix which should also be the name of one of the filters provided. See the "filters" and "Mode" arguments. argument.
# - Prot_FC_is_log: Are the fold changes already log transformed? If not, they will be.
# - FillGaps: Only for "regulated" more. On the X-axis we plot average Z-scored logFCs. Often we do not have valid logFCs for the best proteins. In that case, the Z score will be under-estimated (best case) or even not calculated. If this is TRUE, we will replace missing logFC values by the highest (esp. lower) logFC in the column, if the protein is only found in numerator (resp. denominator) samples.
# - FillGaps_Smpls: Named list of same length as "filters", each item should contain two named vectors, "Numerator" and "Denominator", of sample names used for the logFC calculations.
# - FillGaps_Expr_root: Root of Expression column names.
# - FillGaps_Expr_is_log: Are expression values log?
# - Mode: One of "dataset" or "regulated". In "dataset" Mode, we are comparing the parent database ("DB")  to the observed dataset and using the "Prot_FC_root" argument as the full name of the fold change column. In "regulated" Mode, we use filters to select regulated protein groups and each will come with its own fold change column. See the "filters" and "Prot_FC_root" arguments.
# - filters: A named list of filters used in "regulated" Mode. Each filter should bear a specific name (e.g. "Treatment 1" or "Whole experiment"; if not, defaults will be provided) and correspond to a vector of row indices in the protein groups file to select. The names will also be used as suffix to add to the fold change root (see the "Prot_FC_root" argument) so cannot be random.
# - ref.filters: Optional (default = NA). If provided, then we do not use the whole Prot/DB table provided (depending on Mode) but first apply this filter to it.
# - show: Default = TRUE. Set to FALSE to not print the graphs.
# - save: Should the plot be saved? Default = FALSE. Set it to a vector of acceptable file extensions to save to the corresponding file format.
# - title: It is possible to set a specific title. Failing that a default one will be provided.
# - title.root: Default = "Bubble_plot_". Is added at the beginning of the file names when saving.
# - bars: Also create barplot? Default = TRUE
# - bars_title: Specific title for bar plot. Failing that a default one will be provided.
# - bars_title.root: Default = "Bar_plot_". Is added at the beginning of the file names when saving.
# - subfolder: Name of the sub-folder where graphs are to be saved. Default = ""
# - subfolderpertype: If set to TRUE (default), will save each type of graph in a dedicated sub-folder.
# - P_adjust: Should P-values plotted be the adjusted versions? (uses p.adjust with the "fdr" method)? Default = FALSE
# - GO_FDR: FDR values to consider. Default = c(0.1, 0.2, 0.3)
# - True_Zscore: Default = FALSE. If FALSE, a discretized Z score analog will be plotted on the X axis. If TRUE, the real Z score will be used.
# - repel: Repel overlapping labels? One of "none", "upper" (default: i.e. only labels for P-value = 0), "lower" (i.e. only labels with P-value > 0), or "both".
# - cex: Label text size. Default = 2.5
# - lineheight: Label interline size. Default = 1
# - plotly: Should we generate a plotly?
# - plotly_subfolder: Default = NULL (in which case if will save in the same folder as "subfolder". Name of the subfolder where data is to be saved for plotly.
# - MinCount: Minimum count for a Term to be included. Default = 10
# - MinTerms: How many terms should we at least print labels for? Default = 15
# - MaxTerms: Up to how many terms should we print labels for? Default = 30
# - MaxTerms_bar: Up to how many terms should we print for bar plots? Default = 100
# - MaxChar: GO term labels get truncated at this number of characters. Default = 25
# - graph: Save a graph of enriched GO terms? Default = TRUE
# - OffspringCounts: If set to TRUE (default), counts include each term's descendants. Warning, this can slow down the function considerably, unless you provide pre-mapped term offspring using with "GO.terms" argument!
# - GO.mappings: Really useful if OffspringCounts is TRUE, otherwise doesn't hurt: provide lists of Proteins/Genes mapped to GO terms.
# - GO.terms: Really useful if OffspringCounts is TRUE, otherwise doesn't hurt: provide pre-mapped GO_terms.
#             Should be as exhaustive as possible, i.e. all known mappings for all proteins. The function expects a data.frame with the following columns:
#              - "ID": term accession
#              - "Term": term name
#              - "Ontology": optional, "BP", "CC" or "MF"
#              - "Offspring": optional, but heavily recommended for a faster function, the list of all of each term's offspring terms
# - N.reserved: Default = 1. Number of reserved vCPUs the function is not to use. Note that for obvious reasons it will always use at least one.
# - GlobalScales: If TRUE (default), single consistent X and Y scales are applied on all plots, otherwise plot-specific scales are applied.

# Check our parent cluster
source(parSrc, local = FALSE)

# Argument names
allArgs <- c("db_ID_col",
             "True_Zscore",
             "ID_col",
             "ID_col2",
             "db_Gene_col",
             "Prot",
             "filters", # MUST STAY HERE
             "ref.filters", # MUST STAY HERE
             "Prot_FC_root",
             "Prot_FC_is_log",
             "show",
             "title",
             "bars_title",
             "title.root",
             "bars",
             "bars_title.root",
             "subfolder",
             "subfolderpertype",
             "FillGaps",
             "FillGaps_Smpls",
             "FillGaps_Expr_root",
             "FillGaps_Expr_is_log",
             "GO_FDR",
             "Prot_is_Pep",
             "P_adjust",
             "repel",
             "cex",
             "lineheight",
             "plotly",
             "plotly_subfolder",
             "MinCount",
             "MinTerms",
             "MaxTerms",
             "MaxTerms_bar",
             "MaxChar",
             "graph",
             "OffspringCounts",
             "GlobalScales",
             "N.reserved",
             "GO.mappings",
             "GO.terms",
             "save",
             "subfolderpertype",
             "bars_title")
w <- which(allArgs %in% .obj)
if (length(w)) {
  warning(paste0("Source-specific objects ", paste(allArgs[w], collapse = "/"),
                 " overlap with reserved object names in this workflow and should be renamed!"))
}

# Defaults
DB <- db
db_ID_col <- "Protein ID"
Prot_is_Pep <- FALSE
ID_col <- "Leading protein IDs"
ID_col2 <- "id"
db_Gene_col <- "Gene"
Prot_FC_root <- "Mean log2(Ratio) "
Prot_FC_is_log <- TRUE
FillGaps <- FALSE
FillGaps_Expr_is_log <- TRUE
ref.filters <- NA
show <- TRUE
title.root <- "Bubble_plot_"
bars <- TRUE
bars_title.root <- "Bar_plot_"
subfolder <- ""
subfolderpertype <- TRUE
P_adjust <- FALSE
GO_FDR <- c(0.1, 0.2, 0.3)
True_Zscore <- FALSE
repel <- "upper"
cex <- 2.5
lineheight <- 1
plotly <- create_plotly
plotly_subfolder <- NULL
MinCount <- 10
MinTerms <- 15
MaxTerms <- 30
MaxTerms_bar <- 100
MaxChar <- 25
graph <- TRUE
OffspringCounts <- TRUE
GlobalScales <- TRUE
N.reserved <- 1
GO.mappings <- GO_mappings
GO.terms <- GO_terms
save <- c("jpeg", "pdf")
subfolderpertype <- FALSE
bars_title <- title <- ""
#
# Now define other parameter values based on how the old function was called in each instance
if (Mode == "regulated") {
  True_Zscore <- TRUE
  if (dataType %in% c("PG", "Prot")) {
    if (dataType == "Prot") { ID_col <- "Protein" }
    if (scrptType == "withReps") {
      Prot <- GO_enrich.dat[[tstbee]]
      filters <- flt
      ref.filters <- Ref.Filt
      Prot_FC_root <- GO_enrich.FCRt[[tstbee]]
      show <- (bee == "By condition")
      title <- "GO bubble plot_"
      bars_title <- "GO bar plot_"
      title.root <- ttr
      bars_title.root <- btr
      subfolder <- dir
    }
    if (scrptType == "noReps") {
      Prot <- PG
      filters <- fcFilt
      ref.filters <- filt2
      Prot_FC_root <- PG.rat.col
      FillGaps <- TRUE
      FillGaps_Smpls <- FC_Smpls
      FillGaps_Expr_root <- PG.int.col
      FillGaps_Expr_is_log <- TRUE
      GO_FDR <- c(0.1, 0.2, 0.3)
    }
  }
  if (dataType == "modPeptides") {
    Prot_is_Pep <- TRUE
    ID_col <- "Proteins"
    if (scrptType == "withReps") {
      Prot <- temPTM
      filters <- flt
      ref.filters <- Pep.Ref.Filt
      Prot_FC_root <- PTMs_GO_enrich.FCRt[[Ptm]][[tstbee]]
      show <- (bee == "By condition")
      title <- paste0(Ptm, " GO bubble plot_")
      bars_title <- paste0(Ptm, " GO bar plot_")
      title.root <- ttr
      bars_title.root <- btr
      subfolder <- dir
    }
    if (scrptType == "noReps") {
      Prot <- ptmpep
      filters <- PTMs_FC_filt[[ptm]]
      ref.filters <- filt3
      Prot_FC_root <- paste0(PTMs_ratRf[length(PTMs_ratRf)], " - ")
      title <- paste0(ptm, " GO bubble plot_")
      bars_title <- paste0(ptm, " GO bar plot_")
      FillGaps <- TRUE
      FillGaps_Smpls <- FC_Smpls
      FillGaps_Expr_root <- paste0(PTMs_intRf[length(PTMs_intRf)], " - ")
      FillGaps_Expr_is_log <- FALSE
      GO_FDR <- c(0.1, 0.2, 0.3)
    }
  }
}
if (Mode == "dataset") {
  True_Zscore <- FALSE
  if (dataType == "PG") {
    Prot <- PG
    if (scrptType == "withReps") {
      Prot_FC_root <- "Av. log10 abundance"
      subfolder <- "Reg. analysis/GO enrich/Dataset"
      subfolderpertype <- FALSE
    }
    if (scrptType == "noReps") {
      filters <- xprsFilt
      Prot_FC_root <- ref
    }
  }
  if (dataType == "modPeptides") {
    # Placeholder
  }
}

#Prot = GO_enrich.dat[[tstbee]]; Mode = "regulated"; filters = flt; ref.filters = Ref.Filt; Prot_FC_root = GO_enrich.FCRt[[tstbee]]; show = (bee == "By condition"); title.root = paste0("Bubble_plot_", tolower(bee)); bars_title.root = paste0("Bar_plot_", tolower(bee)); save = c("jpeg", "pdf"); return = TRUE; True_Zscore = TRUE; subfolder = dir; subfolderpertype = FALSE
# OR (dataset)
#Prot = PG; Mode = "dataset"; Prot_FC_root = "Av. log10 abundance"; save = c("jpeg", "pdf"); return = TRUE; True_Zscore = TRUE; subfolder = "Reg. analysis/GO enrich/Dataset"; subfolderpertype = FALSE
# OR (modified peptides)
#Prot = temPTM; Mode = "regulated"; ID_col = "Proteins"; filters = flt; ref.filters = Pep.Ref.Filt; Prot_FC_root = PTMs_GO_enrich.FCRt[[ptm]][[tstbee]]; show = (bee == "By condition"); title.root = paste0("Bubble_plot_", tolower(bee)); save = c("jpeg", "pdf"); return = TRUE; True_Zscore = TRUE; subfolder = dir; subfolderpertype = FALSE
# OR (no replicates script:
#      - sample composition analysis
#Prot = PG; Mode = "dataset"; filters = filt; Prot_FC_root = ref; save = c("jpeg", "pdf"); return = TRUE; True_Zscore = FALSE; subfolderpertype = FALSE
#      -  2 samples comparisons analysis
#Prot = PG; Mode = "regulated"; filters = FC_filt; ref.filters = filt[names(FC_filt)]; Prot_FC_root = paste0(rat.col, " - "); FillGaps = TRUE; FillGaps_Smpls = FC_Smpls; FillGaps_Expr_root = PG.int.col; FillGaps_Expr_is_log = TRUE; save = c("jpeg", "pdf"); return = TRUE; GO_FDR = c(0.1, 0.2, 0.3); True_Zscore = TRUE; subfolderpertype = FALSE
# Preliminary admin stuff...
origWD <- getwd()
if (!exists("plotEval")) { plotEval <- function(plot) { ggplotify::as.ggplot(ggplotify::as.grob(plot)) } }
#
if (!exists("title")) { title <- "" }
if ((!is.logical(GlobalScales))||(!GlobalScales %in% c(TRUE, FALSE))) { GlobalScales <- TRUE }
if ((!is.logical(Prot_is_Pep))||(!Prot_is_Pep %in% c(TRUE, FALSE))) { Prot_is_Pep <- FALSE }
ID_col2_found <- ID_col2 %in% colnames(Prot)
if ((bars)&&(!exists("bars_title"))) { bars_title <- "" }
if (MaxChar < 20) { MaxChar <- 20 } # I mean, come on! A GO term is 10 characters and we also want to see its name.
stopifnot(db_ID_col %in% colnames(DB), nrow(Prot) > 0, ncol(Prot) > 0, is.character(subfolder), is.logical(subfolderpertype))
subfolder <- gsub("^/+|/+$", "", subfolder)
if ((!exists("subfolder"))||(is.null(subfolder))||(!nchar(subfolder))) {
  subfolder <- origWD
} else {
  subfolder <- paste0(origWD, "/", # Check for full path - step 2
                      gsub(proteoCraft::topattern(paste0(origWD, "/")), "", # Check for full path - step 1
                           gsub("\\*|\\?|<|>|\\|", "-", # Remove unallowed characters
                                gsub("^/+|/+$", "", # Check for outward slashes
                                     gsub("\\\\", "/", # Convert to forward slashes
                                          subfolder)))))
  if (!dir.exists(subfolder)) {
    tst <- try(dir.create(subfolder, recursive = TRUE), silent = TRUE)
    if (class(tst) == "try-error") {
      warning("Invalid subfolder, ignoring!")
      subfolder <- origWD
    }
  }
}
if ((length(save) > 1)||(save != FALSE)) {
  save <- unique(gsub("^jpg$", "jpeg", gsub("^\\.", "", tolower(save))))
  for (ss in save) {
    if (subfolderpertype) { sfpt <- paste0(subfolder, "/", ss) } else { sfpt <- subfolder }
    if (!dir.exists(sfpt)) { dir.create(sfpt, recursive = TRUE) }
  }
}
if (nchar(title)) {
  title <- paste0(gsub("[- _\\.\\,;]+$", "", title), " - ")
  if (title == " - ") { title <- "" }
}
if (nchar(title.root)) {
  title.root <- paste0(gsub("[- _\\.\\,;]+$", "", title.root), " - ")
  if (title.root == " - ") { title.root <- "" }
}
if ((title.root == "")&&(title != "")) {
  title.root <- title
} else { if (title != "") { title.root <- paste0(title.root, title) } }
if (bars) {
  if (nchar(bars_title)) {
    bars_title <- paste0(gsub("[- _\\.\\,;]+$", "", bars_title), " - ")
    if (bars_title == " - ") { bars_title <- "" }
  }
  if (nchar(bars_title.root)) {
    bars_title.root <- paste0(gsub("[- _\\.\\,;]+$", "", bars_title.root), " - ")
    if (bars_title.root == " - ") { bars_title.root <- "" }
  }
  if ((bars_title.root == "")&&(bars_title != "")) {
    bars_title.root <- bars_title
  } else { if (bars_title != "") { bars_title.root <- paste0(bars_title.root, bars_title) } }
  
}
stopifnot(ID_col %in% colnames(Prot), Mode %in% c("dataset", "regulated"),
          is.numeric(as.numeric(cex)), is.numeric(as.numeric(lineheight)), is.logical(plotly),
          is.integer(as.integer(MinTerms)), is.integer(as.integer(MaxTerms)),
          is.logical(OffspringCounts))
cex <- as.numeric(cex)
lineheight <- as.numeric(lineheight)
MinTerms <- as.integer(MinTerms)
MaxTerms <- as.integer(MaxTerms)
if (MinTerms < 0) {
  warning("Parameter \"MinTerms\" cannot be negative, defaulting to 0!")
  MinTerms <- 0
}
if (MaxTerms < 0) {
  warning("Parameter \"MaxTerms\" cannot be negative, defaulting to 1!")
  MaxTerms <- 0
}
if (MaxTerms < MinTerms) {
  warning("Parameter \"MaxTerms\" cannot be smaller to \"MinTerms\", setting it to \"MinTerms\"!")
  MaxTerms <- MinTerms
}
if (bars) {
  stopifnot(is.integer(as.integer(MaxTerms_bar)))
  MaxTerms_bar <- as.integer(MaxTerms_bar)
  if (MaxTerms_bar < 0) {
    warning("Parameter \"MaxTerms_bar\" cannot be negative, defaulting to 1!")
    MaxTerms_bar <- 0
  }
}
repel <- tolower(repel)
if (!repel %in% c("none", "upper", "lower", "both")) {
  warning("I cannot understand the value of the \"repel\" argument, defaulting to \"upper\"!")
  repel <- "upper"
}
#
w <- (repel %in% c("lower", "both"))+1
textFun <- ggplot2::geom_text
if (repel %in% c("lower", "both")) { textFun <- ggrepel::geom_text_repel }
if (plotly) {
  if ((!exists("plotly_subfolder"))||(is.null(plotly_subfolder))||(!nchar(plotly_subfolder))) {
    plotly_subfolder <- subfolder
  } else {
    plotly_subfolder <- paste0(origWD, "/", # Check for full path - step 2
                               gsub(proteoCraft::topattern(paste0(origWD, "/")), "", # Check for full path - step 1
                                    gsub("\\*|\\?|<|>|\\|", "-", # Remove unallowed characters
                                         gsub("^/+|/+$", "", # Check for outward slashes
                                              gsub("\\\\", "/", # Convert to forward slashes
                                                   plotly_subfolder)))))
    if (!dir.exists(plotly_subfolder)) {
      tst <- try(dir.create(plotly_subfolder, recursive = TRUE), silent = TRUE)
      if (class(tst) == "try-error") {
        warning("Invalid plotly_subfolder, ignoring!")
        plotly_subfolder <- subfolder
      }
    }
  }
  library(htmlwidgets)
  if ((plotly_subfolder == origWD)&&(subfolderpertype)) {
    plotly_subfolder <- paste0(plotly_subfolder, "/html")
    if (!dir.exists(plotly_subfolder)) { dir.create(plotly_subfolder, recursive = TRUE) }
  }
}
#
cat("\n   topGO analysis\n   --------------\n")
# Data wrangling
cat("     Preparing data...\n")
Ont <- c("BP", "CC", "MF")
if (!exists("GO.terms")) {
  cat(" - Creating GO.terms object...\n")
  GO.terms <- data.frame(ID = gsub(".+ \\[|\\]", "", unique(unlist(strsplit(DB$GO, ";")))),
                         Term = unique(unlist(strsplit(DB$GO, ";"))),
                         Mapping = "UniProtKB")
} else {
  stopifnot("ID" %in% colnames(GO.terms),
            "Term" %in% colnames(GO.terms))
  GO.terms$Mapping <- "Input data"
}
#GO.terms <- GO.terms[which(!is.na(GO.terms$ID)),]
NoOnt <- !"Ontology" %in% colnames(GO.terms)
NoOffspr <- (!"Offspring" %in% colnames(GO.terms))&(OffspringCounts)
if (NoOnt) { GO.terms$Ontology <- NA }
if (NoOffspr) { GO.terms$Offspring <- list(NA) }
if (NoOnt + NoOffspr) {
  for (ont in Ont) { #ont <- Ont[1]
    #print(ont)
    if (NoOnt) {
      w <- which(annotate::filterGOByOntology(GO.terms$ID, ont))
      GO.terms$Ontology[w] <- ont
    } else { w <- which(GO.terms$Ontology == ont) }
    if (NoOffspr) {
      Offspr <- get(paste0("GO.db::GO", ont, "OFFSPRING"))
      # Finally I managed to rewrite the tedious topGO code much faster,
      # and without even using parallel!!!
      Offspr <- AnnotationDbi::toTable(Offspr)
      colnames(Offspr) <- c("Offspring", "Parent")
      Offspr <- data.table::as.data.table(Offspr)
      Offspr <- Offspr[, list(Offspring = list(Offspring)), by = list(Parent = Parent)]
      Offspr <- as.data.frame(Offspr)
      GO.terms$Offspring[w] <- Offspr$Offspring[match(GO.terms$ID[w], Offspr$Parent)]
    }
  }
  if (NoOffspr) {
    GO.terms$Offspring <- apply(GO.terms[, c("ID", "Offspring")], 1, function(x) {
      x <- unique(unlist(x))
      return(x[which(!is.na(x))])
    })
  }
}
# Remove deprecated GO.terms (i.e. those which do not have an ontology):
GO.terms <- GO.terms[which(!is.na(GO.terms$Ontology)),]
# Match terms to rows in the proteins table
GO.terms$"Protein table row(s)" <- NA
if (OffspringCounts) {
  tmp1 <- proteoCraft::listMelt(apply(GO.terms[, c("ID", "Offspring")], 1, function(x) { unlist(x) }), 1:nrow(GO.terms))
} else {
  tmp1 <- data.frame(value = GO.terms$ID, L1 = 1:nrow(GO.terms))
}
#tmp1 <- aggregate(tmp1$L1, list(tmp1$value), c)
w <- which(nchar(Prot$`GO-ID`) > 0)
tmp2 <- proteoCraft::listMelt(strsplit(Prot$`GO-ID`[w], ";"), w)
tmp2 <- aggregate(tmp2$L1, list(tmp2$value), function(x) { as.numeric(x) })
tmp1$Rows <- tmp2$x[match(tmp1$value, tmp2$Group.1)]
GO.terms$"Protein table row(s)" <- tmp1$Rows[match(GO.terms$ID, tmp1$value)] # Includes Offspring if OffspringCounts == TRUE
# Important:
# This column must be the basis of all subsequent filter-based protein/protein group/gene counts/IDs columns!!!
#
#
#
# Check that everything is upper case
#DB[[db_ID_col]] <- toupper(DB[[db_ID_col]])
#Prot[[ID_col]] <- toupper(Prot[[ID_col]])
# Map to protein accessions
# (This is used for Protein/Gene columns)
if (!exists("GO.mappings")) {
  cat(" - Creating GO.mappings object...\n")
  temp <- proteoCraft::GO_map(DB, db_ID_col, db_Gene_col, GO.terms)
  GO.mappings <- temp$Mappings
  GO.terms <- temp$GO.terms
}
if (!"Proteins" %in% colnames(GO.terms)) {
  GO.terms$Proteins <- GO.mappings$Protein$Protein[match(GO.terms$ID, GO.mappings$Protein$GO)]
}
GenTst <- (!is.null(db_Gene_col))&&(db_Gene_col %in% colnames(DB))
if ((GenTst)&&(!"Genes" %in% colnames(GO.terms))) { 
  GO.terms$Genes <- GO.mappings$Gene$Gene[match(GO.terms$ID, GO.mappings$Gene$GO)]
}
# Create mappings for parent dataset
if (Mode == "regulated") {
  pahruhnt <- Prot
  pahrkol <- ID_col
}
if (Mode == "dataset") {
  pahruhnt <- DB
  pahrkol <- db_ID_col
}
w <- which(nchar(GO.mappings$Protein$Protein) > 0)
Mappings2 <- proteoCraft::listMelt(strsplit(GO.mappings$Protein$Protein[w], ";"), GO.mappings$Protein$GO[w])
Mappings2 <- Mappings2[which(Mappings2$value %in% unlist(strsplit(pahruhnt[[pahrkol]], ";"))),]
Mappings2 <- data.table::data.table(L1 = Mappings2$L1, value = Mappings2$value)
Mappings2 <- Mappings2[, list(GO = list(unique(L1))), by = list(Protein = value)]
Mappings2 <- as.data.frame(Mappings2)
Mappings <- setNames(Mappings2$GO, Mappings2$Protein)
# Note: those mappings should be used in priority over the ones provided in Prot and DB,
# because they (optionally) incorporate indirect mappings to proteins only annotated with offspring terms.
#
# Check filters
# NB:
# - These apply to rows of Prot
# - They are optional in "dataset" Mode, obligatory in "regulated" Mode
defltFilt <- FALSE
if ((Mode == "regulated")||((Mode == "dataset")&&(exists("filters")))) {
  stopifnot(class(filters) == "list", length(filters) > 0)
  if (is.null(names(filters))) { names(filters) <- paste0("Filter ", 1:length(filters)) }
} else {
  # Provide default for "dataset" Mode
  defltFilt <- TRUE
  filters <- setNames(list(1:nrow(Prot)), "Observed dataset")
}
# Reference filters (optional in both Modes)
# Those apply to rows of the parent table: Prot in "regulated" Mode, DB in "dataset" Mode!!!
# Defaults also change how we deal with parent counts:
MultRef <- (length(ref.filters) > 1)||(!is.na(ref.filters))
if (MultRef) {
  # Process/check if provided...
  stopifnot(class(ref.filters) == "list", length(ref.filters) > 0, length(ref.filters) == length(filters))
  if (is.null(names(ref.filters))) {
    names(ref.filters) <- names(filters)
  } else {
    stopifnot(sum(sort(names(ref.filters)) != sort(names(filters))) == 0) # Both filter types must have the same names!
  }
  ref.filters <- ref.filters[names(filters)] # Check order (since we call them by index, not name!)
} else {
  # ... otherwise provide defaults
  ref.filters <- setNames(lapply(names(filters), function(x) {
    1:nrow(pahruhnt)
  }), names(filters))
}
# Exclude contaminants
if ("Potential contaminant" %in% colnames(pahruhnt)) {
  wCnt <- which(pahruhnt$"Potential contaminant" == "+")
  if (length(wCnt)) {
    filters <- setNames(lapply(filters, function(x) { x[which(!x %in% wCnt)] }), names(filters))
    ref.filters <- setNames(lapply(ref.filters, function(x) { x[which(!x %in% wCnt)] }), names(ref.filters))
  }
}
# Create new filters applicable to Mappings (list of GO IDs per protein) from user-provided filters
mapFilters <- setNames(lapply(filters, function(x) {
  names(Mappings)[which(names(Mappings) %in% unique(unlist(strsplit(pahruhnt[unlist(x), pahrkol], ";"))))]
}), names(filters))
#vapply(filters, length, 1)
#vapply(mapFilters, length, 1)
ref.mapFilters <- setNames(lapply(ref.filters, function(x) {
  names(Mappings)[which(names(Mappings) %in% unique(unlist(strsplit(pahruhnt[unlist(x), pahrkol], ";"))))]
}), names(ref.filters))
#vapply(ref.filters, length, 1)
#vapply(ref.mapFilters, length, 1)
#
GO_plots <- GO_FDR_thresholds <- GO_plot_ly <- list()
kount <- 0
scrange <- c(1, 30)
grphs <- c()
require(topGO)
FishTst <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
KKol1 <- paste0("Counts - parent ", c("dataset", "database")[match(Mode, c("regulated", "dataset"))])
ProtKol1 <- paste0("Proteins - parent ", c("dataset", "database")[match(Mode, c("regulated", "dataset"))])
PGKol1 <- paste0(c("PG", "Peptide")[Prot_is_Pep+1], " IDs - parent ", c("dataset", "database")[match(Mode, c("regulated", "dataset"))])
if (!MultRef) {
  GenKol1 <- paste0("Genes - parent ", c("dataset", "database")[match(Mode, c("regulated", "dataset"))])
}
# Process and test filters to generate plots input data
wFltL <- which(vapply(filters, length, 1) > 0)
if (length(wFltL)) {
  mapFilters <- mapFilters[wFltL]
  ref.mapFilters <- ref.mapFilters[wFltL]
  tmpGO <- GO.terms[, c("ID", "Term")]
  exports <- list("mapFilters", "ref.mapFilters", "KKol1", "ProtKol1", "PGKol1", "MultRef", "GenTst", "tmpGO", "Mappings", "Ont", "Mode",
                  "Prot_is_Pep", "MinCount", "FishTst")
  if (!MultRef) { exports <- append(exports, "GenKol1") }
  parallel::clusterExport(parClust, exports, envir = environment())
  #parallel::clusterExport(parClust, "FishTst", envir = environment())
  Fisher0 <- function(n1) { #n1 <- names(mapFilters)[1]
    n2 <- gsub("___", " ", n1)
    #if (length(mapFilters) > 1) { cat(paste0("Filter ", n2, "\n")) }
    # Column names
    kkol1 <- KKol1
    if (MultRef) { kkol1 <- paste0("Counts - ", n1 , " (parent)") }
    kkol2 <- paste0("Counts - ", n1)
    protkol2 <- paste0("Proteins - ", n1)
    pgkol2 <- paste0(c("PG", "Peptide")[Prot_is_Pep+1], " IDs - ", n1)
    if (MultRef) {
      protkol1 <- paste0(protkol2, " (parent)")
      pgkol1 <- paste0(pgkol2, " (parent)")
    } else {
      protkol1 <- ProtKol1
      pgkol1 <- PGKol1
    }
    if (GenTst) {
      genkol2 <- paste0("Genes - ", n1)
      if (MultRef) { genkol1 <- paste0(genkol2, " (parent)") } else { genkol1 <- GenKol1 }
    }
    # Reference proteins list (from reference filter)
    allProt <- setNames(1:length(Mappings), names(Mappings))
    allProt <- allProt[ref.mapFilters[[n1]]]
    # topGO
    GOdata <- setNames(lapply(Ont, function(ont) { #ont <- Ont[1]
      cat(paste0(" - ", ont, " terms\n"))
      new("topGOdata",
          description = paste0(n2, " - ", ont),
          ontology = ont,
          allGenes = allProt,
          geneSel = function(x) { x %in% allProt[mapFilters[[n1]]] }, # Apply filter
          nodeSize = MinCount,
          annot = topGO::annFUN.gene2GO,
          gene2GO = Mappings)
    }), Ont)
    cat("Performing Fisher exact test:\n")
    resultFisher <- setNames(lapply(Ont, function(ont) { #ont <- Ont[1]
      cat(paste0(" - ", ont, " terms\n"))
      try(topGO::getSigGroups(GOdata[[ont]], FishTst), silent = TRUE)
    }), Ont)
    Wh1 <- suppressWarnings(which(vapply(resultFisher, function(x) { !"try-error" %in% class(x) }, TRUE)))
    if (!length(Wh1)) {
      msg <- "Analysis failed, investigate!"
      if (Mode == "regulated") { paste0("Filter ", n2, ": ", msg) }
      warning(msg)
      return(list(Outcome = FALSE, Message = msg))
    } else {
      Wh2 <- which(vapply(Ont[Wh1], function(ont) {
        resultFisher[[ont]]@geneData[["Significant"]] > 0
      }, TRUE))
      Wh1 <- Wh1[Wh2]
      if (!length(Wh1)) {
        msg <- "0 significant GO Terms, skipping!"
        if (Mode == "regulated") { paste0("Filter ", n2, ": ", msg) }
        warning(msg)
        res <- list(Outcome = FALSE,
                    Message = msg)
      } else {
        #cat(" -> Success...\n")
        if (length(Wh1) < length(Ont)) {
          msg <- Ont[-Wh1]
          if (length(msg) == 1) { msg <- paste0("Term ", msg) } else {
            msg <- paste0("Terms ", paste(msg[1:(length(msg)-1)], collapse = ", "), " and ", msg[length(msg)])
          }
          msg <- paste0("Analysis failed for ", msg)
          if (Mode == "regulated") { msg <- paste0("Filter ", n2, ": ", msg) }
          warning(msg)
        }
        GO_tbl <- lapply(Ont[Wh1], function(ont) {
          data.frame(Ontology = ont,
                     ID = names(resultFisher[[ont]]@score),
                     Pvalue = resultFisher[[ont]]@score)
        })
        GO_tbl <- reshape::melt.list(GO_tbl, id.vars = c("Ontology", "ID"))
        GO_tbl$Ontology <- Ont[Wh1][GO_tbl$L1]
        GO_tbl$L1 <- NULL
        GO_tbl$variable <- NULL
        colnames(GO_tbl)[which(colnames(GO_tbl) == "value")] <- "Pvalue"
        GO_tbl$Term <- tmpGO$Term[match(GO_tbl$ID, tmpGO$ID)]
        GO_tbl$Mapping <- NA
        res <- list(Outcome = TRUE,
                    Output = GO_tbl,
                    kkol1 = kkol1,
                    kkol2 = kkol2,
                    protkol1 = protkol1,
                    protkol2 = protkol2,
                    pgkol1 = pgkol1,
                    pgkol2 = pgkol2,
                    Wh1 = Wh1)
        if (GenTst) {
          res$genkol1 <- genkol1
          res$genkol2 <- genkol2
        }
      }
      return(res)
    }
  }
  #environment(Fisher0) <- .GlobalEnv # Only needed if code run as function!
  cat("     Running Fisher tests...\n")
  invisible(parallel::clusterCall(parClust, function() { library(topGO); return() }))
  tst <- try({
    GO_tbls <- setNames(parallel::parLapply(parClust, names(mapFilters), Fisher0), names(mapFilters))
  })
  if ("try-error" %in% class(tst)) {
    GO_tbls <- setNames(lapply(names(mapFilters), Fisher0), names(mapFilters))
  }
  #vapply(GO_tbls, function(x) { x$Outcome }, TRUE)
  GO_tbls <- GO_tbls[which(vapply(GO_tbls, function(x) { x$Outcome }, TRUE))]
  if (length(GO_tbls)) {
    # Define filter functions
    f0 <- function(x, filt) {
      if (length(x)) { x <- sum(x %in% filt) } else { x <- 0 }
      return(as.numeric(x))
    }
    #environment(f0) <- .GlobalEnv # Only needed if code run as function!
    f1 <- function(x, filt, ids) {
      x <- x[which(x %in% filt)]
      if (length(x)) { x <- paste(sort(unique(unlist(ids[x]))), collapse = ";") } else { x <- "" }
      return(x)
    }
    #environment(f1) <- .GlobalEnv # Only needed if code run as function!
    #
    tmpGO.terms <- lapply(GO_tbls, function(x) { #x <- GO_tbls[[1]]
      GO_tbl <- x$Output
      GO_tbl[which(is.na(GO_tbl$Term)),]
    })
    tmpGO.terms <- plyr::rbind.fill(tmpGO.terms)
    tmpGO.terms <- tmpGO.terms[which(!tmpGO.terms$ID %in% GO.terms$ID),]
    if ((OffspringCounts)&&(nrow(tmpGO.terms))) { # It only makes sense to use these extra terms if OffspringCounts is TRUE
      #cat("Processing additional terms of interest identified by topGO...\n")
      tmpGO.terms$Pvalue <- NULL
      tmpGO.terms$Term <- NULL
      tmpGO.terms$Mapping <- NULL
      k <- colnames(tmpGO.terms)
      tmpGO.terms <- do.call(paste, c(tmpGO.terms, sep = "___BLEH___"))
      tmpGO.terms <- unique(tmpGO.terms)
      tmpGO.terms <- as.data.frame(t(sapply(strsplit(tmpGO.terms, "___BLEH___"), unlist)))
      colnames(tmpGO.terms) <- k
      tmpGO.terms$Mapping <- "topGO"
      tmpGO.terms$Term <- NA
      tmpGO.terms$Offspring <- list(NA)
      for (ont in Ont) { #ont <- Ont[1]
        wo <- which(annotate::filterGOByOntology(tmpGO.terms$ID, ont))
        if (length(wo)) {
          Offspr <- get(paste0("GO", ont, "OFFSPRING"))
          Offspr <- toTable(Offspr)
          colnames(Offspr) <- c("Offspring", "Parent")
          Offspr <- data.table::as.data.table(Offspr)
          Offspr <- Offspr[, list(Offspring = list(Offspring)), by = list(Parent = Parent)]
          Offspr <- as.data.frame(Offspr)
          tmpGO.terms$Offspring[wo] <- Offspr$Offspring[match(tmpGO.terms$ID[wo], Offspr$Parent)]
          tmpGO.terms$Term[wo] <- annotate::getGOTerm(tmpGO.terms$ID[wo])[[ont]]
        }
      }
      tmpGO.terms <- tmpGO.terms[which(vapply(tmpGO.terms$Offspring, length, 1) > 0),]
      tmpGO.terms$Term <- apply(tmpGO.terms[, c("Term","ID")], 1, function(x) { paste0(unlist(x[[1]]), " [", x[[2]], "]") })
      tmpGO.terms$"Protein table row(s)" <- NA
      #sum(!tmpGO.terms$ID %in% unlist(GO.terms$Offspring)) # Those new terms are all OffSpring terms of existing ones!
      #w <- which(vapply(GO.terms$Offspring, length, 1) > 0)
      tmp1 <- proteoCraft::listMelt(GO.terms$Offspring, 1:nrow(GO.terms), c("ID", "Row")) # Here I can use only Offspring, the IDs in tmpGO.terms are not in GO.terms
      tmp1 <- tmp1[which(tmp1$ID %in% tmpGO.terms$ID),]
      tmp1$"Protein table row(s)" <- GO.terms$"Protein table row(s)"[tmp1$Row]
      tmpGO.terms$"Protein table row(s)" <- tmp1$`Protein table row(s)`[match(tmpGO.terms$ID, tmp1$ID)]
      #
      for (kl in 1:(1+GenTst)) {
        nm <- c("Protein", "Gene")[kl]
        wh <- which(vapply(tmpGO.terms$"Protein table row(s)", length, 1) > 0)
        if (length(wh)) {
          tmp1 <- proteoCraft::listMelt(tmpGO.terms$"Protein table row(s)"[wh], wh)
          tmp1$value2 <- Prot[match(tmp1$value, 1:nrow(Prot)), ID_col]
          tmp1 <- proteoCraft::listMelt(strsplit(tmp1$value2, ";"), tmp1$L1)
          if (kl == 2) {
            tmp1$value <- DB[match(tmp1$value, DB[[db_ID_col]]), db_Gene_col]
          }
          tmp1 <- aggregate(tmp1$value, list(tmp1$L1), function(x) { paste(sort(unique(x)), collapse = ";") })
          tmpGO.terms[[paste0(nm, "s")]] <- tmp1$x[match(1:nrow(tmpGO.terms), tmp1$Group.1)] 
        }
      }
      k <- colnames(GO.terms)[which(!colnames(GO.terms) %in% colnames(tmpGO.terms))]
      if (length(k)) { tmpGO.terms[, k] <- NA }
      GO.terms <- rbind(GO.terms, tmpGO.terms)
      #cat("ok...\n")
    }
    # This part I think is best vectorised locally within the loop, rather than vectorizing the whole loop.
    # (It is more efficient to use many threads occasionally in the vertical (GO terms) direction than systematically in the horizontal (individual filters) direction,
    # with there being only a few filters in most experiments).
    GO_tbls2 <- list()
    for (n1 in names(GO_tbls)) { #n1 <- names(GO_tbls)[1]
      cat(paste0("      - Filter = ", proteoCraft::cleanNms(n1), "\n"))
      kount <- kount + 1
      GO_tbl <- GO_tbls[[n1]]$Output
      kkol1 <- GO_tbls[[n1]]$kkol1
      kkol2 <- GO_tbls[[n1]]$kkol2
      protkol1 <- GO_tbls[[n1]]$protkol1
      protkol2 <- GO_tbls[[n1]]$protkol2
      pgkol1 <- GO_tbls[[n1]]$pgkol1
      pgkol2 <- GO_tbls[[n1]]$pgkol2
      if (GenTst) {
        genkol1 <- GO_tbls[[n1]]$genkol1
        genkol2 <- GO_tbls[[n1]]$genkol2
      }
      Wh1 <- GO_tbls[[n1]]$Wh1
      names(GO_tbls[[n1]])
      w <- which((is.na(GO_tbl$Term)|is.na(GO_tbl$Mapping))&(GO_tbl$ID %in% GO.terms$ID))
      GO_tbl[w, c("Term", "Mapping")] <- GO.terms[match(GO_tbl$ID[w], GO.terms$ID), c("Term", "Mapping")]
      GO_tbl$Rows <- GO.terms$`Protein table row(s)`[match(GO_tbl$ID, GO.terms$ID)]
      # Counts, Proteins and Genes
      #cat("Getting counts of proteins, protein groups, genes...\n")
      ## Counts
      GO_tbl[[kkol1]] <- GO_tbl$Count <- 0
      # No need to distinguish between OffspringCounts TRUE or FALSE here, or below,
      # since the Mappings themselves reflect it:
      fcflt <- filters[[n1]]
      rfflt <- unique(unlist(c(fcflt, ref.filters[[n1]]))) # In case we are missing somethin'
      #
      # Some silly shenanigans to deal with random cluster corruption:
      a <- parallel::parSapply(parClust, GO.terms$`Protein table row(s)`, f0, filt = rfflt)
      if (!length(a) == nrow(GO.terms)) {
        stopCluster(parClust)
        source(parSrc)
        a <- parallel::parSapply(parClust, GO.terms$`Protein table row(s)`, f0, filt = rfflt)
      }
      if (!length(a) == nrow(GO.terms)) { stop() }
      GO.terms[[kkol1]] <- a
      rm(a)
      #
      a <- parallel::parSapply(parClust, GO.terms$`Protein table row(s)`, f0, filt = fcflt)
      if (!length(a) == nrow(GO.terms)) {
        stopCluster(parClust)
        source(parSrc)
        a <- parallel::parSapply(parClust, GO.terms$`Protein table row(s)`, f0, filt = fcflt)
      }
      if (!length(a) == nrow(GO.terms)) { stop() }
      GO.terms[[kkol2]] <- a
      rm(a)
      #
      GO_tbl[, c(kkol1, "Count")] <- GO.terms[match(GO_tbl$ID, GO.terms$ID), c(kkol1, kkol2)]
      stopifnot(max(GO_tbl[[kkol1]] - GO_tbl$Count) >= 0)
      ## Proteins
      tmpProt <- strsplit(Prot[[ID_col]], ";")
      #parallel::clusterExport(parClust, "tmpProt", envir = environment())
      if (ID_col2_found) {
        tmpProtID <- Prot[[ID_col2]]
        #parallel::clusterExport(parClust, "tmpProtID", envir = environment())
      }
      if ((MultRef)||(kount == 1)) {
        a <- parallel::parSapply(parClust, GO.terms$`Protein table row(s)`, f1, filt = rfflt, ids = tmpProt)
        if (!length(a) == nrow(GO.terms)) {
          stopCluster(parClust)
          source(parSrc)
          a <- parallel::parSapply(parClust, GO.terms$`Protein table row(s)`, f1, filt = rfflt, ids = tmpProt)
        }
        if (!length(a) == nrow(GO.terms)) { stop() }
        GO.terms[[protkol1]] <- a
        rm(a)
        #
        if (ID_col2_found) {
          a <- parallel::parSapply(parClust, GO.terms$`Protein table row(s)`, f1, filt = rfflt, ids = tmpProtID)
          if (!length(a) == nrow(GO.terms)) {
            stopCluster(parClust)
            source(parSrc)
            a <- parallel::parSapply(parClust, GO.terms$`Protein table row(s)`, f1, filt = rfflt, ids = tmpProt)
          }
          if (!length(a) == nrow(GO.terms)) { stop() }
          GO.terms[[pgkol1]] <- a
          rm(a)
          #
        }
      }
      #
      a <- parallel::parSapply(parClust, GO.terms$`Protein table row(s)`, f1, filt = fcflt, ids = tmpProt)
      if (!length(a) == nrow(GO.terms)) {
        stopCluster(parClust)
        source(parSrc)
        a <- parallel::parSapply(parClust, GO.terms$`Protein table row(s)`, f1, filt = fcflt, ids = tmpProt)
      }
      if (!length(a) == nrow(GO.terms)) { stop() }
      GO.terms[[protkol2]] <- a
      rm(a)
      #
      tst <- GO.terms[, c(protkol1, protkol2)]
      tst[[protkol1]] <- strsplit(tst[[protkol1]], ";")
      tst[[protkol2]] <- strsplit(tst[[protkol2]], ";")
      tst2 <- apply(tst, 1, function(x) { sum(!x[[2]] %in% x[[1]]) })
      stopifnot(max(tst2) == 0)
      ## Protein groups
      if (ID_col2_found) {
        a <- parallel::parSapply(parClust, GO.terms$`Protein table row(s)`, f1, filt = fcflt, ids = tmpProtID)
        if (!length(a) == nrow(GO.terms)) {
          stopCluster(parClust)
          source(parSrc)
          a <- parallel::parSapply(parClust, GO.terms$`Protein table row(s)`, f1, filt = fcflt, ids = tmpProtID)
        }
        if (!length(a) == nrow(GO.terms)) { stop() }
        GO.terms[[pgkol2]] <- a
        rm(a)
        #    
        tst <- GO.terms[, c(pgkol1, pgkol2)]
        tst[[pgkol1]] <- strsplit(tst[[pgkol1]], ";")
        tst[[pgkol2]] <- strsplit(tst[[pgkol2]], ";")
        stopifnot(max(apply(tst, 1, function(x) { sum(!x[[2]] %in% x[[1]]) })) == 0)
      }
      ## Genes
      if (GenTst) {
        tmp2 <- proteoCraft::listMelt(tmpProt, 1:length(tmpProt))
        tmp2$Gene <- DB[match(tmp2$value, DB[[db_ID_col]]), db_Gene_col]
        tmp2 <- tmp2[which(vapply(tmp2$Gene, nchar, 1) > 0),]
        tmp2 <- proteoCraft::listMelt(strsplit(tmp2$Gene, ";"), tmp2$L1)
        tmp2 <- aggregate(tmp2$value, list(as.numeric(tmp2$L1)), list)
        tmpGn <- data.frame(row = 1:nrow(Prot))
        tmpGn <- tmp2$x[match(tmpGn$row, tmp2$Group.1)]
        #parallel::clusterExport(parClust, "tmpGn", envir = environment())
        if ((MultRef)||(kount == 1)) {
          a <- parallel::parSapply(parClust, GO.terms$`Protein table row(s)`, f1, filt = rfflt, ids = tmpGn)
          if (!length(a) == nrow(GO.terms)) {
            stopCluster(parClust)
            source(parSrc)
            a <- parallel::parSapply(parClust, GO.terms$`Protein table row(s)`, f1, filt = rfflt, ids = tmpGn)
          }
          if (!length(a) == nrow(GO.terms)) { stop() }
          GO.terms[[genkol1]] <- a
          rm(a)
          #
        }
        a <- parallel::parSapply(parClust, GO.terms$`Protein table row(s)`, f1, filt = fcflt, ids = tmpGn)
        if (!length(a) == nrow(GO.terms)) {
          stopCluster(parClust)
          source(parSrc)
          a <- parallel::parSapply(parClust, GO.terms$`Protein table row(s)`, f1, filt = fcflt, ids = tmpGn)
        }
        if (!length(a) == nrow(GO.terms)) { stop() }
        GO.terms[[genkol2]] <- a
        rm(a)
        #
        tst <- GO.terms[, c(genkol1, genkol2)]
        tst[[genkol1]] <- strsplit(tst[[genkol1]], ";")
        tst[[genkol2]] <- strsplit(tst[[genkol2]], ";")
        stopifnot(max(apply(tst, 1, function(x) { sum(!x[[2]] %in% x[[1]]) })) == 0)
      }
      #
      #cat("     Preparing plots...\n")
      tmp <- GO_tbl[, c("ID", "Term")]
      tmp$Term[which(is.na(tmp$Term))] <- ""
      GO_tbl$Label <- do.call(paste, c(tmp, sep = "\n"))
      GO_tbl$Label2 <- GO_tbl$Label
      w <- which(nchar(GO_tbl$Label) > MaxChar)
      GO_tbl$Label2[w] <- paste0(substr(GO_tbl$Label2[w], 1, MaxChar-3), "...")
      GO_tbl$Label3 <- GO_tbl$Term
      w <- which(nchar(GO_tbl$Term) > MaxChar)
      GO_tbl$Label3[w] <- paste0(substr(GO_tbl$Label3[w], 1, MaxChar-3), "...")
      GO.terms[[paste0("Pvalue - ", n1)]] <- GO_tbl$Pvalue[match(GO.terms$ID, GO_tbl$ID)]
      wh <- match(GO_tbl$ID, GO.terms$ID)
      if (P_adjust) {
        GO_tbl$"adj. Pvalue" <- p.adjust(GO_tbl$"Pvalue", "BH")
        GO.terms[[paste0("adj. Pvalue - ", n1)]] <- NA
        GO.terms[wh, paste0("adj. Pvalue - ", n1)] <- GO_tbl$"adj. Pvalue"
        thresh <- data.frame(Threshold = GO_FDR)
      } else {
        GO.terms[, paste0("Significance - ", n1, " ", GO_FDR*100, "%")] <- NA
        thresh <- as.data.frame(sapply(Ont[Wh1], function(x) { GO_FDR*NA }))
        thresh$FDR <- GO_FDR
        for (ont in Ont[Wh1]) { #ont <- Ont[Wh1][1]
          w <- which(GO.terms$Ontology[wh] == ont)
          if (length(w)) {
            f <- proteoCraft::FDR(data = GO.terms[wh[w],], pvalue_col = paste0("Pvalue - ", n1),
                            returns = c(TRUE, TRUE), fdr = GO_FDR)
            GO.terms[wh[w], paste0("Significance - ", n1, " ", GO_FDR*100, "%")] <- f$`Significance vector`
            thresh[[ont]] <- f$Thresholds
          }
        }
        thresh <- reshape::melt.data.frame(thresh, id.vars = "FDR")
        colnames(thresh) <- gsub("^value$", "Threshold", gsub("^variable$", "Ontology", colnames(thresh)))
        GO_FDR_thresholds[[n1]] <- thresh[,c("FDR", "Ontology", "Threshold")]
      }
      GO_tbl$Y <- -log10(GO_tbl[[c("Pvalue", "adj. Pvalue")[P_adjust+1]]]) # I am choosing to call it Y because I am artificially editing some values (+Inf -> Ymax + 1) below
      GO.terms[wh, paste0("Pvalue - ", n1)] <- GO_tbl$Pvalue
      if (P_adjust) { GO.terms[wh, paste0("adj. Pvalue - ", n1)] <- GO_tbl$"adj. Pvalue" }
      thresh$Colour <- colorRampPalette(c("red", "gold"))(length(GO_FDR))
      thresh$"-log10(Threshold)" <- -log10(thresh$Threshold)
      thresh$Label <- paste0(thresh$FDR*100, "% FDR threshold")
      Ymax <- suppressWarnings(max(c(proteoCraft::is.all.good(c(GO_tbl$Y, thresh$`-log10(Threshold)`, 3)))))
      winf <- which(is.infinite(GO_tbl$Y)&(GO_tbl$Y > 0))
      if (length(winf)) { GO_tbl$Y[winf] <- Ymax + 1 }
      # Calculate logFC
      if ((Mode == "regulated")||((Mode == "dataset")&&(!defltFilt))) {
        Prot_FC_root <- paste0(gsub(" (- )?$", "", Prot_FC_root), " - ")
        FCkol <- paste0(Prot_FC_root, n1)
        if (!FCkol %in% colnames(Prot)) {
          if (Mode == "dataset") { FCkol <- gsub(" (- )?$", "", Prot_FC_root) }
        }
      } else {
        FCkol <- Prot_FC_root
      }
      FCkol <- FCkol[which(FCkol %in% colnames(Prot))]
      if (!length(FCkol)) {
        warning(paste0("No fold change column found for this filter (expected name: \"", FCkol, "\"), replacing it with a dummy column!"))
        if (Prot_FC_is_log) { Prot[[FCkol]] <- 0 } else { Prot[[FCkol]] <- 1 }
      }
      w <- which(is.na(Prot[[FCkol]]))
      if ((Mode == "regulated")&&(FillGaps)&&(exists("FillGaps_Expr_root"))&&(!is.null(FillGaps_Expr_root))&&(length(w))) {
        FillGaps_Expr_root <- paste0(gsub(" (- )?$", "", FillGaps_Expr_root), " - ")
        Intkol0 <- paste0(FillGaps_Expr_root, FillGaps_Smpls[[n1]]$Denominator)
        Intkol1 <- paste0(FillGaps_Expr_root, FillGaps_Smpls[[n1]]$Numerator)
        if (!sum(!c(Intkol1, Intkol0) %in% colnames(Prot))) {
          tmp0 <- Prot[w, Intkol0, drop = FALSE]
          tmp1 <- Prot[w, Intkol1, drop = FALSE]
          if (FillGaps_Expr_is_log) {
            tmp0 <- 10^tmp0
            tmp1 <- 10^tmp1
          }
          tmp0 <- apply(tmp0, 1, function(x) { sum(proteoCraft::is.all.good(x, 2)) }) > 0
          tmp1 <- apply(tmp1, 1, function(x) { sum(proteoCraft::is.all.good(x, 2)) }) > 0
          w0 <- which(tmp0&!tmp1)
          w1 <- which(!tmp0&tmp1)
          wB <- which(tmp0&tmp1) # This should always be empty!
          if (length(wB)) { warning("Invalid ratio, yet both sample groups seem to have valid values? Investigate!") }
          Mn <- min(proteoCraft::is.all.good(Prot[[FCkol]]))
          Mx <- max(proteoCraft::is.all.good(Prot[[FCkol]]))
          Prot[[FCkol]][w[w0]] <- Mn
          Prot[[FCkol]][w[w1]] <- Mx
        }
      }
      if (!Prot_FC_is_log) { Prot[[FCkol]] <- log2(Prot[[FCkol]]) }
      GO_tbl$logFC <- vapply(GO_tbl$Rows, function(x) {
        proteoCraft::log_ratio_av(Prot[x, FCkol])
      }, 1)
      GO.terms[[paste0("logFC - ", n1)]] <- NA
      GO.terms[wh, paste0("logFC - ", n1)] <- as.numeric(GO_tbl$logFC)
      # Calculate Z score
      if (True_Zscore) {
        sd <- sd(proteoCraft::is.all.good(GO_tbl$logFC))
        m <- mean(proteoCraft::is.all.good(GO_tbl$logFC))
        GO_tbl$"Z-score" <- (GO_tbl$logFC - m)/sd
        GO.terms[[paste0("Z-score - ", n1)]] <- NA
        GO.terms[wh, paste0("Z-score - ", n1)] <- GO_tbl$"Z-score"
      } else {
        # "Z score" analog as used in package GOplot (see https://cran.r-project.org/web/packages/GOplot/vignettes/GOplot_vignette.html)
        m <- mean(proteoCraft::is.all.good(Prot[[FCkol]]))
        GO_tbl$"Z-score" <- vapply(GO_tbl$Rows, function(x) {
          x <- proteoCraft::is.all.good(Prot[x, FCkol])
          l <- length(x)
          if (l) {
            x <- proteoCraft::is.all.good(x)-m
            x <- x[which(x != 0)]
            x <- ifelse(x > 0, 1, -1)
            x <- sum(x)/sqrt(l)
          } else { x <- NA }
          return(x)
        }, 1)
        GO.terms[[paste0("(N_Up - N_Down)/sqrt(Tot.) - ", n1)]] <- NA
        GO.terms[wh, paste0("(N_Up - N_Down)/sqrt(Tot.) - ", n1)] <- GO_tbl$"Z-score"
      }
      Xmin <- suppressWarnings(min(proteoCraft::is.all.good(c(GO_tbl$"Z-score", -1))))
      Xmax <- suppressWarnings(max(proteoCraft::is.all.good(c(GO_tbl$"Z-score", 1))))
      Xbreadth <- Xmax-Xmin
      Xmin <- Xmin-Xbreadth*0.05
      Xmax <- Xmax+Xbreadth*0.05
      GO_tbl$test <- FALSE
      for (ont in Ont[Wh1]) { #ont <- Ont[Wh1][1]
        w <- which(GO_tbl$Ontology == ont)
        if (P_adjust) { m <- -log10(max(GO_FDR)) } else {
          m <- thresh$"-log10(Threshold)"[which((thresh$Ontology == ont)&(thresh$FDR == max(GO_FDR)))]
        }
        GO_tbl$test[w] <- GO_tbl$Y[w] >= m
        tmp <- GO_tbl[w,]
        ord <- c(1:nrow(tmp))[order(tmp$Y, decreasing = TRUE)]
        if (MinTerms) { tmp$test[which(1:nrow(tmp) %in% ord[1:MinTerms])] <- TRUE }
        if (MaxTerms) { tmp$test[which(1:nrow(tmp) %in% ord[(MaxTerms+1):length(ord)])] <- FALSE }
        GO_tbl[w,] <- tmp
        # aggregate(test, list(test), length)
      }
      w <- which(proteoCraft::is.all.good(GO_tbl$"Z-score", 2)&(proteoCraft::is.all.good(GO_tbl$Y, 2))&(GO_tbl$Count >= MinCount))
      #wN <- which(!proteoCraft::is.all.good(GO_tbl$"Z-score", 2)|(!proteoCraft::is.all.good(GO_tbl$Y, 2))|(GO_tbl$Count < MinCount))
      if (length(w)) {
        GO_tbls2[[n1]] <- list(Success = TRUE,
                               Data = GO_tbl[w,],
                               Thresholds = thresh,
                               Xscales = c("Min" = Xmin, "Max" = Xmax, "Breadth" = Xbreadth),
                               Ymax = Ymax,
                               Infinites = winf,
                               Wh1 = Wh1)
        GO_tbls[[n1]]$Output <- GO_tbl
      } else {
        n2 <- gsub("___", " ", n1)
        msg <- paste0("Not enough terms with valid data for filter ", n2, " -> nothing to plot!")
        warning(msg)
      }
    }
    if (length(GO_tbls2)) {
      if (GlobalScales) {
        Xxtr <- max(c(vapply(names(GO_tbls), function(nm) {
          suppressWarnings(max(abs(proteoCraft::is.all.good(GO_tbls[[nm]]$Data$"Z-score"*1.05))))
        }, 1), 1))
        Xbreadth <- 2*Xxtr
        Xmin <- -Xxtr
        Xmax <- Xxtr
        Ymax <- max(c(vapply(names(GO_tbls), function(nm) {
          suppressWarnings(max(proteoCraft::is.all.good(c(GO_tbls[[nm]]$Data$Y, GO_tbls[[nm]]$Thresholds$"-log10(Threshold)"))))
        }, 1), 3))
      } else {
        Xxtr <- Xbreadth <- Xmin <- Xmax <- Ymax <- ""
      }
      exports <- list("GlobalScales", "Xxtr", "Xbreadth", "Xmin", "Xmax", "Ymax", "GO_tbls2", "GO_plots", "title.root", "P_adjust", "plotly", "show", "grphs",
                      "save", "origWD", "subfolder", "subfolderpertype", "bars", "graph", "True_Zscore", "scrange", "textFun", "cex", "lineheight", "repel",
                      "plotly_subfolder", "MaxTerms", "MaxTerms_bar", "MaxChar",
                      "Ont", "title", "title.root", "bars_title", "bars_title.root", "GO_FDR", "plotEval")
      parallel::clusterExport(parClust, exports, envir = environment())
      plotsF0 <- function(n1) { #n1 <- names(GO_tbls2)[1] #n1 <- names(GO_tbls2)[2]
        GOplts <- list()
        n2 <- gsub("___", " ", n1)
        GO_tbl <- GO_tbls2[[n1]]$Data
        winf <- GO_tbls2[[n1]]$Infinites
        thresh <- GO_tbls2[[n1]]$Thresholds
        Wh1 <- GO_tbls2[[n1]]$Wh1
        if (!GlobalScales) {
          Xxtr <- max(abs(c(GO_tbls2[[n1]]$Xscales["Min"], GO_tbls2[[n1]]$Xscales["Max"], 1)))*1.05
          Xbreadth <- 2*Xxtr
          Xmin <- -Xxtr
          Xmax <- Xxtr
          Ymax <- GO_tbls2[[n1]]$Ymax
        }
        sub1 <- GO_tbl[which(GO_tbl$test),]
        sub2 <- sub1[which(sub1$Y == 0),]
        sub1 <- sub1[which(sub1$Y > 0),]
        dotTtl <- paste0(title.root, n2)
        aes <- data.frame(x = "\`Z-score\`", y = "Y", size = "Count", alpha = "Y", colour = "Ontology")
        #non.aes <- data.frame(alpha = 0.1)
        if (plotly) {
          aes2 <- aes
          aes2$text1 <- "Label"
          aes2$text2 <- "Count"
          aes2$text3 <- "Y"
          aes2$text4 <- "Mapping"
        }
        pluses <- c(paste0("ggplot2::ylab(\"-log10(", c("", "adj. ")[P_adjust+1], "Pvalue)\")"),
                    "ggplot2::scale_radius(range = scrange, guide = \"none\")",
                    "ggplot2::ggtitle(dotTtl)", "ggplot2::facet_wrap(~Ontology)", "ggplot2::theme_bw()",
                    c(paste0("ggplot2::xlab(\"", c("(N(Up) - N(Down))/sqrt(Total)",
                                                   "Z-score")[True_Zscore+1], "\")")),
                    paste0("ggplot2::xlim(", Xmin, ", ", Xmax, ")"),
                    paste0("ggplot2::ylim(0, ", Ymax, ")"))
        aes <- paste(sapply(1:ncol(aes), function(x) { paste(colnames(aes)[x], aes[x], sep = " = ") }), collapse = ", ")
        #non.aes <- paste(sapply(1:ncol(non.aes), function(x) { paste(colnames(non.aes)[x], non.aes[x], sep = " = ") }), collapse = ", ")
        #non.aes <- gsub("^dummy = NA, ", "", non.aes)
        pluses <- paste(pluses, collapse = " + ")
        plot.txt <- paste0("plot <- ggplot2::ggplot(GO_tbl) + ggplot2::geom_point(shape = 16, ggplot2::aes(", aes, "))")
        if ((exists("non.aes", inherits = FALSE))&&(length(non.aes))) {
          plot.txt <- paste0(gsub("\\)$", ", ", plot.txt), non.aes, ")")
        }
        if ((exists("pluses", inherits = FALSE))&&(length(pluses))) {
          plot.txt <- paste0(plot.txt, " + ", pluses)
        }
        #cat(plot.txt)
        #options(warn = -1)
        suppressWarnings(eval(parse(text = plot.txt)))
        plot <- plot +
          ggplot2::geom_hline(data = thresh, ggplot2::aes(yintercept = `-log10(Threshold)`, colour = Colour)) +
          ggplot2::geom_text(data = thresh, ggplot2::aes(label = Label, y = `-log10(Threshold)`+Ymax*0.01, colour = Colour),
                             x = Xmin+Xbreadth*0.2, hjust = 0, cex = 1.8) +
          ggplot2::guides(colour = "none", alpha = "none")
        if (nrow(sub1)) {
          plot <- plot + textFun(data = sub1, ggplot2::aes(label = Label2, x = `Z-score`, y = Y, alpha = Y),
                                 cex = cex, lineheight = lineheight)
        }
        if (nrow(sub2)) {
          plot <- plot + textFun(data = sub2, ggplot2::aes(label = Label2, x = `Z-score`, y = Y),
                                 cex = cex, lineheight = lineheight)
        }
        if (length(winf)) { plot <- plot + ggplot2::geom_hline(yintercept = Ymax, colour = "black", linetype = "dotted") }
        if (show) { proteoCraft::poplot(plot, 12, 20) }        
        GOplts[[paste0("GO bubble plot - ", n1)]] <- plotEval(plot)
        nm <- gsub("/|:|\\*|\\?|<|>|\\|", "-", dotTtl)
        if (nchar(nm) > 98) { nm <- substr(nm, 1, 98) }
        if (nm %in% grphs) {
          fixkount <- 0
          while (nm %in% grphs) {
            fixkount <- fixkount + 1
            if (fixkount == 100) {
              stop("Really? Really?!?! You really have 100 similarly named conditions with very long names?!?! If you tried to break this function then you succeeded with (bad) style!")
            }
            nm <- paste0(substr(nm, 1, 93), "...", c("0", "")[(nchar(fixkount) > 1)+1], fixkount)
          }
        }
        grphs <- c(grphs, nm)
        if ((length(save) > 1)||(save != FALSE)) {
          for (sv in save) {
            if (subfolderpertype) { sfpt <- paste0(subfolder, "/", sv) } else { sfpt <- subfolder }
            if (!dir.exists(sfpt)) { dir.create(sfpt, recursive = TRUE) }
            suppressMessages({
              if (sv %in% c("jpeg", "tiff", "png", "bmp")) { #Note: tiff does not seem to work currently!
                ggplot2::ggsave(paste0(sfpt, "/", nm, ".", sv), plot, dpi = 300, width = 10, height = 10, units = "in")
              } else {
                ggplot2::ggsave(paste0(sfpt, "/",nm, ".", sv), plot)
              }
            })
          }
        }
        if (plotly) {
          aes2 <- paste(sapply(1:ncol(aes2), function(x) { paste(colnames(aes2)[x], aes2[x], sep = " = ") }), collapse = ", ")
          plot.txt2 <- paste0("plot2 <- ggplot2::ggplot(GO_tbl) + ggplot2::geom_point(ggplot2::aes(", aes2, "), ", # non.aes, 
                              ") + ggplot2::guides(colour = \"none\") + ",
                              pluses)
          suppressWarnings(eval(parse(text = plot.txt2)))
          if (length(winf)) {
            plot2 <- plot2 +
              ggplot2::geom_hline(yintercept = Ymax, colour = "black", linetype = "dotted")
          }
          plot_ly <- plotly::ggplotly(plot2, tooltip = c("text1", "text2", "text3", "text4"))
          setwd(plotly_subfolder)
          htmlwidgets::saveWidget(plot_ly, paste0(nm, ".html"), selfcontained = TRUE)
          setwd(origWD)
        }
        if (bars) {
          col_lim <- 3
          GO_tbl2 <- GO_tbl
          GO_tbl2$X <- NULL
          GO_tbl3 <- data.frame()
          wK <- which(GO_tbl2$Count > 0)
          if (length(wK)) {
            wOnt <- unique(GO_tbl2$Ontology[wK])
            GO_tbl3 <- lapply(wOnt, function(ont) { #ont <- wOnt[1]
              tmp <- GO_tbl2[wK,][which(GO_tbl2$Ontology[wK] == ont),]
              if (nrow(tmp)) {
                tmp <- tmp[order(tmp[[paste0(c("", "adj. ")[P_adjust+1], "Pvalue")]], decreasing = FALSE),]
                tmp <- tmp[1:min(c(nrow(tmp), MaxTerms_bar)),]
                tmp <- tmp[order(tmp$`Z-score`, decreasing = FALSE),]
                tmp$"Z-score*" <- tmp$`Z-score`
                tmp$"Z-score*"[which(tmp$"Z-score*" > col_lim)] <- col_lim
                tmp$"Z-score*"[which(tmp$"Z-score*" < -col_lim)] <- -col_lim
                tmp$X <- 1:nrow(tmp)
                tmp[, c("Label", "Label2", "Label3")] <- GO_tbl[match(tmp$ID, GO_tbl$ID), c("Label", "Label2", "Label3")]
                return(tmp)
              } else {
                return()
              }
            })
            GO_tbl3 <- plyr::rbind.fill(GO_tbl3)
          }
          if (nrow(GO_tbl3)) {
            GO_tbl3$Ontology <- factor(GO_tbl3$Ontology, levels = Ont)
            barTtl <- paste0(bars_title.root, n2)
            if (P_adjust) {
              thr <- lapply(unique(GO_tbl3$Ontology), function(x) {
                data.frame(Y = -log10(GO_FDR),
                           Colour = colorRampPalette(c("red", "gold"))(length(GO_FDR)),
                           Label = paste0(GO_FDR*100, "%"))
              })
              thr <- proteoCraft::listMelt(thr, unique(GO_tbl3$Ontology))
              colnames(thr)[which(colnames(thr) == "L1")] <- "Ontology"
              colnames(thr)[which(colnames(thr) == "value")] <- "Y"
              thr$variable <- NULL
            } else {
              thr <- thresh
              colnames(thr)[which(colnames(thr) == "-log10(Threshold)")] <- "Y"
            }
            xmx <- suppressWarnings(max(GO_tbl3$X))+1
            zMsg <- c("", paste0("\n(truncated at ", col_lim, ")"))[(max(abs(GO_tbl3$`Z-score`)) > col_lim) + 1]
            barplot_txt <- paste0("PLOT <- ggplot2::ggplot(GO_tbl3) +
                  ggplot2::geom_bar(stat = \"identity\", ggplot2::aes(x = X, y = Y, fill = `Z-score*`INSERT1)) +
                  ggplot2::facet_grid(Ontology~., switch = \"y\") +
                  ggplot2::xlim(-1, ", xmx, ") +
                  ggplot2::ylim(0, ", Ymax, ") +
                  ggplot2::labs(fill = \"", c("(N(Up) - N(Down))/sqrt(Total)",
                                              "Z-score")[True_Zscore+1], zMsg, "\") +
                  ggplot2::scale_fill_gradient2(low = \"green\", mid = \"grey\", high = \"red\", midpoint = 0, limits = c(-col_lim, col_lim)) +
                  ggplot2::scale_y_continuous(expand = c(0,0)) +
                  ggplot2::ylab(\"-log10(", c("", "adj. ")[P_adjust+1], "Pvalue)\") + ggplot2::ggtitle(barTtl) + ggplot2::theme_bw() +
                  ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
                        axis.text.x = ggplot2::element_blank(),
                        axis.ticks = ggplot2::element_blank(),
                        axis.title.x = ggplot2::element_blank(),
                        strip.text.y = ggplot2::element_text(angle = 90)) +
                  ggplot2::coord_cartesian(clip = \"off\")")
            text_tmp <- gsub("INSERT1", "", gsub("^PLOT", "barplot1", barplot_txt))
            #cat(text_tmp)
            suppressWarnings(eval(parse(text = text_tmp)))
            barplot1 <- barplot1 +
              ggplot2::geom_text(ggplot2::aes(label = Label3, x = X), y = -Ymax/20, hjust = 1, angle = 60, cex = 3) +
              ggplot2::theme(panel.spacing = ggplot2::unit(9, "lines"),
                             plot.margin = ggplot2::unit(c(1, 1, 10, 1), "lines"))
            barplot1 <- barplot1 +
              ggplot2::geom_hline(data = thr, ggplot2::aes(yintercept = Y, colour = Colour), linetype = "dashed") +
              ggplot2::geom_text(data = thr, ggplot2::aes(label = Label, y = Y-Ymax*0.02, colour = Colour), x = -1, hjust = 0, cex = 2.5) +
              ggplot2::guides(colour = "none")
            if (show) { proteoCraft::poplot(barplot1, 12, 20) }
            GOplts[[paste0("GO bar plot - ", n1)]] <- plotEval(barplot1)
            barnm <- gsub("/|:|\\*|\\?|<|>|\\|", "-", barTtl)
            if (nchar(barnm) > 98) { barnm <- substr(barnm, 1, 98) }
            if (barnm %in% grphs) {
              fixkount <- 0
              while (barnm %in% grphs) {
                fixkount <- fixkount + 1
                if (fixkount == 100) {
                  stop("Really? Really?!?! You really have 100 similarly named conditions with very long names?!?! If you tried to break this function then you succeeded with (bad) style!")
                }
                barnm <- paste0(substr(barnm, 1, 93), "...", c("0", "")[(nchar(fixkount) > 1)+1], fixkount)
              }
            }
            grphs <- c(grphs, barnm)    
            if ((length(save) > 1)||(save != FALSE)) {
              for (sv in save) {
                if (subfolderpertype) { sfpt <- paste0(subfolder, "/", sv) } else { sfpt <- subfolder }
                if (!dir.exists(sfpt)) { dir.create(sfpt, recursive = TRUE) }
                setwd(sfpt)
                suppressMessages({
                  if (sv %in% c("jpeg", "tiff", "png", "bmp")) { #Note: tiff does not seem to work currently!
                    ggplot2::ggsave(paste0(sfpt, "/",barnm, ".", sv), barplot1, dpi = 300, width = 10, height = 10, units = "in")
                  } else {
                    ggplot2::ggsave(paste0(sfpt, "/",barnm, ".", sv), barplot1)
                  }
                })
              }
            }
            if (plotly) {
              text_tmp <- gsub("INSERT1",
                               ", text1 = Label, text2 = Count, text3 = Y, text4 = Mapping",
                               gsub("^PLOT", "barplot2", barplot_txt))
              #cat(text_tmp)
              suppressWarnings(eval(parse(text = text_tmp)))
              barplot_ly <- plotly::ggplotly(barplot2, tooltip = c("text1", "text2", "text3", "text4"))
              setwd(plotly_subfolder)
              htmlwidgets::saveWidget(barplot_ly, paste0(barnm, ".html"), selfcontained = TRUE)
              setwd(origWD)
            }
          }
        }
        if (graph) {
          setwd(subfolder)
          for (ont in Ont[Wh1]) {
            tst <- try(topGO::printGraph(GOdata[[ont]], resultFisher[[ont]], firstSigNodes = 5, useInfo = 'all',
                                         fn.prefix = paste0(n2, "_-_", ont, "_-"), pdfSW = TRUE), silent = TRUE)
            fl <- paste0(n2, "_-_", ont, "_-_classic_5_all.pdf")
            #system(paste0("open \"", fl, "\""))
            if (!"try-error" %in% class(tst)) {
              file.copy(fl, paste0("GO_", ont,"_network - ", n2, ".pdf"), overwrite = TRUE)
            }
            if (file.exists(fl)) { unlink(fl, force = TRUE) }
          }
          setwd(origWD)
        }
        res <- list(GO_plots = GOplts)
        return(res)
      }
      #environment(plotsF0) <- .GlobalEnv # Only needed if code run as function!
      cat("     Drawing plots...\n")
      tst <- try({ tstPlots <- setNames(parallel::parLapply(parClust, names(GO_tbls2), plotsF0),
                                        names(GO_tbls2)) }, silent = TRUE)
      if ("try-error" %in% class(tst)) {
        tstPlots <- setNames(lapply(names(GO_tbls2), plotsF0),
                             names(GO_tbls2))
      }
      lapply(names(GO_tbls2), function(flt) { #flt <- names(GO_tbls2)[1]
        nms <- names(tstPlots[[flt]]$GO_plots)
        GO_plots[nms] <<- tstPlots[[flt]]$GO_plots[nms]
      })
    }
  }
}
cat("     Done!\n")
setwd(origWD)
if (kount) {
  KK <- grep("^Counts - ", colnames(GO.terms), value = TRUE)
  wK <- which(rowSums(GO.terms[, KK, drop = FALSE]) > 0)
  goRES <- list(GO_terms = GO.terms[wK,],
                All_GO_terms = GO.terms,
                GO_plots = GO_plots)
  if (!P_adjust) { goRES$GO_FDR_thresholds <- GO_FDR_thresholds }
  if (plotly) { goRES[["GO_plotly"]] <- GO_plot_ly }
} else { goRES <- NA }
