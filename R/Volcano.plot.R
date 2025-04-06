#' Volcano.plot
#'
#' @description
#' A function to create volcano plots.
#' 
#' @param Prot The protein/protein groups table. It will be added columns saying whether a protein group is regulated or not for each test performed (i.e. each "aggregate")
#' @param mode One of "standard", "custom", or "curved". 
#' @param X.root Root of the names for X-axis log-scale values. If missing, defaults are provided in mode = "custom" in the parameters file.
#' @param Y.root Root of the names for Y-axis values (usually -log10(Pvalues)). If missing, default are provided in mode = "custom" in the parameters file.
#' @param X.normalized Is X normalized? Assumed to be TRUE by default.
#' @param experiments.map The experiments map.
#' @param aggregate.map The aggregate map.
#' @param aggregate.list The named list of aggregates.
#' @param aggregate.name The name of the aggregate to use.
#' @param parameters The experiment's parameters file.
#' @param save Should the plot(s) be saved? Default = FALSE. Set it to a vector of acceptable file extensions to save to the corresponding file format.
#' @param labels How to decide which point labels will be printed. Values are:\cr
#'  - "proteins", in which case parameters MaxLabels and MaxLabels_priority will also be used;\cr
#'  - "thresholds", in which case proteins above thresholds will be labelled;
#'  - "FDR", in which case only points that pass Benjamini-Hochberg FDR (or are in the list of proteins of interest) will be displayed.\cr
#'  - The default, "both" (default), which means "proteins + thresholds".
#' @param MaxLabels Maximum number of labels to print. Does not apply to matches to the list provided by "proteins". Default = 50 
#' @param MaxLabels_priority If not printing all labels, will we favour more extreme ratios (default = "X") or more significant values ("Y")?
#' @param proteins Optional, the character list of protein IDs to be highlighted.
#' @param proteins_split Default = FALSE. If set to TRUE, and if there are protein matching the provided list (see "proteins" argument), then these will be highlighted in a separate graph.
#' @param Proteins.col Default = "Protein IDs"; used for identifying matches to the "proteins" argument.
#' @param Ref.Ratio.values Used if the default parameter values for Ratio thresholds are to be replaced with dynamic, reference-based ones. Also requires a few arguments below.
#' @param Ref.Ratio.method Defines how vertical ratios are computed.\cr
#'  - "obs2": current default, will sort observations within sample groups and identify a threshold such that a proportion "ratios.FDR" of most extreme values are selected.\cr
#'  - "obs1": does almost the same as the above, but only for control-to-control or control-to-average-control ratios within the ratios group.\cr
#'  - "SD": same as the previous, but will calculate the threshold using qnorm(ratios.FDR, median(Ref.Ratio.values), sd(Ref.Ratio.values)), i.e. assuming normality.\cr
#' @param ratios.FDR The proportion of reference ratios to be accepted, if Ref.Ratio.values is provided. Default = 0.01.
#' @param FDR.thresh The significance thresholds calculated for each acceptable FDR value.
#' @param FDR.root The root of the columns with Significance for different FDRs. Default = "Significant-FDR="
#' @param arbitrary.lines Default = NULL. Data.frame of arbitrary lines to add to the graph. Must have 5 columns: "slope", "yintercept", "xintercept" (if vertical line, else NA), "colour", "label".
#' @param arbitrary.thresh Default = NULL, ignored if mode = FDR. Similar to the above, but will define only horizontal lines. The lowest will determine the threshold below which no labels are displayed (except when using the "proteins" argument).
#' @param return Should we return a column with significance levels? Default = FALSE.
#' @param return.plot Should we return the plot itself? Default = FALSE.
#' @param show.labels Default = FALSE. If set to TRUE, the popped-up plot will be the one with labels.
#' @param title It is possible to set a specific title. Failing that a default one will be provided based on the value of the mode argument.
#' @param title.root Optional argument to add a character string at the beginning of the file names when saving.
#' @param subfolder Name of the sub-folder where graphs are to be saved. Default = "Reg. analysis/t-tests"
#' @param subfolderpertype If set to TRUE (default), will save each type of graph in a dedicated sub-folder.
#' @param Symmetrical Are we also interested in the left half (down-regulated) or only the right one (up-regulated)? If missing and a value for "parameters" is provided, the function will attempt to infer it from parameters' "Type", otherwise it will default to TRUE.
#' @param Alpha Mapping of each dot's alpha level to a column. Default = 1 (no alpha variations)
#' @param Alpha.identity Do we need to fit the Alphas into a 0:1 scale? Default = FALSE
#' @param Alpha.labels Should labels also be affected by alpha? Default = FALSE
#' @param Alpha.max For when mapping alpha to a column. The highest allowed alpha. Default = 1
#' @param Alpha.min For when mapping alpha to a column. The weakest allowed alpha. Default = 0.01
#' @param Size Mapping of each dot's size to a column. Default = 1 (same size for all)
#' @param Size.max For when mapping size to a column. The largest allowed dot size. Default = 2
#' @param Size.min For when mapping size to a column. The smallest allowed dot size. Default = 0.01
#' @param cex Label text size. Default = 2
#' @param lineheight Label interline size. Default = 1
#' @param plotly Should we generate a plotly? TRUE or FALSE (default); if TRUE, the link will be automatically returned.
#' @param plotly_local Save plotly files as local html? Default = TRUE. If set to FALSE, a valid plotly username and API key will have to be provided.
#' @param plotly_username Plotly username, to be provided if plotly_local = FALSE
#' @param plotly_API_key Plotly API key, to be provided if plotly_local = FALSE
#' @param plotly_subfolder Default = "". Name of the plotly server sub-folder where plotly graphs are to be saved. Does not apply if plotly_local = TRUE
#' @param plotly_sharing One of "public", "private" or "secret". Note that which ones will work depends on the type of plotly subscription you have. Does not apply if plotly_local = TRUE
#' @param plotly_labels Labels to be displayed in the plotly tooltip. A character corresponding to valid names of Prot columns. Naming this vector provides a way to override the original column names: what will be displayed int he tooltip will be the name in this vector, not the column name itself.
#' @param Xlim Set fixed limits for the X axis.
#' @param Ylim Set fixed limits for the Y axis.
#' @param Contaminants Default = FALSE. If TRUE, contaminants will be plotted too.
#' @param N.clust A limit on the number of vCPUs to use. If left as NULL (default), uses the number of available clusters - 1, to a minimum of 1.
#' @param N.reserved Default = 1. Number of reserved vCPUs the function is not to use. Note that for obvious reasons it will always use at least one.
#' @param cl Already have a cluster handy? Why spend time making a new one, which on top of that may invalidate the old one. Just pass it along!
#' @param reg.root Classically, this function would output a "Regulated - ..." column per sample group, based on applying the thresholds (either provided directly - significance is read from the table using "FDR.root", and FC threshold are either applied directly using "arbitrary.thresh" or calculated based on "Ref.Ratio.values"). Providing this input allows the function to bypass this step and apply colors based on a pre-determined decision.
#'
#' @examples                 
#' PG <- Volcano.plot(Prot = PG, mode = "custom", experiments.map = Exp.map,
#'                    aggregate.map = Aggregate.map,
#'                    aggregate.name = Volcano.plots.Aggregate.Level$aggregate,
#'                    parameters = Param,
#'                    save = c("jpeg", "pdf"),
#'                    labels = "both",
#'                    proteins = prot.list,
#'                    return = TRUE,
#'                    title.root = "FDR-type ")
#' 
#' @export
#' 

Volcano.plot <- function(Prot,
                         mode,
                         X.root,
                         Y.root,
                         X.normalized = TRUE,
                         experiments.map,
                         aggregate.map,
                         aggregate.name,
                         aggregate.list,
                         parameters,
                         save = FALSE,
                         labels = "both",
                         MaxLabels = 50,
                         MaxLabels_priority = "X",
                         Ref.Ratio.values = NULL,
                         Ref.Ratio.method = "obs2",
                         ratios.FDR = 0.01,
                         FDR.thresh,
                         FDR.root = "Significant-FDR=",
                         arbitrary.lines = NULL,
                         arbitrary.thresh = NULL,
                         proteins = NULL,
                         proteins_split = FALSE,
                         Proteins.col = "Protein IDs",
                         return = FALSE,
                         return.plot = FALSE,
                         show.labels = "",
                         title = "",
                         title.root = "",
                         subfolder = "Reg. analysis/t-tests",
                         subfolderpertype = TRUE,
                         Symmetrical,
                         Alpha = 1,
                         Alpha.identity = FALSE,
                         Alpha.labels = FALSE,
                         Alpha.max = 1,
                         Alpha.min = 0.01,
                         Size = 1,
                         Size.max = 3,
                         Size.min = 0.01,
                         cex = 2,
                         lineheight = 1,
                         Xlim = Inf,
                         Ylim = Inf,
                         plotly = FALSE,
                         plotly_local = TRUE,
                         plotly_username,
                         plotly_API_key,
                         plotly_subfolder = "",
                         plotly_sharing = "secret",
                         plotly_labels,
                         Contaminants = FALSE,
                         N.clust,
                         N.reserved = 1,
                         cl,
                         reg.root # DO NOT give this one a default non-null value to avoid uncontrolled behavious when rerunning the function!
                         ) {
  TESTING <- FALSE
  #proteoCraft::DefArg(proteoCraft::Volcano.plot); Symmetrical <- TRUE; TESTING <- TRUE; cl <- parClust
  #
  #Prot = PG; mode = "custom"; experiments.map = Exp.map; X.root = paste0("Mean ", Prot.Rat.Root); Y.root = pvalue.col[which(pvalue.use)]; aggregate.map = Aggregate.map; aggregate.name = Volcano.plots.Aggregate.Level$aggregate; aggregate.list = Aggregate.list;parameters = Param; save = c("jpeg", "pdf"); labels = "FDR"; Ref.Ratio.values = Ref.Ratios; ratios.FDR = as.numeric(Param$Ratios.Contamination.Rates); FDR.thresh = FDR.thresholds; arbitrary.lines = arbitrary.thr; proteins = prot.list;  proteins_split = protsplit; return = TRUE;  return.plot = TRUE; title.root = "FDR-type ";  subfolder = "Reg. analysis/t-tests"; subfolderpertype = FALSE; Alpha = "Rel. log10(Peptides count)"; Size = "Av. log10 abundance";  Size.max = 2; plotly = create_plotly; plotly_local = create_plotly_local; FDR.thresh = FDR.thresholds
  # OR (F-test, proteins)
  #Prot = tmpData; mode = "custom"; experiments.map = contr; X.root = paste0("Mean ", ratRef); Y.root = paste0(F_Root, " - "); aggregate.map = aggr_dummy; aggregate.list = aggr_list_dummy; aggregate.name = "Contrast"; parameters = Param; save = c("jpeg", "pdf"); FDR.root = "mod. F-test Significant-FDR="; Ref.Ratio.values = refRat_F; ratios.FDR = as.numeric(Param$Ratios.Contamination.Rates); FDR.thresh = F_FDR.thresholds; arbitrary.lines = arbitrary.thr; proteins = prot.list; proteins_split = protsplit; Proteins.col = idCol; return = TRUE; return.plot = TRUE; title = "F-test volcano plot "; subfolder = ohDeer; subfolderpertype = FALSE; Symmetrical = TRUE; Alpha = Alpha; Size = "Rel. av. log10 abundance"; Size.max = 2; plotly = create_plotly; plotly_local = create_plotly_local; plotly_labels = plotlyLab; Ref.Ratio.method = paste0("obs", RefRat_Mode); reg.root = regRoot_F
  # OR (modified peptides)
  #Prot = ptmpep; mode = "custom"; experiments.map = Exp.map; X.root = paste0("Mean ", ptms.ratios.ref[length(ptms.ratios.ref)]); Y.root = pvalue.col[which(pvalue.use)]; aggregate.map = Aggregate.map; aggregate.list = Aggregate.list; aggregate.name = VPAL$aggregate; parameters = P; save = c("jpeg", "pdf"); labels = "FDR"; Ref.Ratio.values = PTMs_ref.ratios[[ptm]]; ratios.FDR = as.numeric(Param$Ratios.Contamination.Rates); FDR.thresh = PTMs_FDR.thresholds[[ptm]]; arbitrary.lines = arbitrary.thr; proteins = prot.list; Proteins.col = "Name"; return = TRUE; return.plot = TRUE; title = paste0(Ptm, " volcano plot_"); subfolder = dir[2]; subfolderpertype = FALSE; Symmetrical = TwoSided; Size = "Rel. av. log10 abundance"; Size.max = 2; plotly = create_plotly; plotly_local = create_plotly_local; plotly_labels = c(PepLabKol, paste0(Ptm, "-site"))
  # OR (SAINTexpress)
  #Prot = allSAINTs;Proteins.col = "Protein";mode = "custom";experiments.map = Exp.map;X.root = fcRt;Y.root = fdrRt;aggregate.map = Aggregate.map;aggregate.name = VPAL$aggregate;aggregate.list = Aggregate.list;parameters = Parma;save = c("jpeg", "pdf");labels = "thresholds";Ref.Ratio.values = Ref.Ratios;ratios.FDR = as.numeric(Param$Ratios.Contamination.Rates);arbitrary.lines = ArbThr;proteins = prot.list;proteins_split = protsplit;return = TRUE;return.plot = TRUE;title = "SAINTexpress volcano plot_";subfolder = "Reg. analysis/SAINTexpress";subfolderpertype = FALSE;Symmetrical = TwoSided;Alpha = "Av. log10 abundance";Size = "Av. log10 abundance";Size.max = 2;plotly = create_plotly;plotly_local = create_plotly_local;plotly_labels = labKol
  # OR (SSDs)
  #Prot = temp;mode = "custom";experiments.map = Exp.map;X.root = paste0("Mean ", SSD.Root);Y.root = SSD.Pval.Root;aggregate.map = Aggregate.map;aggregate.name = SubCellFracAggr2$aggregate;aggregate.list = Aggregate.list; parameters = Param;save = c("jpeg", "pdf"); labels = "FDR"; Ref.Ratio.values = RefSSDs;ratios.FDR = as.numeric(Param$Ratios.Contamination.Rates);FDR.thresh = SSD.thresh;FDR.root = "Signif. SSDs-FDR=";arbitrary.lines = arbitrary.thr;proteins = prot.list; proteins_split = protsplit;return = TRUE; return.plot = TRUE;title = "Sum of Squared Distance volcano plot_";subfolder = "Reg. analysis/Localisation";subfolderpertype = FALSE; Symmetrical = FALSE;Alpha = "Rel. log10(Peptides count)";Size = "Av. log10 abundance"; Size.max = 2;plotly = create_plotly; plotly_local = create_plotly_local; X.normalized = FALSE
  origWD <- getwd()
  if (TESTING) {
    # Note:
    # This is not a perfect alternative to missing but will work in most cases, unless x matches a function imported by a package 
    misFun <- function(x) { return(!exists(deparse(substitute(x)))) }
  } else { misFun <- missing }
  if (return.plot) {
    plotEval <- function(plot) { ggplotify::as.ggplot(ggplotify::as.grob(plot)) }
  }
  # Create cluster (some steps are slow otherwise)
  cleanUp <- FALSE
  if (misFun(cl)) {
    dc <- parallel::detectCores()
    if (misFun(N.reserved)) { N.reserved <- 1 }
    if (misFun(N.clust)) {
      N.clust <- max(c(dc-N.reserved, 1))
    } else {
      if (N.clust > max(c(dc-N.reserved, 1))) {
        warning("More cores specified than allowed, I will ignore the specified number! You should always leave at least one free for other processes, see the \"N.reserved\" argument.")
        N.clust <- max(c(dc-N.reserved, 1))
      }
    }
    cl <- parallel::makeCluster(N.clust, type = "SOCK")
    cleanUp <- TRUE
  }
  #
  RES <- NA
  #
  if (labels == "both") { labels <- c("proteins", "thresholds") }
  stopifnot(nrow(Prot) > 0, ncol(Prot) > 0, is.character(subfolder), is.logical(subfolderpertype))
  #
  subfolder <- normalizePath(gsub("^/|/$", "", gsub("\\\\", "/", subfolder)),
                             winslash = "/", mustWork = FALSE)
  if (subfolder %in% c("", paste0(LETTERS, ":"))) { subfolder <- getwd() }
  if (!dir.exists(subfolder)) {
    tst <- try(dir.create(subfolder), silent = TRUE)
    if ("try-error" %in% class(tst)) { subfolder <- getwd() } 
  }
  #
  if ((length(save) > 1)||((!is.logical(save))||(save == TRUE))) {
    if ((length(save) == 1)&&(is.logical(save))&&(save == TRUE)) {
      saveExt <- "jpeg" # default
    }
    saveExt <- unique(gsub("^jpg$", "jpeg", gsub("^\\.", "", tolower(save))))
    save <- TRUE
  } else { save <- FALSE}
  #
  if ((!is.integer(MaxLabels))||(MaxLabels < 0)) { MaxLabels <- 100 }
  if (!MaxLabels_priority %in% c("X", "Y")) { MaxLabels_priority <- "X" }
  if (misFun(Symmetrical)) {
    # Provide defaults
    if ((!misFun(parameters))&&("Two.sided" %in% colnames(parameters))&&(is.logical(parameters$Two.sided))) {
      Symmetrical <- parameters$Two.sided
    } else {
      if ((!misFun(parameters))&&("Type" %in% colnames(parameters))) {
        Symmetrical <- !gsub(" |_|-|\\.", "", toupper(parameters$Type)) %in% c("IP", "IMMUNOPRECIPITATION", "BIOID", "PULLDOWN")
      } else { Symmetrical <- TRUE }
    }
  }
  upColRg <- c("blue", "green")
  downColRg <- c("red", "orange")
  if ("FDR" %in% labels) {
    if (misFun(FDR.thresh)) { stop("\"FDR.thresh\" must be provided if labels = \"FDR\"!") }
    f <- grep(proteoCraft::topattern(FDR.root), colnames(Prot), value = TRUE)
    FDR.values <- as.numeric(unique(gsub("%( - .+)?$", "", gsub(proteoCraft::topattern(FDR.root), "", f))))
    g <- as.numeric(gsub("^Threshold-FDR=|%( - .+)?$", "", names(FDR.thresh)))
    tst <- unique(grepl(" - ", names(FDR.thresh)))
    FDR_table <- data.frame(FDR = g,
                            Thresholds = FDR.thresh)
    if ((length(tst) == 1)&&(tst)) {
      FDR_table$Sample <-  gsub("^Threshold-FDR=[1-9][0-9]*\\.*[0-9]*% - ", "", names(FDR.thresh))
    }
    #Automatically correct to a percentage format
    #if (sum(FDR.values < 1) == length(FDR.values)) { FDR.values <- FDR.values*100 } # Removed, this was causing issues.
    #if (sum(g <= 1) == length(g)) { g <- g*100 }
    # In addition, FDR values are now ALWAYS expressed as percentages anyway!
    FDR_table <- FDR_table[order(FDR_table$FDR, decreasing = TRUE),]
    FDR.values <- sort(FDR.values, decreasing = TRUE)
    fdr.col.up <- grDevices::colorRampPalette(upColRg)(length(FDR.values))
    FDR_table$fdr.col.up <- fdr.col.up[sapply(FDR_table$FDR, function(x) { which(FDR.values == x )})]
    fdr.col.down <- grDevices::colorRampPalette(downColRg)(length(FDR.values))
    FDR_table$fdr.col.down <- fdr.col.down[sapply(FDR_table$FDR, function(x) { which(FDR.values == x )})]
    fdr.col.line <- rainbow(n = length(FDR.values), start = 2/6, end = 1/6, v = 0.75)
    FDR_table$fdr.col.line <- fdr.col.line[sapply(FDR_table$FDR, function(x) { which(FDR.values == x )})]
  }
  if (!Contaminants) {
    kontkol <- c("Contaminant", "Potential contaminant")
    kontkol <- kontkol[which(kontkol %in% colnames(Prot))]
    if (length(kontkol) != 1) {
      warning("I could not identify the contaminants column, did you already filter out contaminants?")
      Contaminants <- TRUE
    }
  }
  stopifnot(!misFun(aggregate.name), show.labels %in% c("", TRUE, FALSE),
            mode %in% c("standard", "custom", "curved"))
  if (plotly) {
    volcPlotly <- list()
    if (!plotly_local) {
      volcPlotly2 <- list()
      ArgTst1 <- !misFun(parameters)
      ArgTst2 <- !misFun(plotly_username)
      ArgTst3 <- !misFun(plotly_API_key)
      if ((ArgTst1)&&(!ArgTst2)) { plotly_username <- Param$Plotly_user_name }
      if ((ArgTst1)&&(!ArgTst3)) { plotly_API_key <- Param$Plotly_API_key }
      if (!plotly_sharing %in% c("public", "private", "secret")) {
        warning("\"plotly_sharing\" should be one of \"public\", \"private\" or \"secret\", defaulting to \"secret\"!")
        plotly_sharing <- "secret"
      }
      if (sum(c(!ArgTst2, !ArgTst3))) {
        warning("\"plotly\" is set to TRUE but plotly login details are missing, skipping!")
        plotly <- FALSE
      } else {
        if (nchar(plotly_subfolder)) {
          plotly_subfolder <- gsub(":|\\*|\\?|<|>|\\||/", "-", plotly_subfolder)
          plotly_subfolder <- paste0(gsub("/+$", "", plotly_subfolder), "/")
        }
        Sys.setenv("plotly_username" = plotly_username)
        Sys.setenv("plotly_api_key" = plotly_API_key)
      }
    }
  }
  if (nchar(title.root)) {
    if (!substr(title.root, nchar(title.root), nchar(title.root)) %in% c(" ", ".", "_", ".")) {
      title.root <- paste0(title.root, " ")
    }
  }
  if (title == "") {
    if (mode %in% c("standard", "custom", "curved")) { title <- "Ratios volcano plot " }
  }
  plotNms <- c()
  if (show.labels == "") {
    if (mode == "standard") { show.labels <- FALSE }
    if (mode %in%  c("standard", "custom", "curved")) { show.labels <- TRUE }
  }
  if (misFun(X.root)) {
    a <- as.data.frame(strsplit(unlist(strsplit(parameters$Plot.metrics, "; *")), ": *"))
    X.root <- paste0(gsub("[\\. -_]+$", "", a[2, which(a[1,] == "X")]), " - ")
  }
  dfltPM <- TRUE
  if (mode == "custom") {
    if (toupper(parameters$Plot.threshold.values) != "CURVED") {
      a1 <- set_colnames(as.data.frame(proteoCraft::Isapply(strsplit(unlist(strsplit(parameters$Plot.metrics, "; *")), ": *"), unlist)),
                         c("Axis", "Name"))
      if (misFun(Y.root)) {
        Y.root <- paste0(gsub("[\\. -_]+$", "", a$Name[which(a1$Axis == "Y")]), " - ")
      }
      if (parameters$Plot.threshold.metrics != "") {
        dfltPM <- FALSE
        Plot.metrics <- as.data.frame(strsplit(unlist(strsplit(parameters$Plot.threshold.metrics, "; *")), ": *"))
        Plot.metrics <- as.data.frame(t(Plot.metrics)) 
        rownames(Plot.metrics) <- NULL
        colnames(Plot.metrics) <- c("Levels", "Axis")
        Plot.metrics$Root <- c(X.root, Y.root)[match(Plot.metrics$Axis, a1$Axis)]
        a2 <- set_colnames(as.data.frame(t(sapply(strsplit(unlist(strsplit(parameters$Plot.threshold.values, "; *")), ": *"), unlist))),
                           c("Direction", "Text.value"))
        Plot.metrics$Text.value <- a2$Text.value[match(Plot.metrics$Levels, a2$Direction)]
        Plot.metrics$Value <- sapply(Plot.metrics$Text.value, function(x) { eval(parse(text = x)) })
        a3 <- set_colnames(as.data.frame(t(sapply(strsplit(unlist(strsplit(parameters$Plot.threshold.tests, "; *")), ": *"), unlist))),
                           c("Direction", "Test"))
        Plot.metrics$Test <- a3$Test[match(Plot.metrics$Levels, a3$Direction)]
        a4 <- set_colnames(as.data.frame(t(sapply(strsplit(unlist(strsplit(parameters$Plot.threshold.colours, "; *")), ": *"), unlist))),
                           c("Direction", "Colour"))
        Plot.metrics$Colour <- a4$Colour[match(Plot.metrics$Levels, a4$Direction)]
        a5 <- unlist(strsplit(parameters$Plot.areas.colours, "; *"))
        a5 <- as.data.frame(strsplit(a5, "[_:] *"))
        Plot.colours <- as.data.frame(sapply(Plot.metrics$Levels[which(Plot.metrics$Axis == "X")], function(x) {
          sapply(Plot.metrics$Levels[which(Plot.metrics$Axis == "Y")], function(y) {
            a5[3, which((a5[1,] == x)&(a5[2,] == y))]
          })
        }))
      }
    }
  }
  if (dfltPM) {
    # Vertical (ratio) thresholds
    Plot.metrics <- data.frame(Levels = c("down", "up"),
                               Axis = "X",
                               Name = X.root,
                               Test = c("<=", ">="),
                               Colour = c("red", "red"),
                               Text.value = c("log2(1/2)", "log2(2)"))
    Plot.metrics$Value <- sapply(Plot.metrics$Text.value, function(x) { eval(parse(text = x)) })
    # NB: These defaults may be overwritten by the "Ref.Ratio.values" argument for each plot.
    # Horizontal (P-value) thresholds
    if ((mode == "standard")&&(!misFun(arbitrary.thresh))) {
      tmpmetr <- data.frame(Levels = arbitrary.thresh$label,
                            Axis = "Y",
                            Name = "Significance",
                            Value = arbitrary.thresh$yintercept,
                            Test = ">=",
                            Colour = arbitrary.thresh$colour)
    } else {
      tmpmetr <- data.frame(Levels = c("strict", "loose"),
                            Axis = "Y",
                            Name = Y.root,
                            Value = c(-log10(0.01), -log10(0.05)),
                            Test = ">=",
                            Colour = c("orange", "gold"))
    }
    tmpmetr$Text.value <- as.character(tmpmetr$Value)
    Plot.metrics <- rbind(Plot.metrics, tmpmetr)
    Plot.metrics <- Plot.metrics[which(Plot.metrics$Levels != "Dummy"),]
    w <- which(Plot.metrics$Axis == "Y")
    Plot.colours <- data.frame(up = grDevices::colorRampPalette(downColRg)(length(w)),
                               down = grDevices::colorRampPalette(upColRg)(length(w)))
    rownames(Plot.colours) <- Plot.metrics$Levels[w]
  }
  if (misFun(X.root)) { stop("Argument \"X.root\" is missing, investigate!") }
  if (misFun(Y.root)) { stop("Argument \"Y.root\" is missing, investigate!") }
  #X <- gsub("\\.$", "", X.root)
  #Y <- gsub("\\.$", "", Y.root)
  B <- aggregate.name
  A <- aggregate.list[[B]]
  if (nchar(B) == 3) { B <- aggregate.map$Characteristics[[which(aggregate.map$Aggregate.Name == B)]] }
  if (length(Symmetrical) == 1) { Symmetrical <- rep(Symmetrical, length(A)) }
  if (length(Symmetrical) != length(A)) {
    warning("This length of the \"Symmetrical\" argument is incorrect, check it!")
    Symmetrical <- rep(Symmetrical, length(A))
  }
  names(Symmetrical) <- A
  A <- A[which(A %in% experiments.map[[B]])]
  Symmetrical <- Symmetrical[A]
  xKols <- setNames(paste0(X.root, A), A)
  yKols <- setNames(paste0(Y.root, A), A)
  PorQ <- rev(unlist(strsplit(toupper(gsub("-?value.*", "", Y.root, ignore.case = "")), "")))[1]
  if (!PorQ %in% c("P", "Q")) { PorQ <- "P" }
  if ((proteins_split)&&((is.null(proteins))||(length(proteins) == 0))) {
    proteins_split <- FALSE
  }
  if ((length(labels) == 1)&&(labels == "proteins")&&(proteins_split)) {
    warning("Argument \"proteins_split\" will be ignored since argument \"labels\" is set to \"proteins\".")
    proteins_split <- FALSE
  } # No splitting if we are already only labeling proteins in list.
  if (misFun(plotly_labels)) { # Create those even if plotly off!
    plotly_labels <- setNames(c("Labels", Proteins.col, "Genes"),
                              c("Name(s)", "Accession(s)", "Gene(s)")) # (column "Labels" gets created later, before we use those columns: should not break!)
    plotly_labels <- plotly_labels[which(plotly_labels %in% colnames(Prot))]
  } else {
    wlabkol <- which(plotly_labels %in% colnames(Prot))
    if ((plotly)&&(length(wlabkol) < length(plotly_labels))) {
      warning("Some plotly_labels columns are missing, we won't be able to display them in the tooltip!")
      plotly_labels <- plotly_labels[wlabkol]
    }
    if (is.null(names(plotly_labels))) { names(plotly_labels) <- plotly_labels }
    w <- which(is.na(names(plotly_labels)))
    if (length(w)) { names(plotly_labels)[w] <- plotly_labels[w] }
    w <- which(names(plotly_labels) == "")
    if (length(w)) { names(plotly_labels)[w] <- plotly_labels[w] }
  }
  #
  # Default colors
  myColors <- setNames(c("black", "black", "purple", c("brown", "firebrick1")[proteins_split+1]),
                       c("non significant", "too small FC", "target", "protein in list"))
  #
  # Filter samples for available valid data
  tstTbl <- data.frame(Sample = A,
                       logFC = xKols %in% colnames(Prot),
                       PVal = yKols %in% colnames(Prot))
  tstTbl$All_OK <- tstTbl$logFC & tstTbl$PVal
  if (!misFun(reg.root)) {
    labels <- c("regulated",
                labels[which(labels == c("proteins"))])
    regKols <- setNames(paste0(reg.root, A), A)
    tstTbl$Reg <- regKols %in% colnames(Prot)
    tstTbl$All_OK <- tstTbl$All_OK & tstTbl$PVal
    #
    # Define colors for existing levels
    tmp <- unique(unlist(Prot[, regKols[which(tstTbl$Reg)]]))
    up_Nms <- grep("^up", tmp, value = TRUE)
    upSp <- grep("^Specific", tmp, value = TRUE)
    dwn_Nms <- grep("^down", tmp, value = TRUE)
    downSp <- grep("^Anti-specific", tmp, value = TRUE)
    if (length(up_Nms)) {
      up_Nms <- data.frame(Up = up_Nms)
      if (sum(c("FDR", "regulated") %in% labels)) {
        up_Nms$Val <- as.numeric(gsub(".* FDR = |%$", "", up_Nms$Up))
      } else {
        up_Nms$Val <- as.numeric(gsub("^[^0-9]+|[^0-9]+$", "", up_Nms$Up))
      }
      if (is.numeric(up_Nms$Val)) { up_Nms <- up_Nms$Up[order(up_Nms$Val, decreasing = TRUE)] }
      if (length(upSp)) {
        # This category could be plotted if we have imputed!
        up_Nms <- c(up_Nms, "Specific")
        # We will simplify this category, but not in the original Prot input, since it will be returned as output!
      }
      lUp <- length(up_Nms)
      upCol <- grDevices::colorRampPalette(upColRg)(lUp)
      myColors[up_Nms] <- upCol
    }
    if (length(dwn_Nms)) {
      dwn_Nms <- data.frame(Down = dwn_Nms)
      if (sum(c("FDR", "regulated") %in% labels)) {
        dwn_Nms$Val <- as.numeric(gsub(".* FDR = |%$", "", dwn_Nms$Down))
      } else {
        dwn_Nms$Val <- as.numeric(gsub("^[^0-9]+|[^0-9]+$", "", dwn_Nms$Down))
      }
      if (is.numeric(dwn_Nms$Val)) { dwn_Nms <- dwn_Nms$Down[order(dwn_Nms$Val, decreasing = TRUE)] }
      if (length(downSp)) {
        # This category could be plotted if we have imputed!
        dwn_Nms <- c(dwn_Nms, "Anti-specific")
        # We will simplify this category, but not in the original Prot input, since it will be returned as output!
      }
      lDown <- length(dwn_Nms)
      downCol <- grDevices::colorRampPalette(downColRg)(lDown)
      myColors[dwn_Nms] <- downCol
    }
  }
  wOK <- which(tstTbl$All_OK)
  if (length(wOK)) {
    Wych <- setNames(lapply(wOK, function(x) {
      wych <- (proteoCraft::is.all.good(as.numeric(Prot[[xKols[x]]]), 2))&(proteoCraft::is.all.good(as.numeric(Prot[[yKols[x]]]), 2))
      if (!Contaminants) { wych <- wych & (Prot[[kontkol]] != "+") }
      return(which(wych))
    }), A[wOK])
    tstTbl$All_OK[wOK] <- sapply(Wych, length) > 0
  }
  wOK <- which(tstTbl$All_OK)
  if (!length(wOK)) { stop("Not a single sample group with valid statistics detected!\nCheck inputs or statistical analysis") }
  A <- A[wOK]
  xKols <- xKols[A]
  yKols <- yKols[A]
  if (!misFun(reg.root)) { regKols <- regKols[A] }
  Symmetrical <- Symmetrical[A]
  Wych <- Wych[A]
  #
  # Define global x-limits
  tmp <- Prot[, xKols, drop = FALSE]
  xlim <- proteoCraft::is.all.good(as.numeric(unlist(tmp)))
  if (X.normalized) { xlim <- max(abs(xlim))*c(-1, 1) } else { xlim <- c(min(xlim), max(xlim)) }
  xspan <- xlim[2]-xlim[1]
  xlim <- xlim + xspan*c(-0.05, 0.05)
  #
  # Define global y-limits
  tmp <- Prot[, yKols, drop = FALSE]
  ylim <- proteoCraft::is.all.good(as.numeric(unlist(tmp)))
  ylim <- max(abs(ylim))*1.1
  #
  # Labels
  Prot$Labels <- Prot[[parameters$Plot.labels]]
  weech <- which(Prot$Labels != "")
  colchar <- 25
  clusterExport(cl, "colchar", envir = environment())
  f0 <- function(x) {#x <- strsplit(Prot$Labels[weech], "  ?")[1]
    x <- unlist(x)
    if (length(x) > 8) { x <- c(x[1:8], "...") } # This is because:
    # a) sentences of too many words can cause problems with subsequent combinatorial steps, slowing down the script or even causing it to fail
    # and
    #b) do you really expect a protein label more than 8 words long to be helpful?
    # Usually this will not be an issue with Uniprot but can be for other databases, e.g. TAIR, where protein names contain additional information for unknown proteins.
    nch <- sapply(1:length(x), function(y) { nchar(paste(x[1:y], collapse = " ")) })
    nc <- min(ceiling((sum(nchar(x)) + length(x) - 1)/colchar), length(x))
    if (nc > 1) {
      tstbrk <- cbind(0, gtools::combinations(length(x)-1, nc-1), length(x))
      tstbrk <- apply(tstbrk, 1, function(y) {
        sapply(1:nc, function(z) { paste(x[(y[z]+1):y[z+1]], collapse = " ") })
      })
      tstsd <- apply(tstbrk, 2, function(y) { sum(nchar(y)^2) })
      label <- paste(tstbrk[, which(tstsd == min(tstsd))[1]], collapse = "\n")
    } else { label <- paste(x, collapse = " ") }
    return(label)
  }
  environment(f0) <- .GlobalEnv
  if (length(weech)) {
    Prot$Labels[weech] <- parSapply(cl, strsplit(Prot$Labels[weech], "  ?"), f0)
  }
  weech <- which(Prot$Labels == "")
  id.col <- unique(c(Proteins.col, "Common Names", "Common Name (short)", "Names", "Name", "Genes", "Gene", "Code"))
  id.col <- c(id.col, tolower(id.col), toupper(id.col))
  id.col <- id.col[which(id.col %in% colnames(Prot))]
  if (length(weech)) {
    kount <- 0
    while ((length(weech))&&(kount < length(id.col))) {
      kount <- kount+1
      Prot$Labels[weech] <- sapply(Prot[weech, id.col[kount]], function(x) {
        if (x == "") { res <- "" } else {
          x <- unlist(strsplit(x, ";"))
          if (length(x) == 1) { res <- x } else { res <- paste0(x[1], ";...") }
        }
        return(res)
      })
      weech <- which(Prot$Labels == "")
    }
    if (length(weech)) { Prot$Labels[weech] <- paste0("Default_label_", c(1:length(weech))) }
  }
  #
  plotMetr.lst <- list()
  Plots <- list(Unlabelled = list())
  if (show.labels) { Plots$Labelled <- list() }
  #
  for (i in A) { #i <- A[1]
    symm <- Symmetrical[i]
    prot_split <- proteins_split
    e <- paste0(c(X.root, Y.root), i)
    e <- e[which(e %in% colnames(Prot))]
    plot.metrics <- Plot.metrics
    plot.metrics$Name <- ""
    plot.colours <- Plot.colours
    if (!symm) {
      plot.metrics <- plot.metrics[which(plot.metrics$Levels != "down"),]
      plot.colours <- plot.colours[, which(colnames(plot.colours) != "down"), drop = FALSE]
    }
    temp <- Prot[, "Labels", drop = FALSE]
    if (!Contaminants) { temp[[kontkol]] <- Prot[[kontkol]] }
    if ("Genes" %in% colnames(Prot)) { temp$Genes <- Prot$Genes }
    if (length(id.col)) { temp[, id.col] <- Prot[, id.col] }
    temp[, plotly_labels] <- Prot[, plotly_labels] # Even if plotly off please!
    temp$X <- Prot[[e[1]]]
    temp$Y <- Prot[[e[2]]]
    temp <- temp[Wych[[i]],]
    if ("FDR" %in% labels) {
      fdr_table <- FDR_table
      if ("Sample" %in% colnames(fdr_table)) {
        fdr_table <- fdr_table[which(fdr_table$Sample == i),]
      }
      fdr.values <- fdr_table$FDR
      f1 <- paste0(FDR.root, c(fdr.values, fdr.values/100), "%")
      if ("Sample" %in% colnames(fdr_table)) {
        f1 <- paste0(f1, " - ", i)
      }
      f1 <- f1[which(f1 %in% colnames(Prot))]
      f2 <- paste0("FDR=", fdr.values, "%")
      temp[, f2] <- Prot[Wych[[i]], f1]
      if (length(fdr.values) > 1) {
        test <- apply(temp[, rev(f2)], 1, function(x) {
          c(paste0("FDR ", sort(fdr.values), "%"), "")[which(c(x, "+") == "+")[1]]
        })
      } else {
        test <- sapply(temp[[f2]], function(x) {
          c(paste0("FDR ", sort(fdr.values), "%"), "")[which(c(x, "+") == "+")[1]]
        })
      }
      temp$FDR <- test
      for (f in f2) { temp[[f]] <- NULL }
    }
    if (!is.numeric(Alpha)) {
      if (!Alpha %in% colnames(Prot)) {
        warning("Alpha levels are not numeric or mapped to a column and will be ignored!")
        Alpha <- 1
      } else {
        temp$Alpha <- Prot[Wych[[i]], Alpha]
        if (!Alpha.identity) {
          Alpha2 <- paste0("Alpha mapped to: ", Alpha)
          if ((!is.numeric(Alpha.min))||(Alpha.min < 0)) { Alpha.min <- 0 }
          if ((!is.numeric(Alpha.max))||(Alpha.max > 1)) { Alpha.max <- 1 }
          temp$Alpha[which((temp$Alpha > 0)&(is.infinite(temp$Alpha)))] <- max(proteoCraft::is.all.good(temp$Alpha))
          temp$Alpha[which(!proteoCraft::is.all.good(temp$Alpha, 2))] <- min(proteoCraft::is.all.good(temp$Alpha))
          temp$Alpha <- Alpha.min+(temp$Alpha-min(temp$Alpha))*(Alpha.max-Alpha.min)/(max(temp$Alpha)-min(temp$Alpha))
        }
      }
    } else {
      if ((Alpha > 1)||(Alpha < 0)) {
        warning("The Alpha parameter is not between 0 and 1 and will be ignored!")
        Alpha <- 1
      }
    }
    if (!is.numeric(Size)) {
      if (!Size %in% colnames(Prot)) {
        warning("Sizes are not numeric or mapped to a column and will be ignored!")
        Size <- 1
      } else {
        temp$Size <- Prot[Wych[[i]], Size]
        if (!is.numeric(Size.min)) { Size.min <- 0.01 }
        if (!is.numeric(Size.max)) { Size.max <- 3 }
        #if (!is.numeric(Size.min)) { Size.min <- 0 }
        #if (!is.numeric(Size.max)) { Size.max <- ceiling(max(proteoCraft::is.all.good(temp$Size))) }
        temp$Size[which((temp$Size > 0)&(is.infinite(temp$Size)))] <- max(proteoCraft::is.all.good(temp$Size))
        temp$Size[which(!proteoCraft::is.all.good(temp$Size, 2))] <- min(proteoCraft::is.all.good(temp$Size))
        temp$Size <- Size.min+(temp$Size-min(temp$Size))*(Size.max-Size.min)/(max(temp$Size)-min(temp$Size))
      }
    }
    if ("Target" %in% colnames(experiments.map)) {
      target <- unique(unlist(strsplit(experiments.map$Target[which(experiments.map[[aggregate.name]] == i)], ";")))
      target <- target[which(!target %in% c("", "NA", NA))]
      target <- unique(gsub("^CON_+", "", target))
    }
    i2 <- proteoCraft::cleanNms(i)
    if (!misFun(reg.root)) {
      rgKol <- regKols[i]
      temp$Colour <- Prot[Wych[[i]], rgKol]
      # Important to simplify the specific/anti-specific categories:
      temp$Colour <- gsub("^Specific", "Specific", temp$Colour)
      temp$Colour <- gsub("^Anti-specific", "Anti-specific", temp$Colour)
    } else {
      temp$Colour <- "non significant" # Default Colour values
    }
    w.u <- which(plot.metrics$Levels == "up")
    if (symm) { w.d <- which(plot.metrics$Levels == "down") }
    if (!is.null(Ref.Ratio.values)) {
      x <- data.frame(value = sort(unlist(Ref.Ratio.values[[i]])))
      if (!symm) { x <- x[which(x$value > 0), , drop = FALSE] }
      if (X.normalized) {
        mx <- 0
      } else {
        mx <- median(x$value)  
      }
      R.thresh <- c("Upper" = NA)
      if (Ref.Ratio.method == "SD") {
        offset <- c(median(proteoCraft::is.all.good(temp$X)), 0)[X.normalized + 1]
        sdx <- sd(x$value)
        R.thresh[["Upper"]] <- qnorm(1-ratios.FDR, m = mx, sd = sdx)
        if (symm) { R.thresh[["Lower"]] <- qnorm(ratios.FDR, m = mx, sd = sdx) }
      }
      if (Ref.Ratio.method %in% paste0("obs", c("1", "2"))) {
        # Below, this is not exactly the same thing here as doing quantile(x$value, c(ratios.FDR/2, 1-ratios.FDR/2))
        # but close enough, and better suited to the language used:
        # We are not removing the lower and upper ratios.FDR/2 quantiles, but instead removing the ratios.FDR quantile in absolute value
        # (i.e. the ratios.FDR proportion is applied when considering both tails together, the tail may individually contribute a different proportion).
        # The implication is that we should get symmetrical thresholds - consistent with e.g. limma decideTests applying an absolute lfc threshold
        # (I think this makes perfect sense).
        d <- abs(x$value-mx)
        d <- sort(d, decreasing = TRUE)
        if (ratios.FDR > 0) {
          d <- d[floor(length(d)*ratios.FDR)]
        } else {
          d <- max(d)+0.000001
        }
        R.thresh[["Upper"]] <- d+mx
        if (symm) { R.thresh[["Lower"]] <- -(d+mx) }
      }
      R.thresh.label <- paste0(signif(R.thresh, 3), " = ", c("upper", "lower")[1:(symm+1)], "-tail of ", ratios.FDR*100, "% ",
                               c("ctrl.", "within sample group")[match(Ref.Ratio.method, paste0("obs", c("1", "2")))],
                               " ratios")
      plot.metrics$Name[w.u] <- R.thresh.label[1]
      plot.metrics$Text.value[w.u] <- plot.metrics$Value[w.u] <- R.thresh[1]
      xlim[2] <- max(c(xlim[2], plot.metrics$Value[w.u]*1.1), na.rm = TRUE)  
      if (symm) {
        plot.metrics$Name[w.d] <- R.thresh.label[2]
        plot.metrics$Text.value[w.d] <- plot.metrics$Value[w.d] <- R.thresh[2]
        xlim[1] <- min(c(xlim[1], plot.metrics$Value[w.d]*1.1), na.rm = TRUE)
      }
    }
    if (misFun(reg.root)) {
      if ("FDR" %in% labels) {
        for (f3 in fdr.values) {
          temp$Colour[which(temp$FDR == paste0("FDR ", f3, "%"))] <- "too small FC"
          temp$Colour[which((proteoCraft::logical.op(temp$X,
                                               plot.metrics$Test[w.u],
                                               plot.metrics$Value[w.u]))
                            &(temp$FDR == paste0("FDR ", f3, "%")))] <- paste0("up, FDR = ", f3, "%")
          if (symm) {
            temp$Colour[which((proteoCraft::logical.op(temp$X,
                                                 plot.metrics$Test[w.d],
                                                 plot.metrics$Value[w.d]))&(temp$FDR == paste0("FDR ", f3, "%")))] <- paste0("down, FDR = ", f3, "%")
          }
        }
        myColors[c(paste0("up, FDR = ", fdr_table$FDR, "%"),
                   paste0("down, FDR = ", fdr_table$FDR, "%"))] <- c(fdr_table$fdr.col.up, fdr_table$fdr.col.down)
      } else {
        if (mode == "curved") {
          # Placeholder.
          # This needs rewriting with a real, non arbitrary curved function (see the SAM method)
          temp$Threshold <- curved.threshold(temp$X)
          if (!symm) { temp$Threshold[which(temp$X) < 0] <- NA }
          temp$Colour[which((!is.na(temp$Threshold))&(temp$X > 0)&(temp$Y > temp$Threshold))] <- "up"
          temp$Colour[which((!is.na(temp$Threshold))&(temp$X < 0)&(temp$Y > temp$Threshold))] <- "down"
          myColors[c("up", "down")] <- c("red","green")
        } else {
          test.u <- proteoCraft::logical.op(temp$X,
                                      plot.metrics$Test[w.u],
                                      plot.metrics$Value[w.u])
          if (symm) { test.d <- proteoCraft::logical.op(temp$X,
                                                  plot.metrics$Test[w.d],
                                                  signif(plot.metrics$Value[w.d], 3)) }
          if ((mode == "standard")&&(!misFun(arbitrary.thresh))) {
            fdr <- as.numeric(gsub("% FDR$", "", plot.metrics$Levels[which(plot.metrics$Axis == "Y")]))/100
            fdr <- sort(fdr)
            up_Nms <- setNames(sapply(fdr, function(x) {
              w <- which(plot.metrics$Levels == paste0(x*100, "% FDR"))
              paste0("Ratio ", plot.metrics$Test[w.u], " ", signif(plot.metrics$Value[w.u], 3),
                     ",\n  -log10(p) ", plot.metrics$Test[w], " ", signif(plot.metrics$Value[w], 3))
            }), paste0(fdr*100, "% FDR"))
            test.f <- as.data.frame(sapply(rev(fdr), function(x) { #x <- rev(fdr)[1]
              w <- which(plot.metrics$Levels == paste0(x*100, "% FDR"))
              return(proteoCraft::logical.op(temp$Y,
                                       plot.metrics$Test[w],
                                       plot.metrics$Value[w]))
            }))
            colnames(test.f) <- paste0(rev(fdr)*100, "% FDR")
            for (f in colnames(test.f)) {
              temp$Colour[which(test.f[[f]])] <- "too small FC"
              temp$Colour[which((test.u)&(test.f[[f]]))] <- up_Nms[[f]]
            }
            myColors[up_Nms] <- plot.colours$up[match(names(up_Nms), rownames(plot.colours))]
            if (symm) {
              dwn_Nms <- setNames(sapply(fdr, function(x) {
                w <- which(plot.metrics$Levels == paste0(x*100, "% FDR"))
                paste0("Ratio ", plot.metrics$Test[w.d], " ", signif(plot.metrics$Value[w.d], 3),
                       ",\n  -log10(p) ", plot.metrics$Test[w], " ", signif(plot.metrics$Value[w], 3))
              }), paste0(fdr*100, "% FDR"))
              for (f in colnames(test.f)) { temp$Colour[which((test.d)&(test.f[[f]]))] <- dwn_Nms[[f]] }
              myColors[dwn_Nms] <- plot.colours$down[match(names(dwn_Nms), rownames(plot.colours))]
            }
          } else {
            w.s <- which(plot.metrics$Levels == "strict")
            w.l <- which(plot.metrics$Levels == "loose")
            up_Nms_l <- paste0("Ratio ", plot.metrics$Test[w.u], " ", signif(plot.metrics$Value[w.u], 3),
                               ",\n  -log10(", PorQ, ") ", plot.metrics$Test[w.s], " ", signif(plot.metrics$Value[w.s], 3))
            up_Nms_s <- paste0("Ratio ", plot.metrics$Test[w.u], " ", signif(plot.metrics$Value[w.u], 3),
                               ",\n  -log10(", PorQ, ") ", plot.metrics$Test[w.l], " ", signif(plot.metrics$Value[w.l], 3))
            up_Nms <- c(up_Nms_s, up_Nms_l)
            test.s <- proteoCraft::logical.op(temp$Y,
                                        plot.metrics$Test[w.s],
                                        plot.metrics$Value[w.s])
            test.l <- proteoCraft::logical.op(temp$Y,
                                        plot.metrics$Test[w.l],
                                        plot.metrics$Value[w.l])
            temp$Colour[which(test.l)] <- "too small FC"
            temp$Colour[which((test.u)&(test.l))] <- up_Nms_l
            temp$Colour[which((test.u)&(test.s))] <- up_Nms_s
            myColors[c(up_Nms_s, up_Nms_l)] <- plot.colours$up[match(c("strict", "loose"), rownames(plot.colours))]
            if (symm) {
              dwn_Nms_l <- paste0("Ratio ", plot.metrics$Test[w.d], " ", signif(plot.metrics$Value[w.d], 3),
                                  ",\n  -log10(", PorQ, ") ", plot.metrics$Test[w.l], " ", signif(plot.metrics$Value[w.l], 3))
              dwn_Nms_s <- paste0("Ratio ", plot.metrics$Test[w.d], " ", signif(plot.metrics$Value[w.d], 3),
                                  ",\n  -log10(", PorQ, ") ", plot.metrics$Test[w.s], " ", plot.metrics$Text.value[w.s])
              dwn_Nms <- c(dwn_Nms_s, dwn_Nms_l)
              temp$Colour[which((test.d)&(test.l))] <- dwn_Nms_l
              temp$Colour[which((test.d)&(test.s))] <- dwn_Nms_s
              myColors[c(dwn_Nms_s, dwn_Nms_l)] <- plot.colours$down[match(c("strict", "loose"), rownames(plot.colours))]
            }
          }
        }
      }
    }
    #
    temp$Colour <- factor(temp$Colour, levels = names(myColors))
    colScale <- ggplot2::scale_colour_manual(name = "colour", values = myColors)
    #if ((!is.numeric(Alpha))&&(!is.numeric(Size))) {
    #if (!is.numeric(Size)) {
    #  fillScale <- ggplot2::scale_fill_manual(name = "fill", values = myColors, guide = FALSE)
    #}
    if (misFun(reg.root)) {
      Prot[[paste0("Regulated - ", i)]] <- ""
      Prot[Wych[[i]], paste0("Regulated - ", i)] <- as.character(temp$Colour)
    }
    if (!is.null(proteins)) {
      proteins <- gsub("^CON_+", "", proteins)
      prots <- proteins[which(proteins %in% unique(unlist(strsplit(temp[[Proteins.col]], ";"))))]
      temp$Colour2 <- "not in list"
      temp$Labels2 <- ""
      myColors2 <- setNames(c("lightgrey", "purple", "brown"),
                            c("not in list", "target", "protein in list"))
      colScale2 <- ggplot2::scale_colour_manual(name = "colour", values = myColors2)
      #if (!is.numeric(Size)) { fillScale2 <- ggplot2::scale_fill_manual(name = "fill", values = myColors2, guide = FALSE) }
      if (length(prots)) {
        w <- proteoCraft::grsep2(prots, temp[[Proteins.col]])
        if ((!length(w))&&(prot_split)) {
          warning("No matches for proteins in list, no split plot will be created.")
          prot_split <- FALSE
        }
        if (!prot_split) { temp$Colour[w] <- "protein in list" } else {
          temp$Colour2[w] <- "protein in list"
          temp$Labels2[w] <- temp$Labels[w]
        }
      } else {
        if (prot_split) {
          warning("No matches for proteins in list, no split plot will be created.")
          prot_split <- FALSE
        }
      }
    }
    if ("Target" %in% colnames(experiments.map)) { # Will overwrite "protein in list" tag with "target" tag where relevant
      w <- proteoCraft::grsep2(target, temp[[Proteins.col]])
      temp$Colour[w] <- "target" # The "target" tag should be in all versions of the graph.
      if (prot_split) {
        temp$Colour2[w] <- "target"
        temp$Labels2[w] <- temp$Labels[w]
      }
    }
    noLbl <- 1:nrow(temp)
    if (sum(c("thresholds", "FDR", "regulated") %in% labels)) {
      noLbl <- which(temp$Colour %in% c("non significant", "too small FC"))
    }
    if ("proteins" %in% labels) {
      w <- which(temp$Colour %in% c("protein in list", "target"))
      noLbl <- noLbl[which(!noLbl %in% w)]
    }
    temp$Labels[noLbl] <- ""
    if ("FDR" %in% labels) {
      w1 <- which(temp$Colour %in% c(paste0("up, FDR = ", fdr.values, "%"), paste0("down, FDR = ", fdr.values, "%")))
    } else {
      if (toupper(parameters$Plot.threshold.values) == "CURVED") { w1 <- which(temp$Colour %in% c("up", "down")) } else {
        if (symm) { w1 <- which(temp$Colour %in% c(up_Nms, dwn_Nms)) } else { w1 <- which(temp$Colour %in% up_Nms) } }
    }
    w0 <- which(temp$Colour %in% c("non significant", "too small FC"))
    if ((prot_split)&&("Colour2" %in% colnames(temp))) {
      w2 <- which((temp$Colour2 == "protein in list"))
      w1 <- w1[which(!w1 %in% w2)] # Here it could be that we have overlap of w1 and w2 => we don't want that!
    } else { w2 <- which((temp$Colour == "protein in list")) }
    w3 <- which(temp$Colour == "target")
    temp <- rbind(temp[w0,], temp[w1,], temp[w2,], temp[w3,])
    #temp <- temp[which(proteoCraft::is.all.good(temp$Y, 2)),]
    ttl <- paste0(title, i2)
    #Ylab <- paste0("-log10(", PorQ, "value)")
    #if (grepl("adj\\.|adjusted", tolower(Y.root))) { Ylab <- paste0("-log10(adjusted ", PorQ, "value)") }
    Xlab <- gsub(" - $", "", X.root)
    Ylab <- gsub(" - $", "", Y.root)
    test <- c(is.numeric(Size), is.numeric(Alpha))
    temp$"P-value" <- 10^(-temp$Y)
    pL_lbs <- c(plotly_labels, "X", "Y")
    pL_lbs_nms <- c(names(plotly_labels), Xlab, Ylab)
    temp$plotly_labels <- apply(temp[, pL_lbs], 1, function(x) { #x <- temp[1, kol]
      paste0(pL_lbs_nms, ": ", gsub("\n", " ", x), collapse = "<br>")
    })
    aes <- data.frame(x = "X", y = "Y", text = "plotly_labels")
    non.aes <- data.frame(shape = 16)
    if (is.finite("Xlim")) { xlim <- Xlim }
    if (is.finite("Ylim")) { ylim <- Ylim }
    xlim <- c(min(c(xlim, plot.metrics$Value[which(plot.metrics$Axis == "X")]*1.1)),
              max(c(xlim, plot.metrics$Value[which(plot.metrics$Axis == "X")]*1.1)))
    ylim <- max(c(ylim, plot.metrics$Value[which(plot.metrics$Axis == "Y")]*1.1))
    pluses <- c("ggplot2::ggtitle(ttl, subtitle = paste0(\"Plotted: \", nrow(temp),\" data points.\"))",
                "ggplot2::xlab(Xlab)", "ggplot2::ylab(Ylab)", "ggplot2::theme_bw()",
                "ggplot2::xlim(xlim[1], xlim[2])", "ggplot2::ylim(0, ylim)" , "ggplot2::guides(___)")
    if (length(unique(temp$Colour)) > 1) {
      aes$colour <- "Colour"
      pluses <- c(pluses, "colScale")
    } else { non.aes$colour <- "\"black\"" }
    if (!test[1]) {
      aes$size <- "Size"
      #aes$fill <- "Colour"
      #non.aes$shape <- 21
      #aes$colour <- NULL
      #non.aes$colour <- NULL
      #fillScale <- ggplot2::scale_fill_manual(name = "colour", values = myColors)
      pluses <- c(pluses, "ggplot2::scale_size_identity(Size, guide = \"legend\")"#, "fillScale")
      )
    } else { non.aes$size <- "Size" }
    pluses <- c(pluses, "ggplot2::scale_y_continuous(expand = c(0, 0))")
    if (!test[2]) {
      aes$alpha <- "Alpha"
      if (Alpha.identity) { pluses <- c(pluses, "ggplot2::scale_alpha_identity(Alpha)") } else {
        pluses <- gsub("^ggplot2::guides\\(", "ggplot2::guides(alpha = guide_legend(title = Alpha), ", pluses) 
      }
    } else { non.aes$alpha <- "Alpha" }
    pluses <- gsub(", ___\\)", ")", pluses)
    pluses <- pluses[which(pluses != "ggplot2::guides(___)")]
    test1 <- sapply(aes[1,], function(x) { which(colnames(temp) == x) })
    test2 <- sapply(test1, function(x) { "numeric" %in% class(temp[[x]]) })
    test1 <- test1[which(test2)]
    test <- apply(temp, 1, function(x) {
      length(proteoCraft::is.all.good(as.numeric(unlist(x[test1]))))
    }) == length(test1)
    temp <- temp[which(test),]
    aes[grep(" ", aes)] <- paste0("\"", aes[grep(" ", aes)], "\"")
    aes <- paste(sapply(1:ncol(aes), function(x) {
      paste(colnames(aes)[x], aes[x], sep = " = ")
    }), collapse = ", ")
    non.aes <- paste(sapply(1:ncol(non.aes), function(x) {
      paste(colnames(non.aes)[x], non.aes[x], sep = " = ")
    }), collapse = ", ")
    non.aes <- gsub("^dummy = NA, ", "", non.aes)
    pluses <- paste(pluses, collapse = " + ")
    plot.txt <- paste0("plot <- ggplot2::ggplot(temp) + ggplot2::geom_point(ggplot2::aes(", aes, "), ", non.aes, ") + ", pluses)
    if (prot_split) {
      plot.txt_2 <- gsub("^plot <- ", "plot_prot <- ",
                         gsub("Colour", "Colour2",
                              gsub("colScale", "colScale2",
                                   #gsub("fillScale", "fillScale2",
                                   plot.txt#)
                              )))
    }
    #cat(plot.txt)
    suppressWarnings(eval(parse(text = plot.txt)))
    #
    plotMetr.lst[[i]] <- plot.metrics
    #suppressWarnings(eval(parse(text = plot.txt), envir = globalenv())) # for testing only
    if (prot_split) {
      suppressWarnings(eval(parse(text = plot.txt_2)))
    } # This is the plot for proteins of interest
    if (mode == "curved") {
      offset <- c(median(proteoCraft::is.all.good(temp$X)), 0)[X.normalized + 1]
      plot <- plot +
        ggplot2::stat_function(fun = function(x) {
          proteoCraft::curved.threshold(x) + offset
        }, xlim = xlim)
      if (prot_split) {
        plot_prot <- plot_prot +
          ggplot2::stat_function(fun = function(x) {
            proteoCraft::curved.threshold(x) + offset
          }, xlim = xlim)
      }
    } else {
      for (l in 1:nrow(plot.metrics)) {
        if (plot.metrics$Axis[l] == "X") {
          if ((plot.metrics$Levels[l] == "up")||(symm)) {
            offset <- xspan/50*c(-1,1)[which(c("down", "up") == plot.metrics$Levels[l])]
            plot <- plot +
              ggplot2::geom_vline(xintercept = plot.metrics$Value[l],
                                  color = plot.metrics$Colour[l]) +
              ggplot2::annotate("text", x = plot.metrics$Value[l]+offset,
                                y = ylim*0.99, label = plot.metrics$Name[l],
                                color = plot.metrics$Colour[l], cex = cex*1.2,
                                hjust = 1, angle = 90)
            if (prot_split) {
              plot_prot <- plot_prot +
                ggplot2::geom_vline(xintercept = plot.metrics$Value[l],
                                    color = plot.metrics$Colour[l]) + 
                ggplot2::annotate("text", x = plot.metrics$Value[l]+offset,
                                  y = ylim*0.99, label = plot.metrics$Name[l],
                                  color = plot.metrics$Colour[l], cex = cex*1.2,
                                  hjust = 1, angle = 90)
            }
          }
        }
        if (("FDR" %in% labels)&&(plot.metrics$Axis[l] == "Y")) {
          plot <- plot +
            ggplot2::geom_hline(yintercept = plot.metrics$Value[l],
                                color = plot.metrics$Colour[l])
          if (prot_split) {
            plot_prot <- plot_prot +
              ggplot2::geom_hline(yintercept = plot.metrics$Value[l],
                                  color = plot.metrics$Colour[l])
          }
        }
      }
    }
    if (("FDR" %in% labels)&&(nrow(fdr_table))) {
      plot <- plot +
        ggplot2::geom_hline(data = fdr_table,
                            ggplot2::aes(yintercept = -log10(Thresholds),
                                         colour = I(fdr.col.line)),
                            show.legend = FALSE) + 
        ggplot2::geom_text(data = fdr_table, x = xlim[1]+xspan*0.025,
                           ggplot2::aes(y = -log10(Thresholds)+ylim*0.01,
                                        label = paste0("FDR = ", FDR, "%"),
                                        colour = I(fdr.col.line)),
                           cex = cex * 1.2, hjust = 0,
                           show.legend = FALSE)
      if (prot_split) {
        plot_prot <- plot_prot +
          ggplot2::geom_hline(data = fdr_table,
                              ggplot2::aes(yintercept = -log10(Thresholds),
                                           colour = I(fdr.col.line)),
                              show.legend = FALSE) + 
          ggplot2::geom_text(data = fdr_table, x = xlim[1]+xspan*0.025,
                             ggplot2::aes(y = -log10(Thresholds)+ylim*0.01,
                                          label = paste0("FDR = ", FDR, "%"),
                                          colour = I(fdr.col.line)),
                             cex = cex * 1.2, hjust = 0,
                             show.legend = FALSE)
      }
    }
    if (!is.null(arbitrary.lines)) {
      # (This part could be rewritten without a loop using "data = ...")
      for (j in 1:nrow(arbitrary.lines)) {
        if (is.na(arbitrary.lines$xintercept[j])) {
          if (arbitrary.lines$slope[j] != 0) {
            # Ablines are not really interesting for now so I am ignoring labels.
            plot <- plot +
              ggplot2::geom_abline(intercept = arbitrary.lines$yintercept[j],
                                   slope = arbitrary.lines$slope[j],
                                   colour = arbitrary.lines$colour[j], linetype = 3)
            if (prot_split) {
              plot_prot <- plot_prot +
                ggplot2::geom_abline(intercept = arbitrary.lines$yintercept[j],
                                     slope = arbitrary.lines$slope[j],
                                     colour = arbitrary.lines$colour[j], linetype = 3)
            }
          } else {
            plot <- plot +
              ggplot2::geom_hline(yintercept = arbitrary.lines$yintercept[j],
                                  colour = arbitrary.lines$colour[j], linetype = 3) +
              ggplot2::annotate("text",
                                x = xlim[1]+xspan*0.025,
                                y = arbitrary.lines$yintercept[j] + ylim*0.01,
                                label = arbitrary.lines$label[j],
                                colour = arbitrary.lines$colour[j],
                                cex = cex * 1.2, hjust = 0)
            if (prot_split) {
              plot_prot <- plot_prot +
                ggplot2::geom_hline(yintercept = arbitrary.lines$yintercept[j],
                                    colour = arbitrary.lines$colour[j],
                                    linetype = 3) +
                ggplot2::annotate("text", x = xlim[1]+xspan*0.025,
                                  y = arbitrary.lines$yintercept[j] + ylim*0.01,
                                  label = arbitrary.lines$label[j],
                                  colour = arbitrary.lines$colour[j],
                                  cex = cex * 1.2, hjust = 0)
            }
          }
        } else {
          S <- sign(arbitrary.lines$xintercept)
          plot <- plot +
            ggplot2::geom_vline(xintercept = arbitrary.lines$xintercept[j],
                                colour = arbitrary.lines$colour[j], linetype = 3) +
            ggplot2::annotate("text", x = arbitrary.lines$xintercept[j]+xspan*0.025*S,
                              y = ylim*0.5, label = arbitrary.lines$label[j],
                              colour = arbitrary.lines$colour[j],
                              cex = cex * 1.2, hjust = 0)
          if (prot_split) {
            plot_prot <- plot_prot +
              ggplot2::geom_vline(xintercept = arbitrary.lines$xintercept[j],
                                  colour = arbitrary.lines$colour[j], linetype = 3) +
              ggplot2::annotate("text", x = arbitrary.lines$xintercept[j]+xspan*0.025*S,
                                y = ylim*0.5, label = arbitrary.lines$label[j],
                                colour = arbitrary.lines$colour[j],
                                cex = cex * 1.2, hjust = 0)
          }
        }
      }
    }
    if ((!misFun(arbitrary.thresh))&&(!is.null(arbitrary.thresh))) {
      if (labels != "FDR") {
        for (j in 1:nrow(arbitrary.thresh)) {
          if (!is.na(arbitrary.thresh$xintercept[j])) {
            warning("Significance thresholds are horizontal, \"xintercept\" will be ignored!")
            arbitrary.thresh$xintercept[j] <- NA
          }
          if (!arbitrary.thresh$slope[j] %in% c(NA, 0)) {
            warning("Significance thresholds are horizontal, \"slope\" will be ignored!")
            arbitrary.thresh$xintercept[j] <- 0
          }
          plot <- plot +
            ggplot2::geom_hline(yintercept = arbitrary.thresh$yintercept[j],
                                colour = arbitrary.thresh$colour[j], linetype = 3) +
            ggplot2::annotate("text", x = xlim[1]+xspan*0.025,
                              y = arbitrary.thresh$yintercept[j] + ylim*0.01,
                              label = arbitrary.thresh$label[j],
                              colour = arbitrary.thresh$colour[j],
                              cex = cex * 1.2, hjust = 0)
          if (prot_split) {
            plot_prot <- plot_prot +
              ggplot2::geom_hline(yintercept = arbitrary.thresh$yintercept[j],
                                  colour = arbitrary.thresh$colour[j], linetype = 3) +
              ggplot2::annotate("text", x = xlim[1]+xspan*0.025,
                                y = arbitrary.thresh$yintercept[j] + ylim*0.01,
                                label = arbitrary.thresh$label[j],
                                colour = arbitrary.thresh$colour[j],
                                cex = cex * 1.2, hjust = 0)
          }
        } 
      }
    }
    #
    Plots$Unlabelled[[ttl]] <- plotEval(plot)
    if (show.labels) {
      plot2 <- plot
      W <- which(temp$Labels != "")
      if (length(W)) {
        lab <- temp[W,]
        W2 <- which(!lab$Colour %in% c("protein in list", "target"))
        W3 <- which(lab$Colour %in% c("protein in list", "target"))
        if (length(W2) > MaxLabels) {
          lab2 <- lab[W2,]
          Ord <- order(abs(lab2[[MaxLabels_priority]]), decreasing = TRUE) # ("abs" only needed for X, but makes no difference for Y => more concise code)
          lab2 <- lab2[Ord[1:MaxLabels],]
          lab <- rbind(lab2, lab[W3,])
        }
        if ((!is.numeric(Alpha))&&(Alpha.labels)) {
          plot2 <- plot2 +
            ggrepel::geom_text_repel(data = lab,
                                     ggplot2::aes(label = Labels, x = X, y = Y,
                                                  colour = Colour, alpha = Alpha),
                                     force = 4, cex = cex, lineheight = lineheight,
                                     show.legend = FALSE)
        } else {
          plot2 <- plot2 +
            ggrepel::geom_text_repel(data = lab,
                                     ggplot2::aes(label = Labels, x = X, y = Y,
                                                  colour = Colour), alpha = 1,
                                     force = 4, cex = cex, lineheight = lineheight,
                                     show.legend = FALSE)
        }
        Plots$Labelled[[ttl]] <- plotEval(plot2)
      }
      proteoCraft::poplot(plot2)
      if (prot_split) {
        plot_prot2 <- plot_prot
        W2 <- which(temp$Labels2 != "")
        if (length(W2)) {
          lab2 <- temp[W2,]
          if ((!is.numeric(Alpha))&&(Alpha.labels)) {
            plot_prot2 <- plot_prot2 +
              ggrepel::geom_text_repel(data = lab2,
                                       ggplot2::aes(label = Labels2, x = X, y = Y,
                                                    colour = Colour2, alpha = Alpha),
                                       force = 4, cex = cex, lineheight = lineheight,
                                       show.legend = FALSE)
          } else {
            plot_prot2 <- plot_prot2 +
              ggrepel::geom_text_repel(data = lab2,
                                       ggplot2::aes(label = Labels2, x = X, y = Y,
                                                    colour = Colour2), alpha = 1,
                                       force = 4, cex = cex, lineheight = lineheight,
                                       show.legend = FALSE)
          }
          Plots$"Proteins in list - labelled"[[ttl]] <- plotEval(plot_prot2)
        }
      }
    } else { proteoCraft::poplot(plot) }
    if (prot_split) {
      Plots$"Proteins in list - unlabelled"[[ttl]] <- plotEval(plot_prot)
    }
    tr <- gsub("/|:|\\*|\\?|<|>|\\||/", "-", title.root)
    tt <- gsub("/|:|\\*|\\?|<|>|\\||/", "-", ttl)
    nm <- paste0(tr, tt)
    if (nchar(nm) > 98) { nm <- substr(nm, 1, 98) }
    if (nm %in% plotNms) {
      fixkount <- 0
      while (nm %in% plotNms) {
        fixkount <- fixkount + 1
        if (fixkount == 100) {
          stop("Really? Really?!?! You really have 100 similarly named conditions with very long names?!?! If you tried to break this function then you succeeded with (bad) style!")
        }
        nm <- paste0(substr(nm, 1, 93), "...", c("0", "")[(nchar(fixkount) > 1)+1],
                     fixkount)
      }
    }
    plotNms <- c(plotNms, nm)
    # if (save) {
    #   wLb <- length(which(temp$Labels != ""))
    #   prtSplt <- ((prot_split)&&(length(which(temp$Labels2 != ""))))
    #   for (sv in saveExt) {
    #     if (subfolderpertype) { sfpt <- paste0(subfolder, "/", sv) } else { sfpt <- subfolder }
    #     if (!dir.exists(sfpt)) { dir.create(sfpt, recursive = TRUE) }
    #     ggplot2::ggsave(paste0(sfpt, "/", nm, ".", sv), plot,
    #                     dpi = 300, width = 10, height = 10, units = "in")
    #     if (prot_split) {
    #       ggplot2::ggsave(paste0(sfpt, "/", nm, "_list.", sv), plot_prot,
    #                       dpi = 300, width = 10, height = 10, units = "in")
    #     }
    #     if (show.labels) {
    #       if (wLb) {
    #         ggplot2::ggsave(paste0(sfpt, "/", nm, "_tags.", sv), plot2,
    #                         dpi = 300, width = 10, height = 10, units = "in")
    #       }
    #       if (prtSplt) {
    #         ggplot2::ggsave(paste0(sfpt, "/", nm, "_list_tags.", sv), plot_prot2,
    #                         dpi = 300, width = 10, height = 10, units = "in")
    #       }
    #     }
    #   }
    # }
    if (plotly) {
      plot_ly <- plotly::ggplotly(plot, tooltip = "text")
      volcPlotly[[ttl]] <- list(Ttl = ttl,
                                Plot = plot_ly)
      if (prot_split) {
        plot_ly_list <- plotly::ggplotly(plot_prot, tooltip = "text")
        volcPlotly[[paste0(ttl, "_list")]] <- list(Ttl = paste0(ttl, "_list"),
                                                   Plot = plot_ly_list)
      }
      if (!plotly_local) {
        volcPlotly2[[ttl]] <- plotly::api_create(plot_ly,
                                                 paste0(plotly_subfolder, tr, tt),
                                                 "overwrite",
                                                 plotly_sharing)
        if (prot_split) {
          volcPlotly2[[paste0(ttl, "_list")]] <- plotly::api_create(plot_ly_list,
                                                                    paste0(plotly_subfolder, tr, tt, "_list"),
                                                                    "overwrite",
                                                                    plotly_sharing)
        }
      }
    }
  }
  if (save) {
    if (subfolderpertype) { sfpt <- paste0(subfolder, "/", saveExt[1]) } else { sfpt <- subfolder }
    if (!dir.exists(sfpt)) { dir.create(sfpt, recursive = TRUE) }
    plotsLst <- lst1 <- setNames(lapply(names(Plots$Unlabelled), function(nm) {
      list(Path = sfpt,
           Ttl = nm,
           Plot = Plots$Unlabelled[[nm]],
           Ext = saveExt[1])
    }), paste0(names(Plots$Unlabelled), "_noLabel"))
    if ((show.labels)&&(length(Plots$Labelled))) {
      lst2 <- setNames(lapply(names(Plots$Labelled), function(nm) {
        list(Path = sfpt,
             Ttl = paste0(nm, "_tags"),
             Plot = Plots$Labelled[[nm]],
             Ext = saveExt[1])
      }), paste0(names(Plots$Labelled), "_Label"))
      plotsLst <- c(lst1, lst2)
    }
    if ((prot_split)&&(length(Plots$"Proteins in list - unlabelled"))) {
      lst3 <- setNames(lapply(names(Plots$"Proteins in list - unlabelled"), function(nm) {
        list(Path = sfpt,
             Ttl = paste0(nm, "_list"),
             Plot = Plots$"Proteins in list - unlabelled"[[nm]],
             Ext = saveExt[1])
      }), paste0(names(Plots$"Proteins in list - unlabelled"), "_list_noLabel"))
      plotsLst <- c(plotsLst, lst3)
      if ((show.labels)&&(length(Plots$"Proteins in list - labelled"))) {
        lst4 <- setNames(lapply(names(Plots$"Proteins in list - labelled"), function(nm) {
          list(Path = sfpt,
               Ttl = paste0(nm, "_list_tags"),
               Plot = Plots$"Proteins in list - labelled"[[nm]],
               Ext = saveExt[1])
        }), paste0(names(Plots$"Proteins in list - labelled"), "_list_Label"))
        plotsLst <- c(plotsLst, lst3, lst4)
      }
    }
    l <- length(saveExt)
    if (l > 1) {
      lst0 <- plotsLst
      for (i in 2:l) {
        if (subfolderpertype) { sfpt <- paste0(subfolder, "/", saveExt[i]) } else { sfpt <- subfolder }
        if (!dir.exists(sfpt)) { dir.create(sfpt, recursive = TRUE) }
        lst <- lapply(lst0, function(x) {
          x$Path <- sfpt
          x$Ext <- saveExt[i]
          return(x)
        })
        plotsLst <- c(plotsLst, lst)
      }
    }
    f0 <- function(x) {
      require(ggrepel) # Even though the plot is evaluated, unless I do this the ggrepel labels are absent from the saved plot!
      ggplot2::ggsave(paste0(x$Path, "/", x$Ttl, ".", x$Ext), x$Plot,
                      dpi = 300, width = 10, height = 10, units = "in")
    }
    environment(f0) <- .GlobalEnv
    tst <- parLapply(cl, plotsLst, f0)
  }
  if (plotly&&plotly_local) {
    if (subfolderpertype) { sfpt <- paste0(subfolder, "/html") } else { sfpt <- subfolder }
    if (!dir.exists(sfpt)) { dir.create(sfpt, recursive = TRUE) }
    setwd(sfpt)
    clusterExport(cl, "sfpt", envir = environment())
    f0 <- function(x) {
      htmlwidgets::saveWidget(x$Plot, paste0(sfpt, "/", x$Ttl, ".html"),
                              selfcontained = TRUE)
    }
    environment(f0) <- .GlobalEnv
    tst <- parLapply(cl, volcPlotly, f0)
    setwd(origWD)
  }
  thrsh <- list(Absolute = plotMetr.lst)
  if ("FDR" %in% labels) { thrsh$FDR <- FDR_table }
  RES <- list(Thresholds = thrsh)
  if (return) { RES$Protein_groups_file <- Prot }
  if (return.plot) { RES$Plots <- Plots }
  if ((plotly)&&(!plotly_local)) { RES$"Plotly plots" <- volcPlotly2 }
  #
  setwd(origWD)
  if (cleanUp) { parallel::stopCluster(cl) }
  return(RES)
}
