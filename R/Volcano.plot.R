#' Volcano.plot
#'
#' @description
#' A function to create volcano plots.
#' 
#' @param Prot The protein/protein groups table. It will be added columns saying whether a protein group is regulated or not for each test performed (i.e. each "aggregate")
#' @param mode One of "standard" or "custom". Ignored if SAM = TRUE. 
#' @param X.root Root of the names for X-axis log-scale values. If missing, defaults are provided in mode = "custom" in the parameters file.
#' @param Y.root Root of the names for Y-axis values (usually -log10(Pvalues)). If missing, default are provided in mode = "custom" in the parameters file.
#' @param X.normalized Logical. Is X normalized? Assumed to be TRUE by default.
#' @param experiments.map The data.frame defining the structure of the dataset.
#' @param aggregate.map The aggregate map.
#' @param aggregate.list The named list of aggregates.
#' @param aggregate.name The name of the aggregate to use.
#' @param parameters The experiment's parameters file.
#' @param save Logical or character. Should the plot(s) be saved? Default = FALSE. Set it to a vector of ggplot2-compatible file extensions to save to the corresponding format.
#' @param labels How to decide which point labels will be printed. Values are:\cr
#'  - "proteins", in which case parameters MaxLabels and MaxLabels_priority will also be used;\cr
#'  - "thresholds", in which case proteins above thresholds will be labelled;\cr
#'  - "FDR", in which case only points that pass Benjamini-Hochberg FDR (or are in the list of proteins of interest) will be displayed.\cr
#'  - The default, "both" (default), which means "proteins + thresholds".
#' @param MaxLabels Maximum number of labels to print. Does not apply to matches to the list provided by "proteins". Default = 50 
#' @param MaxLabels_priority If not printing all labels, will we favour more extreme ratios (default = "X") or more significant values ("Y")?
#' @param proteins Optional, the character list of protein IDs to be highlighted.
#' @param proteins_split Logical. Default = FALSE. If TRUE, and if there are protein matching the provided list (see "proteins" argument), then these will be highlighted in a separate graph.
#' @param IDs.col Default = "Protein IDs".
#' @param Proteins.col Default = "Protein IDs"; used for identifying matches to the "proteins" argument. IF missing, IDs.col will be used.
#' @param Ref.Ratio.values Used if the default parameter values for Ratio thresholds are to be replaced with dynamic, reference-based ones. Also requires a few arguments below.
#' @param Ref.Ratio.method Defines how vertical ratios are computed.\cr
#'  - "obs2": current default, will sort observations within sample groups and identify a threshold such that a "ratios.FDR" proportion of most extreme values are selected.\cr
#'  - "obs1": does almost the same as the above, but only for control-to-control or control-to-average-control ratios within the ratios group.\cr
#'  - "SD": same as the previous, but will calculate the threshold using qnorm(ratios.FDR, median(Ref.Ratio.values), sd(Ref.Ratio.values)), i.e. assuming normality.\cr
#' @param ratios.FDR The proportion of reference ratios to be accepted, if Ref.Ratio.values is provided. Default = 0.01.
#' @param FDR.thresh The significance thresholds calculated for each acceptable FDR value.
#' @param FDR.root The root of the columns with Significance for different FDRs. Default = "Significant-FDR="
#' @param arbitrary.lines Default = NULL. Data.frame of arbitrary lines to add to the graph. Must have 5 columns: "slope", "yintercept", "xintercept" (if vertical line, else NA), "colour", "label".
#' @param arbitrary.thresh Default = NULL, ignored if mode = FDR. Similar to the above, but will define only horizontal lines. The lowest will determine the threshold below which no labels are displayed (except when using the "proteins" argument).
#' @param return Logical. Should we return a column with significance levels? Default = FALSE.
#' @param return.plot Logical. Should we return the plot itself? Default = FALSE.
#' @param show.labels Currently unused (deprecated).
#' @param title It is possible to set a specific title. Failing that a default one will be provided based on the value of the mode argument.
#' @param title.root Optional argument to add a character string at the beginning of the file names when saving.
#' @param subfolder Name of the sub-folder where graphs are to be saved. Default = "Reg. analysis/t-tests"
#' @param subfolderpertype Logical. If TRUE (default), will save each type of graph in a dedicated sub-folder.
#' @param Symmetrical Are we also interested in the left half (down-regulated) or only the right one (up-regulated)? May be used to override the symmetry provided by contrasts. If it and contrasts are missing, and "parameters" is provided, the function will attempt to infer it from parameter "Type".
#' @param Alpha Either a single numeric value, or a length 1 character: the name of a numeric column in Prot mapping each dot's alpha level to a column. Default = 1 (no alpha variations)
#' @param Alpha.identity Logical. Do we need to fit the Alphas into a 0:1 scale? Default = FALSE
#' @param Alpha.labels Logical. Should labels also be affected by alpha? Default = FALSE
#' @param Alpha.max For when mapping alpha to a column. The highest allowed alpha. Default = 1
#' @param Alpha.min For when mapping alpha to a column. The weakest allowed alpha. Default = 0.01
#' @param Size Either a single numeric value, or a length 1 character: the name of a numeric column in Prot. Default = 1 (same size for all)
#' @param Size.max For when mapping size to a column. The largest allowed dot size. Default = 2
#' @param Size.min For when mapping size to a column. The smallest allowed dot size. Default = 0.01
#' @param cex Label text size. Default = 2
#' @param lineheight Label interline size. Default = 1
#' @param plotly Logical. Should we generate a plotly? FALSE by default; if TRUE, the link will be automatically returned.
#' @param plotly_local Deprecated (now always TRUE)
#' @param plotly_username Deprecated
#' @param plotly_API_key Deprecated
#' @param plotly_subfolder Default = "". Name of the plotly server sub-folder where plotly graphs are to be saved.
#' @param plotly_sharing Deprecated
#' @param plotly_labels Labels to be displayed in the plotly tooltip. A character corresponding to valid names of Prot columns. Naming this vector provides a way to override the original column names: what will be displayed in the tooltip will be the name in this vector, not the column name itself.
#' @param Xlim Set fixed limits for the X axis.
#' @param Ylim Set fixed limits for the Y axis.
#' @param Contaminants Logical. If FALSE (default), contaminants will be filtered out prior to drawing the plot(s).
#' @param N.clust A limit on the number of vCPUs to use. If left as NULL (default), uses the number of available clusters - 1, to a minimum of 1.
#' @param N.reserved Default = 1. Number of reserved vCPUs the function is not to use. Note that for obvious reasons it will always use at least one.
#' @param cl Already have a cluster handy? Why spend time making a new one, which on top of that may invalidate the old one. Just pass it along!
#' @param reg.root Classically, this function would output a "Regulated - ..." column per sample group, based on applying the thresholds (either provided directly - significance is read from the table using "FDR.root", and FC threshold are either applied directly using "arbitrary.thresh" or calculated based on "Ref.Ratio.values"). Providing this input allows the function to bypass this step and apply colors based on a pre-determined decision (used for the F-test currently). Automatically set to "Regulated - " if SAM is TRUE.
#' @param SAM Logical. If TRUE, then a curved threshold is applied based on the SAM analysis. FALSE by default.
#' @param curved_Thresh Curved threshold parameters, only used if mode = "curved".\cr A named list, with one item for each (same names as the ones appended to X.root and Y.root in Prot column names).\cr Each item should be a named list of S0 (single value), Si (single value) and d (one value per FDR level) numerics.
#' @param saveData Logical. If TRUE, a file containing the processed long-format data used to create each plot will be saved locally. Default = FALSE.
#' @param contrasts If provided, supersedes experiments.map for the purpose of detecting samples (through the "Contrast" column), and if containing the "Up-only" column also defines whether plots are two-sided or not.
#' @param X.root_ind Optional root of the column names of individual expression values, which can then be displayed in the tooltip.
#' @param X.root_ind_base Log base of individual expression values. Default = 10L
#'
#' @details
#' This monster of a function can draw volcano plots, but can also be used for the decision on regulation.
#' 
#' @returns
#' A list with at least one element ("Thresholds"). If return == TRUE, the Protein_groups_file - which may contain additional Regulated columns. If return.plot == TRUE, the evaluated ggplots. If plotly == TRUE, the plotly plots.
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
                         MaxLabels = 50L,
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
                         IDs.col = "Protein IDs",
                         Proteins.col = "Protein IDs",
                         return = FALSE,
                         return.plot = FALSE,
                         title = "",
                         title.root = "",
                         subfolder = "Reg. analysis/t-tests",
                         subfolderpertype = TRUE,
                         Symmetrical,
                         Alpha = 1L,
                         Alpha.identity = FALSE,
                         Alpha.labels = FALSE,
                         Alpha.max = 1L,
                         Alpha.min = 0.01,
                         Size = 1L,
                         Size.max = 3L,
                         Size.min = 0.01,
                         cex = 2L,
                         lineheight = 1L,
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
                         N.reserved = 1L,
                         cl,
                         reg.root, # DO NOT give this one a default non-null value to avoid uncontrolled behaviour when rerunning the function!
                         SAM = FALSE,
                         curved_Thresh,
                         saveData = FALSE,
                         contrasts,
                         X.root_ind,
                         X.root_ind_base = 10L,
                         show.labels = TRUE) {
  TESTING <- FALSE
  #
  #DefArg(Volcano.plot); TESTING <- TRUE; cl <- parClust
  #invisible(lapply(names(volcPlot_args2), \(x) { assign(x, volcPlot_args2[[x]], envir = .GlobalEnv); return() }))
  origWD <- getwd()
  misFun <- if (TESTING) {
    # Note:
    # This is not a perfect alternative to missing but will work in most cases, unless x matches a function imported by a package 
    \(x) { return(!exists(deparse(substitute(x)))) }
  } else { missing }
  # Check logicals
  if ((!is.logical(proteins_split)) || (length(proteins_split) != 1L) || is.na(proteins_split)) {
    proteins_split <- FALSE
  }
  if ((!is.logical(return)) || (length(return) != 1L) || is.na(return)) {
    return <- FALSE
  }
  if ((!is.logical(return.plot)) || (length(return.plot) != 1L) || is.na(return.plot)) {
    return.plot <- FALSE
  }
  if ((!is.logical(subfolderpertype)) || (length(subfolderpertype) != 1L) || is.na(subfolderpertype)) {
    subfolderpertype <- TRUE
  }
  if ((!is.logical(Alpha.identity)) || (length(Alpha.identity) != 1L) || is.na(Alpha.identity)) {
    Alpha.identity <- FALSE
  }
  if ((!is.logical(Alpha.labels)) || (length(Alpha.labels) != 1L) || is.na(Alpha.labels)) {
    Alpha.labels <- FALSE
  }
  if ((!is.logical(plotly)) || (length(plotly) != 1L) || is.na(plotly)) {
    plotly <- TRUE
  }
  plotly_local <- TRUE
  if ((!is.logical(Contaminants)) || (length(Contaminants) != 1L) || is.na(Contaminants)) {
    Contaminants <- FALSE
  }
  if ((!is.logical(SAM)) || (length(SAM) != 1L) || is.na(SAM)) {
    SAM <- FALSE
  }
  if ((!is.logical(saveData)) || (length(saveData) != 1L) || is.na(saveData)) {
    saveData <- FALSE
  }
  #
  # (Don't use as default the value in parameters$Plot.metrics: they are deprecated)
  if (misFun(X.root)) { stop("Argument \"X.root\" is missing, investigate!") }
  if (misFun(Y.root)) { stop("Argument \"Y.root\" is missing, investigate!") }
  #
  # Create cluster
  stopCl <- FALSE
  if (is.null(cl) || (!inherits(cl, "cluster"))) {
    dc <- parallel::detectCores()
    if (misFun(N.reserved)) { N.reserved <- 1L }
    nMax <- max(c(dc - N.reserved, 1L))
    if (misFun(N.clust)) { N.clust <- nMax } else {
      if (N.clust > nMax) {
        warning("More cores specified than allowed, I will ignore the specified number! You should always leave at least one free for other processes, see the \"N.reserved\" argument.")
        N.clust <- nMax
      }
    }
    cl <- parallel::makeCluster(N.clust, type = "SOCK")
    stopCl <- TRUE
  }
  N.clust <- length(cl)
  #
  RES <- NA
  #
  if (labels == "both") { labels <- c("proteins", "thresholds") }
  stopifnot(nrow(Prot) > 0L,
            ncol(Prot) > 0L,
            is.character(subfolder))
  #
  subfolder <- normalizePath(gsub("^/|/$", "", gsub("\\\\", "/", subfolder)),
                             winslash = "/", mustWork = FALSE)
  if (subfolder %in% c("", paste0(LETTERS, ":"))) { subfolder <- getwd() }
  if (!dir.exists(subfolder)) {
    tst <- try(dir.create(subfolder), silent = TRUE)
    if (inherits(tst, "try-error")) { subfolder <- getwd() } 
  }
  #
  if ((length(save) > 1L) || ((!is.logical(save)) || (save == TRUE))) {
    if ((length(save) == 1L) && is.logical(save) && (save == TRUE)) {
      saveExt <- "jpeg" # default
    }
    saveExt <- unique(gsub("^jpg$", "jpeg", gsub("^\\.", "", tolower(save))))
    save <- TRUE
  } else { save <- FALSE}
  #
  if ((!is.integer(MaxLabels)) || (MaxLabels < 0L)) { MaxLabels <- 100L }
  if (!MaxLabels_priority %in% c("X", "Y")) { MaxLabels_priority <- "X" }
  #
  upColRg <- c("blue", "green")
  downColRg <- c("red", "orange")
  regProvided <- ((!misFun(reg.root)) && (!is.null(reg.root)))
  if (SAM) {
    mode <- "curved"
    reg.root <- "Regulated - "
    regProvided <- TRUE
  }
  if (regProvided) {
    labels <- c("regulated",
                labels[which(labels == "proteins")],
                "FDR")
  }
  useFDRtbl <- FALSE
  if (("FDR" %in% labels) || (mode == "curved")) {
    if (regProvided) {
      xKols <- grep(topattern(X.root), colnames(Prot), value = TRUE)
      yKols <- grep(topattern(Y.root), colnames(Prot), value = TRUE)
      rgKols <- grep(topattern(reg.root), colnames(Prot), value = TRUE)
      xTst <- gsub(topattern(X.root), "", xKols)
      yTst <- gsub(topattern(Y.root), "", yKols)
      rgTst <- gsub(topattern(reg.root), "", rgKols)
      tstDF <- data.frame(Group = unique(c(xTst, yTst, rgTst)))
      tstDF$X <- xKols[match(tstDF$Group, xTst)]
      tstDF$Y <- yKols[match(tstDF$Group, yTst)]
      tstDF$Reg <- rgKols[match(tstDF$Group, rgTst)]
      tstDF <- tstDF[which(apply(tstDF[, c("X", "Y", "Reg")], 1L, \(x) { sum(is.na(x)) }) == 0L),]
      stopifnot(nrow(tstDF) > 0L)
      FDR_table <- lapply(1L:nrow(tstDF), \(x) { #x <- 1L
        fdr_table <- Prot[grep("^((up)|(down)), FDR = ", Prot[[tstDF$Reg[[x]]]]),
                          c(tstDF$X[[x]], tstDF$Y[[x]], tstDF$Reg[[x]])]
        if (!nrow(fdr_table)) { return() }
        fdr_table$FDRs <- as.numeric(gsub("^((up)|(down)), FDR = |%$", "", fdr_table[[tstDF$Reg[[x]]]]))
        fdr_table <- aggregate(fdr_table[[tstDF$Y[[x]]]], list(fdr_table$FDRs), min)
        colnames(fdr_table) <- c("FDR", "Thresholds")
        fdr_table$Thresholds <- 10L^-fdr_table$Thresholds # Because of how they are stored
        fdr_table$Sample <- tstDF$Group[[x]]
        return(fdr_table)
      })
      FDR_table <- do.call(rbind, FDR_table)
    } else {
      if (misFun(FDR.thresh)) { stop("\"FDR.thresh\" must be provided if labels = \"FDR\"!") }
      g <- as.numeric(gsub("^Threshold-FDR=|%( - .+)?$", "", names(FDR.thresh)))
      tst <- unique(grepl(" - ", names(FDR.thresh)))
      FDR_table <- data.frame(FDR = g,
                              Thresholds = FDR.thresh)
      if ((length(tst) == 1L) && tst) {
        FDR_table$Sample <-  gsub("^Threshold-FDR=[1-9][0-9]*\\.*[0-9]*% - ", "", names(FDR.thresh))
      }
    }
    useFDRtbl <- (!is.null(FDR_table)) && is.data.frame(FDR_table) && nrow(FDR_table)
  }
  if (useFDRtbl) {
    rownames(FDR_table) <- paste0("Threshold-FDR=", FDR_table$FDR, "% - ", FDR_table$Sample)
    f <- grep(topattern(FDR.root), colnames(Prot), value = TRUE)
    FDR_table <- FDR_table[order(FDR_table$FDR, decreasing = TRUE),]
    FDR.values <- as.numeric(unique(gsub("%( - .+)?$", "", gsub(topattern(FDR.root), "", f))))
    fdr.col.up <- grDevices::colorRampPalette(upColRg)(length(FDR.values))
    fdr.col.down <- grDevices::colorRampPalette(downColRg)(length(FDR.values))
    fdr.col.line <- rainbow(n = length(FDR.values), start = 2/6, end = 1/6, v = 0.75)
    m <- match(FDR_table$FDR, FDR.values)
    FDR_table$fdr.col.up <- fdr.col.up[m]
    FDR_table$fdr.col.down <- fdr.col.down[m]
    FDR_table$fdr.col.line <- fdr.col.line[m]
  }
  if (!Contaminants) {
    kontkol <- intersect(c("Contaminant", "Potential contaminant"),
                         colnames(Prot))
    if (length(kontkol) != 1L) {
      warning("I could not identify the contaminants column, did you already filter out contaminants?")
      Contaminants <- TRUE
    }
  }
  stopifnot(!misFun(aggregate.name),
            mode %in% c("standard", "custom", "curved"))
  if (plotly) {
    volcPlotly <- list()
  }
  if (nchar(title.root)) {
    if (!substr(title.root, nchar(title.root), nchar(title.root)) %in% c(" ", ".", "_", ".")) {
      title.root <- paste0(title.root, " ")
    }
  }
  if (title == "") {
    if (mode %in% c("standard", "custom", "curved")) { title <- "Volcano plot " }
  }
  plotNms <- c()
  dfltPM <- TRUE
  #
  if ((!misFun(parameters)) &&
      nchar(parameters$Plot.metrics) &&
      nchar(parameters$Plot.threshold.metrics) &&
      nchar(parameters$Plot.threshold.values) &&
      nchar(parameters$Plot.threshold.tests) &&
      nchar(parameters$Plot.threshold.colours) &&
      nchar(parameters$Plot.areas.colours)) {
    dfltPM <- try({
      a1 <- set_colnames(as.data.frame(Isapply(strsplit(unlist(strsplit(parameters$Plot.metrics, "; *")), ": *"), unlist)),
                         c("Axis", "Name"))
      Plot.metrics <- as.data.frame(strsplit(unlist(strsplit(parameters$Plot.threshold.metrics, "; *")), ": *"))
      Plot.metrics <- as.data.frame(t(Plot.metrics)) 
      rownames(Plot.metrics) <- NULL
      colnames(Plot.metrics) <- c("Levels", "Axis")
      Plot.metrics$Root <- c(X.root, Y.root)[match(Plot.metrics$Axis, a1$Axis)]
      a2 <- set_colnames(as.data.frame(t(sapply(strsplit(unlist(strsplit(parameters$Plot.threshold.values, "; *")), ": *"), unlist))),
                         c("Direction", "Text.value"))
      Plot.metrics$Text.value <- a2$Text.value[match(Plot.metrics$Levels, a2$Direction)]
      Plot.metrics$Value <- vapply(Plot.metrics$Text.value, \(x) { eval(parse(text = x)) }, 1)
      a3 <- set_colnames(as.data.frame(t(sapply(strsplit(unlist(strsplit(parameters$Plot.threshold.tests, "; *")), ": *"), unlist))),
                         c("Direction", "Test"))
      Plot.metrics$Test <- a3$Test[match(Plot.metrics$Levels, a3$Direction)]
      a4 <- set_colnames(as.data.frame(t(sapply(strsplit(unlist(strsplit(parameters$Plot.threshold.colours, "; *")), ": *"), unlist))),
                         c("Direction", "Colour"))
      Plot.metrics$Colour <- a4$Colour[match(Plot.metrics$Levels, a4$Direction)]
      a5 <- unlist(strsplit(parameters$Plot.areas.colours, "; *"))
      a5 <- as.data.frame(strsplit(a5, "[_:] *"))
      Plot.colours <- as.data.frame(sapply(Plot.metrics$Levels[which(Plot.metrics$Axis == "X")], \(x) {
        sapply(Plot.metrics$Levels[which(Plot.metrics$Axis == "Y")], \(y) {
          a5[3L, which((a5[1L,] == x) & (a5[2L,] == y))]
        })
      }))
    }, silent = TRUE)
    dfltPM <- inherits(dfltPM, "try-error")
  }
  if (dfltPM) { # Provide defaults in case the parameters above are not provided
    # Vertical (ratio) thresholds
    Plot.metrics <- data.frame(Levels = c("down", "up"),
                               Axis = "X",
                               Name = X.root,
                               Test = c("<=", ">="),
                               Colour = c("red", "red"),
                               Text.value = c("log2(1/2)", "log2(2)"))
    Plot.metrics$Value <- vapply(Plot.metrics$Text.value, \(x) { eval(parse(text = x)) }, 1)
    # NB: These defaults may be overwritten by the "Ref.Ratio.values" argument for each plot.
    # Horizontal (P-value) thresholds
    tmpmetr <- if ((mode == "standard") && (!misFun(arbitrary.thresh))) {
      data.frame(Levels = arbitrary.thresh$label,
                 Axis = "Y",
                 Name = "Significance",
                 Value = arbitrary.thresh$yintercept,
                 Test = ">=",
                 Colour = arbitrary.thresh$colour)
    } else {
      data.frame(Levels = c("strict", "loose"),
                 Axis = "Y",
                 Name = setNames(Y.root, NULL), # (names to NULL to suppress a useless warning)
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
  #
  #X <- gsub("\\.$", "", X.root)
  #Y <- gsub("\\.$", "", Y.root)
  B <- aggregate.name
  A <- aggregate.list[[B]]
  if (nchar(B) == 3L) { B <- aggregate.map$Characteristics[[which(aggregate.map$Aggregate.Name == B)]] }
  #
  useContrasts <- !misFun(contrasts) && is.data.frame(contrasts) && !is.null(contrasts)
  #
  A <- if (useContrasts) {
    contrasts$Contrast
  } else {
    intersect(A, experiments.map[[B]])
  }
  lA <- length(A)
  #
  if (!misFun(Symmetrical)) {
    if (length(Symmetrical) == 1L) { Symmetrical <- rep(Symmetrical, lA) }
    if (length(Symmetrical) != lA) {
      warning("The length of the \"Symmetrical\" argument is incorrect, check it!")
      Symmetrical <- rep(Symmetrical, lA)
    }
  } else {
    if (useContrasts && ("Up-only" %in% colnames(contrasts))) {
      Symmetrical <- !contrasts$"Up-only"
    } else {
      # Provide defaults
      Symmetrical <- if ((!misFun(parameters)) && ("Two.sided" %in% colnames(parameters)) && is.logical(parameters$Two.sided)) {
        parameters$Two.sided
      } else {
        if ((!misFun(parameters)) && ("Type" %in% colnames(parameters))) {
          !gsub(" |_|-|\\.", "", toupper(parameters$Type)) %in% c("IP", "IMMUNOPRECIPITATION", "BIOID", "PULLDOWN")
        } else { TRUE }
      }
      Symmetrical <- rep(Symmetrical, lA)
    }
  }
  names(Symmetrical) <- A
  #
  xKols <- setNames(paste0(X.root, A), A)
  yKols <- setNames(paste0(Y.root, A), A)
  PorQ <- rev(unlist(strsplit(toupper(gsub("-value.*", "", Y.root, ignore.case = TRUE)), "")))[1L]
  if (!PorQ %in% c("P", "Q")) { PorQ <- "P" }
  if ((length(labels) == 1L) && (labels == "proteins") && proteins_split) {
    warning("Argument \"proteins_split\" will be ignored since argument \"labels\" is set to \"proteins\".")
    proteins_split <- FALSE
  } # No splitting if we are already only labeling proteins in list.
  if (misFun(plotly_labels)) { # Create those even if plotly off!
    plotly_labels <- setNames(c("Labels", IDs.col, "Genes"),
                              c("Name(s)", "Accession(s)", "Gene(s)")) # (column "Labels" gets created later, before we use those columns: should not break!)
    plotly_labels <- intersect(plotly_labels, colnames(Prot))
  } else {
    wlabkol <- which(plotly_labels %in% colnames(Prot))
    if (plotly && (length(wlabkol) < length(plotly_labels))) {
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
  # List of proteins of interest
  useProtList <- (!misFun(proteins)) & (!is.null(proteins)) & length(proteins)
  if (misFun(proteins_split)) { proteins_split <- FALSE}
  if (useProtList) {
    if (misFun(Proteins.col) || (!Proteins.col %in% colnames(Prot))) { Proteins.col <- IDs.col }
    proteins <- gsub("^CON_+", "", proteins)
    Prot[[Proteins.col]] <- gsub(";CON_+", ";", gsub("^CON_+", "", Prot[[Proteins.col]]))
    proteins <- intersect(proteins, unique(unlist(strsplit(Prot[[Proteins.col]], ";"))))
    useProtList <- length(proteins) > 0L
    if (useProtList) {
      Prot$Found_in_List <- FALSE
      wLst <- grsep2(proteins, Prot[[Proteins.col]])
      if (length(wLst)) {
        Prot$Found_in_List[wLst] <- TRUE
      } else {
        useProtList <- FALSE
        msg <- "No matches for proteins in list!"
        if (proteins_split) {
          msg <- paste0(msg, " (no split plot will be created)")
          proteins_split <- FALSE
        }
      }
    } else { proteins_split <- FALSE }
  } else { proteins_split <- FALSE }
  #
  # Default colors
  myColors <- setNames(c("lightgrey", "lightgrey", "purple", "brown"),
                       c("non significant", "too small FC", "target", "protein in list"))
  if (useFDRtbl) {
    myColors[c(paste0("up, FDR = ", FDR_table$FDR, "%"),
               paste0("down, FDR = ", FDR_table$FDR, "%"))] <- c(FDR_table$fdr.col.up, FDR_table$fdr.col.down)
  }
  if (proteins_split) {
    myColors2 <- setNames(c("lightgrey", "purple", "brown"),
                          c("not in list", "target", "protein in list"))
    colScale2 <- ggplot2::scale_colour_manual(name = "colour", values = myColors2)
  }
  #
  # Filter samples for available valid data
  tstTbl <- data.frame(Sample = A,
                       logFC = xKols %in% colnames(Prot),
                       PVal = yKols %in% colnames(Prot))
  tstTbl$All_OK <- tstTbl$logFC & tstTbl$PVal
  if (regProvided) {
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
      up_Nms$Val <- if (sum(c("FDR", "regulated") %in% labels)) {
        as.numeric(gsub(".* FDR = |%$", "", up_Nms$Up))
      } else {
        as.numeric(gsub("^[^0-9]+|[^0-9]+$", "", up_Nms$Up))
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
      dwn_Nms$Val <- if (sum(c("FDR", "regulated") %in% labels)) {
        as.numeric(gsub(".* FDR = |%$", "", dwn_Nms$Down))
      } else {
        as.numeric(gsub("^[^0-9]+|[^0-9]+$", "", dwn_Nms$Down))
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
    Wych <- setNames(lapply(wOK, \(x) {
      wych <- is.finite(as.numeric(Prot[[xKols[x]]])) & is.finite(as.numeric(Prot[[yKols[x]]]))
      if (!Contaminants) { wych <- wych & (Prot[[kontkol]] != "+") }
      return(which(wych))
    }), A[wOK])
    tstTbl$All_OK[wOK] <- lengths(Wych) > 0L
  }
  wOK <- which(tstTbl$All_OK)
  if (!length(wOK)) { stop("Not a single sample group with valid statistics detected!\nCheck inputs or statistical analysis") }
  A <- A[wOK]
  xKols <- xKols[A]
  yKols <- yKols[A]
  if (regProvided) { regKols <- regKols[A] }
  Symmetrical <- Symmetrical[A]
  Wych <- Wych[A]
  #
  # Define global x-limits
  tmp <- Prot[, xKols, drop = FALSE]
  xlim <- as.numeric(unlist(tmp))
  xlim <- xlim[which(is.finite(xlim))]
  xlim <- if (X.normalized) {
    max(abs(xlim))*c(-1, 1)
  } else {
    c(min(xlim), max(xlim))
  }
  xspan <- xlim[2L]-xlim[1L]
  xlim <- xlim + xspan*c(-0.05, 0.05)
  #
  # Define global y-limits
  tmp <- Prot[, yKols, drop = FALSE]
  ylim <- as.numeric(unlist(tmp))
  ylim <- ylim[which(is.finite(ylim))]
  ylim <- max(abs(ylim))*1.1
  #
  # Labels
  Prot$Labels <- Prot[[parameters$Plot.labels]]
  weech <- which(Prot$Labels != "")
  if (length(weech)) {
    tmp <- strsplit(Prot$Labels[weech], "  ?")
    colchar <- 25L
    #parallel::clusterExport(cl, list("tmp", "colchar"), envir = environment()) # Waste of time
    F0 <- .bind_worker(.labelEdit_worker,
                       list(tmp = tmp,
                            colchar = colchar))
    Prot$Labels[weech] <- parallel::parSapply(cl,
                                              tmp,
                                              F0,
                                              colchar = colchar)
  }
  weech <- which(Prot$Labels == "")
  id.col <- unique(c(IDs.col, "Common Names", "Common Name (short)", "Names", "Name", "Genes", "Gene", "Code"))
  id.col <- c(id.col, tolower(id.col), toupper(id.col))
  id.col <- intersect(id.col, colnames(Prot))
  if (length(weech)) {
    kount <- 0L
    while (length(weech) && (kount < length(id.col))) {
      kount <- kount+1L
      Prot$Labels[weech] <- vapply(Prot[weech, id.col[kount]], \(x) {
        if (x == "") { res <- "" } else {
          x <- unlist(strsplit(x, ";"))
          res <- if (length(x) == 1L) { x } else { paste0(x[1L], ";...") }
        }
        return(res)
      }, "")
      weech <- which(Prot$Labels == "")
    }
    if (length(weech)) { Prot$Labels[weech] <- paste0("Default_label_", c(1L:length(weech))) }
  }
  #
  plotMetr.lst <- list()
  Plots <- list(Simple = list(),
                Unlabelled = list(),
                Labelled = list())
  #
  # Update useProtList
  if (useProtList) {
    tst <- vapply(A, \(i) { sum(Prot[Wych[[i]], "Found_in_List"]) }, 1L)
    if (max(tst) == 0L) {
      useProtList <- FALSE
      msg <- "No data to plot for proteins in list!"
      if (proteins_split) {
        msg <- paste0(msg, " (no split plot will be created)")
        proteins_split <- FALSE
      }
    }
  }
  #
  baseCols <- "Labels"
  #
  Size <- Size[1L]
  if (!is.numeric(Size)) {
    if (!Size %in% colnames(Prot)) {
      warning("Argument Size is not a numeric value or column name and will be ignored!")
      Size <- 1
    } else {
      Prot$dot_size <- Prot[[Size]]
      if (!is.numeric(Size.min)) { Size.min <- 0.01 }
      if (!is.numeric(Size.max)) { Size.max <- 3L }
      #if (!is.numeric(Size.min)) { Size.min <- 0L }
      #if (!is.numeric(Size.max)) { Size.max <- ceiling(max(Prot$dot_size[which(is.finite(Prot$dot_size))])) }
      wF <- which(is.finite(Prot$dot_size))
      Prot$dot_size[which((Prot$dot_size > 0L) & is.infinite(Prot$dot_size))] <- max(Prot$dot_size[wF])
      Prot$dot_size[which(!is.finite(Prot$dot_size))] <- min(Prot$dot_size[wF])
      Prot$dot_size <- Size.min + (Prot$dot_size - min(Prot$dot_size))*(Size.max-Size.min)/(max(Prot$dot_size)-min(Prot$dot_size))
      baseCols <- union(baseCols, "dot_size")
    }
  } else {
    if (is.na(Size) || (Size < 0)) {
      warning("The Size parameter is invalid and will be ignored!")
      Alpha <- 1L
    }
  }
  #
  Alpha <- Alpha[1L]
  if (!is.numeric(Alpha)) {
    if (!Alpha %in% colnames(Prot)) {
      warning("Argument Alpha is not a numeric value or column name and will be ignored!")
      Alpha <- 1
    } else {
      Prot$dot_alpha <- Prot[[Alpha]]
      if (!Alpha.identity) {
        Alpha2 <- paste0("Alpha mapped to: ", Alpha)
        if ((!is.numeric(Alpha.min)) || (Alpha.min < 0L)) { Alpha.min <- 0L }
        if ((!is.numeric(Alpha.max)) || (Alpha.max > 1L)) { Alpha.max <- 1L }
        wF <- which(is.finite(Prot$dot_alpha))
        Prot$dot_alpha[which((Prot$dot_alpha > 0L) & is.infinite(Prot$dot_alpha))] <- max(Prot$dot_alpha[wF])
        Prot$dot_alpha[which(!is.finite(Prot$dot_alpha))] <- min(Prot$dot_alpha[wF])
        Prot$dot_alpha <- Alpha.min+(Prot$dot_alpha-min(Prot$dot_alpha))*(Alpha.max-Alpha.min)/(max(Prot$dot_alpha)-min(Prot$dot_alpha))
        baseCols <- union(baseCols, "dot_alpha")
      }
    }
  } else {
    if (is.na(Alpha) || (Alpha > 1) || (Alpha < 0)) {
      warning("The Alpha parameter is invalid and will be ignored!")
      Alpha <- 1
    }
  }
  #
  for (i in A) { #i <- A[1L]
    i2 <- if (useContrasts) { i } else {
      cleanNms3(i,
                experiments.map = experiments.map,
                aggregate.map = aggregate.map,
                aggregate.list = aggregate.list)
    }
    cat(" ->", i2, "\n")
    symm <- Symmetrical[i]
    xKol <- paste0(X.root, i)
    yKol <- paste0(Y.root, i)
    e <- c(xKol, yKol)
    #
    # Below: should never occur, we already checked and filtered -> candidate code for deletion:
    #e <- intersect(e, colnames(Prot))
    #stopifnot(length(e) == 2L)
    #
    plot.metrics <- Plot.metrics
    plot.metrics$Name <- ""
    w <- which(plot.metrics$Axis == "X")
    plot.metrics$Name[w] <- paste0("log2FC = ", as.character(plot.metrics$Value[w], 3L))
    w <- which(plot.metrics$Axis == "Y")
    plot.metrics$Name[w] <- paste0("P-value = ", as.character(signif(10L^-plot.metrics$Value[w], 3L)))
    plot.colours <- Plot.colours
    if (!symm) {
      plot.metrics <- plot.metrics[which(plot.metrics$Levels != "down"),]
      plot.colours <- plot.colours[, which(colnames(plot.colours) != "down"), drop = FALSE]
    }
    temp <- Prot[, baseCols, drop = FALSE]
    prot_split <- FALSE
    if (useProtList) {
      temp$"Found_in_List" <- Prot$"Found_in_List"
      prot_split <- sum(temp$"Found_in_List")
    }
    if (!Contaminants) { temp[[kontkol]] <- Prot[[kontkol]] }
    if ("Genes" %in% colnames(Prot)) { temp$Genes <- Prot$Genes }
    if (length(id.col)) { temp[, id.col] <- Prot[, id.col] }
    temp[, plotly_labels] <- Prot[, plotly_labels] # Even if plotly off please!
    #
    # Important: as.numeric() here avoids issues where some downstream code doesn't work if values are stored e.g. as matrix 
    temp$X <- as.numeric(Prot[[e[1L]]])
    temp$Y <- as.numeric(Prot[[e[2L]]])
    # This is so critical that to be on the safe side duplicate checks are in place where the issue was identified downstream.
    #
    temp <- temp[Wych[[i]],]
    #
    showHeatMap <- FALSE
    if (!misFun(X.root_ind)) {
      indSmpls <- list()
      if (!misFun(contrasts)) {
        m <- match(i, contrasts$Contrast)
        A_ <- contrasts$A_samples[[m]]
        B_ <- contrasts$B_samples[[m]]
        indSmpls[[contrasts$A[m]]] <- A_
        indSmpls[[contrasts$B[m]]] <- B_
        if (contrasts$isDouble[m]) {
          C_ <- contrasts$C_samples[[m]]
          D_ <- contrasts$D_samples[[m]]
          indSmpls[[contrasts$C[m]]] <- C_
          indSmpls[[contrasts$D[m]]] <- D_
        }
        showHeatMap <- TRUE
      } else {
        smplGrpNm <- parameters$Volcano.plots.Aggregate.Level
        smplGrpNm <- if (nchar(smplGrpNm) == 3L) {
          aggregate.map$Characteristics[[match(smplGrpNm, aggregate.map$Aggregate.Name)]]
        } else {
          gsub(";", "", smplGrpNm)
        }
        compGrpNm <- parameters$Ratios.Groups
        compGrpNm <- if (nchar(compGrpNm) == 3L) {
          aggregate.map$Characteristics[[match(compGrpNm, aggregate.map$Aggregate.Name)]]
        } else {
          gsub(";", "", compGrpNm)
        }
        m1 <- which(experiments.map[[smplGrpNm]] == i)
        A <- unique(experiments.map[m1, smplGrpNm])
        A_ <- experiments.map$Ref.Sample.Aggregate[[m1]]
        grp <- unique(experiments.map[m1, compGrpNm])
        m0 <- which((experiments.map[[compGrpNm]] == grp) & (experiments.map$Reference))
        B <- unique(experiments.map[m0, smplGrpNm])
        B_ <- experiments.map$Ref.Sample.Aggregate[[m0]]
        if (length(A_) && length(B_)) {
          showHeatMap <- TRUE
          indSmpls[[A]] <- A_
          indSmpls[[B]] <- B_
        }
      }
    }
    if (showHeatMap) {
      # For now we always Z-score:
      # When showing individual values, it makes sense to show them Z-scored to highlight changes,
      # otherwise significant values with very small logFC will look unchanged and that will be awkward.
      Zscore_heatMap <- TRUE
      # This could be controlled by a parameter in the future though.
      indCol <- paste0(X.root_ind, unlist(indSmpls))
      indDat <- Prot[Wych[[i]], paste0(X.root_ind, unlist(indSmpls))]
      if (Zscore_heatMap) {
        m <- rowMeans(indDat, na.rm = TRUE)
        sd <- apply(indDat, 1L, sd, na.rm = TRUE)
        indDat <- sweep(sweep(indDat, 1L, m, "-"), 1L, sd, "/")
      }
      indRng <- range(indDat, na.rm = TRUE)
      if (exists("ColScaleList")) {
        lowColor <- ColScaleList$`Individual Expr`[1L]
        highColor <- ColScaleList$`Individual Expr`[2L]
      } else {
        lowColor <- "red"
        highColor <- "green"
      }
      indPal <- scales::col_numeric(palette = c(lowColor, highColor),
                                    domain = indRng)
      indColDat <- matrix(paste0("<span style=\"color:", indPal(as.matrix(indDat)), "\">&#9632;</span>"), ncol = ncol(indDat))
      colnames(indColDat) <- unlist(indSmpls)
      indColDat <- as.data.frame(indColDat)
      tmpHtMp <- setNames(lapply(names(indSmpls), \(grp) { #grp <- names(indSmpls)[1L]
        do.call(paste,  c(indColDat[, unlist(indSmpls[[grp]])], sep = ""))
      }), names(indSmpls))
      tmpHtMp <- setNames(lapply(names(indSmpls), \(grp) { #grp <- names(indSmpls)[1L]
        paste("&nbsp;<span>", grp, ":</span>", tmpHtMp[[grp]])
      }), names(indSmpls))
      tmpHtMp <- as.data.frame(do.call(cbind, tmpHtMp))
      temp$HeatMap <- do.call(paste, c(tmpHtMp, sep = "<br>"))
    }
    #
    if (useFDRtbl) {
      fdr_table <- FDR_table
      if ("Sample" %in% colnames(fdr_table)) {
        fdr_table <- fdr_table[which(fdr_table$Sample == i),]
      }
      fdr.values <- fdr_table$FDR
      if (!regProvided) {
        f1 <- paste0(FDR.root, c(fdr.values, fdr.values/100L), "%")
        if ("Sample" %in% colnames(fdr_table)) {
          f1 <- paste0(f1, " - ", i)
        }
        f1 <- intersect(f1, colnames(Prot))
        f2 <- paste0("FDR=", sort(fdr.values, decreasing = TRUE), "%")
        temp[, f2] <- Prot[Wych[[i]], f1]
        test <- apply(temp[, rev(f2), drop = FALSE], 1L, \(x) {
          c(paste0("FDR ", sort(fdr.values), "%"), "")[which(c(x, "+") == "+")[1L]]
        })
        #unique(test)
        temp$FDR <- test
        temp <- temp[, setdiff(colnames(temp), f2)]
      }
    }
    #
    use_target <- FALSE
    if ("Target" %in% colnames(experiments.map)) { # i-specific!!!
      w <- if (useContrasts) {
        which(experiments.map$Ref.Sample.Aggregate %in% contrasts$A_samples[[which(contrasts$Contrast == i)]])
      } else {
        which(experiments.map[[aggregate.name]] == i)
      }
      target <- unique(unlist(strsplit(experiments.map$Target[w], ";")))
      target <- setdiff(target, c("", "NA", NA))
      target <- unique(gsub("^CON_+", "", target))
      use_target <- length(target) > 0L
    }
    if (regProvided) {
      rgKol <- regKols[i]
      temp$Colour <- Prot[Wych[[i]], rgKol]
      # Important to simplify the specific/anti-specific categories:
      temp$Colour <- gsub("^Specific.*", "Specific", temp$Colour)
      temp$Colour <- gsub("^Anti-specific.*", "Anti-specific", temp$Colour)
    } else {
      temp$Colour <- "non significant" # Default Colour values
    }
    w.u <- which(plot.metrics$Levels == "up")
    if (symm) { w.d <- which(plot.metrics$Levels == "down") }
    if ((!misFun(Ref.Ratio.values)) && (!is.null(Ref.Ratio.values))) {
      x <- data.frame(value = sort(unlist(Ref.Ratio.values[[i]])))
      if (!symm) { x <- x[which(x$value > 0L), , drop = FALSE] }
      mx <- if (X.normalized) { 0 } else { median(x$value) }
      R.thresh <- c("Upper" = NA)
      if (Ref.Ratio.method == "SD") {
        offset <- c(median(temp$X[which(is.finite(temp$X))]), 0)[X.normalized + 1L]
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
        d <- if (ratios.FDR > 0) {
          d[floor(length(d)*ratios.FDR)]
        } else {
          max(d)+0.000001
        }
        R.thresh[["Upper"]] <- d+mx
        if (symm) { R.thresh[["Lower"]] <- -(d+mx) }
      }
      R.thresh.label <- paste0(signif(R.thresh, 3L), " = ", c("upper", "lower")[1L:(symm+1L)], "-tail of ", ratios.FDR*100, "% ",
                               c("ctrl.", "within sample group")[match(Ref.Ratio.method, paste0("obs", c("1", "2")))],
                               " ratios")
      plot.metrics$Name[w.u] <- R.thresh.label[1L]
      plot.metrics$Text.value[w.u] <- plot.metrics$Value[w.u] <- R.thresh[1L]
      xlim[2L] <- max(c(xlim[2L], plot.metrics$Value[w.u]*1.1), na.rm = TRUE)  
      if (symm) {
        plot.metrics$Name[w.d] <- R.thresh.label[2L]
        plot.metrics$Text.value[w.d] <- plot.metrics$Value[w.d] <- R.thresh[2L]
        xlim[1L] <- min(c(xlim[1L], plot.metrics$Value[w.d]*1.1), na.rm = TRUE)
      }
    }
    if (mode == "curved") {
      if (!i %in% names(curved_Thresh)) {
        warning(paste0(i, " not found among the names of curved_Thresh!"))
      } else {
        #curved_Thresh <- SAM_thresh
        dec <- curved_Thresh[[i]]$decision
        mKol <- rev(colnames(dec))[1L]
        dec <- dec[match(dec[[mKol]], Prot[[mKol]]),]
        dec <- dec[Wych[[i]],]
        samS0 <- curved_Thresh[[i]]$S0
        samDF <- curved_Thresh[[i]]$degFr
        samD <- curved_Thresh[[i]]$d
        samD <- samD[which(!is.na(samD$D)),]
        samD <- samD[order(samD$FDR, decreasing = TRUE),]
        samMd <- c(median(temp$X, na.rm = TRUE), 0)[X.normalized + 1L]
      }
    }
    # Decision
    if (!regProvided) {
      if (useFDRtbl) {
        for (f3 in fdr.values) { #f3 <- fdr.values[1L]
          tstDF <- data.frame(FDR = (temp$FDR == paste0("FDR ", f3, "%")),
                              Up = logical.op(as.numeric(temp$X),
                                              plot.metrics$Test[w.u],
                                              plot.metrics$Value[w.u]))
          temp$Colour[which(tstDF$FDR)] <- "too small FC"
          w_u <- which(tstDF$FDR & tstDF$Up)
          #print(length(w_u))
          temp$Colour[w_u] <- paste0("up, FDR = ", f3, "%")
          if (symm) {
            tstDF$Down <- logical.op(temp$X,
                                     plot.metrics$Test[w.d],
                                     plot.metrics$Value[w.d])
            w_d <- which(tstDF$FDR & tstDF$Down)
            temp$Colour[w_d] <- paste0("down, FDR = ", f3, "%")
          }
        }
      } else {
        if (useFDRtbl) {
          for (f3 in samD$FDR) { #f3 <- samD$FDR[1L]
            wA <- which(dec[[paste0(f3, "FDR")]] == "+")
            w <- which((FDR_table$FDR == f3*100) & (FDR_table$Sample == i))
            if (length(wA)) {
              wU <- wA[which(temp$X[wA] > samMd)]
              temp$Colour[wU] <- FDR_table$fdr.col.up[w]
              if (symm) {
                #setdiff(wA, c(wU, wD))
                wD <- wA[which(temp$X[wA] < samMd)]
                temp$Colour[wU] <- FDR_table$fdr.col.down[w]
              }
            }
          }
        } else {
          test.u <- logical.op(temp$X,
                               plot.metrics$Test[w.u],
                               plot.metrics$Value[w.u])
          if (symm) {
            test.d <- logical.op(temp$X,
                                 plot.metrics$Test[w.d],
                                 signif(plot.metrics$Value[w.d], 3L))
          }
          if ((mode == "standard") && (!misFun(arbitrary.thresh))) {
            fdr <- as.numeric(gsub("% FDR$", "", plot.metrics$Levels[which(plot.metrics$Axis == "Y")]))/100
            fdr <- sort(fdr)
            up_Nms <- setNames(vapply(fdr, \(x) {
              w <- which(plot.metrics$Levels == paste0(x*100, "% FDR"))
              paste0("Ratio ", plot.metrics$Test[w.u], " ", signif(plot.metrics$Value[w.u], 3L),
                     ",\n  -log10(p) ", plot.metrics$Test[w], " ", signif(plot.metrics$Value[w], 3L))
            }, ""), paste0(fdr*100, "% FDR"))
            test.f <- as.data.frame(sapply(rev(fdr), \(x) { #x <- rev(fdr)[1L]
              w <- which(plot.metrics$Levels == paste0(x*100, "% FDR"))
              return(logical.op(temp$Y,
                                plot.metrics$Test[w],
                                plot.metrics$Value[w]))
            }))
            colnames(test.f) <- paste0(rev(fdr)*100, "% FDR")
            for (f in colnames(test.f)) {
              temp$Colour[which(test.f[[f]])] <- "too small FC"
              temp$Colour[which(test.u & test.f[[f]])] <- up_Nms[[f]]
            }
            myColors[up_Nms] <- plot.colours$up[match(names(up_Nms), rownames(plot.colours))]
            if (symm) {
              dwn_Nms <- setNames(vapply(fdr, \(x) {
                w <- which(plot.metrics$Levels == paste0(x*100, "% FDR"))
                paste0("Ratio ", plot.metrics$Test[w.d], " ", signif(plot.metrics$Value[w.d], 3L),
                       ",\n  -log10(p) ", plot.metrics$Test[w], " ", signif(plot.metrics$Value[w], 3L))
              }, ""), paste0(fdr*100, "% FDR"))
              for (f in colnames(test.f)) { temp$Colour[which(test.d & test.f[[f]])] <- dwn_Nms[[f]] }
              myColors[dwn_Nms] <- plot.colours$down[match(names(dwn_Nms), rownames(plot.colours))]
            }
          } else {
            w.s <- which(plot.metrics$Levels == "strict")
            w.l <- which(plot.metrics$Levels == "loose")
            up_Nms_l <- paste0("Ratio ", plot.metrics$Test[w.u], " ", signif(plot.metrics$Value[w.u], 3L),
                               ",\n  -log10(", PorQ, ") ", plot.metrics$Test[w.s], " ", signif(plot.metrics$Value[w.s], 3L))
            up_Nms_s <- paste0("Ratio ", plot.metrics$Test[w.u], " ", signif(plot.metrics$Value[w.u], 3L),
                               ",\n  -log10(", PorQ, ") ", plot.metrics$Test[w.l], " ", signif(plot.metrics$Value[w.l], 3L))
            up_Nms <- c(up_Nms_s, up_Nms_l)
            test.s <- logical.op(temp$Y,
                                 plot.metrics$Test[w.s],
                                 plot.metrics$Value[w.s])
            test.l <- logical.op(temp$Y,
                                 plot.metrics$Test[w.l],
                                 plot.metrics$Value[w.l])
            temp$Colour[which(test.l)] <- "too small FC"
            temp$Colour[which(test.u & test.l)] <- up_Nms_l
            temp$Colour[which(test.u & test.s)] <- up_Nms_s
            myColors[c(up_Nms_s, up_Nms_l)] <- plot.colours$up[match(c("strict", "loose"), rownames(plot.colours))]
            if (symm) {
              dwn_Nms_l <- paste0("Ratio ", plot.metrics$Test[w.d], " ", signif(plot.metrics$Value[w.d], 3L),
                                  ",\n  -log10(", PorQ, ") ", plot.metrics$Test[w.l], " ", signif(plot.metrics$Value[w.l], 3L))
              dwn_Nms_s <- paste0("Ratio ", plot.metrics$Test[w.d], " ", signif(plot.metrics$Value[w.d], 3L),
                                  ",\n  -log10(", PorQ, ") ", plot.metrics$Test[w.s], " ", plot.metrics$Text.value[w.s])
              dwn_Nms <- c(dwn_Nms_s, dwn_Nms_l)
              temp$Colour[which(test.d & test.l)] <- dwn_Nms_l
              temp$Colour[which(test.d & test.s)] <- dwn_Nms_s
              myColors[c(dwn_Nms_s, dwn_Nms_l)] <- plot.colours$down[match(c("strict", "loose"), rownames(plot.colours))]
            }
          }
        }
      }
    }
    #
    temp$Colour <- factor(temp$Colour, levels = names(myColors))
    colScale <- ggplot2::scale_colour_manual(name = "colour", values = myColors)
    #if ((!is.numeric(Alpha)) && (!is.numeric(Size))) {
    #if (!is.numeric(Size)) {
    #  fillScale <- ggplot2::scale_fill_manual(name = "fill", values = myColors, guide = FALSE)
    #}
    if (!regProvided) {
      Prot[[paste0("Regulated - ", i)]] <- ""
      Prot[Wych[[i]], paste0("Regulated - ", i)] <- as.character(temp$Colour)
    }
    # List of proteins of interest
    if (useProtList) {
      wLst <- which(temp$"Found_in_List")
      temp$Colour[wLst] <- "protein in list"
      if (prot_split) {
        temp$Colour2 <- "not in list"
        temp$Labels2 <- ""
        temp$Colour2[wLst] <- "protein in list"
        temp$Labels2[wLst] <- temp$Labels[wLst]
      }
    }
    #if (!is.numeric(Size)) { fillScale2 <- ggplot2::scale_fill_manual(name = "fill", values = myColors2, guide = FALSE) }
    #
    # Target
    if (use_target) { # Will overwrite "protein in list" tag with "target" tag where relevant
      if (!useProtList) {
        if (misFun(Proteins.col) || (!Proteins.col %in% colnames(Prot))) { Proteins.col <- IDs.col }
        Prot[[Proteins.col]] <- gsub(";CON_+", ";", gsub("^CON_+", "", Prot[[Proteins.col]]))
      }
      w <- grsep2(target, temp[[Proteins.col]])
      temp$Colour[w] <- "target" # The "target" tag should be in all versions of the graph.
      if (prot_split) {
        temp$Colour2[w] <- "target"
        temp$Labels2[w] <- temp$Labels[w]
      }
    }
    noLbl <- 1L:nrow(temp)
    if (sum(c("thresholds", "FDR", "regulated") %in% labels)) {
      noLbl <- which(temp$Colour %in% c("non significant", "too small FC"))
    }
    if ("proteins" %in% labels) {
      w <- which(temp$Colour %in% c("protein in list", "target"))
      noLbl <- setdiff(noLbl, w)
    }
    temp$Labels[noLbl] <- ""
    w1 <- if (useFDRtbl || (mode == "curved")) {
      which(temp$Colour %in% c(paste0("up, FDR = ", fdr.values, "%"), paste0("down, FDR = ", fdr.values, "%")))
    } else {
      if (symm) { which(temp$Colour %in% c(up_Nms, dwn_Nms)) } else { which(temp$Colour %in% up_Nms) }
    }
    w2 <- integer()
    if (prot_split) {
      w2 <- if ("Colour2" %in% colnames(temp)) { which(temp$Colour2 == "protein in list") } else {
        which(temp$Colour2 == "protein in list")
      }
      w1 <- setdiff(w1, w2) # Here it could be that we have overlap of w1 and w2 => we don't want that!
    }
    w3 <- which(temp$Colour == "target")
    w0 <- setdiff(1L:nrow(temp), c(w1, w2, w3))
    temp <- rbind(temp[w0,], temp[w1,], temp[w2,], temp[w3,])
    #temp <- temp[which(is.finite(temp$Y)),]
    ttl <- paste0(title, i2)
    subTtl <- paste0("Plotted: ", nrow(temp), " data points")
    #
    Xlab <- gsub(" - $", "", X.root)
    Ylab <- gsub(" - $", "", Y.root)
    test <- c(is.numeric(Size), is.numeric(Alpha))
    temp$"P-value" <- 10L^(-temp$Y)
    pL_lbs <- c(plotly_labels, "X", "Y")
    pL_lbs_nms <- c(names(plotly_labels), Xlab, Ylab)
    temp$plotly_labels <- apply(temp[, pL_lbs, drop = FALSE], 1L, \(x) { #x <- temp[1, kol]
      paste0(pL_lbs_nms, ": ", gsub("\n", " ", x), collapse = "<br>")
    })
    if (showHeatMap) {
      temp$plotly_labels <- do.call(paste, c(temp[, c("plotly_labels", "HeatMap")],
                         sep = paste0("<br>", c("", "Z-scored ")[Zscore_heatMap+1L], "log",
                                      X.root_ind_base, " expr. values:<br>")))
    }
    if (saveData) {
      temp$Table_labels <- apply(temp[, pL_lbs], 1L, \(x) { #x <- temp[1L, kol]
        paste0(pL_lbs_nms, ": ", gsub("\n", " ", x), collapse = "<br>")
      })
    }
    #
    #
    if (is.finite("Xlim")) { xlim <- Xlim }
    if (is.finite("Ylim")) { ylim <- Ylim }
    xFrPar <- plot.metrics$Value[which(plot.metrics$Axis == "X")]*1.1
    xlim <- c(min(c(xlim, xFrPar)),
              max(c(xlim, xFrPar)))
    ylim <- max(c(ylim, plot.metrics$Value[which(plot.metrics$Axis == "Y")]*1.1))
    xspan <- xlim[2L]-xlim[1L] # update
    #
    xTxtOffset <- xspan*0.01
    xTxtPos <- xlim[1L]+xTxtOffset
    yTxtOffset <- ylim*0.01
    #
    # Start building plots
    # (The code here has been extensively re-written to avoid the old, clunky text mining-based approach.
    # Not all parameters could be tested so there may be some hiccups.)
    guides_ <- list()
    aes_ <- list(x = rlang::expr(.data[[!!"X"]]),
                 y = rlang::expr(.data[[!!"Y"]]))
    args_gg_ <- list(shape = 16L)
    pluses_gg_ <- list(title = ggplot2::ggtitle(ttl,
                                                subtitle = subTtl),
                       xlab = ggplot2::xlab(Xlab),
                       ylab = ggplot2::ylab(Ylab),
                       theme_bw = ggplot2::theme_bw())
    if (plotly) {
      args_ly_ <- list(data = temp,
                       type = "scatter",
                       mode = "markers",
                       x = ~X,
                       y = ~Y,
                       text = ~plotly_labels,
                       hoverinfo = "text")
      layout_ly_ <- list(title = list(text = ttl,
                                      x = 0.1,
                                      xanchor = "left",
                                      y = 1,
                                      subtitle = list(text = subTtl,
                                                      lineposition = "under")),
                         xaxis = list(title = Xlab, minallowed = xlim[1L], maxallowed = xlim[2L]),
                         yaxis = list(title = Ylab, minallowed = 0L, maxallowed = ylim))
      traces_ly_ <- list()
    }
    pluses_gg_basic_ <- pluses_gg_
    # Aesthetics which can be mapped per plot
    if (length(unique(temp$Colour)) > 1L) {
      aes_$color <- rlang::expr(.data[[!!"Colour"]])
      pluses_gg_$colScale <- colScale
      guides_$color <- ggplot2::guide_legend(title = "Category")
      if (plotly) {
        args_ly_$color <- I(myColors[as.character(temp$Colour)])
      }
    } else {
      args_gg_$color <- "darkgrey"
      if (plotly) {
        args_ly_$color <- "darkgrey"
      }
    }
    # Aesthetics which (currently) can be mapped once (not per plot)
    test <- c(Size = !is.numeric(Size) | (length(Size) != 1L),
              Alpha = !is.numeric(Alpha) | (length(Alpha) != 1L))
    if (test["Size"]) {
      aes_$size <- rlang::expr(.data[[!!"dot_size"]])
      pluses_gg_$scale_size <- ggplot2::scale_size_identity(Size, guide = "legend")
      if (plotly) {
        args_ly_$size <- I(temp$dot_size*10) # *10 because it seems the dots are smaller in plotly
      }
    } else {
      args_gg_$size <- Size
      if (plotly) {
        args_ly_$size <- "darkgrey"
      }
    }
    if (test[2L]) {
      aes_$alpha <- rlang::expr(.data[[!!"dot_alpha"]])
      if (Alpha.identity) { pluses_gg_$scale_alpha <- ggplot2::scale_alpha_identity(Alpha) } else {
        guides_$alpha <- ggplot2::guide_legend(title = Alpha)
      }
      if (plotly) {
        if (!"marker" %in% names(args_ly_)) { args_ly_$marker <- list() }
        args_ly_$marker$opacity = ~dot_alpha
      }
    } else {
      args_gg_$alpha <- Alpha
      if (plotly) {
        args_ly_$alpha <- 1
      }
    }
    pluses_gg_$guides <- do.call(ggplot2::guides, guides_)
    #
    # Labels (ggplot version only)
    wLabs <- which(temp$Labels != "")
    plotLabels <- length(wLabs)
    if (plotLabels) {
      labels_dat <- temp[wLabs,]
      kolkol <- "Colour"
      labels_col <- "Labels"
      W2 <- which(!labels_dat[[kolkol]] %in% c("protein in list", "target"))
      W3 <- which(labels_dat[[kolkol]] %in% c("protein in list", "target"))
      if (length(W2) > MaxLabels) {
        labels_dat2 <- labels_dat[W2,]
        Ord <- order(abs(labels_dat2[[MaxLabels_priority]]), decreasing = TRUE) # ("abs" only needed for X, but makes no difference for Y => more concise code)
        labels_dat2 <- labels_dat2[Ord[1L:MaxLabels],]
        labels_dat <- rbind(labels_dat2, labels_dat[W3,])
      }
      labels_aes_ <- list(x = aes_$x,
                          y = aes_$y,
                          label = rlang::expr(.data[[!!labels_col]]),
                          colour = aes_$color)
      labels_args_ <- list(data = labels_dat,
                           force = 4L,
                           cex = cex,
                           lineheight = lineheight,
                           show.legend = FALSE)
      if (("alpha" %in% names(aes_)) && Alpha.labels) {
        labels_aes_ <- aes_$alpha
      } else {
        labels_args_$alpha <- 1L
      }
      labels_args_$mapping <- do.call(ggplot2::aes, labels_aes_)
      labels_layer <- do.call(ggrepel::geom_text_repel, labels_args_)
    }
    #
    # Significance thresholds
    if ((mode == "curved") && (i %in% names(curved_Thresh))) {
      # Option 1: curved SAM thresholds
      #
      # See "Uses and Misuses of the Fudge Factor in Quantitative Discovery Proteomics", Gianetto et al., Proteomics 2016
      # https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/epdf/10.1002/pmic.201600132
      # Note that in formula (9), "FC" is actually a logFC - it is clearly meant to be negative in some cases and expected to be normally distributed.
      # Also see their supplementary materials, with R code.
      SAM_thresh <- \(x, s0, conf, df) {
        -log10(2*(1 - pt(conf*(1 + s0/(abs(x)/conf - s0)), df)))
      }
      ta <- data.frame(a = c(0.01, 0.05),
                       Colour = c("blue4", "blue"))
      ta$Ta <- vapply(ta$a, \(x) { qt(1-x, samDF) }, 1)
      xlim[1L] <- min(c(-ta$Ta*samS0*2, xlim[1L]))
      xlim[2L] <- max(c(ta$Ta*samS0*2, xlim[2L]))
      pluses_gg_$curve <- list(curv_A1 = ggplot2::annotate("text", x = xTxtPos, y = ylim*0.98,
                                                           label = paste0("SAM: s0 = ", signif(samS0, 5L)),
                                                           vjust = 1L, hjust = 0L, size = 3.5),
                               curv_A2 = ggplot2::annotate("text", x = xlim[2L], y = -log10(2*(1 - pt(ta$Ta[1L], samDF)))*1.05,
                                                           label = paste0(100*(1 - ta$a[1L]), "% conf. lev."),
                                                           color = ta$Colour[1L], hjust = 1L, size = 3.5),
                               curv_A3 = ggplot2::annotate("text", x = xlim[2L], y = -log10(2*(1 - pt(ta$Ta[2L], samDF)))*1.05,
                                                           label = paste0(100*(1 - ta$a[2L]), "% conf. lev."),
                                                           color = ta$Colour[2L], hjust = 1L, size = 3.5))
      myCurves <- list(curv_F1 = list(range = c(ta$Ta[1L]*samS0, xlim[2L]), i = 1L),
                       curv_F2 = list(range = c(xlim[1L], -ta$Ta[1L]*samS0), i = 1L),
                       curv_F3 = list(range = c(ta$Ta[2L]*samS0, xlim[2L]), i = 2L),
                       curv_F4 = list(range = c(xlim[1L], -ta$Ta[2L]*samS0), i = 2L))
      curvNms <- names(myCurves)
      for (crvnm in curvNms) {
        i <- myCurves[[crvnm]]$i
        pluses_gg_$curve[[curvNms]] <- ggplot2::stat_function(fun = \(x) { SAM_thresh(x, samS0, ta$Ta[i], samDF) },
                                                              color = ta$Colour[i],
                                                              xlim = myCurves[[crvnm]]$range,
                                                              linetype = "dotted")
      }
      if (plotly) {
        xCurveRanges2 <- setNames(lapply(curvNms, \(crvnm) {
          x <- myCurves[[crvnm]]$range
          c(0L:100L/100)*(x[[2L]]-x[[1L]])+x[[1L]]
        }), curvNms)
        for (crvnm in curvNms) {
          i <- myCurves[[crvnm]]$i
          traces_ly_[[crvnm]] <- list(Trace = "plotly::add_segments",
                                      x = xCurveRanges2[[x]],
                                      y = SAM_thresh(x, samS0, ta$Ta[i], samDF),
                                      type = "scatter",
                                      mode = "lines",
                                      line = list(color = ta$Colour[i],
                                                  dash = "dot"),
                                      name = crvnm,
                                      inherit = FALSE)
        }
      }
    } else {
      # Option 2: straight horizontal/vertical thresholds
      wVL <-  which( { if (symm) { plot.metrics$Axis == "X" } else { plot.metrics$Levels == "up" } } & (plot.metrics$Value != 0L) )
      wHL <- which(plot.metrics$Axis == "Y")
      if (length(c(wVL, wHL))) {
        pluses_gg_$newScale_1 <- ggnewscale::new_scale_color()
        pluses_gg_$newScale_2 <- ggplot2::scale_color_identity()
      }
      if (length(wVL)) {
        pluses_gg_$vlines <- list()
        plot.metrics$hOffset <- plot.metrics$Value + xTxtOffset*c(-1, 1)[match(plot.metrics$Levels, c("down", "up"))]
        pluses_gg_$vlines$L <- ggplot2::geom_vline(data = plot.metrics[wVL,], linetype = "dotted",
                                                   ggplot2::aes(xintercept = Value, color = Colour))
        pluses_gg_$vlines$A <- ggplot2::geom_text(data = plot.metrics[wVL,],
                                                  ggplot2::aes(label = Name, x = hOffset, color = Colour),
                                                  y = ylim*0.8, cex = cex*1.2, hjust = 1L, angle = 90)
        if (plotly) {
          for (ii in seq_along(wVL)) { # Apparently colors cannot be vectorized...
            traces_ly_$vlines[[paste0("L", ii)]] <- list(Trace = "plotly::add_segments",
                                                         x = plot.metrics$Value[wVL[ii]], xend = plot.metrics$Value[wVL[ii]],
                                                         y = 0, yend = ylim,
                                                         type = "scatter", mode = "lines",
                                                         line = list(color = plot.metrics$Colour[wVL[ii]], dash = "dot"),
                                                         inherit = FALSE, text = plot.metrics$Name[wVL[ii]], hoverinfo = "text",
                                                         showlegend = FALSE)
          }
        }
      }
      if (length(wHL)) {
        pluses_gg_$hlines <- list()
        pluses_gg_$hlines$L <- ggplot2::geom_hline(data = plot.metrics[wHL,], linetype = "dotted",
                                                   ggplot2::aes(yintercept = Value, color = Colour))
        pluses_gg_$hlines$A <- ggplot2::geom_text(data = plot.metrics[wHL,], x = -xTxtPos, cex = cex*1.2,
                                                  hjust = 1L, show.legend = FALSE,
                                                  ggplot2::aes(label = Name, y = Value + yTxtOffset,
                                                               color = Colour))
        if (plotly) {
          for (ii in seq_along(wHL)) { # Apparently colors cannot be vectorized...
            traces_ly_$hlines[[paste0("L", ii)]] <- list(Trace = "plotly::add_segments",
                                                         x = xlim[1L], xend = xlim[2L],
                                                         y = plot.metrics$Value[wHL[ii]], yend = plot.metrics$Value[wHL[ii]],
                                                         type = "scatter", mode = "lines",
                                                         line = list(color = plot.metrics$Colour[wHL[ii]], dash = "dot"),
                                                         inherit = FALSE, text = plot.metrics$Name[wHL[ii]], hoverinfo = "text",
                                                         showlegend = FALSE)
          }
        }
      }
    }
    if (useFDRtbl && nrow(fdr_table)) {
      pluses_gg_$fdr <- list(L = ggplot2::geom_hline(data = fdr_table, linetype = "dotted",
                                                     ggplot2::aes(yintercept = -log10(Thresholds), colour = fdr.col.line),
                                                     show.legend = FALSE),
                             T = ggplot2::geom_text(data = fdr_table, x = xTxtPos, cex = cex*1.2,
                                                    hjust = 0L, show.legend = FALSE,
                                                    ggplot2::aes(label = paste0("FDR = ", FDR, "%"),
                                                                 y = -log10(Thresholds) + yTxtOffset,
                                                                 colour = fdr.col.line)))
      if (plotly) {
        for (ii in 1L:nrow(fdr_table)) { # Apparently colors cannot be vectorized...
          traces_ly_[[paste0("F", ii)]] <- list(Trace = "plotly::add_segments",
                                                x = xlim[1L], xend = xlim[2L],
                                                y = -log10(fdr_table$Thresholds[[ii]]), yend = -log10(fdr_table$Thresholds[[ii]]),
                                                type = "scatter", mode = "lines",
                                                line = list(color = fdr_table$fdr.col.line[[ii]], dash = "dot"),
                                                inherit = FALSE, text = paste0("FDR = ", fdr_table$FDR[[ii]], "%"), hoverinfo = "text",
                                                showlegend = FALSE)
        }
      }
    }
    # Arbitrary lines
    if ((!misFun(arbitrary.lines)) && (!is.null(arbitrary.lines)) && is.data.frame(arbitrary.lines) && nrow(arbitrary.lines)) {
      w1 <- which(is.na(arbitrary.lines$xintercept) & arbitrary.lines$slope != 0L)
      w2 <- which(is.na(arbitrary.lines$xintercept) & arbitrary.lines$slope == 0L)
      w3 <- which(!is.na(arbitrary.lines$xintercept))
      if (length(w1)) {
        pluses_gg_$Arbi1 <- ggplot2::geom_abline(data = arbitrary.lines[w1,], linetype = "dotted",
                                                 ggplot2::aes(intercept = yintercept, slope = slope,
                                                              colour = colour))
      }
      if (length(w2)) {
        pluses_gg_$Arbi2l <- ggplot2::geom_hline(data = arbitrary.lines[w2,], linetype = "dotted",
                                                 ggplot2::aes(yintercept = yintercept, colour = colour)) 
        pluses_gg_$Arbi2a <-  ggplot2::geom_text(data = arbitrary.lines[w2,], cex = cex * 1.2, hjust = 1L, x = -xTxtPos,
                                                 ggplot2::aes(y = yintercept + yTxtOffset, label = label, colour = colour))
      }
      if (length(w3)) {
        arbitrary.lines$Sign <- sign(arbitrary.lines$xintercept)
        pluses_gg_$Arbi3l <- ggplot2::geom_vline(data = arbitrary.lines[w3,], linetype = "dotted",
                                                 ggplot2::aes(xintercept = xintercept,  colour = colour))
        pluses_gg_$Arbi3a <- ggplot2::geom_text(data = arbitrary.lines[w3,], cex = cex * 1.2, hjust = 0L, y = ylim * 0.5,
                                                ggplot2::aes(x = xintercept + xTxtOffset * Sign,
                                                             label = label, colour = colour))
      }
      if (plotly) {
        w12 <- c(w1, w2)
        if (length(w12)) {
          for (ii in seq_along(w12)) { # Apparently colors cannot be vectorized...
            jj <- w12[ii]
            traces_ly_[[paste0("Arbi12_", ii)]] <- list(Trace = "plotly::add_segments",
                                                        x = xlim[1L], xend = xlim[2L],
                                                        y = arbitrary.lines$slope[[jj]] * xlim[1L] + arbitrary.lines$yintercept[[jj]],
                                                        yend = arbitrary.lines$slope[[jj]] * xlim[2L] + arbitrary.lines$yintercept[[jj]],
                                                        type = "scatter", mode = "lines",
                                                        line = list(color = arbitrary.lines$colour[[jj]], dash = "dot"),
                                                        text = arbitrary.lines$label[[jj]], hoverinfo = "text",
                                                        inherit = FALSE, showlegend = FALSE)
          }
        }
        if (length(w3)) {
          for (ii in seq_along(w3)) { # Apparently colors cannot be vectorized...
            jj <- w3[ii]
            traces_ly_[[paste0("Arbi3_", ii)]] <- list(Trace = "plotly::add_segments",
                                                       x = arbitrary.lines$xintercept[[jj]], xend = arbitrary.lines$xintercept[[jj]],
                                                       y = 0, yend = ylim, type = "scatter", mode = "lines",
                                                       line = list(color = arbitrary.lines$colour[[jj]], dash = "dot"),
                                                       text = arbitrary.lines$label[[jj]], hoverinfo = "text",
                                                       inherit = FALSE, showlegend = FALSE)
          }
        }
      }
    }
    # Arbitrary thresholds
    if ((!misFun(arbitrary.thresh)) && (!is.null(arbitrary.thresh)) && is.data.frame(arbitrary.thresh) && nrow(arbitrary.thresh)) {
      w <- which(is.na(arbitrary.thresh$xintercept) | (arbitrary.thresh$slope %in% c(NA, 0)))
      if (length(w)) {
        warning("Significance thresholds are horizontal, \"xintercept\" will be ignored!")
        arbitrary.thresh$xintercept[w] <- NA_real_
        arbitrary.thresh$slope[w] <- 0
      }
      pluses_gg_$Arbi4l <- ggplot2::geom_hline(data = arbitrary.thresh, linetype = "dotted",
                                               ggplot2::aes(yintercept = yintercept, colour = colour)) 
      pluses_gg_$Arbi4a <-  ggplot2::geom_text(data = arbitrary.thresh, cex = cex * 1.2, hjust = 0L, x = xTxtPos,
                                               ggplot2::aes(y = yintercept + yTxtOffset, label = label, colour = colour))
      if (plotly) {
        for (ii in 1L:nrow(arbitrary.thresh)) { # Apparently colors cannot be vectorized...
          traces_ly_[[paste0("Arbi4_", ii)]] <- list(Trace = "plotly::add_segments",
                                                     x = arbitrary.thresh$xintercept[[ii]], xend = arbitrary.thresh$xintercept[[ii]],
                                                     y = 0, yend = ylim, type = "scatter", mode = "lines",
                                                     line = list(color = arbitrary.thresh$colour[[ii]], dash = "dot"),
                                                     text = arbitrary.lines$label[[jj]], hoverinfo = "text",
                                                     inherit = FALSE, showlegend = FALSE)
        }
      }
    }
    #
    # Final ggplot2 layers
    pluses_gg_basic_$theme1 <- ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                                              panel.grid.minor = ggplot2::element_blank())
    pluses_gg_$scale_x <- pluses_gg_basic_$scale_x <- ggplot2::scale_x_continuous(limits = xlim,
                                                                                  expand = c(0L, 0L))
    pluses_gg_$scale_y <- pluses_gg_basic_$scale_y <- ggplot2::scale_y_continuous(limits = c(0, ylim),
                                                                                  expand = c(0L, 0L))
    #
    # Evaluate plots
    args_gg_$mapping <- do.call(ggplot2::aes, aes_)
    main_layer <- do.call(ggplot2::geom_point, args_gg_)
    plot <- simPlot <- ggplot(temp) + main_layer
    if (plotLabels) {
      plot_labels <- plot + labels_layer
    }
    for (ii in seq_along(pluses_gg_)) { #ii <- 1L #ii <- ii+1L
      if (inherits(pluses_gg_[[ii]], "list") && (!inherits(pluses_gg_[[ii]], "gg"))) {
        for (jj in seq_along(pluses_gg_[[ii]])) { #jj <- 1L #jj <- jj+1L
          plot <- plot + pluses_gg_[[ii]][[jj]]
          if (plotLabels) { plot_labels <- plot_labels + pluses_gg_[[ii]][[jj]] }
        }
      } else {
        plot <- plot + pluses_gg_[[ii]]
        if (plotLabels) { plot_labels <- plot_labels + pluses_gg_[[ii]] }
      }
    }
    for (ii in seq_along(pluses_gg_basic_)) { #ii <- 1L #ii <- ii+1L
      if (inherits(pluses_gg_basic_[[ii]], "list") && (!inherits(pluses_gg_basic_[[ii]], "gg"))) {
        for (jj in seq_along(pluses_gg_basic_[[ii]])) { #jj <- 1L #jj <- jj+1L
          simPlot <- simPlot + pluses_gg_basic_[[ii]][[jj]]
        }
      } else {
        simPlot <- simPlot + pluses_gg_[[ii]]
      }
    }
    #poplot(simPlot)
    #poplot(plot)
    #poplot(plot_labels)
    Plots$Simple[[ttl]] <- plotEval(simPlot)
    Plots$Unlabelled[[ttl]] <- plotEval(plot)
    if (plotLabels) { Plots$Labelled[[ttl]] <- plotEval(plot_labels) }
    if (prot_split) {
      # This we do just for ggplot2, not with plotly (there will be a better solution)
      args_gg_prot <- args_gg_
      aes_prot <- aes_
      aes_prot$color <- rlang::expr(.data[[!!"Colour2"]])
      pluses_gg_prot <- pluses_gg_
      pluses_gg_prot$colScale <- colScale2
      args_gg_prot <- args_gg_
      args_gg_prot$mapping <- do.call(ggplot2::aes, aes_prot)
      main_layer_prot <- do.call(ggplot2::geom_point, args_gg_prot)
      plot_prot <- ggplot(temp) + main_layer_prot
      if (plotLabels) {
        labels_dat_prot <- labels_dat[which(labels_dat$Colour2 %in% c("protein in list", "target")),]
        labels_args_prot_ <- labels_args_
        labels_args_prot_$data <- labels_dat_prot
        labels_args_prot_$mapping$colour <- main_layer_prot$mapping$colour
        labels_layer_prot <- do.call(ggrepel::geom_text_repel, labels_args_prot_)
        plot_prot_labels <- plot_prot + labels_layer_prot
      }
      for (ii in seq_along(pluses_gg_prot)) { #ii <- 1L #ii <- ii+1L
        if (inherits(pluses_gg_prot[[ii]], "list") && (!inherits(pluses_gg_prot[[ii]], "gg"))) {
          for (jj in seq_along(pluses_gg_prot[[ii]])) { #jj <- 1L #jj <- jj+1L
            plot_prot <- plot_prot + pluses_gg_prot[[ii]][[jj]]
            if (plotLabels) { plot_prot_labels <- plot_prot_labels + pluses_gg_prot[[ii]][[jj]] }
          }
        } else {
          plot_prot <- plot_prot + pluses_gg_prot[[ii]]
          if (plotLabels) { plot_prot_labels <- plot_prot_labels + pluses_gg_prot[[ii]] }
        }
      }
      #poplot(plot_prot)
      #poplot(plot_prot_labels)
      Plots$"Proteins in list - unlabelled"[[ttl]] <- plotEval(plot_prot)
      if (plotLabels) { Plots$"Proteins in list - labelled"[[ttl]] <- plotEval(plot_prot_labels ) }
      #
      Symb <- c("circle", "star")[temp$Found_in_List+1L]
      args_ly_$symbol <- I(Symb)
    }
    #
    poplot(simPlot)
    # Make valid file name
    tr <- gsub("/|:|\\*|\\?|<|>|\\||/", "-", title.root)
    tt <- gsub("/|:|\\*|\\?|<|>|\\||/", "-", ttl)
    nm <- paste0(tr, tt)
    if (nchar(nm) > 98L) { nm <- substr(nm, 1L, 98L) }
    if (nm %in% plotNms) {
      fixkount <- 0L
      while (nm %in% plotNms) {
        fixkount <- fixkount + 1L
        if (fixkount == 100L) {
          stop("Really? Really?!?! You really have 100 similarly named conditions with very long names?!?! If you tried to break this function then you succeeded with (bad) style!")
        }
        nm <- paste0(substr(nm, 1L, 93L), "...", c("0", "")[(nchar(fixkount) > 1L)+1L],
                     fixkount)
      }
    }
    plotNms <- c(plotNms, nm)
    # Save data for users to replot with their own methods if they feel so inclined
    if (saveData) {
      flPth <- paste0(subfolder, "/", nm, "_dat.csv")
      tmpDat <- data.frame("ID" = gsub("\n|<br>", " | ", temp$Table_labels),
                           "log2FC" = temp$X,
                           "P-value" = 10L^-temp$Y,
                           "-log10(P-value)" = temp$Y,
                           "adj. P-value" = p.adjust(10L^-temp$Y, "BH"),
                           check.names = FALSE)
      data.table::fwrite(tmpDat, flPth, sep = ",", row.names = FALSE, na = "NA")
    }
    if (plotly) {
      plotLy <- do.call(plotly::plot_ly,
                        args_ly_)
      for (ii in seq_along(traces_ly_)) { #ii <- 1L #ii <- ii+1L
        if ("Trace" %in% names(traces_ly_[[ii]])) {
          tmp <- traces_ly_[[ii]]
          tmp$p <- plotLy
          fun <- eval(parse(text = tmp$Trace))
          tmp$Trace <- NULL
          plotLy <- do.call(fun, tmp)
        } else {
          for (jj in seq_along(traces_ly_[[ii]])) { #jj <- 1L #jj <- jj+1L
            tmp <- traces_ly_[[ii]][[jj]]
            tmp$p <- plotLy
            fun <- eval(parse(text = tmp$Trace))
            tmp$Trace <- NULL
            plotLy <- do.call(fun, tmp)
          }
        }
      }
      #
      layout_ly_$p <- plotLy
      plotLy <- do.call(plotly::layout, layout_ly_)
      # Remove selection tools
      plotLy <- plotly::config(plotLy,
                               modeBarButtonsToRemove = c("select2d", "lasso2d"))
      # Build
      plotLy <- plotly::plotly_build(plotLy)
      # Add search java-script
      searchCol <- intersect(setdiff(plotly_labels, "PEP"), colnames(temp))
      searchDat <- as.matrix(temp[, searchCol, drop = FALSE])
      plotLy <- .plotlySearch(plotLy, searchDat, 0.1)
      #
      volcPlotly[[ttl]] <- list(Ttl = ttl,
                                Plot = plotLy)
      
    }
    plotMetr.lst[[i]] <- plot.metrics
  }
  if (save) {
    cat(" ---> Saving ggplots...\n")
    # Move me out of the function!!!
    sfpt <- if (subfolderpertype) { paste0(subfolder, "/", saveExt[1L]) } else { subfolder }
    if (!dir.exists(sfpt)) { dir.create(sfpt, recursive = TRUE) }
    lst0 <- setNames(lapply(names(Plots$Simple), \(nm) {
      list(Path = sfpt,
           Ttl = paste0("_", nm),
           Plot = Plots$Simple[[nm]],
           Ext = saveExt[1L])
    }), paste0("_", names(Plots$Simple)))
    lst1 <- setNames(lapply(names(Plots$Unlabelled), \(nm) {
      list(Path = sfpt,
           Ttl = nm,
           Plot = Plots$Unlabelled[[nm]],
           Ext = saveExt[1L])
    }), paste0(names(Plots$Unlabelled), "_noLabel"))
    plotsLst <- c(lst0, lst1)
    if (length(Plots$Labelled)) {
      lst2 <- setNames(lapply(names(Plots$Labelled), \(nm) {
        list(Path = sfpt,
             Ttl = paste0(nm, "_tags"),
             Plot = Plots$Labelled[[nm]],
             Ext = saveExt[1L])
      }), paste0(names(Plots$Labelled), "_Label"))
      plotsLst <- c(plotsLst, lst2)
    }
    if (prot_split && length(Plots$"Proteins in list - unlabelled")) {
      lst3 <- setNames(lapply(names(Plots$"Proteins in list - unlabelled"), \(nm) {
        list(Path = sfpt,
             Ttl = paste0(nm, "_list"),
             Plot = Plots$"Proteins in list - unlabelled"[[nm]],
             Ext = saveExt[1L])
      }), paste0(names(Plots$"Proteins in list - unlabelled"), "_list_noLabel"))
      plotsLst <- c(plotsLst, lst3)
      if (length(Plots$"Proteins in list - labelled")) {
        lst4 <- setNames(lapply(names(Plots$"Proteins in list - labelled"), \(nm) {
          list(Path = sfpt,
               Ttl = paste0(nm, "_list_tags"),
               Plot = Plots$"Proteins in list - labelled"[[nm]],
               Ext = saveExt[1L])
        }), paste0(names(Plots$"Proteins in list - labelled"), "_list_Label"))
        plotsLst <- c(plotsLst, lst4)
      }
    }
    l <- length(saveExt)
    if (l > 1L) {
      lst0 <- plotsLst
      for (j in 2L:l) {
        sfpt <- if (subfolderpertype) { paste0(subfolder, "/", saveExt[j]) } else { subfolder }
        if (!dir.exists(sfpt)) { dir.create(sfpt, recursive = TRUE) }
        lst <- lapply(lst0, \(x) {
          x$Path <- sfpt
          x$Ext <- saveExt[j]
          return(x)
        })
        plotsLst <- c(plotsLst, lst)
      }
    }
    #parallel::clusterExport(cl, "plotsLst", envir = environment()) #Slows down
    F0 <- .bind_worker(.ggSave_worker,
                       list(plotsLst = plotsLst))
    tst <- parallel::parLapply(cl,
                               plotsLst,
                               F0)
  }
  #
  thrsh <- list(Absolute = plotMetr.lst)
  if (useFDRtbl) { thrsh$FDR <- FDR_table }
  RES <- list(Thresholds = thrsh)
  if (return) { RES$Protein_groups_file <- Prot }
  if (return.plot) { RES$Plots <- Plots }
  if (plotly) {
    RES$"Plotly plots" <- volcPlotly
  }
  #
  setwd(origWD)
  #
  if (stopCl) { parallel::stopCluster(cl) }
  return(RES)
}
