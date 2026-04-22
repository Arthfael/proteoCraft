#' Time.Profile
#'
#' @description 
#' A function to plot the time profile of protein(s).
#' The function assumes that expression column names follow the following patttern:
#' "root.cond.tp" (i.e. a root string then the condition then the time point).
#' If there are error columns, these should follow a similar patttern with their own root.
#' 
#' @param df The protein/protein groups file, which must include expression data (and SE)
#' @param prot Optional: the protein ID to be plotted. If left at its default value of "", all protein groups will be plotted on one graph.
#' @param prot.col The name of the Protein IDs columns. Default = "Protein IDs"
#' @param exprs.root The root of the names of expression/ratio columns. Default = "Ratio_Mean log2_"
#' @param exprs.name Optional. What should we call the expression dimension? Default = "Ratio"
#' @param error.root Optional argument (default = FALSE). The root of the names of error columns, if any are provided. It is probably not a good idea to provide it when plotting thousands of proteins.
#' @param conds The conditions to plot.
#' @param tp The time points to plot. This should be a named numeric (the name is what the function expects to find in the expression and error column names).
#' @param save Should the plot be saved? default = FALSE. Set it to a valid ggsave name to save it to the required format, e.g. "Time profiles.jpeg".
#' @param title Optional title (default = "").
#' @param labels default = FALSE. If set to a column name, will display these as protein labels to the left of each profile. Ignored when only plotting a single protein ID.
#' 
#' @examples
#' # To plot a single protein's profile:
#' Time.Profiledf(df = PG, prot = "Q23098", prot.col = "Protein IDs",
#'                          exprs.root = "Ratio_Mean log2_ ", error.root = "SE.",
#'                          conds = Cond, tp = Time.Pts, save = FALSE)
#' # To plot the profiles of all the proteins in the data frame: 
#' Time.Profile(df = PG, prot.col = "Protein IDs",
#'                          exprs.root = "Ratio_Mean log2_", error.root = FALSE,
#'                          conds = Cond, tp = Time.Pts, save = FALSE)
#' @export

Time.Profile <- function(df,
                         prot = "",
                         prot.col = "Protein IDs",
                         exprs.root = "Ratio_Mean log2_",
                         exprs.name = "Ratio",
                         error.root = FALSE,
                         conds,
                         tp,
                         save = FALSE,
                         title = "",
                         labels = FALSE) {
  a1 <- as.character(sapply(conds, \(x) {
    vapply(names(tp), \(y) { paste0(exprs.root, x, ".", y) }, "")
  }))
  if (!is.character(error.root)) { error.root <- as.character(error.root) }
  plotRibbon <- FALSE
  if (error.root != "FALSE") {
    a2 <- as.character(sapply(conds, \(x) {
      vapply(names(tp), \(y) { paste0(error.root, x, ".", y) }, "")
    }))
    plotRibbon <- a2 %in% colnames(df)
    if (!plotRibbon) {
      warning("The \"error.root\" parameter provided cannot be mapped to any columns! Errors margins will not be plotted.")
    }
  }
  if (prot != "") {
    test <- grep(prot, PG[[prot.col]])
    for (i in test) {
      temp <- data.frame(Expression = as.numeric(df[i, a1]))
      temp$Condition <- as.character(sapply(conds, \(x) { rep(x, length(tp)) }))
      temp$Time.point <- as.numeric(sapply(1L:length(conds), \(x) { tp }))
      plot <- ggplot2::ggplot(temp) +
        ggplot2::geom_line(data = temp,
                           ggplot2::aes(x = Time.point, y = Expression, group = Condition,
                                        colour = Condition)) +
        ggplot2::ggtitle(paste0("Temporal profile: ", prot)) + ggplot2::theme_bw() +
        ggplot2::labs(y = exprs.name) + ggplot2::ggtitle(title)
      ggplot2::labs(y = exprs.name)
      if (plotRibbon) {
        temp$Error <- as.data.frame(as.numeric(df[i,a2]))
        temp$min <- temp$Expression - temp$Error
        temp$max <- temp$Expression + temp$Error
        plot <- plot +
          ggplot2::geom_ribbon(data = temp,
                               ggplot2::aes(x = Time.point, ymin = min, ymax = max,
                                            group = Condition, fill = Condition), alpha = 0.2)
      }
      poplot(plot)
      if (as.character(save) != "FALSE") {
        save <- gsub("^jpg$", "jpeg", gsub("^\\.", "", save))
        if (save == "pdf") {
          ggplot2::ggsave(paste0("Time profile_", df[i, prot.col], ".pdf"), plot)
        }
        if (save %in% c("jpeg", "tiff", "png", "bmp")) {
          ggplot2::ggsave(paste0("Time profile_", df[i, prot.col], ".", save), plot,
                          dpi = 300L, width = 10L, height = 10L, units = "in")
        }
      }
    }
  } else {
    if (labels != FALSE) {
      temp <- reshape2::melt(df[, c(prot.col, labels, a1)])
      colnames(temp) <- c("IDs", "Tag", "variable", "Expression")
      temp$Tag <- vapply(strsplit(temp$Tag, split = ";"), \(x) {
        x <- unlist(x)
        if (length(x) == 1L) { return(x) }
        return(paste0(x[1L], "..."))
      }, "")
    } else {
      temp <- reshape2::melt(df[, c(prot.col, a1)])
      colnames(temp) <- c("IDs", "variable", "Expression")
    }
    temp$IDs <- vapply(strsplit(temp$IDs, split = ";"), \(x) {
      x <- unlist(x)
      if (length(x) == 1L) { return(x) }
      return(paste0(x[1L], "..."))
    }, "")
    c1 <- gsub(topattern(exprs.root), "", temp$variable)
    temp$Time.point <- gsub(paste(paste0("^", gsub("\\.", "\\\\.", conds), "\\."), collapse = "|"), "", c1)
    temp$Time.point <- tp[as.character(temp$Time.point)]
    temp$Condition <- gsub(paste(paste0("\\.", gsub("\\.", "\\\\.", names(tp)), "$"), collapse = "|"), "", c1)
    temp$X <- do.call(paste, c(temp[c("IDs", "Condition")], sep = "_"))
    plot <- ggplot2::ggplot(temp) +
      ggplot2::geom_line(data = temp,
                ggplot2::aes(x = Time.point, y = Expression, group = IDs, color = IDs)) +
      ggplot2::ggtitle(paste0("Global temporal profile", prot)) +
      ggplot2::facet_wrap(~Condition) + ggplot2::theme_bw() +
      ggplot2::labs(y = exprs.name) + ggplot2::ggtitle(title)
    if (plotRibbon) {
      temp2 <- reshape2::melt(df[, c(prot.col, a2)])
      temp$Error <- temp2$value
      temp$min <- temp$Expression - temp$Error
      temp$max <- temp$Expression + temp$Error
      plot <- plot +
        ggplot2::geom_ribbon(data = temp,
                             ggplot2::aes(x = Time.point, ymin = min, ymax = max,
                                          group = IDs, fill = IDs),
                    alpha = 0.2)
    }
    if (labels != FALSE) {
      x1 <- min(is.all.good(temp$Time.point))
      x2 <- max(is.all.good(temp$Time.point))
      temp2 <- temp[which(temp$Time.point == x1),]
      plot <- plot +
        ggplot2::geom_text(data = temp2,
                           ggplot2::aes(label = Tag, x = Time.point, y = Expression,
                                        group = IDs, color = IDs), hjust = 0)
    }
    plot <- if (plotRibbon) {
      plot +
        ggplot2::guides(group = FALSE, color = FALSE, fill = FALSE)
    } else {
      plot +
        ggplot2::guides(group = FALSE, color = FALSE)
    }
    poplot(plot)
    if (as.character(save) != "FALSE") {
      save <- gsub("^jpg$", "jpeg", gsub("^\\.", "", save))
      if (save == "pdf") {
        ggplot2::ggsave(paste0("Global time profile.pdf"), plot)
      }
      if (save %in% c("jpeg", "tiff", "png", "bmp")) {
        ggplot2::ggsave(paste0("Global time profile.", save), plot,
                        dpi = 300L, width = 10L, height = 10L, units = "in")
      }
    }
  }
}
