#' basic.heatmap
#'
#' @description
#' A function to create a simple heatmap of a data frame.
#' 
#' @param matr The data frame or matrix to plot. If it has rownames or colnames, these will be displayed.
#' @param title The title of the graph.
#' @param size The character size. Default = 2.5, you may want to reduce it for large proteins.
#' @param h.margins Optional: size of the horizontal margins. Default = c(5, 10)
#' @param v.margin Optional: size of the upper margin. Default = 3
#' @param colours Since October 2024, this should be a viridis scale name. The previous approach did not work.
#' @param colours_invert Default = FALSE. If TRUE, the viridis scale of choice is reverted.
#' @param na.colour Optional: What colour should NA values be displayed as?
#' @param print TRUE by default. But can be set to FALSE because sometimes you just want to save a plot for later without printing it.
#' @param save FALSE by default. Set to format you want pictures saved as ("jpeg", "png", "pdf"...)
#' @param folder Where to save the graph. Default = getwd()
#' @param return FALSE by default. Set to TRUE to return the created plot.
#' @param print_values If TRUE (default), values are printed on the heatmap.
#' @param normRows Default = FALSE. If TRUE, the data is re-normalized by row.
#' @param isLog Default = FALSE. Is the data log-transformed?
#' @param hTree If TRUE (default), a horizontal (x) dendrogram will be plotted.
#' @param vTree If TRUE (default), a vertical (y) dendrogram will be plotted.
#' @param split Logical (default = FALSE), allows splitting labels over several lines.
#' 
#' @returns
#' If return is TRUE, the evaluated ggplot is returned.
#' 
#' @export

basic.heatmap <- function(matr,
                          title = "Heatmap",
                          subtitle = "",
                          normRows = FALSE,
                          isLog = FALSE,
                          size = 2L,
                          h.margins = c(5L, 10L),
                          v.margin = 10L,
                          colours = "D",
                          colours_invert = FALSE,
                          na.colour = "black",
                          print = TRUE,
                          save = FALSE,
                          folder,
                          return = FALSE,
                          print_values = FALSE,
                          hTree = TRUE,
                          vTree = TRUE,
                          split = FALSE) {
  #DefArg(basic.heatmap); h.margins <- v.margin <- 2.5
  #matr <- temp2; subtitle <- paste0(tstrt, "\n(", tolower(bee), ")"); colours = c(darkblue = 0, red = 1)
  #matr <- pepSLA1[, exprsCol]; title = "Normalization - before"; save = TRUE
  # lC <- length(colours)
  # if (!lC %in% 2L:3L) { stop("Invalid colour names provided!") }
  # test <- sapply(c(names(colours), na.colour), \(x) {
  #   try(grDevices::col2rgb(x), silent = TRUE)
  # })
  # if (!"matrix" %in% class(test)) { stop("Invalid colour names provided!") }
  # colDF <- as.data.frame(t(test[, 1L:lC]))
  # colnames(colDF) <- c("red", "green", "blue")
  # colDF$RangeVal <- colours
  # if (lC == 3L) { colDF$Range <- c("low", "mid", "high") } else {
  #   colDF <- rbind(colDF, colMeans(colDF))
  #   colDF$Range <- c("low", "high", "mid")
  #   colDF <- colDF[match(c("low", "mid", "high"), colDF$Range),]
  # }
  # colDF$hex <- apply(colDF[, c("red", "green", "blue")], 1, \(x) {
  #   x <- x/255
  #   grDevices::rgb(x[[1L]], x[[2L]], x[[3L]])
  # })
  if (!colours %in% c("A", "B", "C", "D", "E", "F", "G", "H")) { warning("Invalid colour scale name provided, defaulting to \"D\" ( = viridis)!") }
  if (missing(folder)) { folder <- getwd() }
  if (length(h.margins) == 1L) { h.margins <- c(h.margins, 2L) }
  if (normRows) {
    if (isLog) {
      normVect <- log10(rowMeans(10L^matr, na.rm = TRUE))
      matr <- sweep(matr, 1L, normVect, "-")
    } else {
      normVect <- rowMeans(matr, na.rm = TRUE)
      matr <- sweep(matr, 1L, normVect, "/")
    }
  }
  Nrow <- nrow(matr)
  Ncol <- ncol(matr)
  matr2 <- as.matrix(matr)
  h.Marg <- h.margins
  v.Marg <- v.margin
  #h.Padd <- max(nchar(unlist(strsplit(rownames(matr2), "\n"))))/25
  h.Padd <- 0L
  clnms <- colnames(matr2)
  if (split) {
    clnms <- unlist(strsplit(colnames(matr2), "\n"))
  }
  v.Padd <- max(nchar(clnms))/10
  h.Marg[1L] <- h.Marg[1L]+h.Padd
  v.Marg[1L] <- v.Marg[1L]+v.Padd
  # Dendrograms
  if (hTree+vTree) {
    matr3 <- as.matrix(Data_Impute2(matr)$Imputed_data)
    if (hTree) {
      Hclust <- hclust(dist(t(matr3)))
      Hdendro <- as.dendrogram(Hclust)
      Hddata <- ggdendro::dendro_data(Hdendro)
      Hlabs <- ggdendro::label(Hddata)
      # Re-order our matrix based on extracted dendrogram labels
      matr2 <- matr2[, match(Hlabs$label, colnames(matr2))]
      if (is.null(rownames(matr2))) { rownames(matr2) <- paste0("Sample", as.character(1L:Nrow)) }
    }
    if (vTree) {
      Vclust <- hclust(dist(matr3))
      Vdendro <- as.dendrogram(Vclust)
      Vddata <- ggdendro::dendro_data(Vdendro)
      Vlabs <- ggdendro::label(Vddata)
      # Re-order our matrix based on extracted dendrogram labels
      matr2 <- matr2[match(Vlabs$label, rownames(matr2)),]
      if (is.null(colnames(matr2))) { colnames(matr2) <- paste0("Sample", as.character(1L:Ncol)) }
    }
    
  }
  #
  if (split) {
    row_nms <- gsub(" - ", "\n", rownames(matr2))
    col_nms <- gsub(" - ", "\n", colnames(matr2))
    row_nms <- gsub(" vs ", "\nvs ", row_nms)
    col_nms <- gsub(" vs ", "\nvs ", col_nms)
    rownames(matr2) <- gsub(" VS ", "\nVS ", row_nms)
    colnames(matr2) <- gsub(" VS ", "\nVS ", col_nms)
  }
  RowNms <- data.frame(Name = rownames(matr2),
                       X = h.Marg[1L]+Ncol+1,
                       Y = Nrow:1)
  nChar <- max(nchar(RowNms$Name))
  ColNms <- data.frame(Name = colnames(matr2),
                       X = 1L:Ncol + h.Marg[1L],
                       Y = Nrow + 0.75)
  Vmax <- Nrow + 2L + v.Marg
  matr3 <- cbind(expand.grid(dimnames(matr2)),
                 value = as.vector(matr2)) # Complicated, but works nicely and avoids warnings
  colnames(matr3) <- c("Row", "Column", "value")
  lim <- ceiling(max(abs(is.all.good(as.numeric(matr3$value)))))
  matr3$value[which((!is.finite(matr3$value))&(matr3$value < 0))] <- -lim
  matr3$value[which((!is.finite(matr3$value))&(matr3$value > 0))] <- lim
  matr3$value[which(!is.all.good(matr3$value, mode = "logical"))] <- NA
  matr3$Y <- Nrow:1L
  matr3$X <- as.numeric(sapply(1L:Ncol, \(x) { rep(x, Nrow) })) + h.Marg[1L]
  ttlX <- (h.Marg[1L] + 1L)/2
  Title <- data.frame(Title = unlist(strsplit(as.character(title), "\n")))
  Title$X <- ttlX
  Title$Y <- Vmax*0.9 - (0L:(nrow(Title)-1L))*0.6
  # Fill scale
  fillScale <- viridis::scale_fill_viridis(option = colours, discrete = FALSE,
                                           direction = c(1L, -1L)[(colours_invert)+1L])
  # fillScale <- ggplot2::scale_fill_gradient2(low = colDF$hex[1L], mid = colDF$hex[2L], high = colDF$hex[3L],
  #                                           limits = c(min(colDF$RangeVal), max(colDF$RangeVal)),
  #                                           na.value = na.colour)
  #
  if (normRows) { colnames(matr3)[which(colnames(matr3) == "value")] <- "normValue" }
  plot <- ggplot2::ggplot(matr3)
  # Cells
  plot <- if (normRows) {
    plot +
      ggplot2::geom_rect(ggplot2::aes(xmin = X-0.5, xmax = X+0.5, ymin = Y-0.5, ymax = Y+0.5,
                                      fill = normValue))
  } else {
    plot +
      ggplot2::geom_rect(ggplot2::aes(xmin = X-0.5, xmax = X+0.5, ymin = Y-0.5, ymax = Y+0.5,
                                      fill = value))
  }
  plot <- plot +
    # Cells fill scale
    fillScale +
    # Row names
    ggplot2::geom_text(data = RowNms, ggplot2::aes(label = Name, x = X, y = Y),
              hjust = 0, cex = size) +
    # Column names
    ggplot2::geom_text(data = ColNms, ggplot2::aes(label = Name, x = X, y = Y),
              hjust = 0, angle = 60, cex = size) +
    # Title
    ggplot2::geom_text(data = Title, ggplot2::aes(label = Title, x = X, y = Y),
              hjust = 0.5, cex = size*1.5, fontface = "bold") +
    # Theme
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(), 
                   axis.ticks.y = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   axis.line = ggplot2::element_blank())
  if (!missing(subtitle)) {
    # Subtitle
    tmp <- unlist(strsplit(as.character(subtitle), "\n"))
    if (!length(tmp)) { tmp <- "" }
    subTitle <- data.frame(Title = tmp)
    subTitle$X <- ttlX
    subTitle$Y <- Vmax*0.9 - nrow(Title)*0.6 - (0L:(nrow(subTitle)-1L))*0.6
    plot <- plot +
      ggplot2::geom_text(data = subTitle, ggplot2::aes(label = Title, x = X, y = Y),
                         hjust = 0.5, cex = size*1.3)
  }
  if (print_values) {
    plot <- plot + ggplot2::geom_text(ggplot2::aes(label = value, x = X, y = Y), colour = "white", cex = 1.5)
  }
  plot <- plot +
    ggplot2::coord_fixed(ratio = 1L,
                         xlim = c(-10L, Ncol + sum(h.Marg) + 0.5 + nChar/5),
                         ylim = c(0L, Vmax),
                         expand = FALSE)
  #poplot(plot)
  # Dendrograms
  if (hTree) {
    # - Horizontal
    Hdendr <- Hddata$segments
    Hdendr$Type <- c("horizontal", "vertical")[(Hdendr$x == Hdendr$xend)+1L]
    xScl <- summary(c(Hdendr$x, Hdendr$xend))
    yScl <- summary(c(Hdendr$y, Hdendr$yend))
    XScl <- summary(ColNms$X)
    YScl <- summary(ColNms$Y)
    Hdendr$X <- (Hdendr$x-xScl["Min."])/(xScl["Max."]-xScl["Min."])
    Hdendr$X <- Hdendr$X*(XScl["Max."]-XScl["Min."])+XScl["Min."]
    Hdendr$Xend <- (Hdendr$xend-xScl["Min."])/(xScl["Max."]-xScl["Min."])
    Hdendr$Xend <- Hdendr$Xend*(XScl["Max."]-XScl["Min."])+XScl["Min."]
    Hdendr$Y <- (Hdendr$y-yScl["Min."])/(yScl["Max."]-yScl["Min."]) + Nrow + v.Padd + 1.5
    Hdendr$Yend <- (Hdendr$yend-yScl["Min."])/(yScl["Max."]-yScl["Min."]) + Nrow + v.Padd + 1.5
    plot <- plot +
      ggplot2::geom_segment(data = Hdendr, linewidth = 0.5,
                            ggplot2::aes(x = X, y = Y, xend = Xend, yend = Yend))
    #poplot(plot)
  }
  if (vTree) {
    # - Vertical
    Vdendr <- Vddata$segments
    mtch <- match(c("x", "y", "xend", "yend"), colnames(Vdendr))
    colnames(Vdendr)[mtch] <- c("y", "x", "yend",  "xend")
    Vdendr <- Vdendr[, c("x", "y", "xend", "yend")]
    Vdendr$Type <- c("horizontal", "vertical")[(Vdendr$x == Vdendr$xend)+1]
    xScl <- summary(c(Vdendr$x, Vdendr$xend))
    yScl <- summary(c(Vdendr$y, Vdendr$yend))
    XScl <- summary(RowNms$X)
    YScl <- summary(RowNms$Y)
    Vdendr$Y <- (Vdendr$y-yScl["Min."])/(yScl["Max."]-yScl["Min."])
    Vdendr$Y <- Vdendr$Y*(YScl["Max."]-YScl["Min."])+YScl["Min."]
    Vdendr$Yend <- (Vdendr$yend-yScl["Min."])/(yScl["Max."]-yScl["Min."])
    Vdendr$Yend <- Vdendr$Yend*(YScl["Max."]-YScl["Min."])+YScl["Min."]
    Vdendr$X <- h.Marg[1L] - 0.5 - h.Padd - (Vdendr$x-xScl["Min."])/(xScl["Max."]-xScl["Min."])
    Vdendr$Xend <- h.Marg[1L] - 0.5 - h.Padd - (Vdendr$xend-xScl["Min."])/(xScl["Max."]-xScl["Min."])
    # Symmetrize
    Vdendr$Y <- YScl["Max."]+1L-Vdendr$Y
    Vdendr$Yend <- YScl["Max."]+1L-Vdendr$Yend
    plot <- plot +
      ggplot2::geom_segment(data = Vdendr, linewidth = 0.5,
                            ggplot2::aes(x = X, y = Y, xend = Xend, yend = Yend))
    #poplot(plot)
  }
  # Print
  if (print) {
    poplot(plot)
  }
  # Save
  save <- unique(gsub("^jpg$", "jpeg", gsub("^\\.", "", tolower(save))))
  if ((length(save) > 1)||(save != "false")) {
    if (!dir.exists(folder)) { dir.create(folder, recursive = TRUE) }
    folder <- normalizePath(folder, winslash = "/")
    ttl <- title
    if ((!missing(subtitle))&&(nchar(subtitle))) { ttl <- paste0(ttl, " - ", subtitle) }
    ttl <- gsub("\n|/|:|\\*|\\?|<|>|\\|", "-", ttl)
    if (identical(save, "true")) { save <- "jpeg" }
    for (sv in save) {
      ttlSv <- paste0(ttl, ".", sv)
      fl <- paste0(folder, "/", ttlSv)
      #cat(fl, "\n")
      ggplot2::ggsave(fl, plot, dpi = 300L, width = 10L, height = 10L, units = "in")
    }
  }
  # Return
  if (return) {
    plot <- plotEval(plot)
    return(plot)
  }
}
