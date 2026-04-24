# Ranked abudance and profile plots at peptides level
# This chunk has been vastly improved, and the others should be improved on the same model!!! 
source(parSrc, local = FALSE)
clusterExport(parClust, "abbrFun", envir = environment())
plotPepProf %<o% runProfPlots # For now, should come under control of a parameter eventually
if (plotPepProf) {
  #
  # This can be best parallelized thusly:
  # - Create data.frame with one row per combination of sample / quant_type / way to plot
  # - Parallelize over those
  #
  pepQuantTypes %<o% "intensities"
  PltTst <- c("All", "List")[1L:(1L+("+" %in% pep$`In list`))]
  plotTypes <- c("jpeg", "pdf", "html")
  l1 <- length(Exp)
  l2 <- length(pepQuantTypes)
  l3 <- length(PltTst)
  l4 <- length(plotTypes)
  plotsDF <- data.frame(Exp = rep(Exp, l2*l3*l4))
  plotsDF$Type <- as.character(sapply(pepQuantTypes, \(x) { rep(x, l1) }))
  plotsDF$Mode <- as.character(sapply(PltTst, \(x) {  rep(x, l1*l2) }))
  plotsDF$Ext <- as.character(sapply(plotTypes, \(x) {  rep(x, l1*l2*l3) }))
  #
  tmp <- listMelt(strsplit(pep$Proteins, ";"), 1L:nrow(pep))
  m <- match(tmp$value, db$`Protein ID`)
  tmp$name <- unlist(db$`Common Name`[m])
  w <- which(!nchar(tmp$name) == 0L)
  tmp$name[w] <- do.call(paste, c(tmp[w, c("name", "value")], sep = " "))
  w <- which(nchar(tmp$name) == 0L)
  tmp$name[w] <- tmp$value[w]
  tmp3 <- aggregate(tmp$name, list(tmp$L1), paste, collapse = "\n") #This one is faster than with data.table!
  tmp3 <- data.frame(Seq = gsub("^_|_$", "", pep$`Modified sequence`),
                     Name = tmp3$x[match(1L:nrow(pep), tmp3$Group.1)])
  pep$Peptide_ID <- do.call(paste, c(tmp3, sep = "\n"))
  pep$"Peptide ID" <- gsub("\n", " ", pep$Peptide_ID)
  if (tstOrg) {
    tmp$Cont <- db$`Potential contaminant`[m]
    tmp$Org <- db[m, dbOrgKol]
    tmp1 <- aggregate(tmp$Cont, list(tmp$L1), \(x) {
      x <- if ("+" %in% x) { "+" } else { "" }
      return(x)
    })
    tmp2 <- aggregate(tmp$Org, list(tmp$L1), \(x) {
      x <- if ("Contaminant" %in% x) { "Contaminant" } else { paste(unique(x), collapse = ";") }
      return(x)
    })
    pep$`Potential contaminant` <- tmp1$x[match(1L:nrow(pep), tmp1$Group.1)]
    pep$temp <- tmp2$x[match(1L:nrow(pep), tmp2$Group.1)]
  }
  if (prot.list.Cond) {
    pep$"In list" <- ""
    pep$"In list"[grsep2(prot.list, pep$Proteins)] <- "+"
  }
  #
  plotsDF$kolnm <- plotsDF$root <- NA
  w <- which(plotsDF$Type == "intensities")
  plotsDF$kolnm[w] <- "Intensity"
  plotsDF$root[w] <- paste0(plotsDF$kolnm[w], " - ")
  w <- which(is.na(plotsDF$root))
  if (length(w)) { stop("This case is not addressed in the code yet!") }
  plotsDF$kol <- do.call(paste, c(plotsDF[, c("root", "Exp")], sep = ""))
  #
  # Prepare data
  temp_pep <- dfMelt(pep[, unique(plotsDF$kol), drop = FALSE])
  colnames(temp_pep)[which(colnames(temp_pep) == "value")] <- "Y"
  temp_pep$variable <- as.character(temp_pep$variable)
  varkol <- c("Sequence", "Proteins", "Peptide ID", "Peptide_ID", "id")
  if (prot.list.Cond) { varkol <- c(varkol, "In list") }
  temp_pep[, varkol] <- pep[, varkol]
  if (prot.list.Cond) {
    temp_pep$`In list`[which(temp_pep$`In list` == "")] <- "-"
    temp_pep$"In list" <- factor(temp_pep$"In list", levels = c("-", "+"))
  }
  temp_pep$Category <- if (lOrg) { pep$temp[match(temp_pep$id, pep$id)] } else { "-" }
  temp_pep <- temp_pep[which(temp_pep$Y > 0),]
  temp_pep <- temp_pep[which(is.all.good(temp_pep$Y, 2L)),]
  g <- grep("^Intensity - ", temp_pep$variable)
  temp_pep$Y[g] <- log10(temp_pep$Y[g])
  temp_pep$variable[g] <- gsub_Rep("^Intensity - ", "log10(Int.) - ", temp_pep$variable[g])
  g <- grep("^Intensity - ", plotsDF$kol)
  plotsDF$kol[g] <- gsub_Rep("^Intensity - ", "log10(Int.) - ", plotsDF$kol[g])
  plotsDF$root[g] <- "log10(Int.) - "
  temp_pep$Sample <- gsub_Rep(".* - ", "", temp_pep$variable)
  lev <- c("In list", "+", tstOrg2, "-", "Contaminant")
  lev <- lev[which(lev %in% unique(temp_pep$Category))]
  temp_pep$Category <- factor(temp_pep$Category, levels = lev)
  #
  # Serialize data
  tmp <- aggregate(1L:nrow(temp_pep), list(temp_pep$Sample), list)
  sapply(Exp, \(exp) { #exp <- Exp[1L]
    readr::write_rds(temp_pep[tmp$x[[which(tmp$Group.1 == exp)]],], paste0(wd, "/tmp_", exp, ".RDS"))
  })
  #
  MainDir <- paste0(wd, "/Ranked abundance")
  #
  # Scales - re-using from the PG-level source
  #myColorsB <- myColors <- setNames("black", "-")
  #myColors2 <- setNames(c("lightgrey", "brown"), c("-", "+"))
  #lOrg <- length(tstOrg2)
  #if (lOrg) {
    #myColorsB <- myColors <- setNames(grDevices::colorRampPalette(c("blue", "green"))(lOrg), tstOrg2)
    #w <- which(names(myColorsB) != "In list")
    #names(myColorsB)[w] <- gsub("[a-z]+ ", ". ", names(myColorsB)[w])
    #names(myColorsB) <- gsub(";", "\n", names(myColorsB))
    #myColorsB[["potential contaminant"]] <- myColors[["potential contaminant"]] <- "grey"
  #}
  #myColorsB[[c("+", "In list")[(lOrg > 0L)+1L]]] <- myColors[[c("+", "In list")[(lOrg > 0L)+1L]]] <- "red"
  #colScale <- scale_colour_manual(name = "Category", values = myColors)
  #colScaleB <- scale_colour_manual(name = "Category", values = myColorsB)
  #fillScale <- scale_fill_manual(name = "Category", values = myColors)
  #colScale2 <- scale_colour_manual(name = "In list", values = myColors2)
  #fillScale2 <- scale_fill_manual(name = "In list", values = myColors2)
  #
  clusterExport(parClust, list("plotsDF", "tstOrg2", "wd", "MainDir", #"myColorsB",
                               "myColors2", "colScale", #"colScaleB",
                               "colScale2", "fillScale", "fillScale2"), envir = environment())
  cat("Plotting peptide ranked abundance plots\n")
  invisible(clusterCall(parClust, \() {
    library(ggplot2)
    assign("lOrg", length(tstOrg2), .GlobalEnv)
    return()
  }))
  pepSortTst <- parLapply(parClust, 1L:nrow(plotsDF), \(ii) { #ii <- 1L
    exp <- plotsDF$Exp[ii]
    pepQuantType <- plotsDF$Type[ii]
    plotMode <- plotsDF$Mode[ii]
    kolnm <- plotsDF$kolnm[ii]
    root <- plotsDF$root[ii]
    fileType <- plotsDF$Ext[ii]
    #
    SubDir <- paste0(MainDir, "/Peptide ", pepQuantType)
    if (!dir.exists(SubDir)) { dir.create(SubDir, recursive = TRUE) }
    #
    # Load data
    temp <- readr::read_rds(paste0(wd, "/tmp_", exp, ".RDS"))
    nr <- nrow(temp)
    test <- aggregate(temp$Peptide_ID, list(temp$Peptide_ID), length)
    w <- which(test$x > 1L)
    if (length(w)) {
      test <- test[w,]
      for (i in test$Group.1) {
        w <- which(temp$Peptide_ID == i)
        temp$Peptide_ID[w[2L:length(w)]] <- paste0(temp$Peptide_ID[w[2L:length(w)]], "_", 2L:length(w))
        temp$"Peptide ID"[w[2L:length(w)]] <- paste0(temp$"Peptide ID"[w[2L:length(w)]], "_", 2L:length(w))
      }
    }
    #
    temp$Value <- paste0(kolnm, ": ", temp$Y)
    temp <- temp[order(temp$Y, decreasing = TRUE),]
    temp$Peptide_ID <- factor(temp$Peptide_ID, levels = unique(temp$Peptide_ID))
    temp$"Peptide ID" <- factor(temp$"Peptide ID", levels = unique(temp$"Peptide ID"))
    intmin <- floor(min(temp$Y))
    intmax <- ceiling(max(temp$Y))
    intscale <- intmax-intmin
    xmax <- nr*18/15
    ttl <- paste0("Peptide ", gsub(" - $", "", root), " ranked abundance plots - ", exp)
    if (plotMode == "All") {
      catnm <- "Category"
      filt <- 1L:nr
    }
    if (plotMode == "List") {
      ttl <- paste0(ttl, ", proteins of interest")
      catnm <- "In list"
      filt <- which(temp[[catnm]] == "+")
    }
    temp$y <- temp$Y + intmax*0.01
    flPath <- paste0(SubDir, "/", ttl, ".", fileType)
    plot <- ggplot(temp) +
      geom_bar(stat = "identity",
               aes(x = Peptide_ID, y = Y, fill = .data[[catnm]], text = Value)) +
      annotate("text", nr/2, intmax*1.2+intscale*0.04, label = paste0(length(unique(temp$id)), " peptidoforms"),
               hjust = 0.5) + theme_bw() + ggtitle(ttl) + ylab(kolnm) + 
      coord_cartesian(xlim = c(1, xmax),
                      ylim = c(intmin, intmax*1.2+intscale*0.05)) +
      theme(panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_blank(),
            axis.ticks = element_blank(),
            plot.margin = margin(r = 100L))
    if (plotMode == "All") { plot <- plot + colScale + fillScale }
    if (plotMode == "List") { plot <- plot + colScale2 + fillScale2 }
    if (fileType %in% c("jpeg", "pdf")) {
      plot <- plot +
        geom_text(data = temp[filt,], angle = 45, hjust = 0, cex = 3.5, show.legend = FALSE,
                  aes(Peptide_ID, y, colour = .data[[catnm]], label = `Peptide ID`))
      #poplot(plot, 12L, 22L)
      suppressMessages({
        ggsave(filename = flPath, plot, dpi = 600L, width = 30L, height = 10L, units = "in")
      })
    }
    if (fileType == "html") {
      plot_ly <- plotly::ggplotly(plot, tooltip = c("Peptide_ID", "text"))
      setwd(SubDir) # For some reason, unless I do this the default selfcontained = TRUE argument gets ignored and
      # a folder with external resources is created for each html plot!
      tst <- try(htmlwidgets::saveWidget(plotly::partial_bundle(plot_ly), flPath), silent = TRUE)
      if (inherits(tst, "try-error")) { tst <- try(htmlwidgets::saveWidget(plot_ly, flPath), silent = TRUE) }
      setwd(wd)
    }
  })
  unlink(paste0(wd, "/tmp_", Exp, ".RDS"))
  if (length(Exp) > 1L) {
    plotsDF2 <- plotsDF[which((plotsDF$Ext != "html")|(plotsDF$Mode == "List")),]
    # We only draw the html version if we are drawing for proteins in the list
    plotsDF2$Exp <- NULL
    plotsDF2$kol <- NULL
    tst <- do.call(paste, c(plotsDF2, sep = "_______"))
    tst <- aggregate(1L:nrow(plotsDF2), list(tst), min)
    plotsDF2 <- plotsDF2[tst$x,]
    # temp_pep is all the data we need
    MainDir <- paste0(wd, "/Profile plots")
    clusterExport(parClust, list("plotsDF2", "MainDir", "tstOrg2"), envir = environment())
    readr::write_rds(temp_pep, paste0(wd, "/tmp.RDS"))
    cat("Plotting peptide profile plots\n")
    pepProfTst <- parLapply(parClust, 1L:nrow(plotsDF2), \(ii) { #ii <- 1L
      pepQuantType <- plotsDF2$Type[ii]
      plotMode <- plotsDF2$Mode[ii]
      kolnm <- plotsDF2$kolnm[ii]
      root <- plotsDF2$root[ii]
      fileType <- plotsDF2$Ext[ii]
      #
      temp <- readr::read_rds(paste0(wd, "/tmp.RDS"))
      temp$Value <- temp$Y
      #
      SubDir <- paste0(MainDir, "/Peptide ", pepQuantType)
      if (!dir.exists(SubDir)) { dir.create(SubDir, recursive = TRUE) }
      #
      ttl <- paste0("Peptides ", gsub(" - $", "", root), " profiles")
      if (plotMode == "All") {
        lev <- levels(temp$Category)
        w <- which(!lev %in% c("In list", "Contaminant"))
        lev <- c(c("In list", "Contaminant"), lev[w])
        w <- which(!lev %in% c("In list", "Contaminant"))
        lev2 <- lev
        lev2[w] <- sapply(strsplit(as.character(lev[w]), ""), abbrFun)
        lev2[which(lev2 == "Contaminant")] <- "Cont."
        lev2 <- gsub(";", "\n", lev2)
        temp$Category <- as.character(temp$Category)
        temp$Category <- lev2[match(temp$Category, lev)]
        temp$Category <- factor(temp$Category, levels = lev2)
        catnm <- "Category"
        filt <- 1L:nrow(temp)
      }
      if (plotMode == "List") {
        ttl <- paste0(ttl, ", proteins of interest")
        catnm <- "In list"
        filt <- which(temp[[catnm]] == "+")
      }
      m <- match(temp$Category, c("-", "Contaminant", tstOrg2, "+", "In list"))
      temp$LineType <- c(rep("dotted", 2L), rep("dashed", lOrg), rep("solid", 2L))[m]
      #temp$DotSize <- c(rep(0.2, 2L), rep(0.3, lOrg), rep(0.5, 2L))[m]
      temp$DotSize <- 1L
      #temp$Alpha <- c(rep(0.25, 2L), rep(0.5, lOrg), rep(1L, 2L))[m]
      temp$Alpha <- 1L
      Ngl <- c(0, 90)[(length(levels(temp[[catnm]])) > 5L) + 1L]
      wTxt <- which(temp$Sample == rev(levels(temp$Sample))[1L])
      frm <- as.formula(paste0("~`", catnm, "`"))
      flPath <- paste0(SubDir, "/", ttl, ".", fileType)
      if (fileType %in% c("jpeg", "pdf")) {
        dat <- temp
      }
      if (fileType == "html") {
        dat <- temp[which((!is.na(temp$`In list`))&(temp$`In list` == "+")),]
      }
      plot <- ggplot(dat, aes(text1 = Peptide_ID, text2 = Value)) +
        geom_line(aes(x = Sample, y = Y, group = id, color = id), alpha = 0.1, show.legend = FALSE) +
        geom_point(aes(x = Sample, y = Y, size = DotSize, color = id)) +
        ggtitle(ttl) + ylab(kolnm) +
        theme_bw() + facet_grid(frm) +
        scale_size_identity(guide = "none") +
        scale_linetype_identity(guide = "none") +
        theme(legend.position = "none",
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
              strip.text = element_text(face = "bold", size = 8, lineheight = 0.8, angle = Ngl),
              strip.background = element_rect(fill = "lightblue", colour = "black", linewidth = 1L)) +
        scale_alpha_identity(guide = "none") +
        viridis::scale_color_viridis(option = "D")
      #poplot(plot, 12L, 22L)
      if (fileType %in% c("jpeg", "pdf")) {
        plot <- plot +
          geom_text(data = dat[wTxt,], aes(label = `Peptide ID`, x = Sample, y = Y, alpha = Alpha, color = id),
                    hjust = 0, cex = 2L)
        #poplot(plot, 12L, 22L)
        suppressMessages({
          ggsave(filename = flPath, plot, dpi = 600L, width = 30L, height = 10L, units = "in")
        })
      }
      if (fileType == "html") {
        plot_ly <- plotly::ggplotly(plot, tooltip = c("Peptide_ID", "text")) # Super slowwwwww
        setwd(SubDir) # For some reason, unless I do this the default selfcontained = TRUE argument gets ignored and
        # a folder with external resources is created for each html plot!
        tst <- try(htmlwidgets::saveWidget(plotly::partial_bundle(plot_ly), flPath), silent = TRUE)
        if (inherits(tst, "try-error")) { tst <- try(htmlwidgets::saveWidget(plot_ly, flPath), silent = TRUE) }
        #system(paste0("open \"", flPath, "\""))
        setwd(wd)
      }
    })
  }
  unlink(paste0(wd, "/tmp.RDS"))
  pep$temp <- NULL
}
invisible(clusterCall(parClust, \(x) { rm(list = ls());gc() }))
