# Re-normalize peptide ratios
# Legacy, currently unused.
# If planning to re-use it, first update the main shiny app for parameters,
# as it is currently not dealing with the parameters needed here.
#
stop("Currently this is not supported! The scripts exist but need a revision as the code hasn't been checked for long... Also, does it even make sense to do this?")
# # Step 1:
# # Classic normalisation by the median for each column;
# # NB: This does not work the same way as intensities, so here the global scale of ratios need not be conserved:
# # We really want them centred on 0 in log scale.
# Norm.log2.Pep.Ratios <- c()
# NormGrps <- unique(pep$`Normalisation group`)
# g <- grep(topattern(pep.ratios.ref[1L]), colnames(pep), value = TRUE)
# pep[, paste0("norm. ", g)] <- NA
# for (nrmgrp in NormGrps) {
#   w <- which(pep$`Normalisation group` == nrmgrp)
#   for (k in g) {
#     m <- median(is.all.good(pep[w, k]))
#     #m <- mlv(is.all.good(pep[[k]]), method = "Parzen")[1L]
#     pep[w, paste0("norm. ", k)] <- pep[w, k]-m
#     Norm.log2.Pep.Ratios[paste0(nrmgrp, " - ", gsub(topattern(pep.ratios.ref[1L]), "", k))] <- m
#   }
# }
# pep.ratios.ref <- unique(c(pep.ratios.ref, paste0("norm. ", pep.ratios.ref[1L])))
# # Step 2
# # Optional advanced normalisation:
# if (("Adv.Norma.Pep.Ratio" %in% colnames(Param))&&(Param$Adv.Norma.Pep.Ratio != FALSE)) {
#   if (Param$Adv.Norma.Pep.Ratio.Type == "C") {
#     k <- Adv.Norma.Pep.Ratio.Type.Group$column
#     test <- vapply(Adv.Norma.Pep.Ratio.Type.Group$values, \(i) { #i <- Adv.Norma.Pep.Ratio.Type.Group$values[1L]
#       i <- Exp.map$Ref.Sample.Aggregate[which(Exp.map[[k]] == i)]
#       return(length(which(paste0(pep.ratios.ref[length(pep.ratios.ref)], i) %in% colnames(pep))))
#     }, 1L)
#     agg <- Adv.Norma.Pep.Ratio.Type.Group$values[which(test > 1L)]
#     exports <- list("agg", "Adv.Norma.Pep.Ratio.Type.Group", "Exp.map", "pep.ratios.ref", "pep", "Param")
#     clusterExport(parClust, exports, envir = environment())
#     invisible(clusterCall(parClust, \() {
#       library(proteoCraft)
#       return()
#     }))
#     norm_temp <- parSapply(parClust, 0:length(agg), \(i) { #i <- 1L
#       if (i == 0) {
#         kol <- grep(paste0(topattern(pep.ratios.ref[1L]), ".+_REF\\.to\\.REF_"), colnames(pep), value = TRUE)
#       } else {
#         x <- agg[i]
#         j <- setNames(unlist(strsplit(x, "___")),
#                       Adv.Norma.Pep.Ratio.Type.Group$names)
#         temp <- lapply(Adv.Norma.Pep.Ratio.Type.Group$names, \(x) {
#           if (j[[x]] == "NA") {
#             return(which((is.na(Exp.map[[x]]))|(Exp.map[[x]] == j[[x]])))
#           }
#           return(which(Exp.map[[x]] == j[[x]]))
#         })
#         temp2 <- sort(unique(unlist(temp)))
#         test <- vapply(temp2, \(x) { sum(vapply(temp, \(y) { x %in% unlist(y) }, TRUE)) }, 1L)
#         temp2 <- temp2[which(test == length(temp))]
#         temp3 <- Exp.map[temp2,]
#         kol <- paste0(pep.ratios.ref[length(pep.ratios.ref)], temp3$Ref.Sample.Aggregate)
#         w <- which(kol %in% colnames(pep))
#         temp3 <- temp3[w,]
#         kol <- kol[w]
#         if (length(w) <= 1L) {
#           warning(paste0("There are no columns to normalize for samples group ", x))
#         }
#       }
#       if (length(kol) > 1L) {
#         temp2 <- temp <- pep[, c("Modified sequence", kol)]
#         colnames(temp2) <- gsub("^AdvNorm\\.norm\\.", "AdvNorm. ", colnames(temp2))
#         for (nrmgrp in unique(pep$`Normalisation group`)) {
#           w <- which(pep$`Normalisation group` == nrmgrp)
#           tmp <- AdvNorm.IL(temp[w,], "Modified sequence", kol, TRUE, 5L)
#           colnames(tmp) <- gsub("^AdvNorm\\.norm\\.", "AdvNorm. ", colnames(tmp))
#           temp2[w, colnames(tmp)] <- tmp[, colnames(tmp)]
#         }
#         #cat("Advanced peptides ratio normalisation done for samples group", x, "\n")
#         temp2$"Modified sequence" <- NULL
#         return(list(temp2))
#       }
#     })
#   } else {
#     stop("Not implemented yet, I need to update the current AdvNorm function as it takes too long or crashes.")
#   }
#   for (i in seq_along(norm_temp)) {
#     tmp <- norm_temp[[i]]
#     pep[, colnames(tmp)] <- tmp
#   }
#   pep.ratios.ref <- unique(c(pep.ratios.ref, paste0("AdvNorm. ", pep.ratios.ref[1L])))
# }
# if (Param$Norma.Pep.Ratio.show) {
#   for (i in Ratios.Plot.split$values) { #i <- Ratios.Plot.split$values[1L]
#     j <- setNames(unlist(strsplit(i, "___")),
#                   unlist(Aggregate.map$Characteristics[which(Aggregate.map$Aggregate.Name == Ratios.Plot.split$aggregate)]))
#     k <- lapply(names(j), \(x) {
#       if (j[[x]] == "NA") {
#         return(which((is.na(Exp.map[[x]]))|(Exp.map[[x]] == j[[x]])))
#       }
#       return(which(Exp.map[[x]] == j[[x]]))
#     })
#     l <- sort(unique(unlist(k)))
#     test <- vapply(l, \(x) { sum(vapply(k, \(y) { x %in% y }, TRUE)) == length(k) }, 1L)
#     temp <- Exp.map$Ref.Sample.Aggregate[l[which(test)]]
#     a1 <- paste0(pep.ratios.ref[1L], temp)
#     a2 <- paste0(pep.ratios.ref[length(pep.ratios.ref)], temp)
#     a1 <- a1[which(a1 %in% colnames(pep))]
#     a2 <- a2[which(a2 %in% colnames(pep))]
#     if (length(a1)) {
#       temp <- pep[, c("Modified sequence", a1, a2)]
#       temp <- reshape2::melt(temp)
#       temp$Norm <- grepl(topattern(pep.ratios.ref[length(pep.ratios.ref)]), temp$variable)
#       temp$Norm <- c("Original", "Normalised")[temp$Norm + 1L]
#       temp$Norm <- factor(temp$Norm, levels = c("Original", "Normalised"))
#       temp$variable <- as.character(temp$variable)
#       temp$variable[which(temp$Norm == "Original")] <- gsub_Rep(topattern(pep.ratios.ref[1L]), "", temp$variable[which(temp$Norm == "Original")])
#       temp$variable[which(temp$Norm == "Normalised")] <- gsub_Rep(topattern(pep.ratios.ref[length(pep.ratios.ref)]), "", temp$variable[which(temp$Norm == "Normalised")])
#       temp2 <- Isapply(strsplit(temp$variable, "___"), unlist)
#       colnames(temp2) <- unlist(Aggregate.map$Characteristics[which(Aggregate.map$Aggregate.Name == RSA$aggregate)])
#       temp[, colnames(temp2)] <- temp2
#       if (length(Ratios.Plot.wrap$names) > 1L) {
#         temp$Wrap <- do.call(paste, c(temp2[, Ratios.Plot.wrap$names], sep = "_"))
#       } else { temp$Wrap <- temp2[[Ratios.Plot.wrap$names]] }
#       if (length(Ratios.Plot.colour$names) > 1L) {
#         temp$Colour <- do.call(paste, c(temp2[, Ratios.Plot.colour$names], sep = "_"))
#       } else { temp$Colour <- temp2[[Ratios.Plot.colour$names]] }
#       temp$X <- do.call(paste, c(temp[, c("Norm", "Colour")], sep = "_"))
#       temp$X <- factor(temp$X, levels = unique(unlist(vapply(c("Original", "Normalised"), \(x) {
#         paste(x, unique(temp2[[Ratios.Plot.colour$names]]), sep = "_")
#       }, ""))))
#       temp <- temp[which(is.all.good(temp$value, 2L)),]
#       dir <- paste0(wd, "/Workflow control/Peptides/Ratios")
#       if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
#       dirlist <- unique(c(dirlist, dir))
#       ttl <- paste0("Peptide ratios distribution_sample group: ", i)
#       plot <- ggplot(temp) +
#         geom_violin(aes(x = X, y = value, color = Colour, fill = Colour), alpha = 0.25) +
#         geom_boxplot(aes(x = X, y = value, color = Colour, fill = Colour), alpha = 0.5) +
#         theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
#         scale_color_viridis_d(begin = 0.25) +
#         scale_fill_viridis_d(begin = 0.25) +
#         facet_grid(. ~ Norm, scales = "free", space = "free") +
#         ggtitle(ttl)
#       if (length(unique(temp$Wrap)) > 1) { plot <- plot + facet_wrap(~Wrap) }
#       print(plot) # This type of QC plot does not need to pop up, the side panel is fine
#       suppressMessages({
#         ggsave(paste0(dir, "/", gsub(":", "_", ttl), ".jpeg"), plot,
#                dpi = 300, width = 10, height = 10, units = "in")
#         ggsave(paste0(dir, "/", gsub(":", "_", ttl), ".pdf"), plot,
#                dpi = 300, width = 10, height = 10, units = "in")
#       })
#       ReportCalls <- AddPlot2Report(Title = ttl)
#     } else { warning(paste0("Nothing to plot for level ", i)) }
#   }
#   # Also look at Ref-to-Ref ratios:
#   a1 <- grep(paste0(topattern(pep.ratios.ref[1L]), ".+_REF.to.REF_[0-9]+$"), colnames(pep), value = TRUE)
#   a2 <- grep(paste0(topattern(pep.ratios.ref[length(pep.ratios.ref)]), ".+_REF.to.REF_[0-9]+$"), colnames(pep), value = TRUE)
#   if (length(a1)) {
#     temp <- pep[, c("Modified sequence", a1, a2)]
#     temp <- reshape2::melt(temp)
#     temp$Norm <- grepl(topattern(pep.ratios.ref[length(pep.ratios.ref)]), temp$variable)
#     temp$Norm <- c("Original", "Normalised")[temp$Norm + 1L]
#     temp$Norm <- factor(temp$Norm, levels = c("Original", "Normalised"))
#     temp$variable <- as.character(temp$variable)
#     temp$variable[which(temp$Norm == "Original")] <- gsub_Rep(topattern(pep.ratios.ref[1L]), "", temp$variable[which(temp$Norm == "Original")])
#     temp$variable[which(temp$Norm == "Normalised")] <- gsub_Rep(topattern(pep.ratios.ref[length(pep.ratios.ref)]), "", temp$variable[which(temp$Norm == "Normalised")])
#     temp$Ratios.Group <- gsub_Rep("_REF.to.REF_[0-9]+$", "", temp$variable)
#     temp <- temp[which(is.all.good(temp$value, 2L)),]
#     dir <- paste0(wd, "/Workflow control/Peptides/Ratios")
#     if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
#     dirlist <- unique(c(dirlist, dir))
#     ttl <- "Peptide ratios distribution_References-to-References"
#     plot <- ggplot(temp) +
#       geom_violin(aes(x = variable, y = value, color = variable, fill = variable), alpha = 0.25) +
#       geom_boxplot(aes(x = variable, y = value, color = variable, fill = variable), alpha = 0.5) +
#       scale_color_viridis_d(begin = 0.25) +
#       scale_fill_viridis_d(begin = 0.25) +
#       theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
#       facet_wrap(~Ratios.Group, scales = "free") +
#       ggtitle(ttl)
#     print(plot) # This type of QC plot does not need to pop up, the side panel is fine
#     suppressMessages({
#       ggsave(paste0(dir, "/", gsub(":", "_", ttl), ".jpeg"), plot,
#              dpi = 300L, width = 10L, height = 10L, units = "in")
#       ggsave(paste0(dir, "/", gsub(":", "_", ttl), ".pdf"), plot,
#              dpi = 300L, width = 10L, height = 10L, units = "in")
#     })
#     ReportCalls <- AddPlot2Report(Title = ttl)
#     l <- length(DatAnalysisTxt)
#     DatAnalysisTxt[l] <- paste0(DatAnalysisTxt[l], " Peptide ratios were then re-normalized.")
#   } else { warning("Nothing to plot for Reference-to-Reference ratios!") }
# }

