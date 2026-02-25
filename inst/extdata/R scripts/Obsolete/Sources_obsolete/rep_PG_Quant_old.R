# Currently 9 methods are implemented, but the last 3 are just for testing and do not create all columns required for the script to complete
# (These could be added relatively easily if necessary though)
QuantData <- setNames(paste0("quant.data", seq_along(QuantMethods)), QuantMethods)
QuantMethods_all <- FALSE
.obj <- unique(c("QuantData", "QuantMethods_all", .obj))
exprsCol <- paste0("log10(Expr.) - ", RSA$values)
# Weights:
# - Higher for peptides with low intra-sample group CV on average
# - Higher for peptides with low PEP
if ((length(inDirs) == 1)&&(QuantUMS)) {
  weightsInsrt <- "mean DiaNN \"Quantity Quality\"" 
  pep$Weights <- pep$"Quantity Quality"
} else {
  weightsInsrt <- "-log(PEP)/CV" 
  source(parSrc, local = FALSE)
  Kols <- paste0(pep.ref[length(pep.ref)], Exp.map$Ref.Sample.Aggregate)
  Kols <- Kols[which(Kols %in% colnames(pep))]
  tmp <- pep[, Kols]
  clusterExport(parClust, list("Exp.map", "VPAL", "pep.ref", "is.all.good", "tmp"), envir = environment())
  CV <- parSapply(parClust, VPAL$values, function(x) { #x <- VPAL$values[1]
    smpls <- Exp.map$Ref.Sample.Aggregate[which(Exp.map[[VPAL$column]] == x)]
    kols <- paste0(pep.ref[length(pep.ref)], smpls)
    kols <- kols[which(kols %in% colnames(tmp))]
    x <- apply(tmp[, kols, drop = FALSE], 1, function(y) {
      y <- is.all.good(log10(unlist(y)))
      if (length(y)) {
        y <- sd(y)/mean(y)
      } else { y <- NA }
      return(y)
    })
    return(x)
  })
  CV <- rowMeans(CV, na.rm = TRUE)
  pep$Weights <- -log10(pep$PEP)/CV
}
summary(pep$Weights)
m <- max(is.all.good(pep$Weights))
pep$Weights <- pep$Weights/m
pep$Weights[which((is.na(pep$Weights))|(pep$Weights < 0.001))] <- 0.001
if ("Prot.Quant.Mod.Excl.is.strict" %in% colnames(Param)) {
  Mod.Excl.is.strict <- Param$Prot.Quant.Mod.Excl.is.strict
  if (!Mod.Excl.is.strict %in% c(1, 0, TRUE, FALSE)) { Mod.Excl.is.strict <- FALSE }
} else { Mod.Excl.is.strict <- FALSE }
Discard.unmod <- Mod.Excl.is.strict+1
if (Discard.unmod == 1) { Discard.unmod <- as.logical(Discard.unmod) }
.obj <- unique(c("Mod.Excl.is.strict", "Discard.unmod", .obj))
if (!grepl("^Prot\\.Quant", Param$QuantMeth)) {
  stop("NB: currently only methods 1 to 6 provide all necessary columns, not just quantitative columns per sample but also reference columns, ratios, etc... Until those are added they cannot be used by this script. Defaulting to method 3!")
  Param$QuantMeth <- "Prot.Quant.Unique"
}
if ((Param$QuantMeth %in% c("Prot.Quant", "Prot.Quant + weights", "Prot.Quant.Unique", "Prot.Quant.Unique + weights"))||(QuantMethods_all)) {
  # Calculates individual protein group expression and ratio values per sample, excluding some modifications and
  # (usually) requiring at least 2 peptides with different modified sequence.
  # This code has on occasion failed for very large datasets with the following error:
  # "Error in serialize(data, node$con) : error writing to connection"
  # Admittedly, when I write code, I tend to favour complexity over speed.
  # To avoid this issue, I added a bit of code to estimate how many threads to use:
  # - Size of clusters should not exceed N-1 threads (N = the computer's number of vCPUs).
  # - Somehow the LFQ function seems to require less than 7x the size of the pep & PG data frames in RAM
  #   per vCPU, regardless of number of threads.
  #   This is true even after I removed useless variable duplicates from the code!
  #   This is at least true for one project (ADAPTED project E), not sure yet about other datasets.
  #
  # This function is adapted from: https://stackoverflow.com/questions/27788968/how-would-one-check-the-system-memory-available-using-r-on-a-windows-machine
  # The sections for Linux and Mac-OS are untested, guesses based on: https://apple.stackexchange.com/questions/4286/is-there-a-mac-os-x-terminal-version-of-the-free-command-in-linux-systems
  # They are highly likely to need some adjustments before they work, but I have no way to test them!!!
  get_free_ram <- function() {
    sysnm <- tolower(Sys.info()[["sysname"]])
    if (sysnm == "windows") {
      x <- try(system2("wmic", args =  "OS get FreePhysicalMemory /value", stdout = TRUE), silent = TRUE)
      if (!"try-error" %in% class(x)) {
        x <- as.numeric(gsub("^FreePhysicalMemory=|\r$", "", grep("FreePhysicalMemory", x, value = TRUE)))/1000000
      } else {
        x <- try(system2("systeminfo.exe", stdout = TRUE), silent = TRUE)
        if (!"try-error" %in% class(x)) {
          x <- grep("^Total Physical Memory:     ", x, value = TRUE)
          x <- gsub("^Total Physical Memory:     |,", "", x)
          if (grepl(" B$", x)) { x <- as.numeric(gsub(" KB$", "", x))/1024^3 }
          if (grepl(" KB$", x)) { x <- as.numeric(gsub(" KB$", "", x))/1024^2 }
          if (grepl(" MB$", x)) { x <- as.numeric(gsub(" MB$", "", x))/1024 } # The only case I expect
          if (grepl(" GB$", x)) { x <- as.numeric(gsub(" GB$", "", x)) }
          if (grepl(" TB$", x)) { x <- as.numeric(gsub(" TB$", "", x))*1024 }
        } else {
          x <- try(system2("powershell", args = c("Get-WmiObject", "-Class", "Win32_ComputerSystem"),
                           stdout = TRUE), silent = TRUE)
          if (!"try-error" %in% class(x)) {
            x <- grep("^TotalPhysicalMemory : ", x, value = TRUE)
            x <- as.numeric(gsub("^TotalPhysicalMemory : ", "", x))/1024^3
          } else {
            # (If I really cannot detect RAM by any method, I will have to trust in the Omnissiah that I have enough...)
            x <- 10000
          }
        }
      }
    }
    if (sysnm %in% c("mac", "macos")) {
      warning("This code has not been tested and is likely to fail! Check it the first time you run on a Mac-OS system!")
      x <- try(system2("$ vm_stat", stdout = TRUE), silent = TRUE)
      if (!"try-error" %in% class(x)) {
        ps <- as.integer(gsub(".+\\(page size of | bytes\\)", "",x[1]))
        x <- as.numeric(gsub("^Pages free: +|\\.$", "", x[2]))*ps/1000000
      } else {
        x <- 10000 # If I cannot detect RAM, I will have to trust in the Omnissiah that I have enough...
      }
    }
    if (sysnm == "linux") {
      warning("This code has not been tested and is likely to fail! Check it the first time you run on a Linux system!")
      x <- system2("$free", stdout = TRUE)
      x <- as.integer(unlist(strsplit(x[1], " +"))[4])/1000000
    }
    if (!sysnm %in% c("windows", "mac", "macos", "linux")) { stop("Unknown computer system type!") }
    return(x)
  }
  free_ram <- get_free_ram()
  mem_per_thread <- as.numeric(gsub(" Gb$", "", (format(object.size(pep) + object.size(PG), "Gb"))))*3
  N.clust2 <- N.clust
  while (N.clust2*mem_per_thread > free_ram) { N.clust2 <- N.clust2-1 }
  tmpMod2Xclud <- Modifs$`Full name`[match(Mod2Xclud$Mark, Modifs$Mark)]
  if ((Param$QuantMeth == "Prot.Quant")||(QuantMethods_all)) {
    # Classic profiles method, no weights
    DatAnalysisTxt <- gsub("\\.$",
                           paste0(", and quantified using an in-house algorithm which: ",
                                  "i) computes a mean protein-level profile across samples using individual, normalized peptidoform profiles (\"relative quantitation\" step), ",
                                  "ii) following the best-flyer hypothesis, normalizes this profile to the mean intensity level of the most intense peptidoform (\"unscaled absolute quantitation\" step); ",
                                  "only unique and razor peptidoforms were used; ",
                                  tmpMod2Xclud, " peptidoforms and their unmodified counterparts were excluded from the calculations.",
                                  collapse = " "), DatAnalysisTxt)
    source(parSrc, local = FALSE)
    quant.data1 <- Prot.Quant(Prot = PG, Peptide.IDs = Pep4Quant, Pep = pep[Pep2Use,], id = "id",
                              Summary.method = "mean", #Summary.weights = "Weights",
                              Intensity.weights = FALSE,
                              experiments.map = Exp.map, param = Param,
                              Pep.Intens.root = pep.ref[length(pep.ref)],
                              Pep.Ratios.root = pep.ratios.ref[length(pep.ratios.ref)],
                              log.Pep.Intens = FALSE, log.Pep.Ratios = 2,
                              Prot.LFQ.to.log = TRUE, Prot.Ratios.to.log = TRUE,
                              Mods = Mod4Quant, Mods.to.Exclude = Mod2Xclud,
                              Discard.unmod = Discard.unmod,
                              Min.N = 1, Priority = c("int", "rat")[(LabelType %in% c("SILAC"))+1],
                              Refs_Mode = RefRat_Mode,
                              cl = parClust)
    #colnames(quant.data1)
    saveFun(quant.data1, file = "quant.data1.RData")
    #loadFun("quant.data1.RData")
  }
  if ((Param$QuantMeth == "Prot.Quant + weights")||(QuantMethods_all)) {
    # Classic profiles method, weights
    DatAnalysisTxt <- gsub("\\.$",
                           paste0(", and quantified using an in-house algorithm which: ",
                                  "i) computes a protein group-level profile across samples as the weighted mean of individual, normalized peptidoform profiles (\"relative quantitation\" step, weights = ",
                                  weightsInsrt,
                                  ", where PEP and CV = the peptidoform's Posterior Error Probability and Coefficient of Variation, resp.), ",
                                  "ii) following the best-flyer hypothesis, normalizes this profile to the mean intensity level of the most intense peptidoform (\"unscaled absolute quantitation\" step); ",
                                  "only unique and razor peptidoforms were used; ",
                                  tmpMod2Xclud, "peptidoforms and their unmodified counterparts were excluded from the calculations.",
                                  collapse = " "), DatAnalysisTxt)
    source(parSrc, local = FALSE)
    quant.data2 <- Prot.Quant(Prot = PG, Peptide.IDs = Pep4Quant, Pep = pep[Pep2Use,], id = "id",
                              Summary.method = "weighted.mean", Summary.weights = "Weights",
                              Intensity.weights = FALSE,
                              experiments.map = Exp.map, param = Param,
                              Pep.Intens.root = pep.ref[length(pep.ref)],
                              Pep.Ratios.root = pep.ratios.ref[length(pep.ratios.ref)],
                              log.Pep.Intens = FALSE, log.Pep.Ratios = 2,
                              Prot.LFQ.to.log = TRUE, Prot.Ratios.to.log = TRUE,
                              Mods = Mod4Quant, Mods.to.Exclude = Mod2Xclud,
                              Discard.unmod = Discard.unmod,
                              Min.N = 1, Priority = c("int", "rat")[(LabelType %in% c("SILAC"))+1],
                              Refs_Mode = RefRat_Mode,
                              cl = parClust)
    #colnames(quant.data2)
    saveFun(quant.data2, file = "quant.data2.RData")
    #loadFun("quant.data2.RData")
  }
  if ((Param$QuantMeth == "Prot.Quant.Unique")||(QuantMethods_all)) {
    # Profiles method, prefer unique, no weights
    DatAnalysisTxt <- gsub("\\.$",
                           paste0(", and quantified using an in-house algorithm which: ",
                                  "i) computes a mean protein-level profile across samples using individual, normalized peptidoform profiles (\"relative quantitation\" step), ",
                                  "ii) following the best-flyer hypothesis, normalizes this profile to the mean intensity level of the most intense peptidoform (\"unscaled absolute quantitation\" step); ",
                                  "for protein groups with at least 3 unique peptidoforms, only unique ones were used, otherwise razor peptidoforms were also included; ",
                                  tmpMod2Xclud, "peptidoforms and their unmodified counterparts were excluded from the calculations.",
                                  collapse = " "), DatAnalysisTxt)
    source(parSrc, local = FALSE)
    quant.data3 <- Prot.Quant(Prot = PG, Mode = "PreferUnique",
                              Peptide.IDs = "Razor peptide IDs", Unique.peptide.IDs = "Unique peptide IDs",
                              Pep = pep[Pep2Use,], id = "id",
                              Summary.method = "mean", #Summary.weights = "Weights",
                              Intensity.weights = FALSE,
                              experiments.map = Exp.map, param = Param,
                              Pep.Intens.root = pep.ref[length(pep.ref)],
                              Pep.Ratios.root = pep.ratios.ref[length(pep.ratios.ref)],
                              log.Pep.Intens = FALSE, log.Pep.Ratios = 2,
                              Prot.LFQ.to.log = TRUE, Prot.Ratios.to.log = TRUE,
                              Mods = Mod4Quant, Mods.to.Exclude = Mod2Xclud,
                              Discard.unmod = Discard.unmod,
                              Priority = c("int", "rat")[(LabelType %in% c("SILAC"))+1],
                              Refs_Mode = RefRat_Mode,
                              cl = parClust)
    #colnames(quant.data3)
    saveFun(quant.data3, file = "quant.data3.RData")
    #loadFun("quant.data3.RData")
  }
  if ((Param$QuantMeth == "Prot.Quant.Unique + weights")||(QuantMethods_all)) {
    # Profiles method, prefer unique, weights
    DatAnalysisTxt <- gsub("\\.$",
                           paste0(", and quantified using an in-house algorithm which: ",
                                  "i) computes a protein-level profile across samples as the weighted mean of individual, normalized peptidoform profiles (\"relative quantitation\" step, weights = ",
                                  weightsInsrt,
                                  ", where PEP and CV = the peptidoform's Posterior Error Probability and Coefficient of Variation, resp.), ",
                                  "ii) following the best-flyer hypothesis, normalizes this profile to the mean intensity level of the most intense peptidoform (\"unscaled absolute quantitation\" step); ",
                                  "for protein groups with at least 3 unique peptidoforms, only unique ones were used, otherwise razor peptidoforms were also included; ",
                                  tmpMod2Xclud, "peptidoforms and their unmodified counterparts were excluded from the calculations.",
                                  collapse = " "), DatAnalysisTxt)
    source(parSrc, local = FALSE)
    quant.data4 <- Prot.Quant(Prot = PG, Mode = "PreferUnique",
                              Peptide.IDs = "Razor peptide IDs", Unique.peptide.IDs = "Unique peptide IDs",
                              Pep = pep[Pep2Use,], id = "id",
                              Summary.method = "weighted.mean", Summary.weights = "Weights",
                              Intensity.weights = FALSE,
                              experiments.map = Exp.map, param = Param,
                              Pep.Intens.root = pep.ref[length(pep.ref)],
                              Pep.Ratios.root = pep.ratios.ref[length(pep.ratios.ref)],
                              log.Pep.Intens = FALSE, log.Pep.Ratios = 2,
                              Prot.LFQ.to.log = TRUE, Prot.Ratios.to.log = TRUE,
                              Mods = Mod4Quant, Mods.to.Exclude = Mod2Xclud,
                              Discard.unmod = Discard.unmod,
                              Min.N = 1, Priority = c("int", "rat")[(LabelType %in% c("SILAC"))+1],
                              Refs_Mode = RefRat_Mode,
                              cl = parClust)
    #colnames(quant.data4)
    saveFun(quant.data4, file = "quant.data4.RData")
    #loadFun("quant.data4.RData")
  }
}
if ((Param$QuantMeth == "Prot.Quant2 + weights")||(QuantMethods_all)) {
  # Averaging method, weights
  DatAnalysisTxt <- gsub("\\.$",
                         paste0(", and quantified using an in-house algorithm which: ",
                                "i) computes protein group expression values as the mean of peptidoform intensity values (\"relative quantitation\" step), then ",
                                # "i) computes protein group expression values as the weighted mean of peptidoform intensity values (\"relative quantitation\" step, weights = ",
                                weightsInsrt,
                                ", where PEP and CV = the peptidoform's Posterior Error Probability and Coefficient of Variation, resp.), then ",
                                "ii) following the best-flyer hypothesis, normalizes this profile to the mean intensity level of the most intense peptidoform (\"unscaled absolute quantitation\" step); ",
                                "only unique and razor peptidoforms were used; ",
                                tmpMod2Xclud, " peptidoforms and their unmodified counterparts were excluded from the calculations.",
                                collapse = " "), DatAnalysisTxt)
  quant.data5 <- Prot.Quant2(Prot = PG, Peptide.IDs = Pep4Quant, Pep = pep[Pep2Use,], id = "id",
                             Summary.method = "weighted.mean", Summary.weights = "Weights",
                             experiments.map = Exp.map, param = Param,
                             Pep.Intens.root = pep.ref[length(pep.ref)],
                             Pep.Ratios.root = pep.ratios.ref[length(pep.ratios.ref)],
                             log.Pep.Intens = FALSE, log.Pep.Ratios = 2,
                             Prot.LFQ.to.log = TRUE, Prot.Ratios.to.log = TRUE,
                             Mods = Mod4Quant, Mods.to.Exclude = Mod2Xclud,
                             Discard.unmod = Discard.unmod,
                             Min.N = 1,
                             Priority = c("int", "rat")[(LabelType %in% c("SILAC"))+1])
  #colnames(quant.data5)
  saveFun(quant.data5, file = "quant.data5.RData")
  #loadFun("quant.data5.RData")
}
if ((Param$QuantMeth == "Prot.Quant2")||(QuantMethods_all)) {
  # Averaging method, no weights
  DatAnalysisTxt <- gsub("\\.$",
                         paste0(", and quantified using an in-house algorithm which: ",
                                "i) computes protein group expression values as the mean of peptidoform intensity values (\"relative quantitation\" step), then ",
                                "ii) following the best-flyer hypothesis, normalizes this profile to the mean intensity level of the most intense peptidoform (\"unscaled absolute quantitation\" step); ",
                                "only unique and razor peptidoforms were used; ",
                                tmpMod2Xclud, "peptidoforms and their unmodified counterparts were excluded from the calculations.",
                                collapse = " "), DatAnalysisTxt)
  quant.data6 <- Prot.Quant2(Prot = PG, Peptide.IDs = Pep4Quant, Pep = pep[Pep2Use,], id = "id",
                             Summary.method = "mean", #Summary.weights = "Weights",
                             experiments.map = Exp.map, param = Param,
                             Pep.Intens.root = pep.ref[length(pep.ref)],
                             Pep.Ratios.root = pep.ratios.ref[length(pep.ratios.ref)],
                             log.Pep.Intens = FALSE, log.Pep.Ratios = 2,
                             Prot.LFQ.to.log = TRUE, Prot.Ratios.to.log = TRUE,
                             Mods = Mod4Quant, Mods.to.Exclude = Mod2Xclud,
                             Discard.unmod = Discard.unmod,
                             Min.N = 1,
                             Priority = c("int", "rat")[(LabelType %in% c("SILAC"))+1]
  )
  #colnames(quant.data6)
  saveFun(quant.data6, file = "quant.data6.RData")
  #loadFun("quant.data6.RData")
}
if ((Param$QuantMeth == "IQ_MaxLFQ")||(QuantMethods_all)) {
  # MaxLFQ
  if (!require(iq, quietly = TRUE)) { pak::pak("iq") }
  library(iq)
  kol <- grep(topattern(pep.ref[length(pep.ref)]), colnames(pep), value = TRUE)
  kol <- grep("\\.REF$", kol, value = TRUE, invert = TRUE)
  Pep2Use2 <- Pep2Use
  if (nrow(Mod2Xclud)) {
    pat <- paste(apply(Mod2Xclud, 1, function(x) {
      paste0(x[[2]], "(", x[[1]], ")", collapse = "|")
    }), collapse = "|")
    if (!Mod.Excl.is.strict) {
      Pep2Use2 <- Pep2Use2[grep(pat, pep$"Modified sequence"[Pep2Use2], invert = TRUE)]
    } else {
      seq <- pep$Sequence[Pep2Use2][grep(pat, pep$"Modified sequence"[Pep2Use2])]
      Pep2Use2 <- Pep2Use2[which(!pep$Sequence[Pep2Use2] %in% seq)]
    }
  }
  # That part rewritten to bypass the very slow iq::create_protein_list() 
  tmpPGlist <- listMelt(strsplit(PG$"Peptide IDs", ";"),
                        1:nrow(PG),
                        c("id", "Row"))
  tmpPGlist$id <- as.integer(tmpPGlist$id)
  tmpPGlist <- tmpPGlist[which(tmpPGlist$id %in% Pep$id[Pep2Use2]),]
  tmpPep2 <- tmpPep[, c("id", kol)]
  tmpPGlist[, kol] <- tmpPep2[match(tmpPGlist$id, tmpPep2$id), kol]
  colnames(tmpPGlist) <- gsub(topattern(pep.ref[length(pep.ref)]), "", colnames(tmpPGlist))
  kol2 <- gsub(topattern(pep.ref[length(pep.ref)]), "", kol)
  tmpPGlist <- as.data.table(tmpPGlist[, c("Row", kol2)])
  tmpPGlist <- split(tmpPGlist[, ..kol2], tmpPGlist$Row)
  names(tmpPGlist) <- Prot$`Leading protein IDs`[as.integer(names(tmpPGlist))]
  #
  # This part below is horrendously SLOWWWWWWWW...
  quant.data7 <- iq::create_protein_table(tmpPGlist)
  
  quant.data7 <- create_protTbl(tmpPGlist, cl = parClust)
  
  
  # One more reason to use in house quant...
  # I need to look into what that function is doing...
  # Also to put a warning in the parameters app if this becomes an option...
  #
  quant.data7 <- quant.data7$estimate/log2(10) # Convert to log10
  quant.data7 <- as.data.frame(quant.data7)
  colnames(quant.data7) <- paste0("log10(Expr.) - ", colnames(quant.data7))
  quant.data7$"Leading protein IDs" <- names(tmpPGlist)
  w <- which(!PG$"Leading protein IDs" %in% quant.data7$"Leading protein IDs")
  temp <- as.data.frame(matrix(rep(NA, length(w)*length(exprsCol)), ncol = length(exprsCol)))
  colnames(temp) <- exprsCol
  temp$"Leading protein IDs" <- PG$"Leading protein IDs"[w]
  quant.data7 <- rbind(quant.data7, temp)
  quant.data7 <- quant.data7[match(PG$"Leading protein IDs", quant.data7$"Leading protein IDs"),]
  #colnames(quant.data7)
  saveFun(quant.data7, file = "quant.data7.RData")
  #loadFun("quant.data7.RData")
  DatAnalysisTxt <- gsub("\\.$",
                         paste0(", and quantified using the MaxLFQ algorithm as implemented in the iq package. ",
                                "Only unique and razor peptidoforms were used; ",
                                tmpMod2Xclud, "peptidoforms and their unmodified counterparts were excluded from the calculations.",
                                collapse = " "), DatAnalysisTxt)
}
if ((Param$QuantMeth == "Top3")||(QuantMethods_all)) {
  # Top3
  DatAnalysisTxt <- gsub("\\.$",
                         paste0(", and quantified using a local implementation of the Top3 algorithm, allowing for a value to be calculated also if only 1 or 2 peptidoforms were present (the resulting systematic bias between proteins with 1, 2 or 3 peptidoforms is then corrected for); ",
                                "only unique and razor peptidoforms were used; ",
                                tmpMod2Xclud, "peptidoforms and their unmodified counterparts were excluded from the calculations.",
                                collapse = " "), DatAnalysisTxt)
  quant.data8 <- TopN(3, PG, "Razor peptide IDs", pep[Pep2Use,],
                      Pep.Intens.Nms = paste0(pep.ref[length(pep.ref)], RSA$values),
                      log.Pep.Intens = FALSE, Mods = Mod4Quant, Out.Norm = FALSE, corr = "global")
  colnames(quant.data8) <- paste0("log10(Expr.) - ", RSA$values)
  saveFun(quant.data8, file = "quant.data8.RData")
  #loadFun("quant.data8.RData")
}
if ((Param$QuantMeth == "Top1")||(QuantMethods_all)) {
  # Top1
  DatAnalysisTxt <- gsub("\\.$",
                         paste0(", and quantified using the highest peptidoform intensity value in each sample. ",
                                "Only unique and razor peptidoforms were used; ",
                                tmpMod2Xclud, "peptidoforms and their unmodified counterparts were excluded from the calculations.",
                                collapse = " "), DatAnalysisTxt)
  quant.data9 <- TopN(1, Prot = PG, Pep = pep[Pep2Use,],
                      Pep.Intens.Nms = paste0(pep.ref[length(pep.ref)], RSA$values),
                      log.Pep.Intens = FALSE, Mods = Mod4Quant, Out.Norm = FALSE, corr = "global")
  colnames(quant.data9) <- paste0("log10(Expr.) - ", RSA$values)
  saveFun(quant.data9, file = "quant.data9.RData")
  #loadFun("quant.data9.RData")
}
if (QuantMethods_all) {
  # Comparison of LFQ profiles and summation approaches:
  cv <- CV <- CV_SD <- list()
  MaxPep <- 5
  # Create average intensity normalisation across methods
  # We want to look at the variability between methods,
  # but need to correct for differences in estimates of absolute abundance,
  # as we are really interested in precision of sample-to-sample variability
  NormFact1 <- sapply(seq_along(QuantMethods), function(i) {
    tmp <- get(QuantData[i])
    tmp$"Leading protein IDs" <- PG$"Leading protein IDs"
    tmp2 <- tmp[, exprsCol]
    tmp2 <- rowMeans(tmp2, na.rm = TRUE)
    w <- which(PG$"Leading protein IDs" %in% tmp$"Leading protein IDs")
    res <- rep(NA, nrow(PG))
    res[w] <- tmp2[match(PG$"Leading protein IDs"[w], tmp$`Leading protein IDs`)]
    return(res)
  })
  NormFact2 <- rowMeans(NormFact1, na.rm = TRUE)
  for (i in seq_along(QuantMethods)) {
    tmp1 <- get(QuantData[i])
    tmp1$"Leading protein IDs" <- PG$"Leading protein IDs"
    tmp1 <- tmp1[, c(exprsCol, "Leading protein IDs")]
    tmp2 <- sapply(VPAL$values, function(x) {
      kol <- paste0("log10(Expr.) - ", Exp.map$Ref.Sample.Aggregate[which(Exp.map[[VPAL$column]] == x)])
      x <- apply(tmp1[ , kol, drop = FALSE], 1, function(y) {
        y <- is.all.good(y)
        sd(y)/mean(y)
      })
    })
    tmp2 <- data.frame(CV = rowMeans(tmp2, na.rm = TRUE))
    colnames(tmp1) <- gsub(topattern("log10(Expr.) - "), "", colnames(tmp1))
    tmp1[, c("id", "Razor + unique peptides")] <- PG[match(tmp1$"Leading protein IDs", PG$"Leading protein IDs"),
                                                     c("id", "Peptide counts (razor+unique)")]
    tmp1$"Razor + unique peptides" <- as.integer(gsub(";.*", "", tmp1$"Razor + unique peptides"))
    tmp2$"Razor + unique peptides" <- tmp1$"Razor + unique peptides"
    tmp2$`Razor + unique peptides`[which(tmp2$`Razor + unique peptides` >= MaxPep)] <- paste0(MaxPep, " or more") # No need to look at too rare cases
    tmp2$`Razor + unique peptides` <- as.factor(tmp2$`Razor + unique peptides`)
    cv[[QuantMethods[i]]] <- tmp2[which(is.all.good(tmp2$CV, 2)),]
    CV[[QuantMethods[i]]] <- mean(cv[[QuantMethods[i]]]$CV)
    tmp3 <- tmp1
    m <- match(tmp3$`Leading protein IDs`, PG$`Leading protein IDs`)
    tmp3[, gsub(topattern("log10(Expr.) - "), "", exprsCol)] <- sweep(tmp3[, gsub(topattern("log10(Expr.) - "), "", exprsCol)], 1, NormFact1[m,i], "-")+NormFact2[m]
    tmp1 <- reshape2::melt(tmp1, id.vars = c("id", "Razor + unique peptides", "Leading protein IDs"))
    tmp3 <- reshape2::melt(tmp3, id.vars = c("id", "Razor + unique peptides", "Leading protein IDs"))
    colnames(tmp1) <- c("id", "Razor + unique peptides", "Leading protein IDs", "Sample", QuantMethods[i])
    tmp1[[paste0(QuantMethods[i], " - renorm")]] <- tmp3$value
    tmp1$Sample <- as.character(tmp1$Sample)
    tmp1[, RSA$names] <- Isapply(strsplit(tmp1$Sample, "___"), unlist)
    tmp1$Sample <- cleanNms(tmp1$Sample)
    if (i == 1) { tmp <- tmp1 } else {
      tmp[[QuantMethods[i]]] <- NA
      w <- which(tmp$"Leading protein IDs" %in% tmp1$"Leading protein IDs")
      tmp[w, QuantMethods[i]] <- tmp1[match(tmp$"Leading protein IDs"[w], tmp1$"Leading protein IDs"), QuantMethods[i]] 
    }
  }
  #plot <- ggplot(tmp) + geom_scattermore(aes(x = Prot.Quant, y = Prot.Quant2, colour = Sample), size = 0.1, alpha = 0.01) +
  #  coord_fixed() + theme_bw() + facet_wrap(~Sample)
  comb <- gtools::combinations(length(QuantMethods), 2, QuantMethods)
  dir <- paste0(wd, "/Workflow control/Protein groups/Quantitative methods")
  if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE) }
  dirlist <- unique(c(dirlist, dir))
  for (i in 1:nrow(comb)) { #i <- 1
    tmp1 <- tmp
    tmp1$X <- tmp1[[comb[i, 1]]]
    tmp1$Y <- tmp1[[comb[i, 2]]]
    tmp1 <- tmp1[which((is.all.good(tmp1$X, 2))&(is.all.good(tmp1$Y, 2))),]
    tmp1$`Razor + unique peptides`[which(tmp1$`Razor + unique peptides` >= MaxPep)] <- paste0(MaxPep, " or more") # No need to look at too rare cases
    tmp1$`Razor + unique peptides` <- as.factor(tmp1$`Razor + unique peptides`)
    tmp2 <- tmp1
    tmp2$X <- tmp2[[paste0(comb[i, 1], " - renorm")]]
    tmp2$Y <- tmp2[[paste0(comb[i, 2], " - renorm")]]
    tst <- vapply(RSA$names, function(x) { length(get(substr(x, 1, 3))) }, 1)
    tst <- RSA$names[which(tst > 1)]
    form <- as.formula(paste0(tst[1], "~", paste(tst[2:length(tst)], collapse = "+")))
    ttl1 <- paste0("LFQ method comparison: ", comb[i, 1], " VS ", comb[i, 2])
    plot1 <- ggplot(tmp1) +
      geom_scattermore(aes(x = X, y = Y, colour = `Razor + unique peptides`),
                       shape = 16, size = 0.1, alpha = 0.1) +
      scale_color_viridis_d(begin = 0.25) +
      coord_fixed() + theme_bw() + facet_grid(form) + xlab(comb[i, 1]) + ylab(comb[i, 2]) + ggtitle(ttl1)
    #poplot(plot1, 12, 20)
    ttl1a <- gsub("\\:", " -", ttl1)
    suppressMessages({
      ggsave(paste0(dir, "/", ttl1a, ".jpeg"), plot1, dpi = 300)
      ggsave(paste0(dir, "/", ttl1a, ".pdf"), plot1, dpi = 300)
    })
    ReportCalls <- AddPlot2Report(plot1, Title = ttl1a)
    ttl2 <- paste0("LFQ method comparison (renormalized): ", comb[i, 1], " VS ", comb[i, 2])
    plot2 <- ggplot(tmp1) +
      geom_scattermore(aes(x = X, y = Y, colour = `Razor + unique peptides`),
                       shape = 16, size = 0.1, alpha = 0.1) +
      scale_color_viridis_d(begin = 0.25) +
      coord_fixed() + theme_bw() + facet_grid(form) + xlab(comb[i, 1]) + ylab(comb[i, 2]) + ggtitle(ttl2)
    print(plot2) # This type of QC plot does not need to pop up, the side panel is fine
    ttl2a <- gsub("\\:", " -", ttl1)
    suppressMessages({
      ggsave(paste0(dir, "/", ttl2a, ".jpeg"), plot2, dpi = 300)
      ggsave(paste0(dir, "/", ttl2a, ".pdf"), plot2, dpi = 300)
    })
    ReportCalls <- AddPlot2Report(plot2, Title = ttl2a)
  }
  for (meth in QuantMethods) {
    cat(paste0(meth, ": ", nrow(cv[[meth]]), " quantified protein groups, mean intra-sample-groups CV = ", round(CV[[meth]]*100, 2), "%\n"))
  }
  tst <- reshape2::melt(cv)
  tst$variable <- NULL
  colnames(tst) <- c("Razor + unique peptides", "CV", "Method")
  ttl <- "Distribution of intra-sample groups coefficients of variation"
  plot <- ggplot(tst) +
    geom_histogram(aes(x = CV, fill = Method), bins = 100) + theme_bw() + ggtitle(ttl) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_viridis_d(begin = 0.25) +
    facet_grid(Method~`Razor + unique peptides`) +
    theme(strip.text.y.right = element_text(angle = 0))
  print(plot) # This type of QC plot does not need to pop up, the side panel is fine
  ttla <- gsub(":", " -", ttl)
  suppressMessages({
    ggsave(paste0(dir, "/", ttla, ".jpeg"), plot, dpi = 300)
    ggsave(paste0(dir, "/", ttla, ".pdf"), plot, dpi = 300)
  })
  ReportCalls <- AddPlot2Report(Title = ttla)
}
ReportCalls <- AddTxt2Report(paste0("Protein groups quantitation done using method: ", Param$QuantMeth))
if (!exists(QuantData[match(Param$QuantMeth, QuantMethods)])) { # (for when rerunning script without rerunning quant methods)
  load(paste0(QuantData[match(Param$QuantMeth, QuantMethods)], ".RData"))
}
quant.data <- get(QuantData[match(Param$QuantMeth, QuantMethods)])
Prot.Expr.Root %<o% c(Original = "log10(Expr.) - ")
Prot.Rat.Root %<o% c(Original = "log2(Ratio) - ")
stopifnot(length(grep(topattern(Prot.Expr.Root), colnames(quant.data))) > 0,
          length(grep(topattern(Prot.Rat.Root), colnames(quant.data))) > 0)
#write.csv(quant.data, file = "Quantitative data.csv", row.names = FALSE)
.obj <- unique(c("quant.data", "Prot.Expr.Root", "Prot.Rat.Root", .obj))
DatAnalysisTxt <- paste0(DatAnalysisTxt, " Estimated expression values were log10-converted...")
