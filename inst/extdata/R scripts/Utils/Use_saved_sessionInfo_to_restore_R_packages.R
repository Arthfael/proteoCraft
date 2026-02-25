#####################################################
#                                                   #
#  A script to reload the correct package versions  #
#      using our saved output of sessionInfo()      #
#                                                   #
#       This is a BAD solution, but it works.       #
#                                                   #
#                 Use renv instead!                 #
#                                                   #
#####################################################
# If the output of sessionInfo() is all you have, then this will get you some way... but maybe not all the way!
# Assumes all packages are from cran or bioconductor (not github)
# Other repositories are currently not supported
# There will be some exceptions which we will just make sure we have installed/loaded from the beginning, regardless of version!
install.packages("pak", repos = sprintf("https://r-lib.github.io/p/pak/stable/%s/%s/%s", .Platform$pkgType, R.Version()$os, R.Version()$arch))
pak::pkg_install("rstudioapi")

# These packages are either base or exceptions which this code needs to work and for which we will use whichever version we have:
except <- c("parallel", "snow", "rstudioapi", "BiocManager", "pak")
# These include packages which this script itself uses, and some RCpp-related packages which fail because of some gcc-related issues
# (it seems to forget to release some handles?)
sessInfFl <- rstudioapi::selectFile()
#system(paste0("open \"", sessInfFl, "\""))
sessInf <- readLines(sessInfFl)

currRVers <- paste(sessionInfo()$R.version[c("major", "minor")], collapse = ".")
reloadRVers <- gsub("^R +version +| +\\(.*", "", sessInf[1])
if (currRVers != reloadRVers) {
  stop(paste0("RStudio is currently using R version ", currRVers, " but the session you are trying to reload used ", reloadRVers, "!\n",
              "This script cannot change RStudio's version during a session, so please go to 'Tools > Global Options...', change the version, restart the session then start again."))
}

# Available packages
avail <- as.data.frame(available.packages(type = "both"))
avail$Origin <- "CRAN"
avail2 <- as.data.frame(available.packages(repos = BiocManager::repositories()[1], type = "both"))
avail2$Origin <- "Bioc"
avail <- rbind(avail, avail2)
avail$Dep <- gsub(" *\\([^\\)]+\\)", "", avail$Depends)
avail$Dep <- gsub("^R$", "", gsub("^R, |^R,\n", "", avail$Dep))
avail$Dep <- strsplit(avail$Dep, ",[ \n]")
avail$Dep <- sapply(avail$Dep, function(x) { x[which(!is.na(x))] })
avail$nDep <- sapply(avail$Dep, length)
except <- unique(c(except, unlist(avail$Dep[match(except, avail$Package)])))

for (pck in except) {
  if (!require(pck, character.only = TRUE, quietly = TRUE)) {
    a <- try(install.packages(pck, dependencies = TRUE), silent = TRUE)
  }
}
for (pck in except) {
  library(pck, character.only = TRUE)
}

pcksFull <- sessInf[which(sessInf == "other attached packages:"):length(sessInf)]
pcksFull <- grep("\\[[0-9]+\\]", pcksFull, value = TRUE)
pcksFull <- unlist(strsplit(gsub(" *\\[[0-9]+\\] *", "", pcksFull), " +"))
proteoCraft_Vers <- grep("^proteoCraft_?", pcksFull, value = TRUE)
aRmel_Vers <- grep("^aRmel_?", pcksFull, value = TRUE) # Old name
stopifnot(length(c(proteoCraft_Vers, aRmel_Vers)) > 0)
pcksFull <- grep("^((proteoCraft)|(aRmel))_?", pcksFull, value = TRUE, invert = TRUE)
if (exists("except")) {
  for (pck in except) {
    pcksFull <- grep(paste0("^", pck, "_?"), pcksFull, value = TRUE, invert = TRUE)
  }
}
pckDF <- data.frame(Package = gsub("_.*", "", pcksFull),
                    version = NA)

w <- grep("_", pcksFull)
pckDF$version[w] <- gsub("[^_]+_", "", pcksFull[w])

# Skip packages for which the correct version is already installed
inst <- as.data.frame(installed.packages())
tst1 <- do.call(paste, c(pckDF[, c("Package", "version")]))
tst2 <- do.call(paste, c(inst[, c("Package", "Version")]))
w <- which(!tst1 %in% tst2)

if (length(w)) {
  pckDF <- pckDF[w,]
  # Re-order
  w1 <- which(avail$Package %in% pckDF$Package)
  pckDF <- pckDF[c(match(avail$Package[w1], pckDF$Package),
                   which(!pckDF$Package %in% avail$Package)),]
  m <- match(pckDF$Package, avail$Package)
  pckDF$Depends <- avail$Depends[m]
  pckDF$Dep <- avail$Dep[m]
  pckDF$nDep <- avail$nDep[m]
  pckDF <- pckDF[order(pckDF$nDep, decreasing = FALSE),]
  pckDF$Origin <- avail$Origin[m]
  w <- which(pckDF$nDep > 0)
  rg <- 1:nrow(pckDF)
  for (i in w) { #i <- w[1]
    pck <- pckDF$Package[i]
    dep <- pckDF$Dep[[i]]
    w1 <- rev(which(pckDF$Package %in% dep))[1]
    if (!is.na(w1)) {
      w0 <- which(pckDF$Package == pck)
      wRst <- rg[which(!rg %in% c(w1, w0))]
      wAbove <- wRst[which(wRst < w1)]
      wBelow <- wRst[which(wRst > w1)]
      ord <- c(wAbove, w1, w0, wBelow)
      stopifnot(length(ord) == nrow(pckDF))
      pckDF <- pckDF[ord,]
    }
  }
  
  # Make/check cluster
  N.clust <- parallel::detectCores()-1
  a <- 1
  tst <- try(clusterExport(parClust, "a", envir = environment()), silent = TRUE)
  if ("try-error" %in% class(tst)) {
    try(parallel::stopCluster(parClust), silent = TRUE)
    parClust <- parallel::makeCluster(N.clust, type = "SOCK")
  }
  
  # Uninstall packages
  parallel::clusterExport(parClust, list("pckDF", "inst"), envir = environment())
  w <- which(pckDF$Package %in% inst$Package)
  if (length(w)) {
    kount <- 0
    while ((kount < 5)&&(length(w))) {
      for (pck in pckDF$Package[w]) {
        parallel::parLapply(parClust, w, function(i) { #i <- w[1]
          tst <- try(unloadNamespace(pck), silent = TRUE)
        })
      }
      parallel::parLapply(parClust, w, function(i) { #i <- w[1]
        pck <- pckDF$Package[i]
        if (pck %in% inst$Package) {
          tst <- try({
            a <- try({ remove.packages(pck) })
            if ("try-error" %in% class(a)) { pkg_remove(pkg) }
          }, silent = TRUE)
          return(!"try-error" %in% class(tst))
        } else { return(TRUE) }
      })
      inst <- as.data.frame(installed.packages())
      clusterExport(parClust, list("inst"), envir = environment())
      w <- which(pckDF$Package %in% inst$Package)
      kount <- kount + 1
    }
  }
  
  # There may be some which we did not manage to remove, we'll still try to install them just in case
  inst <- as.data.frame(installed.packages())
  Tsts <- setNames(pckDF$Package %in% inst$Package, pckDF$Package)
  w <- which(Tsts)
  tst2 <- sapply(w, function(i) { #i <- w[1]
    pck <- pckDF$Package[i]
    a <- try(library(pck, character.only = TRUE), silent = TRUE)
    a <- !"try-error" %in% class(a)
    if (a) {
      try({
        unloadNamespace(pck) # Can stay here as long as it is within sapply
        # If you ever turn this sapply into parSapply then you will need to move unloadNamespace to a separate sapply(... parLapply...) (see above)
        # so it applies to the whole cluster before we even 
        b <- try({ remove.packages(pck) }, silent = TRUE)
        if ("try-error" %in% class(b)) { pak::pkg_remove(pck) }
      }, silent = TRUE)
    }
    return(a)
  })
  print(length(which(!tst2)))
  Tsts[w[which(!tst2)]] <- FALSE
  
  # Install packages
  inst <- as.data.frame(installed.packages())
  pckDF$goAhead <- sapply(pckDF$Dep, function(x) {
    x <- x[which(x %in% pckDF$Package)]
    sum(!x %in% inst$Package) == 0
  })
  pckDF$Installed <- pckDF$Package %in% inst$Package
  #i <- 1
  kount <- 0
  wNO <- which((pckDF$goAhead)&(!pckDF$Installed))
  if (!exists("wNO2")) { wNO2 <- NA }
  while ((length(wNO))&&(!identical(wNO, wNO2))) {
    #wNO_ord <- wNO <- which(!pckDF$Installed)
    wNO <- which((pckDF$goAhead)&(!pckDF$Installed))
    #pckDF$Package[wNO]
    #if (i == -1) { wNO_ord <- rev(wNO_ord) } # In case package installation is important
    clusterExport(parClust, list("pckDF", "wNO"), envir = environment())
    A <- parallel::parLapply(parClust, wNO#_ord
                             , function(i) { #i <- wNO[1] #i <- wNO[2]
                               pck <- pckDF$Package[i]
                               vers <- pckDF$version[i]
                               orig <- pckDF$Origin[i]
                               retry <- FALSE
                               if (is.na(orig)) {
                                 orig <- "CRAN"
                                 retry <- TRUE
                               }
                               dep <- pckDF$Dep[[i]]
                               dep <- dep[which(!is.na(dep))]
                               dep <- dep[which(dep %in% pckDF$Package)]
                               inst <- as.data.frame(installed.packages())
                               ok <- sum(!dep %in% inst$Package) == 0
                               if (ok) {
                                 retry <- !((orig == "CRAN")&(!is.na(vers)))
                                 if (!retry) {
                                   url <- paste0("https://cran.r-project.org/src/contrib/Archive/", pck, "/", pck, "_", vers, ".tar.gz")
                                   a <- try(install.packages(url, repos = NULL, type = "source"), silent = TRUE)
                                   if ("try-error" %in% class(a)) {
                                     tst <- unlist(strsplit(vers, "\\."))
                                     l <- length(tst)
                                     if (l > 1) {
                                       vers <- paste(tst[1:(l-1)], collapse = ".")
                                       url <- paste0("https://cran.r-project.org/src/contrib/Archive/", pck, "/", pck, "_", vers, ".tar.gz")
                                       a <- try(install.packages(url, repos = NULL, type = "source"), silent = TRUE)
                                     } 
                                   }
                                   inst <- as.data.frame(installed.packages())
                                   retry <- !pck %in% inst$Package
                                 }
                                 if (retry) {
                                   a <- try(pak::pkg_install(pck), silent = TRUE)
                                 }
                                 # if (orig == "CRAN") {
                                 #   if (is.na(vers)) {
                                 #     #a <- try(utils::install.packages(pck, repos = "https://cran.wu.ac.at"), silent = TRUE) # Somehow the explicit mirror is only required within the parallel call
                                 #     a <- try(pak::pkg_install(pck), silent = TRUE)
                                 #   } else {
                                 #     a <- try(devtools::install_version(pck, vers, force = TRUE, upgrade = "never"), silent = TRUE)
                                 #     if ("try-error" %in% class(a)) { a <- try(utils::install.packages(pck, repos = "https://cran.wu.ac.at"), silent = TRUE) }
                                 #     if (("try-error" %in% class(a))&&(retry)) { orig <- "Bioc" }
                                 #     inst <- as.data.frame(installed.packages())
                                 #     if (!pck %in% inst$Package) {
                                 #       # In case we wrongly predicted where we can find a package
                                 #       warning(paste0("Could not install ", pck, " from CRAN, trying Bioconductor..."))
                                 #       orig <- "Bioc"
                                 #     }
                                 #   }
                                 # }
                                 # if (orig == "Bioc") {
                                 #   # Here we will ignore versions as it appears that they are linked to R versions and that we do not have any easy way to get a specific one
                                 #   a <- try(BiocManager::install(pck, update = FALSE), silent = TRUE)
                                 # }
                                 return(a)
                               } else { return(c()) }
                               #inst <- as.data.frame(installed.packages())
                               #return(pck %in% inst$Package)
                             })
    #A
    #
    inst <- as.data.frame(installed.packages())
    pckDF$Installed <- pckDF$Package %in% inst$Package
    pckDF$goAhead <- sapply(pckDF$Dep, function(x) {
      x <- x[which(x %in% pckDF$Package)]
      sum(!x %in% inst$Package) == 0
    })
    wNO2 <- which((pckDF$goAhead)&(!pckDF$Installed))
    solved <- pckDF$Package[wNO]
    solved <- solved[which(!solved %in% pckDF$Package[wNO2])]
    if (length(solved)) {
      print("Resolved the following packages:")
      print(solved)
    }
    print(paste0(sum(!pckDF$Installed), " packages remaining..."))
    #i <- -i
    kount <- kount + 1
    print(paste0("Current loop count = ", kount))
    print("")
  }
  inst <- as.data.frame(installed.packages())
  pckDF$Installed <- pckDF$Package %in% inst$Package
  wNO <- which(!pckDF$Installed)
  names(pckDF$Package)[wNO]
  # Hard fix for uchardet - as always
  pck <- "uchardet"
  m <- match(pck, pckDF$Package)
  if ((!is.na(m))&&(!pckDF$Installed[m])&&(!require(pck, character.only = TRUE, quietly = TRUE))) {
    vers <- pckDF$version[m]
    url <- paste0("https://cran.r-project.org/src/contrib/Archive/uchardet/uchardet_", vers, ".tar.gz")
    destfile <- paste0("uchardet_", vers, ".tar.gz")
    tst <- try(download.file(url, destfile, "curl"), silent = TRUE)
    if ("try-error" %in% class(tst)) { try(download.file(url, destfile, "wget"), silent = TRUE) }
    install.packages(destfile)
    unlink(destfile)
  }
  #
  inst <- as.data.frame(installed.packages())
  pckDF$Installed <- pckDF$Package %in% inst$Package
  wNO <- which(!pckDF$Installed)
  names(pckDF$Package)[wNO]
  #
  # If any remain, you will have to install them manually from github or whichever source you like
  #
  # Final check
  inst <- as.data.frame(installed.packages())
  tst1 <- do.call(paste, c(pckDF[, c("Package", "version")]))
  tst2 <- do.call(paste, c(inst[, c("Package", "Version")]))
  w <- which(!tst1 %in% tst2)
  length(w)
  # Conclusion
  # We were able to install all packages but the version number varies from the original for some
  # So... prefer renv!!!
  #View(pckDF[w, c("Package", "version")]);View(inst[match(pckDF$Package[w], inst$Package), c("Package", "Version")])
} else {
  cat("All packages are installed with the correct version...\nThat was easy ^^ !\n")
}

# Remember to also install proteoCraft...
if (length(proteoCraft_Vers)) {
  pak::pak(paste0("Arthfael/", gsub("_", "@", proteoCraft_Vers), ".tar.gz"))
  library(proteoCraft)
  proteoCraft::Configure()
}
# ... or its older ancestor, aRmel
if (length(aRmel_Vers)) {
  install.packages(paste0("H:/aRmel_package/", aRmel_Vers, ".tar.gz"))
  library(aRmel)
  aRmel::Configure()
}

# Also setwd to expected work directory
wd <- gsub("/Workflow control/?$", "", dirname(sessInfFl))
setwd(wd)

# Now load a backup - which may not be the one - if any - available in wd!
#load_Bckp()

# Find all R scripts currently in wd...
#fls <- list.files(wd, "\\.R$", full.names = TRUE);fls
#fl <- fls[6]
#rstudioapi::documentOpen(fl)
#openwd()
