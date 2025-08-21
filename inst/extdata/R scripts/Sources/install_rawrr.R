# Check that the correct version of rawrr is installed
# Currently we are limiting ourselves to 1.11.14 which we get from Bioconductor

rawrrVers <- rawrrVersions <- c(#"github",
  #"bioc",
  "bioc_1.11.14")
rawrr_tst <- try({
  #unloadNamespace("rawrr")
  #remove.packages("rawrr")
  #
  rawrrVers_inst <- packageVersion("rawrr")
  if ((rawrrVers_inst != "1.11.14")||(!require(rawrr, quietly = TRUE))) {
    # if (length(rawrrVers) > 1) {
    #   rawrrVers <- dlg_list(rawrrVersions, rawrrVersions[3], title = "Select which rawrr version should be installed")$res
    # }
    # if (rawrrVers == "github") {
    #   pak::pak("cpanse/rawrr")
    #   #install.packages('http://fgcz-ms.uzh.ch/~cpanse/rawrr_0.2.1.tar.gz', repo = NULL) # Old address, now on github
    # }
    # if (rawrrVers == "bioc") {
    #   if (!require("BiocManager", quietly = TRUE)) { install.packages("BiocManager") }
    #   BiocManager::install("rawrr", version = "1.11")
    #   BiocManager::install("rawrr")
    # }
    # if (rawrrVers == "bioc_1.11.14") {
    url <- "https://bioconductor.org/packages/3.19/bioc/src/contrib/Archive/rawrr/rawrr_1.11.14.tar.gz"
    require(curl)
    dstfl <- paste0("C:/Users/", Sys.getenv("USERNAME"), "/Downloads/rawrr_1.11.14.tar.gz")
    curl_download(url, dstfl)
    install.packages(dstfl)
    #}
    # 
  }
  yesRawFileReaderLicenseIsAccepted <- function () {
    licenseFile <- file.path(system.file(package = "rawrr"), 
                             "rawrrassembly", "RawFileReaderLicense.txt")
    stopifnot(file.exists(licenseFile))
    eulaFile <- file.path(rawrr::rawrrAssemblyPath(), "eula.txt")
    msg <- c("# By changing the setting below to TRUE you are accepting ", 
             "the Thermo License agreement.")
    if (!file.exists(eulaFile)) {
      file.show(licenseFile)
      response <- "y"
      if (tolower(response) == "y") {
        if (isFALSE(dir.exists(dirname(eulaFile)))) {
          dir.create(dirname(eulaFile), recursive = TRUE)
        }
        fileConn <- file(eulaFile)
        writeLines(paste(msg, paste0("# ", date()), "eula=true", 
                         sep = "\n"), fileConn)
        close(fileConn)
        return(TRUE %in% grepl("eula=true", tolower(readLines(eulaFile))))
      }
    } else { return(TRUE %in% grepl("eula=true", tolower(readLines(eulaFile)))) }
    msg <- "Yes, we accept Thermo's License agreement, get on with it!"
    cat(msg, "\n")
  }
  installRawFileReaderDLLsNoAcpt <- function (sourceUrl = rawrr:::.thermofisherlsmsUrl(), ...) { #sourceUrl = rawrr:::.thermofisherlsmsUrl()
    rawfileReaderDLLsPath <- rawrr::rawrrAssemblyPath()
    if (isTRUE(dir.exists(rawfileReaderDLLsPath))) {
      msg <- sprintf("removing files in directory '%s'", rawfileReaderDLLsPath)
      message(msg)
      file.remove(file.path(rawrr::rawrrAssemblyPath(), list.files(rawrr::rawrrAssemblyPath())))
    }
    if (isFALSE(dir.exists(rawfileReaderDLLsPath))) {
      dir.create(rawfileReaderDLLsPath, recursive = TRUE)
    }
    #
    stopifnot(yesRawFileReaderLicenseIsAccepted())
    .rawfileReaderDLLs <- getAnywhere(.rawfileReaderDLLs)
    .rawfileReaderDLLs <- .rawfileReaderDLLs$objs[[1]]
    rv <- vapply(.rawfileReaderDLLs(), function(dll) { #dll <- .rawfileReaderDLLs()[1]
      destfile <- file.path(rawfileReaderDLLsPath, dll)
      download.file(file.path(paste0(sourceUrl, c("", "Assemblies/")[(rawrrVersions == "github")+1],
                                     dll)), destfile = destfile, 
                    mode = "wb")
    }, 0)
    rv
  }
  installRawFileReaderDLLsNoAcpt(sourceUrl = rawrr:::.thermofisherlsmsUrl())
  #rawrr::installRawFileReaderDLLs(sourceUrl = .thermofisherlsmsUrl())
  #require(rawrr); getAnywhere(".isRawFileReaderLicenseAccepted")
  rawrr::installRawrrExe()
  library(rawrr)
}, silent = TRUE)
if ("try-error" %in% class(rawrr_tst)) {
  warning("Could not install rawrr 1.11.14!")
}
