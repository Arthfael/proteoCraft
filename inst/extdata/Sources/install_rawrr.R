# Check that the correct version of rawrr is installed
# Currently we are limiting ourselves to 1.11.14 which we get from Bioconductor

rawrrVers <- rawrrVersions <- c(#"github",
  #"bioc",
  "bioc_1.11.14") # NB: At least v1.18.0 seems to work too...
rawrr_tst <- try({
  #unloadNamespace("rawrr")
  #remove.packages("rawrr")
  #
  rawrrVers_inst <- packageVersion("rawrr")
  if ((rawrrVers_inst != "1.11.14") || (!requireNamespace("rawrr", quietly = TRUE))) {
    url <- "https://bioconductor.org/packages/3.19/bioc/src/contrib/Archive/rawrr/rawrr_1.11.14.tar.gz"
    dstfl <- paste0("C:/Users/", Sys.getenv("USERNAME"), "/Downloads/rawrr_1.11.14.tar.gz")
    curl::curl_download(url, dstfl)
    install.packages(dstfl)
  }
  yesRawFileReaderLicenseIsAccepted <- \(showLicense = FALSE) {
    licenseFile <- file.path(system.file(package = "rawrr"), 
                             "rawrrassembly", "RawFileReaderLicense.txt")
    stopifnot(file.exists(licenseFile))
    eulaFile <- file.path(rawrr::rawrrAssemblyPath(), "eula.txt")
    msg <- c("# By changing the setting below to TRUE you are accepting ", 
             "the Thermo License agreement.")
    if (!file.exists(eulaFile)) {
      if (showLicense) { file.show(licenseFile) }
      response <- "y"
      if (tolower(response) == "y") {
        if (isFALSE(dir.exists(dirname(eulaFile)))) {
          dir.create(dirname(eulaFile), recursive = TRUE)
        }
        fileConn <- file(eulaFile)
        writeLines(paste(msg, paste0("# ", base::date()), "eula=true", 
                         sep = "\n"), fileConn)
        base::close(fileConn)
        return(TRUE %in% grepl("eula=true", tolower(readr::read_lines(eulaFile))))
      }
    } else { return(TRUE %in% grepl("eula=true", tolower(readr::read_lines(eulaFile)))) }
    msg <- "Yes, we accept Thermo's License agreement, get on with it!"
    cat(msg, "\n")
  }
  installRawFileReaderDLLsNoAcpt <- \(sourceUrl = rawrr:::.thermofisherlsmsUrl(), ...) { #sourceUrl = rawrr:::.thermofisherlsmsUrl()
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
    .rawfileReaderDLLs <- .rawfileReaderDLLs$objs[[1L]]
    rv <- vapply(.rawfileReaderDLLs(), \(dll) { #dll <- .rawfileReaderDLLs()[1L]
      destfile <- file.path(rawfileReaderDLLsPath, dll)
      download.file(file.path(paste0(sourceUrl, c("", "Assemblies/")[(rawrrVersions == "github")+1L],
                                     dll)), destfile = destfile, 
                    mode = "wb")
    }, 0)
    rv
  }
  installRawFileReaderDLLsNoAcpt()
  #rawrr::installRawFileReaderDLLs(sourceUrl = .thermofisherlsmsUrl())
  #requireNamespace("rawrr"); getAnywhere(".isRawFileReaderLicenseAccepted")
  rawrr::installRawrrExe()
}, silent = TRUE)
if (inherits(rawrr_tst, "try-error")) {
  warning("Could not install rawrr 1.11.14!")
}
