options(stringsAsFactors = TRUE)

require(RCurl)

wd <- "L:/"
Sys.sleep(1L) # Required because sometimes R will attempt to set the word directory too quickly and fail.
setwd(wd)

url <- "https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/GetDIALibs"

f <- RCurl::CFILE("bfile.zip", mode = "wb")
RCurl::curlPerform(url = url, writedata = f@ref)
base::close(f)
