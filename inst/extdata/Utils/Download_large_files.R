options(stringsAsFactors = TRUE)

require(RCurl)

wd <- "L:/"
Sys.sleep(1) # Required because sometimes R will attempt to set the word directory too quickly and fail.
setwd(wd)

url <- "https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/GetDIALibs"

f <- CFILE("bfile.zip", mode = "wb")
curlPerform(url = url, writedata = f@ref)
close(f)
