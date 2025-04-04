#
if (!require(R.utils)) { install.packages("R.utils") }
library(R.utils)

# Set up virtual environment using reticulate
#Sys.setenv(RETICULATE_PYTHON = "C:\\Program Files\\")
if (!require(reticulate)) { install.packages("reticulate") }
library(reticulate) # Note! We use reticulate mostly for a few convenient functions,
# but after that we run Python through command line!
# Using reticulate has these issues:
#   - For us at least it fails in some cases where command line works.
#   - One cannot deactivate a virtualenv through reticulate the way it can be done through command line.
#   - Also it does not seem possible with it to pass arguments to a python script (at least easily).
#
myHome <- paste0(shell("echo %homedrive%", intern = TRUE),
                 shell("echo %homepath%", intern = TRUE))
#cat(myHome)
pyEnv <- "Bleh" # Virtual environment
pyPacks <- c("numpy",
             "pandas",
             "argparse",
             "matplotlib",
             "seaborn",
             "tensorflow",
             "keras",
             "scikit-learn",
             "sklearn",
             "scipy") # dependencies
tmp <- list.files("C:/Program files", "python\\.exe$", recursive = TRUE, full.names = TRUE)
tmp <- grep("/Lib/venv/", tmp, value = TRUE, invert = TRUE)
tmp <- normalizePath(tmp)
tmp <- gsub("\\\\Program Files\\\\", "\\\\PROGRA~1\\\\", tmp)
tmp2 <- data.frame(type = "PythonCore", hive = "HLM",
                   install_path = gsub("\\\\python\\.exe$", "\\\\", tmp),
                   executable_path = gsub("\\\\python\\.exe$", "\\\\\\\\python.exe", tmp),
                   arch = "x64")
tmp2$version <- sapply(paste0(tmp2$install_path, "NEWS.txt"), function(x) { #x <- paste0(tmp2$install_path, "NEWS.txt")[1]
  if (file.exists(x)) {
    x <- readLines(x)
    x <- grep("What's New in Python ", x, value = TRUE)[1]
    x <- gsub(" .+", "", gsub("What's New in Python ", "", x))
  } else { x <- NA }
  return(x)
})
pyPaths <- py_versions_windows()
pyPaths <- rbind(pyPaths, tmp2)
pyFound <- FALSE
if (!nrow(pyPaths)) {
  stop("A valid installation of python is required!")
  # Below:
  # Code for installing python from reticulate... which currently fails!
  #vers <- rev(grep("-win32$", install_python(list = TRUE), invert = TRUE, value = TRUE))[1]
  #vers <- rev(grep("-win32$", grep("^3\\.10\\.", install_python(list = TRUE), value = TRUE), invert = TRUE, value = TRUE))[1]
  #vers <- "3.10.5"
  #vers <- "3.10.4"
  #pyPaths <- install_python(vers, force = TRUE)
  #pyFound <- TRUE
} else {
  if (nrow(pyPaths) > 1) {
    w <- c(which(pyPaths$type == "PythonCore"),
           which(pyPaths$type == "Anaconda"),
           which(!pyPaths$type %in% c("PythonCore", "Anaconda")))
    pyPaths <- pyPaths[w,]
    tst <- strsplit(pyPaths$version, "\\.")
    mx <- max(sapply(tst, length))
    tst <- as.data.frame(t(sapply(tst, function(x) {
      l <- length(x)
      if (l < mx) { x <- c(x, rep(NA, mx-l)) }
      return(as.integer(x))
    })))
    pyPaths[, c("V1", "V2", "V3")] <- tst
    # Currently, tensorflow is not supported for Python > 3.10 
    w <- which((pyPaths$V1 < 3)|((pyPaths$V1 == 3)&(pyPaths$V2 <= 10)))
    pyPaths <- pyPaths[w,]
    if (length(w)) { pyPaths <- pyPaths[order(pyPaths$V1, pyPaths$V2, pyPaths$V3, decreasing = TRUE),] }
  }
  if (nrow(pyPaths)) {
    k <- 0
    while ((k == 0)||("try-error" %in% class(tst))) {
      k <- k+1
      pyPath <- pyPaths$executable_path[k]
      tst <- try(use_python(pyPath), silent = TRUE) # Wrapped in try
    }
    pyFound <- k <= nrow(pyPaths)
  }
}
stopifnot(pyFound)
pyPath <- gsub("\\\\+", "\\\\", pyPath)
#
cmd <- paste0("\"", pyPath, "\" -m pip install --upgrade pip")
#cat(cmd)
system(cmd)
cmd <- paste0("\"", pyPath, "\" -m pip install virtualenv")
#cat(cmd)
system(cmd)
for (pyPack in pyPacks) {
  cmd <- paste0("\"", pyPath, "\" -m pip install ", pyPack)
  #cat(cmd)
  system(cmd)
  cat("\n")
}
#
# Remove virtual environment
pyEnvPath <- paste0(myHome, gsub("~/", "\\\\Documents\\\\", virtualenv_root()), "\\", pyEnv)
cmd <- paste0("\"", pyPath, "\" -m virtualenv --python=\"", pyPath,"\" \"", pyEnvPath, "\"")
cat(cmd)
system(cmd)
# Load venv and modules
#cmd <- paste0(myHome, gsub("~/", "\\\\Documents\\\\", virtualenv_root()), "\\", pyEnv, "\\Scripts\\deactivate.bat")
#system(cmd)
cmd <- paste0(myHome, gsub("~/", "\\\\Documents\\\\", virtualenv_root()), "\\", pyEnv, "\\Scripts\\activate")
#cat(cmd)
system(cmd)

# Important!
############
# I added "scikit-learn" to pyPacks, the vector of packages to install,
# even though it is not explicitly loaded in the python script.
# Simply speaking, sklearn won't work without it.
# For why, see https://towardsdatascience.com/fix-modulenotfounderror-sklearn-99db60ed5d

