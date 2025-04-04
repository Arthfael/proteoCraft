##############################################
# Calculate MS instrument time for a project #
##############################################
# This only works if:
# - files are local,
# or
# - files have been copied with a script which keeps date-created information.
# or
# - the data created time stamp has been added in the name
#
# TO DO:
# - Add detection and processing of time stamp with check against date created for inconsistencies
# - Add option to generate method file to extract method length using ScanHeadsman
#
# This lets the user select a folder, then 2 raw files in the folder (all files must be in the same folder)
# The instrument running time of all of these runs will be summed up.
#
if (!require(svDialogs)) { install.packages("svDialogs") }
require(svDialogs)

TargDir <- normalizePath(choose.dir("D:/Data/", "Select Raw files subfolder"), winslash = "/")
setwd(TargDir)
Raws <- grep("\\.raw$", list.files(TargDir, recursive = FALSE), value = TRUE)
if (!length(Raws)) { warning("Not a single Raw file found, skipping!") } else {
  Raws <- paste0(TargDir,"/", Raws)
  Raws <- data.frame(Full = Raws,
                     Name = basename(Raws),
                     Created = file.info(Raws)$ctime,
                     Modified = file.info(Raws)$mtime)
  Raws$RunLength <- difftime(Raws$Modified, Raws$Created, units = "hours")
  Raws <- Raws[order(Raws$Created, decreasing = FALSE),]
  Raw1 <- dlg_list(Raws$Name, 1, title = "Select the 1st file chronologically:")$res
  m1 <- match(Raw1, Raws$Name)
  Raw2 <- dlg_list(Raws$Name[m1:nrow(Raws)], 1, title = "Select the last file chronologically:")$res
  m2 <- match(Raw2, Raws$Name)
  CalcRaws <- m1:m2
  msg <- paste0("\n--> Folder = ", TargDir,
                "\n¯¯¯¯¯¯¯¯¯¯¯¯",
                "\n     - First file: ", Raw1,
                "\n     - Last file: ", Raw2)
  if (m2 > m1+1) {
    XcludedRaws <- dlg_list(Raws$Name[(m1+1):(m2-1)], multiple = TRUE,
                                       title = "Select any samples you wish to exclude from the calculations:")$res
    if (length(XcludedRaws)) {
      mXcl <- match(XcludedRaws, Raws$Name)
      CalcRaws <- CalcRaws[which(!CalcRaws %in% mXcl)]
      msg <- paste0(msg, "\n Excluding: ", paste(XcludedRaws, collapse = "/"))
    }
  }
  msg <- paste0(msg, "\n\n",
                "Instr. time: ", round(sum(Raws$RunLength[CalcRaws]), 1), " h\n",
                "¯¯¯¯¯¯¯¯¯¯¯¯\n",
                "Sequence start to end: ",
                round(difftime(Raws$Modified[m2], Raws$Created[m1], units = "hours"), 1)," h\n",
                "¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯\n")
  cat(msg)
}
#system(paste0("open \"", TargDir, "\""))
