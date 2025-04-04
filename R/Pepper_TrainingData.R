#' Pepper_TrainingData
#'
#' @description
#' A function to take an ev data.frame and convert it to a Pepper-compatible table (which can then be one-hot encoded for the purpose
#' of training a model on it).
#' This is just a wrapper to Pepper_ProcessData
#' 
#' 
#' Note:
#' An issue may be silent fixed modifications (TMT, cysteine alkylations) in MaxQuant. The final model should thus specify
#' what type of data it is suitable for.
#' 
#' @param Ev The evidences (PSMs) data frame.
#' @param Modifs Modifications table as output by DIANN_to_MQ, FP_to_MQ, PD_to_MQ... Should ideally contain a "Mass shift" column, but for MQ workflows this will be usually missing: two additional arguments are provided to deal with those cases: SearchSoft and MQFold.
#' @param experiments.map The experiments map.
#' @param Ev_Sample.col Name of the sample column in Ev
#' @param experiments.map_Sample.col Name of the sample column in experiments.map
#' @param path The path to which to save the file.
#' @param SearchSoft Search software used. Important because Modifs may miss the "Mass shift" column for data originating from some software (MaxQuant). If it is missing, and software is MaxQuant, then providing the next argument (MQFold) can allow parsing directly modifications and calculating the associated mass shift. 
#' @param MQFold Main MaxQuant folder, containing /bin/conf subfolder with modification xml(s). Only used if SearchSoft == "MaxQuant"
#' @param intCol Which intensity column to use? Default = "Intensity"
#' 
#' @import data.table
#' @export

Pepper_TrainingData <- function(Ev,
                                Modifs,
                                experiments.map,
                                Ev_Sample.col,
                                experiments.map_Sample.col,
                                path = "",
                                SearchSoft = SearchSoft,
                                MQFold,
                                intCol) {
  proteoCraft::Pepper_ProcessData(Ev, 
                            Modifs,
                            experiments.map,
                            Ev_Sample.col,
                            experiments.map_Sample.col,
                            path,
                            SearchSoft,
                            MQFold,
                            intCol,
                            filter = TRUE)
}

