#' MQ.load
#'
#' @description 
#' A function to quickly import MaxQuant's main output datatables into the above environment.
#' 
#' @param ev Default = TRUE. If set to FALSE, no evidence file will be created.
#' @param ev.name The name of the evidence object that should be created in the above environment, if "assign" = TRUE.
#' @param ev.file.name The name of the evidence file that should be loaded.
#' @param pep Default = TRUE. If set to FALSE, no peptide file will be created.
#' @param pep.name The name of the peptides object that should be created in the above environment, if "assign" = TRUE.
#' @param pep.file.name The name of the peptides file that should be loaded.
#' @param prot Default = TRUE. If set to FALSE, no protein groups file will be created.
#' @param prot.name The name of the protein groups object that should be created in the above environment, if "assign" = TRUE.
#' @param prot.file.name The name of the protein groups file that should be loaded.
#' @param sub Are the files in a subfolder? Default = FALSE. If set to TRUE, will search in subfolders automatically, and will return an error if it finds more than one or no matches; can instead be set to any string specifying a subfolder.
#' @param assign Should we assign the files (defaut = TRUE) to a variable in the parent frame? Use arguments "ev.name", "pep.name" and "prot.name" to define the variable names.
#' @param return Should we return the files (in case the function is used in the form "object <- MQ.load(...)"); default = FALSE.
#' @param check.names logical (default = FALSE). If TRUE, then the names of the variables in the data frame are checked to ensure that they are syntactically valid variable names. If necessary they are adjusted (by make.names) so that they are, and also to ensure that there are no duplicates. 
#'
#' @examples
#' MQ.load()
#' Result: the environment should now contain ev, pep and prot objects.
#' 
#' @export

MQ.load <- function(ev = TRUE, ev.name = "ev", ev.file.name = "evidence.txt",
                    pep = TRUE, pep.name = "pep", pep.file.name = "peptides.txt",
                    prot = TRUE, prot.name = "PG", prot.file.name = "proteinGroups.txt",
                    sub = FALSE, assign = TRUE, envir = .GlobalEnv, return = FALSE, check.names = FALSE) {
  #proteoCraft::DefArg(proteoCraft::MQ.load)
  if (sum(c(assign, return)) == 0) {
    stop("Since both \"assign\" and \"return\" are set to FALSE, there is nothing for me to do!")
  }
  # When I have time, I can rewrite this so the data.frame below is used to avoid repeating so much code:
  #df <- data.frame(n1 = c("ev", "pep", "prot"),
  #                 n2 = c("Ev", "Pep", "Prot"),
  #                 n3 = c("evidences", "peptides", "protein groups"))
  #
  # Check or process names
  if ((is.logical(sub))&&(sub)) {
    a <- list.files(getwd(), recursive = TRUE, include.dirs = TRUE)
    if (ev) {
      ev.file.name <- grep(paste(paste0(c("^", "/"), gsub("\\.", "\\\\.", ev.file.name)), collapse = "|"), a, value = TRUE)
      l <- length(ev.file.name)
      if (l != 1) {
        ev <- FALSE
        warning(c("I could not find any file with that name, skipping loading evidences!",
               "I can only load a single evidence file at a time, skipping loading evidences!")[which(c(l == 0, l > 1))])
      }
    }
    if (pep) {
      pep.file.name <- grep(paste(paste0(c("^", "/"), gsub("\\.", "\\\\.", pep.file.name)), collapse = "|"), a, value = TRUE)
      
      l <- length(pep.file.name)
      if (l != 1) {
        pep <- FALSE
        warning(c("I could not find any file with that name, skipping loading peptides!",
               "I can only load a single peptides file at a time, skipping loading peptides!")[which(c(l == 0, l > 1))])
      }
    }
    if (prot) {
      prot.file.name <- grep(paste(paste0(c("^", "/"), gsub("\\.", "\\\\.", prot.file.name)), collapse = "|"), a, value = TRUE)
      l <- length(prot.file.name)
      if (l != 1) {
        pep <- FALSE
        warning(c("I could not find any file with that name, skipping loading protein groups!",
               "I can only load a single protein groups file at a time, skipping loading protein groups!")[which(c(l == 0, l > 1))])
      }
    }
  } else {
    if (ev) {
      if ((!is.logical(sub))||(sub)) { ev.file.name <- gsub("//","/", paste(sub, ev.file.name, sep = "/")) }
      if (!file.exists(ev.file.name)) {
        warning("I could not find any file with that name, skipping loading evidences!")
        ev <- FALSE
      }
    }
    if (pep) {
      if ((!is.logical(sub))||(sub)) { pep.file.name <- gsub("//","/", paste(sub, pep.file.name, sep = "/")) }
      if (!file.exists(pep.file.name)) {
        warning("I could not find any file with that name, skipping loading peptides!")
        pep <- FALSE
      }
    }
    if (prot) {
      if ((!is.logical(sub))||(sub)) { prot.file.name <- gsub("//","/",paste(sub, prot.file.name, sep = "/")) }
      if (!file.exists(prot.file.name)) {
        warning("I could not find any file with that name, skipping loading protein groups!")
        prot <- FALSE
      }
    }
  }
  res <- list()
  if (ev) {
    Ev <- data.table::fread(ev.file.name,
                            sep = "\t",
                            integer64 = "numeric",
                            check.names = FALSE,
                            data.table = FALSE)
    if (assign) { assign(ev.name, Ev, pos = parent.frame()) }
    res$evidences <- Ev
  }
  if (pep) {
    Pep <- data.table::fread(pep.file.name,
                             sep = "\t",
                             integer64 = "numeric",
                             check.names = FALSE,
                             data.table = FALSE)
    if (assign) { assign(pep.name, Pep, pos = parent.frame()) }
    res$peptides <- Pep
  }
  if (prot) {
    Prot <- data.table::fread(prot.file.name,
                              sep = "\t",
                              integer64 = "numeric",
                              check.names = FALSE,
                              data.table = FALSE)
    if (assign) { assign(prot.name, Prot, pos = parent.frame()) }
    res$protein.groups <- Prot
  }
  if ((return)&&(sum(c(ev, pep, prot)) > 0)) {
    return(res)
  }
}
