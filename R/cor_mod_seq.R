#' cor_mod_seq
#' 
#' @description
#' A function to change MaxQuant modified sequences to the old format (pre MaxQuant 1.6.7.0)
#' 
#' Returns a list of two objects:
#' - Peptides: the peptides table with fixed modified sequence (if applicable)
#' - PTMs: the modifications table
#' 
#' If MaxQuant is from an older version, the modified sequences will not be modified, but the list of two objects will be created nonetheless.
#' 
#' @param Pep A peptides table with a modified sequence column.
#' @param modseqcol The name of the modified sequence column. Default = "Modified sequence"
#' @param modcol The name of the modifications column. Only used if the file is actually in the old MaxQuant format. Default = "Modifications"
#' 
#' @examples
#' temp <- cor_mod_seq(ev)
#' # Check results:
#' View(temp$Peptides)
#' View(temp$PTMs)
#' # Assign results to discrete objects:
#' ev <- temp$Peptides
#' Modifs <- temp$PTMs
#' 
#' @export

cor_mod_seq <- function(Pep,
                        modseqcol = "Modified sequence",
                        modcol = "Modifications") {
  #proteoCraft::DefArg(proteoCraft::cor_mod_seq)
  #Pep <- ev
  stopifnot(modseqcol %in% colnames(Pep),
            modcol %in% colnames(Pep))
  seq <- Pep[[modseqcol]]
  witch <- grep("\\(", seq)
  if (length(witch)) {
    seq <- seq[witch]
    seq1 <- strsplit(seq, "\\(|\\)")
    taest <- unique(nchar(unlist(sapply(seq1, function(x) { x[(1:floor(length(x)/2))*2] }))))
    taest <- (length(taest) > 1)||(taest != 2)
    if (taest) {
      warning("Modified sequences are in the new MaxQuant format and will be converted to the old (more concise) one!")
      Pep[[paste0(modseqcol, "_verbose")]] <- Pep[[modseqcol]]
      w1 <- which(colnames(Pep) == modseqcol)
      w2 <- which(colnames(Pep) == paste0(modseqcol, "_verbose"))
      w0 <- which(!colnames(Pep) %in% c(modseqcol, paste0(modseqcol, "_verbose")))
      w0l <- w0[which(w0 < w1)]
      w0h <- w0[which(w0 > w1)]
      if (length(w0l) > 0) { temp <- Pep[, colnames(Pep)[c(w0l, w1, w2)]] } else { temp <- Pep[, colnames(Pep)[c(w1, w2)]] }
      if (length(w0h) > 0) { temp[, colnames(Pep[w0h])] <- Pep[, colnames(Pep[w0h])] }
      Pep <- temp
      seq2 <- strsplit(gsub("[^\\(\\)]", "", seq), "")
      seq2 <- sapply(seq2, function(x) { #x <- seq2[[1]]
        tst1 <- sapply(1:length(x), function(y) { length(which(x[1:y] == "(")) })
        tst2 <- sapply(1:length(x), function(y) { length(which(x[1:y] == ")")) })
        ends <- which(tst2-tst1 == 0)
        starts <- 1
        if (length(ends) > 1) {  starts <- c(starts, ends[1:(length(ends)-1)]+1) }
        x[which((x == "("))] <- "_[_"
        x[which((x == ")"))] <- "_]_"
        x[starts] <- "("
        x[ends] <- ")"
        return(x)
      })
      seq1 <- sapply(sapply(1:length(seq1), function(x) { #x <- 1
        l <- length(seq1[[x]])
        r <- rep("", l*2-1)
        r[(1:l)*2-1] <- seq1[[x]]
        r[(1:(l-1))*2] <- seq2[[x]]
        return(r)
      }), paste, collapse = "")
      seq1 <- strsplit(seq1, "\\(|\\)")
      mods <- unique(unlist(sapply(seq1, function(x) { #x <- seq1[1]
        x[(1:floor(length(x)/2))*2]
      })))
      mods <- data.frame(ModName = mods, `Full name` = gsub("_\\[_", "(", gsub("_\\]_", ")", mods)), check.names = FALSE)
      mods$Mark <- tolower(substr(mods$"Full name", 1, 2))
      mods$Type <- "Variable" # MaxQuant only reports variable modifications
      ## Identify affected AA
      mods$AA <- list(rep(c(), nrow(mods)))
      for (i in 1:nrow(mods)) { #i <- 1
        m <- proteoCraft::topattern(paste0("(", mods$"Full name"[i], ")"), start = FALSE)
        w <- grep(m, Pep[[modseqcol]], value = TRUE)
        if (length(w) > 0) {
          w <- sapply(strsplit(w, m), function(x) {
            x <- unlist(x)
            x <- x[1:(length(x)-1)]
            return(unique(unlist(sapply(x, function(y) { rev(unlist(strsplit(y, "")))[1] }))))
          })
          mods$AA[i] <- list(sort(unique(unlist(w))))
        }
      }
    } else {
      message("No need to convert modified sequences as they are already in the old format, skipping.")
      mods <- data.frame(`Full name` = sort(unique(gsub("[0-9]+ ", "", unlist(sapply(strsplit(Pep[which(Pep[[modcol]] != "Unmodified"), modcol], ","), unlist))))),
                         check.names = FALSE)
      mods$Type <- "Variable"
      mods$Mark <- tolower(substr(mods$"Full name", 1, 2))
      ## Identify affected AA
      mods$AA <- list(rep(c(), nrow(mods)))
      for (i in 1:nrow(mods)) {
        m <- proteoCraft::topattern(mods$"Full name"[i], start = FALSE)
        W <- grep(paste0("^", m, "$|^[0-9]+ ", m, "$"), Pep[[modcol]])
        if (length(W) > 0) {
          e <- Pep[W, modseqcol]
          e <- sapply(strsplit(e, paste0("\\(", mods$Mark[i], "\\)")), function(x) {
            x <- unlist(x)
            x <- x[1:(length(x)-1)]
            return(unique(unlist(sapply(x, function(y) { rev(unlist(strsplit(y, "")))[1] }))))
          })
          mods$AA[i] <- list(sort(unique(unlist(e))))
        }
      }
    }
    ## Sometimes, some marks are duplicates, e.g. if you searched for "Acetyl (Protein N-term)" and "Acetyl (K)" together!
    ## We want to fix that so that each modification has a unique mark:
    test <- aggregate(mods$Mark, list(mods$Mark), length)
    W <- which(test$x > 1)
    #Modifs$Mark <- Modifs$"Old mark"
    if (length(W) > 0) {
      mods$"Old mark" <- mods$Mark
      for (i in W) {
        #i <- W[1]
        w <- which(mods$Mark == test$Group.1[i])
        m <- mods[w,]
        m$AA[which(sapply(m$AA, length) == 0)] <- "X"
        if ("Acetyl (Protein N-term)" %in% m$"Full name") { r <- which(m$"Full name" == "Acetyl (Protein N-term)") } else { r <- 1 }
        s <- c(1:nrow(m)); s <- s[which(s != r)]
        test <- apply(m[s, c("AA", "Mark")], 1, function(x) { paste0(tolower(x[[1]]), substr(x[[2]], 1, 1)) })
        w1 <- which(!test %in% mods$Mark)
        m$Mark[s][w1] <- test
        w1 <- which(test %in% mods$Mark)
        if (length(w1) > 0) {
          # not tested
          s <- s[w1]
          test <- sapply(s, list)
          kount <- 1
          char <- c(c(0:9), letters)
          taken <- unique(c(mods$Mark, m$Mark))
          for (j in s) {
            tst <- paste0(tolower(m$AA[s]), char[kount])
            while (((tst) %in% taken)&&(kount < length(char))) {
              kount <- kount+1
              tst <- paste0(tolower(m$AA[s]), char[kount])
            }
            if (kount == length(char)) {
              stop("I am really out of options here, never thought this would go this far! Check the code just in case, this should never happen.")
            } else {
              test[[j]] <- tst
            }
          }
          m$Mark[[s]] <- unlist(test)
        }
        mods[w,] <- m
      }
      if (!taest) {
        ## Now fix modification marks in modified sequence:
        ## - this is for the rare case where we have old style modified sequences but some duplicate marks:
        w <- which(mods$"Old mark" != mods$Mark)
        test <- apply(mods[,c("AA", "Old mark")], 1, function(x) {paste0(x[[1]], "\\(", x[[2]], "\\)")})
        repl <- apply(mods[,c("AA", "Mark")], 1, function(x) {paste0(x[[1]], "(", x[[2]], ")")})
        for (i in w) {
          m <- test[[i]]
          r <- repl[[i]]
          for (j in 1:length(m)) {
            test2 <- sapply(test, function(x) {m[[j]] %in% unlist(x)})
            if (sum(which(test2) != i) == 0) {
              seq <- gsub(m[[j]], r[[j]], seq)
            } else { stop("Don't know yet how to deal with this case!!!") }
          }
        }
      }
    }
    if (taest) {
      ## Now fix modification marks in modified sequence
      ## - this is the general case for new style modified sequences:
      for (i in 1:nrow(mods)) {
        seq <- gsub(proteoCraft::topattern(paste0("(", mods$"Full name"[i], ")"), start = FALSE),
                    paste0("(", mods$Mark[i], ")"), seq)
      }
      mods$ModName <- NULL
    }
    Pep[witch, modseqcol] <- seq
  } else {
    warning("Not a single modified peptide was identified, is this normal?")
    mods <- data.frame(`Full name` = NA,
                       Type = NA,
                       Mark = NA,
                       AA = NA,
                       check.names = FALSE)
    mods <- mods[-1,]
  }
  return(list(Peptides = Pep,
              PTMs = mods))
}


