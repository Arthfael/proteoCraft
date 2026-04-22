#' cor_mod_seq
#' 
#' @description
#' A function to edit recent MaxQuant verbose modified sequences to the old, more succinct two-letters format (pre MaxQuant 1.6.7.0)
#' 
#' If MaxQuant is from an older version, the modified sequences will not be modified, but the list of two objects will be created nonetheless.
#' 
#' @param Pep A peptides table with a modified sequence column.
#' @param modseqcol The name of the modified sequence column. Default = "Modified sequence"
#' @param modcol The name of the modifications column. Only used if the file is actually in the old MaxQuant format. Default = "Modifications"
#' 
#' @returns
#' A list of two objects:
#' - Peptides: the peptides table with (where applicable) modified sequence in 2-letters format
#' - PTMs: a modifications table
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
  #DefArg(cor_mod_seq)
  #Pep <- ev
  stopifnot(modseqcol %in% colnames(Pep),
            modcol %in% colnames(Pep))
  Seq <- Pep[[modseqcol]]
  uSeq <- unique(Seq)
  w_uMd <- grep("\\(", uSeq)
  if (length(w_uMd)) {
    uMdSeq <- uSeq[w_uMd]
    seq1 <- strsplit(uMdSeq, "\\(|\\)")
    taest <- unique(nchar(unlist(lapply(seq1, \(x) { x[(1L:floor(length(x)/2))*2] }))))
    taest <- (length(taest) > 1L)||(taest != 2L)
    if (taest) {
      warning("Modified sequences are in the new MaxQuant format and will be converted to the old (more concise) one!")
      Pep[[paste0(modseqcol, "_verbose")]] <- Pep[[modseqcol]]
      w1 <- match(modseqcol, colnames(Pep))
      w2 <- match(paste0(modseqcol, "_verbose"), colnames(Pep))
      w0 <- which(!colnames(Pep) %in% c(modseqcol, paste0(modseqcol, "_verbose")))
      w0l <- w0[which(w0 < w1)]
      w0h <- w0[which(w0 > w1)]
      tmPep <- if (length(w0l) > 0L) { Pep[, colnames(Pep)[c(w0l, w1, w2)]] } else { Pep[, colnames(Pep)[c(w1, w2)]] }
      if (length(w0h) > 0L) { tmPep[, colnames(Pep[w0h])] <- Pep[, colnames(Pep[w0h])] }
      Pep <- tmPep
      seq2 <- strsplit(gsub("[^\\(\\)]", "", uMdSeq), "")
      seq2 <- lapply(seq2, \(x) { #x <- seq2[[1L]]
        l <- length(x)
        tst1 <- vapply(1L:l, \(y) { length(which(x[1L:y] == "(")) }, 1L)
        tst2 <- vapply(1L:l, \(y) { length(which(x[1L:y] == ")")) }, 1L)
        ends <- which(tst2-tst1 == 0)
        starts <- 1L
        if (length(ends) > 1L) { starts <- c(starts, ends[1L:(length(ends)-1L)]+1L) }
        x[which((x == "("))] <- "_[_"
        x[which((x == ")"))] <- "_]_"
        x[starts] <- "("
        x[ends] <- ")"
        return(x)
      })
      seq1 <- vapply(lapply(1L:length(seq1), \(x) { #x <- 1
        l <- length(seq1[[x]])
        r <- rep("", l*2-1L)
        r[(1L:l)*2L-1L] <- seq1[[x]]
        r[(1L:(l-1L))*2L] <- seq2[[x]]
        return(r)
      }), paste, "", collapse = "")
      seq1 <- strsplit(seq1, "\\(|\\)")
      mods <- unique(unlist(lapply(seq1, \(x) { #x <- seq1[1L]
        x[(1L:floor(length(x)/2))*2L]
      })))
      mods <- data.frame(ModName = mods,
                         `Full name` = gsub("_\\[_", "(", gsub("_\\]_", ")", mods)),
                         check.names = FALSE)
      mods$Mark <- tolower(substr(mods$"Full name", 1L, 2L))
      mods$Type <- "Variable" # MaxQuant only reports variable modifications
      ## Identify affected AA
      nr <- nrow(mods)
      mods$AA <- list(rep(c(), nr))
      for (i in 1L:nr) { #i <- 1
        m <- topattern(paste0("(", mods$"Full name"[i], ")"), start = FALSE)
        w <- grep(m, uMdSeq, value = TRUE)
        if (length(w)) {
          w <- lapply(strsplit(w, m), \(x) {
            x <- unlist(x)
            x <- x[1L:(length(x)-1L)]
            return(unique(unlist(lapply(x, \(y) { rev(unlist(strsplit(y, "")))[1L] }))))
          })
          mods$AA[i] <- list(sort(unique(unlist(w))))
        }
      }
    } else {
      message("No need to convert modified sequences as they are already in the old format, skipping.")
      mods <- data.frame(`Full name` = sort(unique(gsub("[0-9]+ ", "", unlist(strsplit(unique(Pep[which(Pep[[modcol]] != "Unmodified"), modcol]), ","))))),
                         check.names = FALSE)
      mods$Type <- "Variable"
      mods$Mark <- tolower(substr(mods$"Full name", 1L, 2L))
      nr <- nrow(mods)
      ## Identify affected AA
      mods$AA <- list(rep(c(), nr))
      for (i in 1L:nr) {
        m <- topattern(mods$"Full name"[i], start = FALSE)
        e <- grep(paste0("^", m, "$|^[0-9]+ ", m, "$"), unique(Pep[[modcol]]), value = TRUE)
        if (length(e)) {
          e <- lapply(strsplit(e, paste0("\\(", mods$Mark[i], "\\)")), \(x) {
            x <- unlist(x)
            x <- x[1L:(length(x)-1L)]
            return(unique(unlist(lapply(x, \(y) { rev(unlist(strsplit(y, "")))[1L] }))))
          })
          mods$AA[i] <- list(sort(unique(unlist(e))))
        }
      }
    }
    ## Sometimes, some marks are duplicates, e.g. if you searched for "Acetyl (Protein N-term)" and "Acetyl (K)" together!
    ## We want to fix that so that each modification has a unique mark:
    test <- aggregate(mods$Mark, list(mods$Mark), length)
    W <- which(test$x > 1L)
    nuMdSeq <- uMdSeq
    #Modifs$Mark <- Modifs$"Old mark"
    if (length(W)) {
      mods$"Old mark" <- mods$Mark
      for (i in W) {
        #i <- W[1L]
        w <- which(mods$Mark == test$Group.1[i])
        m <- mods[w,]
        m$AA[which(lengths(m$AA) == 0L)] <- "X"
        r <- if ("Acetyl (Protein N-term)" %in% m$"Full name") { which(m$"Full name" == "Acetyl (Protein N-term)") } else { 1L }
        s <- c(1L:nrow(m)); s <- s[which(s != r)]
        test <- apply(m[s, c("AA", "Mark")], 1L, \(x) { paste0(tolower(x[[1L]]), substr(x[[2L]], 1L, 1L)) })
        w1 <- which(!test %in% mods$Mark)
        m$Mark[s][w1] <- test
        w1 <- which(test %in% mods$Mark)
        if (length(w1)) {
          # not tested
          s <- s[w1]
          test <- lapply(s, list)
          kount <- 1L
          char <- c(0L:9L, letters)
          taken <- unique(c(mods$Mark, m$Mark))
          for (j in s) {
            tst <- paste0(tolower(m$AA[s]), char[kount])
            while (((tst) %in% taken)&&(kount < length(char))) {
              kount <- kount+1L
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
        test <- apply(mods[, c("AA", "Old mark")], 1L, \(x) { paste0(x[[1L]], "\\(", x[[2L]], "\\)") })
        repl <- apply(mods[, c("AA", "Mark")], 1L, \(x) { paste0(x[[1L]], "(", x[[2L]], ")") })
        for (i in w) {
          m <- test[[i]]
          r <- repl[[i]]
          for (j in 1L:length(m)) {
            test2 <- vapply(test, \(x) { m[[j]] %in% unlist(x) }, TRUE)
            if (!sum(which(test2) != i)) {
              nuMdSeq <- gsub(m[[j]], r[[j]], nuMdSeq)
            } else { stop("Don't know yet how to deal with this case!!!") }
          }
        }
      }
    }
    if (taest) {
      ## Now fix modification marks in modified sequence
      ## - this is the general case for new style modified sequences:
      for (i in 1L:nrow(mods)) {
        nuMdSeq <- gsub(topattern(paste0("(", mods$"Full name"[i], ")"), start = FALSE),
                        paste0("(", mods$Mark[i], ")"), nuMdSeq)
      }
      mods$ModName <- NULL
    }
    w_Md <- which(Seq %in% uMdSeq)
    Pep[[modseqcol]][w_Md] <- nuMdSeq[match(Seq[w_Md], uMdSeq)]
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


