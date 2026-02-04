#' AdvNorm.IL
#' 
#' @description 
#' This function performs simplified "advanced normalization" on the columns of a data.frame.
#' 
#' @param df The input data frame, containing unique peptide identifiers and expression data.
#' @param ids.col The name of the column in df containing peptide identifiers. Default is "Unique State", a column I usually create which contains modified sequence pasted to charge.
#' @param exprs.col The name of the columns in df which contain expression values.
#' @param exprs.log Default = FALSE, meaning the input data is not log-transformed. Set instead to the base of the log if input data is log-transformed. If set to TRUE, the data is assumed to be log10 base. The data returned will be transformed or not as per the input.
#' @param K A fold factor defining the minimum (1/K) and maximum (*K) range of accepted values for correction factors.
#' 
#' @details
#' This normalizes the data using the Levenberg Marquard procedure (aka. "advanced normalization") to minimize column-wise the sum square difference between different quantitative vectors.
#' The ".IL" suffix in the name stands for "isobaric labeling", because at first this was written with a TMT or iTRAQ reporter intensities data frame in mind.
#' However, this can be used for any data frame of quantitative values.
#' The input data.frame df must contain a column of IDs and 2 or more expression columns.
#' This is slower but is better than median only.
#' 
#' @returns
#' Normalized data. New column names are recycled from the original ones (in exprs.col), with "AdvNorm. " added as a prefix. 
#' 
#' @examples
#' adv.norm.data <- AdvNorm.IL(data, "Unique State", paste0("Reporter intensity ", c(0:9)), FALSE, K = 5)
#' require(ggplot2)
#' require(reshape2)
#' test <- melt(adv.norm.data)
#' test$value <- log10(test$value)
#' plot <- ggplot(test) + geom_density(stat = "density", aes(x = value, colour = variable))
#' grDevices::windows(width = 10, height = 10)
#' print(plot)
#' 
#' @import data.table
#' @export

AdvNorm.IL <- function(df,
                       ids.col = "Unique State",
                       exprs.col,
                       exprs.log = FALSE,
                       K = 5) {
  # NB:
  # The function does not lend itself well to parallelization. It is much easier to parallelize different calls to the function.
  TESTING <- FALSE
  #proteoCraft::DefArg(proteoCraft::AdvNorm.IL); TESTING <- TRUE
  #df = temp[, c("id", g)]; ids.col = "id"; exprs.col = g; exprs.log = TRUE
  #
  if (TESTING) {
    tm1 <<- Sys.time()
  }
  if (length(exprs.col) == 1) { stop("There must be at least 2 columns to normalise the data! Check the \"exprs.col\" argument.") }
  if ((!is.logical(exprs.log))||(exprs.log != FALSE)) {
    if (is.logical(exprs.log)) { exprs.log.base <- 10 } else {
      exprs.log.base <- exprs.log
      exprs.log <- TRUE
    }
  }
  #
  # Convert to linear if required:
  if (exprs.log) { for (klnm in exprs.col) { df[[klnm]] <- exprs.log.base^df[[klnm]] } }
  # Original medians
  ## Global
  M1 <- proteoCraft::is.all.good(unlist(df[, exprs.col]))
  M1 <- median(M1[which(M1 > 0)])
  # Per column
  m1 <- apply(df[, exprs.col], 2, function(x) {
    x <- proteoCraft::is.all.good(unlist(x))
    return(median(x[which(x > 0)]))
  })
  df2 <- df
  df2[, exprs.col] <- sweep(df2[, exprs.col], 2, m1, "/")
  #
  n.exprs <- length(exprs.col) # Number of columns
  # Starting factors
  h.fact <- setNames(lapply(seq_len(n.exprs), function(x) { 1 }),
                     paste0("exprs.", seq_len(n.exprs)))
  #
  # Create main optimization function:
  N <- length(exprs.col)
  comb <- gtools::combinations(N, 2, exprs.col)
  diffLog2 <- function(p, nms, dat, append = TRUE) {
    # There is also proteoCraft::diffLog! Can I replace one by the other?
    p <- c(1, unlist(p))
    #stopifnot(length(p) == length(nms))
    dat <- dat[, c(ids.col, nms)]
    dat[, nms] <- sweep(dat[, nms], 2, p, "/")
    dat2 <- data.table::as.data.table(dat)
    colnames(dat2)[which(colnames(dat2) == ids.col)] <- "Group.1"
    # Aggregate (if necessary) by ID, using sum:
    dat2 <- dat2[, lapply(.SD, sum, na.rm = TRUE), by = Group.1, , .SDcols = nms]
    #dat2 <- as.data.frame(dat2); dat3 <- as.matrix(dat2[, nms])
    dat3 <- as.matrix(dat2[, ..nms]) # As matrix...
    dat3 <- log10(dat3)
    dat3 <- dat3[, comb[,1]] - dat3[, comb[,2]]
    dat3 <- unlist(dat3)
    dat3 <- dat3[which(!is.na(dat3))]
    dat3 <- dat3[which(is.finite(dat3))]
    return(dat3)
  }
  diffLog2_h <- function(...) {
    p <- list(...)
    res <- diffLog2(p, exprs.col, df2)
    return(res)
  }
  LM <- minpack.lm::nls.lm(par = h.fact[2:n.exprs],
                           fn = diffLog2_h,
                           lower = unlist(h.fact[2:n.exprs])/K,
                           upper = unlist(h.fact[2:n.exprs])*K)
  h <- c(1, as.numeric(LM$par))
  h <- h/(prod(h)^(1/n.exprs)) # Re-centre, since the 1st value was kept at 1
  df2[, exprs.col] <- sweep(df2[, exprs.col], 2, h, "/")
  # Median after
  M2 <- proteoCraft::is.all.good(unlist(df2[, exprs.col]))
  M2 <- median(M2[which(M2 > 0)])
  # Re-apply original values scaling:
  df2[, exprs.col] <- df2[, exprs.col]*M1/M2
  # Apply to expression data
  if (exprs.log) { for (i in exprs.col) { df2[[i]] <- log(df2[[i]], exprs.log.base) } }
  w <- which(colnames(df2) %in% exprs.col)
  colnames(df2)[w] <- paste0("AdvNorm.", colnames(df2)[w])
  if (TESTING) {
    tm2 <<- Sys.time()
    print(tm2-tm1)
  } else {
    return(df2)
  }
}
