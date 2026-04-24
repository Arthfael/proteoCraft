# This may evolve into a full function
require(reshape2)
require(gtools)
require(ggplot2)
require(svDialogs)
require(proteoCraft)

RPath <- as.data.frame(library()$results)
RPath <- normalizePath(RPath$LibPath[match("proteoCraft", RPath$Package)], winslash = "/")
libPath <- paste0(RPath, "/proteoCraft")
homePath <- paste0(normalizePath(Sys.getenv("HOME"), winslash = "/"), "/R/proteoCraft")
dfltLocsFl <- paste0(homePath, "/Default_locations.xlsx")
dfltLocs <- openxlsx2::read_xlsx(dfltLocsFl)
searchDir <- dfltLocs$Path[match("Search folder", dfltLocs$Folder)]

# Choose working directory
wd <- rstudioapi::selectDirectory(path = searchDir)
setwd(wd)

#Seq <- "MEVRNMVDYELLKKVVEAPGVSGYEFLGIRDVVIEEIKDYVDEVKVDKLGNVIAHKKGEGPKVMIAAHMDQIGLMVTHIEKNGFLRVAPIGGVDPKTLIAQRFKVWIDKGKFIYGVGASVPPHIQKPEDRKKAPDWDQIFIDIGAESKEEAEDMGVKIGTVITWDGRLERLGKHRFVSIAFDDRIAVYTILEVAKQLKDAKADVYFVATVQEEVGLRGARTSAFGIEPDYGFAIDVTIAADIPGTPEHKQVTHLGKGTAIKIMDRSVICHPTIVRWLEELAKKHEIPYQLEILLGGGTDAGAIHLTKAGVPTGALSVPARYIHSNTEVVDERDVDATVELMTKALENIHELKI"
Seq <- "MGNIFANLFKGLFGKKEMRILMVGLDAAGKTTILYKLKLGEIVTTIPTIGFNVETVEYKNISFTVWDVGGQDKIRPLWRHYFQNTQGLIFVVDSNDRERVNEAREELMRMLAEDELRDAVLLVFANKQDLPNAMNAAEITDKLGLHSLRHRNWYIQATCATSGDGLYEGLDWLSNQLRNQK"

# Load isotopic frequencies table
# Add to this table each time I need a new element
# e.g.:
temp <- c("F", "18.9984031629_1")
if (!temp[1] %in% IsotopeProbs$Atom) {
  L <- length(temp)
  Nc <- ncol(IsotopeProbs)
  if (L < Nc) { temp <- c(temp, rep("", Nc-L)) }
  if (L > Nc) {
    tst <- grep("^N\\+[0-9]+", colnames(IsotopeProbs), value = TRUE)
    tst <- max(as.numeric(gsub("^N\\+", "", tst)))
    for (i in L-Nc) { IsotopeProbs[[paste0("N+", tst+i)]] <- "" }
  }
  IsotopeProbs <- rbind(IsotopeProbs, temp)
  
}

# Option: remove N-terminal methionine
if (grepl("^M", Seq)) {
  msg <- "The sequence starts with a Methionine, should we remove it?"
  RemoveNTermMet <- c(TRUE, FALSE)[match(dlg_message(msg, "yesno")$res, c("yes", "no"))]
  if (RemoveNTermMet) { Seq <- gsub("^M", "", Seq) }
}

# Check amino acid counts
AACounts <- unlist(strsplit(Seq, ""))
AACounts <- aggregate(AACounts, list(AACounts), length)
colnames(AACounts) <- c("AA", "Count")
AACounts <- AACounts[order(AACounts$AA),]
AACounts

# Choose type of mass:
alkTbl <- data.frame(Alk = c("IAM", "NEM", "CAM", "MMTS", "iST_NHS"),
                     Full_name = c("Iodoacetamide", "N-Ethylmaleimide", "Chloroacetamide", "Methyl methanethiosulfonate", "Proprietory alkylating agent from iST-NHS kit"),
                     C = c(2, 6, 2, 1, 6),
                     H = c(3, 7, 3, 2, 11),
                     N = c(1, 1, 1, 0, 1),
                     O = c(1, 2, 1, 0, 1),
                     S = c(0, 0, 0, 1, 0))

Native <- FALSE
Alkylator <- FALSE
nCfree <- 0
if (grepl("C", Seq)) {
  msg <- "Is the protein denatured?"
  Native <- !c(TRUE, FALSE)[match(dlg_message(msg, "yesno")$res, c("yes", "no"))]
  if (Native) {
    nC <- AACounts$Count[match("C", AACounts$AA)]
    if (nC > 1) {
      nIntra2Smax <- floor(nC/2)
      msg <- paste0("There are ", nC, " cysteines in the protein, how many intra-molecular disulfide bonds (1-",
                    nIntra2Smax, ") do you expect?")
      nIntra2S <- as.integer(dlg_input(msg, n2Smax)$res)
      while((!is.integer(nIntra2S))||(!nIntra2S %in% 1:n2Smax)) { nIntra2S <- as.integer(dlg_input(msg, nIntra2Smax)$res) }
    } else { nIntra2S <- 0 }
    nCfree <- nC-2*nIntra2S
    if (nCfree) {
      msg <- paste0("There ", c("is", "are still")[(nC > 1)+1], " " , nCfree, " free cysteine", c("", "s")[(nC > 1)+1],
                    " in the protein, how many inter-molecular disulfide bonds (0-", nCfree, ") do you expect?")
      nInter2S <- as.integer(dlg_input(msg, 0)$res)
      while((!is.integer(nInter2S))||(!nInter2S %in% 0:nCfree)) { nInter2S <- as.integer(dlg_input(msg, 0)$res) }
    } else { nInter2S <- 0 }
    nCfree  <- nC-2*nIntra2S-nInter2S
  } else { nCfree <- nC }
  if (nCfree) {
    msg <- paste0("Did you alkylate the protein?\n - 0: no",
                  paste0("\n - ", 1:nrow(alkTbl), ": ", alkTbl$Full_name, collapse = ""), "\n")
    defAlk <- c(match("iodoacetamide", tolower(alkTbl$Full_name)), 0)[Native+1]
    whichalk <- as.integer(dlg_input(msg, defAlk)$res)
    while ((!is.integer(whichalk))||(!whichalk %in% 0:nrow(alkTbl))) { whichalk <- as.integer(dlg_input(msg, defAlk)$res) }
    if (whichalk > 0) { Alkylator <- alkTbl$Alk[whichalk] } else { Alkylator <- FALSE }
  } else { Alkylator <- FALSE }
}

# Note:
# Up to now I have always used IAA as an abbreviation for Iodoacetamide.
# However, I should use IAM instead, as it seems Iodoacetic acid is often meant under IAA.
# Calculate elemental formula
ProtFormula <- function(seq, collapse = TRUE, freeCyst = nCfree, alkylated = Alkylator, alktbl = alkTbl) {
  seq <- unlist(strsplit(toupper(seq), ""))
  wC1 <- which(seq == "C")
  freeCyst <- as.integer(freeCyst)
  if ((is.na(freeCyst))||(freeCyst < 0)) {
    warning("Ignoring invalid \"freeCyst\" argument!")
    freeCyst <- 0
  }
  if (freeCyst > length(wC1)) {
    warning("Argument \"freeCyst\" cannot exceed the number of cysteines in the sequence!")
    freeCyst <- length(wC1)
  }
  nC <- which(seq == "C")
  L <- length(seq)
  atms <- c("C", "H", "N", "O", "S", "Se")
  atmsMonoMasses <- setNames(as.numeric(gsub("_.*", "", IsotopeProbs$Monoisotopic[match(atms, IsotopeProbs$Atom)])), atms)
  Tbl <- proteoCraft::AA_table
  Tbl$State <- "standard"
  Tbl$Tag <- Tbl$AA
  wC <- which(Tbl$Name == "Cysteine")
  oxC <- as.data.frame(t(c("C", "oxidized Cysteine", "Oxidized", "C(ox)")))
  colnames(oxC) <- c("AA", "Name", "State", "Tag")
  for (atm in atms) {
    oxC[[atm]] <- Tbl[wC, atm]
    if (atm == "H") { oxC[[atm]] <- oxC[[atm]]-1 }
  }
  oxC$"Monoisotopic mass" <- sapply(1:nrow(oxC), function(x) {
    sum(oxC[x, atms]*atmsMonoMasses)
  })
  oxC$"Residue monoisotopic mass" <- oxC$"Monoisotopic mass" - sum(c(2, 1)*atmsMonoMasses[c("H", "O")])
  Tbl <- rbind(Tbl, oxC)
  alk <- (!is.logical(alkylated))&&(!is.na(alkylated))&&(nchar(as.character(alkylated)) > 0)
  if ((nCfree < nC)&&(alk)) {
    warning("Did you really alkylate but did not reduce? Ignoring alkylations.")
    alk <- FALSE
  }
  if ((alk)&&(!alkylated %in% alktbl$Alk)) {
    warning(paste0("\"", alkylated, "\" is not recognized as a valid alkylation and will be ignored!"))
    alk <- FALSE
  }
  if (alk) {
    mAlk <- match(alkylated, alktbl$Alk)
    alkC <- as.data.frame(t(c("C", paste0(alkylated, "-alkylated Cysteine"), "Alkylated", "C(alk)")))
    colnames(alkC) <- c("AA", "Name", "State", "Tag")
    for (atm in atms) {
      alkC[[atm]] <- Tbl[wC, atm]
      if (atm %in% colnames(alkTbl)) { alkC[[atm]] <- alkC[[atm]] + alkTbl[mAlk, atm] }
    }
    alkC$"Monoisotopic mass" <- sapply(1:nrow(alkC), function(x) {
      sum(alkC[x, atms]*atmsMonoMasses)
    })
    alkC$"Residue monoisotopic mass" <- alkC$"Monoisotopic mass" - sum(c(2, 1)*atmsMonoMasses[c("H", "O")])
    Tbl <- rbind(Tbl, alkC)
  }
  w <- which(!seq %in% Tbl$AA)
  if (length(w)) { warning("Bypassing some un-recognised amino acids!") }
  seq <- seq[which(seq %in% Tbl$AA)]
  wC2 <- which(seq == "C")
  if (freeCyst) { seq[wC2[1:freeCyst]] <- "C(ox)" }
  wC2 <- which(seq == "C")
  if (alk) { seq[wC2] <- "C(alk)" }
  seq <- aggregate(seq, list(seq), length)
  colnames(seq) <- c("AA", "Count")
  for (atm in atms) { seq[[atm]] <- sapply(seq$AA, function(aa) { Tbl[match(aa, Tbl$Tag), atm] })*seq$Count }
  seq2 <- data.frame(Element = atms, Count = colSums(seq[, atms]))
  w <- which(atms %in% seq2$Element)
  seq2 <- seq2[match(atms[w], seq2$Element),]
  # Effect of peptide bonds:
  seq2$Count[which(seq2$Element == "O")] <- seq2$Count[which(seq2$Element == "O")]-(L-1)
  seq2$Count[which(seq2$Element == "H")] <- seq2$Count[which(seq2$Element == "H")]-(L-1)*2
  stopifnot(min(seq2$Count) >= 0)
  #
  if (collapse) {
    seq2$Count <- as.character(seq2$Count)
    seq2 <- paste(apply(seq2, 1, paste, collapse = ""), collapse = "")
  }
  return(seq2)
}
Form <- ProtFormula(Seq, collapse = FALSE)
# Check monoisotopic mass
MonoIsoMass <- sum(apply(Form, 1, function(x) {
  as.numeric(gsub("_.+", "", IsotopeProbs$Monoisotopic[match(x[[1]], IsotopeProbs$Atom)]))*as.numeric(x[[2]])
}))
MonoIsoMass

# Parse isotopic frequencies table
ParseElementProbs <- function(IsoProbs, ProbThresh = 0) {
  setNames(lapply(1:nrow(IsoProbs), function(x) {
    x <- unlist(IsoProbs[x, 2:ncol(IsoProbs)])
    w <- which(x != "")
    x <- x[w]
    x <- as.data.frame(strsplit(x, "_"))
    for (i in 1:ncol(x)) { x[[i]] <- as.numeric(x[[i]]) }
    colnames(x) <- colnames(IsoProbs)[w+1]
    x <- x[, which(x[2,] >= ProbThresh), drop = FALSE]
    x <- x[, order(unlist(x[1,]), decreasing = FALSE), drop = FALSE]
    x[2,] <- x[2,]*1/sum(x[2,])
    return(x)
  }), IsoProbs$Atom)
}
AtomProbs <- ParseElementProbs(IsotopeProbs)
atProbs <- AtomProbs[which(names(AtomProbs) %in% Form$Element)]

# Calculate isotopic mass envelope
MaxN <- Inf
ProbThresh <- 10^-8
# For simplicity only considering cases which are more probable than the above threshold
# This allows me to use very inefficient code, because the mathematically correct.
Masses <- setNames(lapply(Form$Element, function(x) { #x <- Form$Element[1]
  AtPs <- atProbs[[x]]
  N <- Form$Count[match(x, Form$Element)]
  K <- max(c(min(c(ncol(AtPs), MaxN)), 1))
  # Mathematically correct code:
  Crazy <- FALSE
  if (Crazy) {
    masses <- combinations(K, N, repeats.allowed = TRUE)
    masses <- as.data.frame(t(apply(masses, 1, function(x) {
      sapply(1:K, function(y) { sum(x == y)} )
    })))
    colnames(masses) <- paste0("N", as.numeric(1:K))
    masses$Mass <- apply(masses[, paste0("N", as.numeric(1:K))], 1, function(x) { sum(unlist(x)*AtPs[1,]) })
    masses <- masses[order(masses$Mass, decreasing = FALSE),]
    masses$Probability <- apply(masses[, paste0("N", as.numeric(1:K))], 1, function(x) {
      #x <- masses[4, paste0("N", as.numeric(1:K))]
      x <- unlist(x)
      Ncum <- N-c(0, cumsum(x)[1:(length(x)-1)])
      prod(choose(Ncum, x)*(AtPs[2,]^x))
    })
    #sum(masses$Probability)
  }
  # ... but this fails for computational reasons.
  # Heuristic:
  masses <- data.frame(N1 = N)
  if (K > 1) {
    for (i in 2:K) { #i <- 2
      masses[[paste0("N", as.character(i))]] <- 0
      masses2 <- masses
      p <- AtPs[2,i]
      k <- 1; p0 <- p
      while (p0 > ProbThresh) { k <- k+1; p0 <- p0*p }
      for (j in 1:k) { #j <- 1
        tmp <- combinations(i-1, j, repeats.allowed = TRUE)
        tmp <- as.data.frame(t(apply(tmp, 1, function(x) { sapply(1:(i-1), function(y) { sum(x == y) }) })))
        tmp <- as.data.frame(cbind(-tmp, j))
        colnames(tmp) <- colnames(masses)
        tmp2 <- data.frame(Row1 = rep(1:nrow(masses), nrow(tmp)))
        tmp2$Row2 <- as.numeric(sapply(1:nrow(tmp), function(x) { rep(x, nrow(masses)) }))
        tmp2 <- masses[tmp2$Row1,]+tmp[tmp2$Row2,]
        tmp2 <- tmp2[which(apply(tmp2, 1, min) >= 0),]
        colnames(tmp2) <- colnames(masses)
        stopifnot(rowSums(tmp2) == N)
        masses2 <- rbind(masses2, tmp2)
        masses2 <- masses2[which(masses2[[paste0("N", as.character(i))]] > 0),]
      }
      masses <- rbind(masses, masses2)
      tst <- apply(masses, 1, paste, collapse = "_")
      tst <- aggregate(1:nrow(masses), list(tst), min)
      masses <- masses[tst$x,]
      stopifnot(rowSums(masses) == N)
    }
    masses$Mass <- apply(masses[, paste0("N", as.numeric(1:K))], 1, function(x) { sum(unlist(x)*AtPs[1,]) })
    masses <- masses[order(masses$Mass, decreasing = FALSE),]
    masses$Probability <- apply(masses[, paste0("N", as.numeric(1:K))], 1, function(x) {
      #x <- masses[4, paste0("N", as.numeric(1:K))]
      x <- unlist(x)
      Ncum <- N-c(0, cumsum(x)[1:(length(x)-1)])
      prod(choose(Ncum, x)*(AtPs[2,]^x))
    })
    masses$Probability <- masses$Probability/sum(masses$Probability)
    #sum(masses$Probability)
    masses <- masses[, c("Mass", "Probability"), drop = FALSE]
  } else {
    mass <- data.frame(Mass = AtPs[1,1], Probability = 1)
  }
}), Form$Element)
Envelope <- Masses[[1]]
if (length(Masses) > 1) {
  for (i in 2:length(Masses)) {
    nr1 <- nrow(Envelope)
    nr2 <- nrow(Masses[[i]])
    temp <- rep(1:nr1, nr2)
    Envelope <- Envelope[temp,]
    for (j in 1:nr2) {
      rg <- (1:nr1)+(nr1*(j-1))
      Envelope$Mass[rg] <- Envelope$Mass[rg]+Masses[[i]]$Mass[j]
      Envelope$Probability[rg] <- Envelope$Probability[rg]*Masses[[i]]$Probability[j]
    }
  }  
}
Envelope$Intensity <- Envelope$Probability*10^6
Res <- 0.01
Round <- -log10(signif(Res, 1))
Envelope$"Rounded_Mass" <- round(Envelope$Mass, Round)
Envelope2 <- aggregate(Envelope$Probability, list(Envelope$"Rounded_Mass"), mean)
colnames(Envelope2) <- c("Mass", "Density")
m <- Envelope2$Mass
m <- unique(c(m-Res/2, m+Res/2))
m <- m[which(!m %in% Envelope2$Mass)]
temp <- data.frame(Mass = m, Density = 0)
Envelope2 <- rbind(Envelope2, temp)
m2 <- Envelope2$Mass[which(Envelope2$Density > 0)]
stopifnot(sum(!(c(m2-Res/2, m2+Res/2) %in% Envelope2$Mass)) == 0)
Envelope2 <- Envelope2[order(Envelope2$Mass, decreasing = FALSE),]
ttl <- "Mass envelope"
plot <- ggplot(Envelope2) +
  geom_path(aes(x = Mass, y = Density), size = 0.1) + theme_bw() + ggtitle(ttl)
poplot(plot, 12, 20)
ggsave(paste0(gsub("/", "", ttl), ".jpeg"), plot, dpi = 300)

# This would be the part to incorporate variable modifications, allowing for a range of incorporation rates.
# Example: incorporate Fluorine (F has a single main isotope. Yes! Why is this not always so simple?)
Envelopes <- Envelope
IncorpFl <- FALSE # Eventually should be rewritten to accept a generic PTM (or multiple ones?)
if (IncorpFl) {
  nW <- nchar(Seq) - nchar(gsub("W", "", Seq))
  temp <- setNames(lapply(0:nW, function(x) {
    env <- Envelope
    env$Mass <- env$Mass + AtomProbs$F[1,1]*x
    env$Fluorines <- x
    return(env)
  }), paste0("F", 0:nW))
  Envelopes <- temp[[paste0("F", 0)]]
  for (i in 1:nW) { Envelopes <- rbind(Envelopes, temp[[paste0("F", i)]]) }
  # Effect of Fluorine incorporation:
  IncRate <- (0:10)/10
  AllEnvs <- data.frame(Row = rep(1:nrow(Envelopes), length(IncRate)))
  AllEnvs$IncRate <- as.numeric(sapply(IncRate, function(x) { rep(x, nrow(Envelopes))}))
  AllEnvs[, colnames(Envelopes)] <- Envelopes[AllEnvs$Row,]
  AllEnvs$Probability <- AllEnvs$Probability*dbinom(AllEnvs$Fluorines, nW, AllEnvs$IncRate)
  AllEnvs$Intensity <- AllEnvs$Probability*10^6
  #
  AllEnvs2 <- aggregate(AllEnvs$Probability, list(AllEnvs$"Rounded_Mass", AllEnvs$Fluorines, AllEnvs$IncRate), sum)
  colnames(AllEnvs2) <- c("Mass", "Fluorines", "Incorporation Rate", "Density")
  AllEnvs2$Fluorines <- factor(AllEnvs2$Fluorines, levels = nW:0)
  AllEnvs2$"Incorporation Rate" <- factor(AllEnvs2$"Incorporation Rate", levels = IncRate)
  ttl <- "Mass envelope with effect of Fluorine incorporation"
  plot <- ggplot(AllEnvs2) +
    geom_bar(stat = "identity", aes(x = Mass, y = Density, fill = Fluorines), width = Res,#position = "stack"
    ) + facet_grid(`Incorporation Rate`~.) + theme_bw() + ggtitle(ttl)
  poplot(plot, 12, 20)
  ggsave(paste0(gsub("/", "", ttl), ".jpeg"), plot, dpi = 300)
} else {
  AllEnvs <- Envelopes
}

#
# We cannot predict charge distribution, so we will simulate it from observed charge states
# The lowest observed charge was 25+ for denatured Tet, 9+ for native ARF1
# The most intense observed charge was 51+ for denatured Tet, 19-21 for native ARF1
msg <- "Enter the most intense observed charge state:"
Zmean <- as.integer(dlg_input(msg, c(50, 20)[Native+1])$res)
while ((!is.integer(Zmean))||(is.na(Zmean))) { Zmean <- dlg_input(msg, c(50, 20)[Native+1])$res }
msg <- "Enter the width of the observed charge states envelope (= Zmax - Zmin):"
Zwidth <- as.integer(dlg_input(msg, c(50, 20)[Native+1])$res)
while ((!is.integer(Zwidth))||(is.na(Zwidth))) { Zwidth <- dlg_input(msg, c(50, 20)[Native+1])$res }
Zmin <- Zmean - ceiling(Zwidth/2)
Zmax <- Zmean + floor(Zwidth/2)
ZPr <- data.frame(Z = Zmin:Zmax, Probability = dnorm(Zmin:Zmax, Zmean, (Zwidth)/6))
ttl <- "Simulated charge distribution"
plot <- ggplot(ZPr) + geom_bar(stat = "identity", aes(x = Z, y = Probability)) + theme_bw() + ggtitle(ttl) 
poplot(plot, 12, 20)
ggsave(paste0(gsub("/", "", ttl), ".jpeg"), plot, dpi = 300)
ZPr <- ZPr[which(ZPr$Probability >= ProbThresh),]
ZPr$Probability <- ZPr$Probability/sum(ZPr$Probability)
MZs <- data.frame(Row = rep(1:nrow(AllEnvs), nrow(ZPr)))
k <- colnames(AllEnvs)
k <- k[which(k != "Row")]
MZs[, k] <- AllEnvs[,k]
MZs$Z <- as.integer(sapply(ZPr$Z, function(x) { rep(x, nrow(AllEnvs)) }))
MZs$"M/Z" <- (MZs$Mass + 1.0078250319*MZs$Z)/MZs$Z
MZs$Intensity <- MZs$Intensity*ZPr$Probability[match(MZs$Z, ZPr$Z)]
MZs$Probability <- MZs$Intensity/(10^6)
MZs <- MZs[which(MZs$Probability > 0),]
min(MZs$Probability)
max(MZs$Probability)

Res <- 0.1
Round <- -log10(signif(Res, 1))
MZs$"Rounded_M/Z" <- round(MZs$"M/Z", Round)
aggrlist <- list(MZs$"Rounded_M/Z")
if (IncorpFl) {
  aggrlist <- list(MZs$"Rounded_M/Z", MZs$Fluorines, MZs$IncRate)
}
MZs2 <- aggregate(MZs$Intensity, aggrlist, mean)
colnames(MZs2)[1] <- "M/Z"
colnames(MZs2)[which(colnames(MZs2) == "x")] <- "Intensity"
MZs2 <- MZs2[which(MZs2$Intensity > 0),]
MZmin <- min(MZs2$`M/Z`)
MZmax <- max(MZs2$`M/Z`)
rg <- ((MZmin/Res):(MZmax/Res))*Res
#MZs2Bckp <- MZs2
#MZs2 <- MZs2Bckp
#print(nrow(MZs2))
wLst <- list("Global" = 1:nrow(MZs2))
if (IncorpFl) {
  wLst <- list()
  for (i in IncRate) { #i <- IncRate[1]
    for (f in 0:nW) { #f <- 0
      wLst[[paste0(i, " ", f)]] <- which((MZs2$`Incorporation Rate` == i)&(MZs2$Fluorines == f))
    }
  }
}
for (wNm in names(wLst)) {
  w <- wLst[[wNm]]
  if (length(w)) {
    m <- MZs2$`M/Z`[w]
    m <- unique(c(m-Res/2, m+Res/2))
    m <- m[which(!m %in% MZs2$`M/Z`[w])]
    temp <- data.frame("M/Z" = m, "Intensity" = 0, check.names = FALSE)
    if (IncorpFl) {
      temp$Fluorines <- f
      temp$"Incorporation Rate" <- i
    }
    MZs2 <- rbind(MZs2, temp)
    w2 <- 1:nrow(MZs2)
    if  (IncorpFl) {
      w2 <- which((MZs2$`Incorporation Rate` == i)&(MZs2$Fluorines == f))
    }
    m2 <- MZs2$`M/Z`[w2][which(MZs2$Intensity[w2] > 0)]
    stopifnot(sum(!(c(m2-Res/2, m2+Res/2) %in% MZs2$`M/Z`[w2])) == 0)
  }
}
ttl <- "M/Z envelope"
if (IncorpFl) {
  MZs2$Fluorines <- factor(MZs2$Fluorines, levels = 0:nW)
  MZs2$"Incorporation Rate" <- factor(MZs2$"Incorporation Rate", levels = IncRate)
  MZs2 <- MZs2[order(MZs2$Fluorines, MZs2$`Incorporation Rate`, decreasing = FALSE),]
  ttl <- paste0(ttl, " with effect of Fluorine incorporation")
}
MZs2 <- MZs2[order(MZs2$`M/Z`, decreasing = FALSE),]
m <- max(MZs2$Intensity)
#w <- which(MZs2$Intensity >= 1)
#plot <- ggplot(MZs2[w,]) +
#  geom_bar(stat = "identity", aes(x = `M/Z`, y = Intensity, fill = Fluorines), width = Res, position = "stack") +
plot <- ggplot(MZs2)
plot <- plot + geom_path(aes(x = `M/Z`, y = Intensity), size = 0.1)
if (IncorpFl) {
  plot <- plot + geom_path(aes(x = `M/Z`, y = Intensity, colour = Fluorines, group = Fluorines), size = 0.1) +
    facet_grid(`Incorporation Rate`~Fluorines)
}
plot <- plot + theme_bw() + ggtitle(ttl)
poplot(plot, 12, 20)
ggsave(paste0(gsub("/", "", ttl), ".jpeg"), plot, dpi = 300)
