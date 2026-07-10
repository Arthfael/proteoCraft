#' Coverage
#'
#' @description
#' This function can be used in 5 different modes:\cr
#' - Coverage: calculates the proportion (from 1 to 0) of parent protein sequence(s) covered by the peptide sequence(s) provided.\cr
#' - Align: plots the protein(s), highlighting amino acids covered by the peptides sequence(s) provided.\cr
#' - Align2L: as "Align", but peptides are displayed individually as bars; accepts PTM-modified sequences and will display PTM positions on the map; in addition, peptide colour can be mapped to the "intensities" argument.\cr
#' - Heat: Amino acid colour reflects the sum of intensities for all peptides observed covering said position. There is in addition one row for each PTM - also with color mapped to intensity.\cr
#' - XML: formats covered sequence as XML to display coverage in Excel.\cr
#' \cr
#' NB:\cr
#' This version incorporates a code chunk which was designed with the help of AIs to very significantly improve the speed of mode Align2 where there is a large number of peptides to plot.\cr
#' The order in which peptides are plotted may differ from before as a result.\cr
#' Display issues from when too many peptides were plotted (and spilled into the next row) have also been addressed (by enforcing a minimum peptide segment width).\cr
#' 
#' @param proteins The (named) protein sequence(s) on which to map the peptides.
#' @param peptides The collection of peptide sequences to map on the protein sequence(s).
#' @param Mode One of:\cr
#'  - "Coverage": outputs only the percentage of each protein sequence covered by the peptide(s);\cr
#'  - "Align": creates simple alignment plots;\cr
#'  - "Align2": creates more advanced alignment plots with individual peptides and PTMs displayed;\cr
#'  - "Heat": amino acid color reflects sum of intensities at said position; PTMs are represented with one row for each;\cr
#'  - "XML": creates xml-formatted strings for writing into Excel tables using openxlsx2.
#' @param new.window Only used if Mode = "Align", "Align2" or "Heat". If set to TRUE (default), will create plots in new window.
#' @param display Only used if Mode = "Align", "Align2" or "Heat". Logical: should we print the plot?
#' @param scale Only used if Mode = "Align", "Align2" or "Heat". Number of amino acids to display per row of the plots (default = 60).
#' @param title Only used if Mode = "Align", "Align2" or "Heat". Default = "Coverage": title of the plot Can be either a single character string (all proteins will be plotted and saved as one graph) or a vector of the same length as "proteins" (each protein will have its graph)
#' @param colour Ignored if Mode = "Coverage". The colour in which covered regions of protein sequences should be shown. Default = "red".
#' @param colscale Only used if Mode = "Align2". Which colour scale do you want to use for peptide intensities? Default = 1 ("darkblue"->"orange"). Other values are 2L: "darkblue"->"red", 3: "darkblue"->"green", 4: "black"->"blue", 5: "black"->"orange", 6: "black"->"red", 7: "black"->"green" and 8: "red"->"green"
#' @param size Only used if Mode = "Align", "Align2" or "Heat". The character size. Default = 2.5, you may want to reduce it for large proteins.
#' @param save Only used if Mode = "Align", "Align2" or "Heat". Preferred alternative to "save.path", though less flexible: vector of ggsave-compatible file extensions, (or .html to save as plotly) to add to the "title" argument to save the graph to. Allows saving as multiple formats through a single function call. Default = "jpeg"; set to FALSE to not save it.
#' @param save.path Only used if Mode = "Align", "Align2" or "Heat". Prefer the "save" argument. If set, must be a character vector with the name of the file(s) in which to save each plot, so must have the same length as \"title\". The name must end in a ggsave-compatible extension. Default = FALSE
#' @param intensities Only used if Mode = "Align2" or "Heat". Values to map to the "colscale" argument so that peptides can be printed with different colours, e.g. mapped to abundance.
#' @param na Only used if Mode = "Align2" or "Heat". Colour for NA values.
#' @param maxInt If provided, the maximum of the intensity scale is not detected from the values in intensities but provided externally. Useful if drawing maps from several intensity vectors and wanting to apply a single scale to all.
#' @param bgcol Only used if Mode = "Align2" or "Heat". Which colour should the background be? Default = "black".
#' @param I_eq_L Should we consider I and L identical? Currently, by default, TRUE for both DIA and DDA: see https://github.com/vdemichev/DiaNN/discussions/1631
#' 
#' @returns
#' The function only returns something in these conditions:
#'  - If mode is "Coverage", the coverage of proteins by peptides.\cr
#'  - If mode is "XML", Excel-compatible xml-formatted sequence coverage strings.\cr
#'  - If mode is "Align2" and a plotly html plot is saved, it is returned.\cr
#' 
#' @examples
#' proteins <- "MNTTDCFIALVQAIREIKALFLSRTTGKMELTLYNGEKKTFYSRPNNHDNCWLNAILQLFRYVEEPFFDWVYSSPENLTLEAIKQLEDLTGLELHEGGPPALVIWNIKHLLHTGIGTASRPSEVCMVDGTDMCLADFHAGIFLKGQEHAVFACVTSNGWYAIDDEDFYPWTPDPSDVLVFVPYDQEPLNGEWKAKVQRKLKGAGQSSPATGSQNQSGNTMHMDIVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYKGPHHHHHHTEYKPTVRLATRDDVPRAVRTLAAAFADYPATRHTVDPDRHIERVTELQELFLTRVGLDIGKVWVADDGPAVAVWTTPESVEAGAVFAEIGPRMAELSGSRLAAQQQMEGLLAPHRPKEPAWFLATVGVSPDHQGKGLGSAVVLPGVEAAERAGVPAFLETSAPRNLPFYERLGFTVTADVEVPEGPRTWCMTRKPGATRVTELLYRMKRAETYCPRPLLAIHPTEARHKQKIVAPVKQTLNFDLLKLAGDVESNPGPFFFSDVRSNFSKLVETINQMQEDMSTKHGPDFNRLVSAFEELAIGVKAIRTGLDEAKPWYKLIKLLSRLSCMAAVAARSKDPVLVAIMLADTGLEILDSTFVVKKISDSLSSLFHVPAPVFSFGAPVLLAGLVKVASSFFRSTPEDLERAEKQLKARDINDIFAILKNGEWLVKLILAIRDWIKAWIASEEKFVTMTDLVPGILEKQRDLNDPSKYKEAKEWLDNARQACLKSGNVHIANLCKVVAPAPSKSRPEPVVVCLRGKSGQGKSFLANVLAQAISTHFTGRIDSVWYCPPDPDHFDGYNQQTVVVMDDLGQNPDGKDFKYFAQMVSTTGFIPPMASLEDKGKPFNSKVIIATTNLYSGFTPRTMVCPDALNRRFHFDIDVSAKDGYKINSKLDIIKALEDTHANPVAMFQYDCALLNGMAVEMKRMQQDMFKPQPPLQNVYQLVQEVIDRVELHEKVSSHPIFKQISIPSQKSVLYFLIEKGQHEAAIEFFEGMVHDSIKEELRPLIQQTSFVKRAFKRLKENFEIVALCLTLLANIVIMIRETRKRQKMVDDAVNEYIEKANITTDDKTLDEAEKSPLETSGASTVGFRERTLPGQKACDDVNSEPAQPVEEQPQAEGPYAGPLERQKPLKVRAKLPQQEGPYAGPMERQKPLKVKAKAPVVKEGPYEGPVKKPVALKVKAKNLIVTESGAPPTDLQKMVMGNTKPVELILDGKTVAICCATGVFGTAYLVPRHLFAEKYDKIMVDGRAMTDSDYRVFEFEIKVKGQDMLSDAALMVLHRGNRVRDITKHFRDTARMKKGTPVVGVINNADVGRLIFSGEALTYKDIVVCMDGDTMPGLFAYRAATKAGYCGGAVLAKDGADTFIVGTHSAGGNGVGYCSCVSRSMLLKMKAHIDPEPHHEGLIVDTRDVEERVHVMRKTKLAPTVAHGVFNPEFGPAALSNKDPRLNEGVVLDEVIFSKHKGDTKMSEEDKALFRRCAADYASRLHSVLGTANAPLSIYEAIKGVDGLDAMEPDTAPGLPWALQGKRRGALIDFENGTVGPEVEAALKLMEKREYKFVCQTFLKDEIRPLEKVRAGKTRIVDVLPVEHILYTRMMIGRFCAQMHSNNGPQIGSAVGCNPDVDWQRFGTHFAQYRNVWDVDYSAFDANHCSDAMNIMFEEVFRTEFGFHPNAEWILKTLVNTEHAYENKRITVGGGMPSGCSATSIINTILNNIYVLYALRRHYEGVELDTYTMISYGDDIVVASDYDLDFEALKPHFKSLGQTITPADKSDKGFVLGHSITDVTFLKRHFHMDYGTGFYKPVMASKTLEAILSFARRGTIQEKLISVAGLAVHSGPDEYRRLFEPFQGLFEIPSYRSLYLRWVNAVCGDA"
#' #  OR...
#' proteins <- c(proteins, "MNTTDCFIALVQAIREIKALFLSRTTGKMELTLYNGEKKTFYSRPNNHDNCWLNAILQLFRYVEEPFFDWVYSSPENLTLEAIKQLEDLTGLELHEGGPPALVIWNIKHLLHTGIGTASRPSEVCMVDGTDMCLADFHAGIFLKGQEHAVFACVTSNGWYAIDDEDFYPWTPDPSDVLVFVPYDQEPLNGEWKAKVQRKLK", "PFFFSDVRSNFSKLVETINQMQEDMSTKHGPDFNRLVSAFEELAIGVKAIRTGLDEAKPWYKLIKLLSRLSCMAAVAARSKDPVLVAIMLADTGLEILDSTFVVKKISDSLSSLFHVPAPVFSFGAPVLLAGLVKVASSFFRSTPEDLERAEKQ")
#' names(proteins) <- c("polyprotein", "Lpro", "2B")
#' peptides <- c("ACDDVNSEPAQPVEEQPQAE","ACDDVNSEPAQPVEEQPQAEGPYAGPLER","ACDDVNSEPAQPVEEQPQAEGPYAGPLERQK","AGVPAFLETSAPR","AHIDPEPHHE","AHIDPEPHHEGLIVDTR","AIRTGLDEAK","AIRTGLDEAKPWYK","ALEDTHANPVAMFQYDCALLNGMAVEMK","ALFLSRTTGK","ALFLSRTTGKMELTLYNGEK","AMTDSDYR","AMTDSDYRVFEFEIK","ANITTDDK","ANITTDDKTLDEAEK","ARDINDIFAILK","AWIASEEK","DEIRPLEK","DIVVCMDGDTMPGLFAYR","DIVVCMDGDTMPGLFAYRAATK","DLNDPSK","DPRLNEGVVLDEVIFSK","EDGNILGHK","EELRPLIQQTSFVK","EGPYEGPVK","EPAWFLATVGVSPDHQGK","ERTLPGQK","EWLDNAR","EWLDNARQACLK","FEGDTLVNR","FEGDTLVNRIELK","FGTHFAQYR","FHFDIDVSAK","FHFDIDVSAKDGYK","FSVSGEGEGDATYGK","FSVSGEGEGDATYGKLTLK","FVCQTFLK","FVTMTDLVPGILEK","FVTMTDLVPGILEKQR","FVTMTDLVPGILEKQRDLNDPSK","GALIDFENGTVGPEVEAALK","GEELFTGVVPILVELDGDVNGHK","GFVLGHSITDVTFLK","GIDFKEDGNILGHK","GLGSAVVLPGVEAAER","GLIVDTRDVEER","GLIVDTRDVEERVHVMRK","GQDMLSDAALMVLHR","GQDMLSDAALMVLHRGNRVRDITK","GQHEAAIEFFEGMVHDSIK","GTPVVGVINNADVGR","GTPVVGVINNADVGRLIFSGEALTYK","GVDGLDAMEPDTAPGLPWALQGK","HGPDFNR","HLLHTGIGTASR","HTVDPDR","HTVDPDRHIER","IMVDGRAMTDSDYRVFEFEIK","INSKLDIIK","ISIPSQK","IVDVLPVEHILYTR","KGTPVVGVINNADVGR","KGTPVVGVINNADVGRLIFSGEALTYK","LAAQQQMEGLLAPHRPK","LAGDVESNPGPFFFSDVR","LAGDVESNPGPFFFSDVRSNFSK","LAPTVAHGVFNPEFGPAALSNK","LAPTVAHGVFNPEFGPAALSNKDPR","LEYNYNSHNVYIMADK","LGFTVTADVEVPEGPR","LHSVLGTANAPLSIYEAIK","LIFSGEALTYK","LISVAGLAVHSGPDEYR","LISVAGLAVHSGPDEYRR","LNEGVVLDEVIFSK","LPQQEGPYAGPMER","LPQQEGPYAGPMERQK","LSCMAAVAAR","LVETINQMQEDMSTK","LVSAFEELAIGVK","MAELSGSR","MELTLYNGEK","MELTLYNGEKK","MVDDAVNEYIEK","MVDDAVNEYIEKANITTDDK","MVMGNTKPVELILDGK","NLIVTESGAPPTDLQK","NLPFYER","PFFFSDVR","PFFFSDVRSNFSK","PLIQQTSFVK","PLLAIHPTEAR","PVELILDGK","QISIPSQK","QRDLNDPSK","QRDLNDPSKYK","RAETYCPRPLLAIHPTEARHK","RFHFDIDVSAK","RFHFDIDVSAKDGYK","RHFHMDYGTGFYK","RHFHMDYGTGFYKPVMASK","RQKMVDDAVNEYIEK","RRGALIDFENGTVGPEVEAALK","SAMPEGYVQER","SAMPEGYVQERTIFFK","SGAPPTDLQK","SGNVHIANLCK","SLGQTITPADK","SLGQTITPADKSDK","SNFSKLVETINQMQEDMSTK","SPLETSGASTVGFR","SPLETSGASTVGFRER","SPLETSGASTVGFRERTLPGQK","SRPEPVVVCLR","SRPEPVVVCLRGK","STPEDLER","STPEDLERAEK","STPEDLERAEKQ","SVLYFLIEK","TEFGFHPNAEWILK","TGLDEAK","TGLDEAKPWYK","TGLDEAKPWYKLIK","TIFFKDDGNYK","TLAAAFADYPATR","TLDEAEK","TLDEAEKSPLETSGASTVGFR","TLDEAEKSPLETSGASTVGFRERTLPGQK","TLEAILSFAR","TLEAILSFARRGTIQEK","TLVNTEHAYENK","TLVNTEHAYENKR","TMVCPDALNR","TRIVDVLPVEHILYTR","VASSFFR","VASSFFRSTPEDLER","VASSFFRSTPEDLERAEK","VASSFFRSTPEDLERAEKQ","VFEFEIK","VGLDIGK","VIIATTNLYSGFTPR","VSSHPIFK","VSSHPIFKQ","VTELLYR","VTELQELFLTR","VVAPAPSK","YFAQMVSTTGFIPPMASLEDK")
#' Coverage(proteins, peptides, "Coverage")
#' Coverage(proteins, peptides, "Align", TRUE, c(TRUE,TRUE,TRUE), 30)
#' Coverage(proteins, peptides, "Align", TRUE, c(TRUE,TRUE,TRUE), 60)
#' Coverage(proteins, peptides, "Align2", TRUE, c(TRUE,TRUE,TRUE), 60)
#' 
#' @export

Coverage <- function(proteins,
                     peptides,
                     Mode = "Coverage",
                     new.window = TRUE, 
                     display = TRUE,
                     scale = 60L,
                     title = "Coverage",
                     colour = "red",
                     colscale = 1L,
                     size = 2.5,
                     save = "jpeg",
                     save.path,
                     intensities = NULL,
                     na = "red",
                     maxInt = NULL,
                     bgcol = "black",
                     I_eq_L = TRUE) {
  TESTING <- FALSE
  #DefArg(Coverage);TESTING = TRUE
  #proteins = P; peptides = tmp$"Modified sequence"; Mode = "Align2"; title = ttl; save = c("jpeg", "pdf"); intensities = tmp$`log10(Intensity)`
  #proteins = x[[1L]]; peptides = x[[2L]]; Mode = "XML"; colour = "green"
  #proteins = sq; peptides = tmpSq$"Modified sequence_verbose"; intensities = tmpSq$Intensity; Mode = "Align2"; save = FALSE
  #proteins = sq; peptides = tmpSq$"Modified sequence_verbose"; intensities = tmpSq$Intensity; Mode = "Heat"; save = FALSE
  #proteins = seq; peptides = p1$Sequence; Mode = "Align2"; title = paste0("Coverage map - ", nm); save = c("jpeg", "pdf"); intensities = p1$Intensity; display = FALSE
  #proteins = P; peptides = tmp$"Modified sequence"; Mode = "Heat"; display = FALSE;title = ttl1b; save = c("jpeg", "pdf"); intensities = tmp$`log10(Intensity)`;maxInt = mxInt; na = "cyan"
  misFun <- if (TESTING) {
    # Note:
    # This is not a perfect alternative to missing but will work in most cases, unless x matches a function imported by a package 
    \(x) { return(!exists(deparse(substitute(x)))) }
  } else { missing }
  #
  if (misFun(I_eq_L) || (!is.logical(I_eq_L)) || is.na(I_eq_L)) {
    #if (exists("isDIA") && is.logical(isDIA) && (!is.na(isDIA))) {
    #  I_eq_L <- !isDIA # NO! Cf. Vadim's reply that at the moment models do not allow discriminating between I and L.
    #} else {
    I_eq_L <- TRUE
    #}
  }
  if (!Mode %in% c("Coverage", "Align", "Align2", "XML", "Heat")) {
    stop("Accepted value for \"Mode\" argument:\n - \"Coverage\": only return percentage of sequence coverage\n - \"Align\": create simple peptide coverage map\n - \"Align2\": create peptide coverage map with individual peptide intensities\n - \"XML\": writes xml formated coverage text for writing into Excel tables.")
  }
  if ((Mode == "Align") && (colour == "black")) {
    warning("Any R-compatible colour but \"black\" is allowed - though I would advise something shiny...\nDefaulting to \"red\".\nYou're welcome!")
    colour <- "red"
  }
  if (is.null(names(proteins))) { names(proteins) <- 1L:length(proteins) }
  namez <- names(proteins)
  if (!length(peptides)) {
    if (Mode == "Coverage") {
      cat("No peptides provided -> coverage = 0%\n")
      return()
    }
    if (Mode %in% c("Align", "Align2")) {
      warning("No peptides to plot!")
      return()
    }
    if (Mode == "XML") {
      return(fmt_txt(paste(proteins, collapse = ";"), color = wb_color(hex = c("grey"))))
    }
  }
  if (Mode != "Coverage") {
    Table <- data.frame(Title = title)
    if (length(title) == 1L) { Table$"Protein(s)" <- list(namez) } else {
      if (length(title) != length(proteins)) {
        if (length(title) == 1L) { title <- paste0(title, " - ", namez) } else {
          warning("\"title\"'s length must be either one of the same length as the protein vector, creating a default title!")
          title <- paste0("Coverage - ", namez)
          Table <- data.frame(Title = gsub(":|\\*|\\?|<|>|\\||/", "-", title))
        }
      }
      Table$"Protein(s)" <- as.list(namez)
    }
    if ((length(save) > 1L) || (as.character(save) != "FALSE")) {
      save.grph <- misFun(save.path) || (length(save.path) != length(title))
      if (save.grph) {
        if ((!misFun(save.path)) && (length(save.path) != length(title))) {
          warning("Argument \"save.path\" will be ignored as its length is not compatible with that of the \"title\" argument!")
        }
        save.grph <- as.character("save") != "FALSE"
        if (save.grph) { Table$Path <- lapply(Table$Title, \(x) { paste0(x, ".", save) }) }
      } else {
        Table$Path <- if (as.character("save") == "FALSE") {
          save.path
        } else {
          #warning("Both the \"save.path\" and \"save\" arguments were provided. I will use \"save.path\" for the file name but \"save\" for the extension(s).")
          lapply(gsub("\\.[^A-Z,a-z,0-9]+$", "", save.path), \(x) { paste0(x, ".", save) })
        }
      }
    } else { save.grph <- FALSE }
  }
  if (misFun(intensities) || is.null(intensities)) { intensities <- rep(1, length(peptides)) } else {
    intensities <- as.numeric(unlist(intensities))
    if (length(intensities) != length(peptides)) {
      stop("The \"intensities\" vector should be the same length as the \"peptides\" vector!")
    }
    w2 <- which(!is.finite(intensities))
    intensities[w2] <- NA_real_
  }
  peptides <- gsub("^_|_$", "", unlist(peptides))
  # Account for the fact that Leucine and Isoleucine are indistinguishable by MS
  peptidesOrig <- peptides
  proteinsOrig <- proteins
  if (I_eq_L) {
    peptides <- gsub("I", "L", peptides)
    proteins <- gsub("I", "L", proteins)
  }
  #
  if (Mode == "Align2") {
    maxInt <- if (misFun(maxInt)) { NA_real_ } else { suppressWarnings(as.numeric(maxInt)) }
    maxInt <- max(c(maxInt, intensities), na.rm = TRUE)
    maxInt <- ceiling(maxInt)
    l <- nchar(as.character(maxInt))
    maxInt2 <- round(maxInt, 1L-l)
    while (maxInt2 < maxInt) { maxInt2 <- maxInt2 + 10L^(l-1L) }
    brks <- sort(unique(c(c(0L:4L)*(maxInt2/4), maxInt)))
    #brks <- brks[which(brks <= maxInt)]
    brkstxt <- as.character(brks)
    w <- which(brkstxt == as.character(maxInt))
    brkstxt[w] <- paste0(" ---> max. int.: ", brkstxt[w])
  }
  peptides <- data.frame(Original.seq = peptidesOrig,
                         Match.seq = peptides,
                         Sequence = gsub("\\([^\\)]+\\)", "", peptides),
                         id = 1L:length(peptides),
                         stripped.Original.seq = gsub("\\([^\\)]+\\)", "", peptidesOrig))
  if (!is.null(intensities)) { peptides$Intensity <- intensities }
  peptides$Length <- nchar(peptides$Sequence)
  seq <- unique(peptides$Sequence)
  myCover <- c()
  pos <- list()
  pos2 <- data.frame(Seq = seq)
  for (P in namez) { #P <- namez[1L]
    protein <- paste0("_", proteins[P], "_")
    proteinOrig <- paste0("_", proteinsOrig[P], "_")
    pos2[[P]] <- lapply(seq, \(sq) { #sq <- seq[1L]
      temp <- unlist(strsplit(protein, sq))
      l <- length(temp)
      res <- if (l > 1L) { cumsum(nchar(temp[1L:l-1L])) + cumsum(c(0L, rep(nchar(sq), l-2L))) } else { NA }
      return(res)
    })
    temp <- setNames(lapply(1L:length(seq), \(x) {
      sq <- seq[x]
      x <- data.frame(start = pos2[[match(sq, pos2$Seq), P]])
      w <- which(!is.na(x$start))
      x <- x[w, , drop = FALSE]
      x$end <- x$start + nchar(sq) - 1L
      x$spread <- apply(x[, c("start", "end")], 1L, \(x) { list(x[[1L]]:x[[2L]])})
      return(x)
    }), seq)
    spread <- unique(unlist(sapply(seq, \(sq) { temp[[sq]]$spread })))
    res2 <- rep(0L, nchar(protein)-2L)
    res2[spread] <- 1L
    #setNames(res2, unlist(strsplit(substr(protein, 2L, nchar(protein)-1L), "")))
    pos[[P]] <- res2
    myCover[P] <- mean(pos[[P]])
  }
  if (Mode == "Coverage") { return(myCover) }
  if (Mode %in% c("Align2", "Heat")) {
    colscales <- rbind(c("darkblue", "orange"),
                       c("darkblue", "red"),
                       c("darkblue", "green"),
                       c("black", "blue"),
                       c("black", "orange"),
                       c("black", "red"),
                       c("black", "green"),
                       c("red", "green"))
    if (!suppressWarnings(as.numeric(colscale[1L])) %in% 1L:nrow(colscales)) {
      warning(paste0("Argument \"colscale\" must be a (single) number between 1 and ", nrow(colscales), ", defaulting to 1!"))
      colscale <- 1L
    }
    colscale <- colscales[colscale,]
  }
  xmList <- plotLyst <- list()
  for (Ttl in unique(title)) { #Ttl <- unique(title)[1L]
    plotCreated <- FALSE
    tbl <- Table[which(Table$Title == Ttl),]
    nmz <- unlist(tbl$Protein)
    prots <- proteins[nmz]
    protsOrig <- proteinsOrig[nmz]
    poZ2 <- pos2[,c("Seq", nmz)]
    align.temp <- list()
    YlabZ <- list()
    for (P in nmz) { #P <- nmz[1L]
      protein <- prots[P]
      proteinOrig <- protsOrig[P]
      tmpAlgn <- data.frame(AA = unlist(strsplit(protein, "")))
      tmpAlgn$N <- 1L:nrow(tmpAlgn)
      tmpAlgn$Amino_acid <- unlist(strsplit(proteinOrig, ""))
      tmpAlgn$Protein <- P
      tmpAlgn$Legend <- c("not found", "found")[pos[[P]] + 1L]#poZ[[P]] + 1L]
      tmpAlgn$X <- ((1L:length(pos[[P]])) - 1L) %% scale + 1L#poZ[[P]])) - 1L) %% scale + 1L
      tmpAlgn$Y <- - vapply(1L:length(pos[[P]]) - 1L, \(x) {
        x - (x %% scale)
      }, 1)/scale#Z[[P]]) - 1L, \(x) {x - (x %% scale)})/scale
      tmpAlgn$Y <- -min(tmpAlgn$Y) + 1L + tmpAlgn$Y
      tmpAlgn$Face <- vapply(pos[[P]], \(x) { c("plain", "bold")[x+1L] }, "")#poZ[[P]], \(x) {c("plain", "bold")[x+1L]})
      tmpAlgn$Hjust <- 0.5
      YlabZ[[P]] <- c(rev(c(c(0L:(max(tmpAlgn$Y, na.rm = TRUE) - 1L)) * scale + 1L)), "")
      covStr <- paste0("Coverage: ", round(100*myCover[P], 1L), "%")
      align.temp[[P]] <- rbind(tmpAlgn,
                               data.frame(AA = paste0("Protein name: ",  P),
                                          N = NA,
                                          Amino_acid = "",
                                          Protein = "",
                                          Legend = "not found",
                                          X = scale/2,
                                          Y = max(tmpAlgn$Y) + 1.2,
                                          Face = "bold.italic",
                                          Hjust = 0.5),
                               data.frame(AA = covStr,
                                          N = NA,
                                          Amino_acid = "",
                                          Protein = "",
                                          Legend = "not found",
                                          X = scale/2,
                                          Y = max(tmpAlgn$Y) + 0.6,
                                          Face = "italic",
                                          Hjust = 0.5))
    }
    if (length(prots) > 1L) {
      offset <- vapply(prots[2L:length(prots)], \(x) { (nchar(x)/scale) }, 1) + 2L
      if (length(offset) > 1L) {
        offset <- vapply(c(1L:length(offset)), \(x) { sum(offset[x:length(offset)]) }, 1)
      }
      for (i in 1L:length(offset)) { align.temp[[i]]$Y <- align.temp[[i]]$Y + offset[i] }
      for (i in 1L:length(offset)) {
        align.temp[[1L]] <- rbind(align.temp[[1L]], align.temp[[i + 1L]])
        YlabZ[[length(offset) + 1L]] <- c(YlabZ[[length(offset) + 1L]], "", YlabZ[[length(offset) + 1L - i]])
      }
      YlabZ <- c(YlabZ[[length(offset) + 1L]], "")
    } else { YlabZ <- c(YlabZ[[1L]], "") }
    align.temp <- align.temp[[1L]]
    align.temp <- align.temp[order(align.temp$Y, decreasing = TRUE),]
    if (Mode == "Align") {
      myColors <- setNames(c("black", colour), c("not found", "found"))
      colScale <- ggplot2::scale_colour_manual(name = "colour", values = myColors)
      Xleft <- -3L
      Xright <- ceiling(max(align.temp$X, na.rm = TRUE)) + 4L
      Ytop <- max(align.temp$Y, na.rm = TRUE) + 3L
      Ybottom <- min(align.temp$Y, na.rm = TRUE) - 3L
      yxRat <- yxRat2 <- 2L
      covPlot <- ggplot2::ggplot(align.temp) + ggplot2::coord_fixed(yxRat) +
        ggplot2::geom_text(ggplot2::aes(label = AA,  x = X, y = Y, colour = Legend, fontface = Face,
                                        hjust = Hjust), cex = size) +
        colScale + ggplot2::scale_x_discrete(breaks = NULL) +
        ggplot2::scale_y_discrete(breaks = NULL, #labels = YlabZ
        ) +
        ggplot2::labs(x = "", y = "Position") +
        ggplot2::xlim(Xleft, Xright) + ggplot2::ylim(Ybottom, Ytop) +
        ggplot2::ggtitle(title) +
        ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank(),
                       panel.background = ggplot2::element_blank(),
                       legend.title = ggplot2::element_blank(),
                       #aspect.ratio = 1.5
        )
      #poplot(covPlot)
      plotCreated <- TRUE
    }
    if (Mode %in% c("Align2", "Heat")) {
      matches <- lapply(nmz, \(n) { #n <- nmz[1L]
        peptides[[paste0("Matches_", n)]] <- lapply(peptides$Sequence, \(x) { #x <- peptides$Sequence[1L]
          x <- poZ2[[match(x, poZ2$Seq), n]]
          return(x[which(!is.na(x))])
        })
        tst <- which(lengths(peptides[[paste0("Matches_", n)]]) > 0L)
        if (!length(tst)) { return() }
        mtchs <- peptides[tst, c(paste0("Matches_", n), "Match.seq", "Sequence", "id", "stripped.Original.seq", "Original.seq")]
        colnames(mtchs) <- c("Match", "Match.seq", "Sequence", "pepid", "stripOrig", "Orig")
        mtchs$Protein <- n
        return(mtchs)
      })
      matches <- plyr::rbind.fill(matches)
      if (is.data.frame(matches) && nrow(matches)) {
        matches$temp <- do.call(paste, c(matches[, c( "Match.seq", "Sequence", "pepid", "Protein", "stripOrig", "Orig")],
                                         sep = "_____"))
        matches <- listMelt(matches$Match, matches$temp)
        colnames(matches)[1L] <- "Match"
        matches[, c( "Match.seq", "Sequence", "pepid", "Protein", "stripOrig", "peptide")] <- Isapply(strsplit(matches$L1, "_____"), unlist)
        matches$L1 <- NULL
        matches$Length <- nchar(matches$Sequence)
        matches$pepid <- as.numeric(matches$pepid)
        matches[, c("X1", "Y1")] <- as.data.frame(t(apply(matches[,c("Protein", "Match")], 1L, \(x) {
          unlist(align.temp[which(align.temp$Protein == as.character(x[[1L]]))[as.numeric(x[[2L]])], c("X", "Y")])
        })))
        matches[, c("X2", "Y2")] <- as.data.frame(t(apply(matches[,c("Length", "X1", "Y1")], 1L, \(x) {
          x2 <- x[[2L]]+x[[1L]]-1L
          y2 <- x[[3L]] - floor(x2/scale)+(x2 %% scale == 0L)
          x2 <- (x2 %% scale) + scale*(x2 %% scale == 0L)
          return(c(x2, y2))
        })))
        matches$Intensity <- if (!is.null(intensities)) { peptides$Intensity[match(matches$pepid, peptides$id)] } else { 1 }
        # Deal with peptides spanning several rows
        matches$Start <- 1L # Not the match in the protein sequence! Used only for peptides spanning several rows!
        matches$End <- matches$Length # As above
        matches$Is_start <- TRUE
        matches$Is_end <- TRUE
        w <- which(matches$Y1 != matches$Y2) # Peptides spanning several rows
        if (length(w)) {
          tmpMatches <- lapply(w, \(i) { #i <- w[1L] #i <- w[5L]
            #print(i)
            n <- matches$Y1[i]-matches$Y2[i] # (we almost never expect this to be more than 1: this would mean a super long peptide!!!)
            nr <- n+1L
            temp <- matches[rep(i, nr),]
            temp$Y1 <- temp$Y2 <- max(temp$Y1)-(0L:n)
            temp$X1[2L:nr] <- 1L
            temp$X2[1L:(nr-1L)] <- scale
            temp$Length <- temp$X2 - temp$X1 + 1L
            temp2 <- as.data.frame(annot_to_tabl(unique(temp$Match.seq),
                                                 Nterm = FALSE,
                                                 Cterm = FALSE))
            colnames(temp2) <- c("A", "B")
            temp2$B[which(temp2$B != "")] <- paste0("(", temp2$B[which(temp2$B != "")], ")")
            temp2$A[which(temp2$A == "_")] <- ""
            rrange <- cumsum(temp$Length)+1L
            temp$Match.seq <- vapply(1L:nr, \(x) {
              paste(apply(temp2[c(1L, rrange+1L)[x]:rrange[x], , drop = FALSE], 1L, paste, collapse = ""), collapse = "")
            }, "")
            temp$Sequence <- gsub("\\([^\\)]+\\)", "", temp$Match.seq)
            temp$End <- cumsum(temp$Length)
            temp$Start <- c(0L, temp$End[1L:(nr-1L)]) + 1L
            temp$stripOrig <- substr(temp$stripOrig, temp$Start, temp$End)
            temp$Is_start[2L:nr] <- FALSE
            temp$Is_end[1L:(nr-1L)] <- FALSE
            return(temp)
          })
          tmpMatches <- do.call(rbind, tmpMatches)
          matches <- rbind(matches[-w,],
                           tmpMatches)
          matches <- matches[order(matches$pepid),]
        }
        matches <- matches[order(matches$Protein, matches$Match, decreasing = FALSE),]
        matches$Y <- matches$Y1
        matches$pepid <- as.character(matches$pepid)
        matches$Match <- as.character(matches$Match)
        matches$ID <- do.call(paste, c(matches[, c("pepid", "Match", "Protein")], sep = "_"))
        matches$Match <- as.integer(matches$Match)
        #
        # ---- Code chunk developed with istaGPT and chatGPT - START ---- #
        Xextent <- 0.3
        if (Mode == "Align2") {
          matches_list <- lapply(split(matches, matches$Protein), \(sub) {
            ## One entry per peptide
            peps <- split(sub, sub$ID)
            pepinfo <- do.call(rbind, lapply(peps, \(p) {
              data.frame(ID = p$ID[1L],
                         start = min((p$Y1 - 1L) * scale + p$X1),
                         length = sum(p$X2 - p$X1 + 1L),
                         stringsAsFactors = FALSE)
            }))
            ## Left-to-right, longest first
            pepinfo <- pepinfo[order(pepinfo$start, -pepinfo$length),]
            ## Each lane contains occupied segments:
            ## row, X1, X2
            lanes <- list()
            pepinfo$offset <- integer(nrow(pepinfo))
            for (i in seq_len(nrow(pepinfo))) {
              pep <- peps[[pepinfo$ID[i]]]
              assigned <- FALSE
              ## Try existing offsets
              for(k in seq_along(lanes)) {
                overlap <- FALSE
                for (j in seq_len(nrow(pep))) {
                  seg <- pep[j,]
                  occ <- lanes[[k]][[as.character(seg$Y1)]]
                  if (!is.null(occ)) {
                    if(any(occ$X1 <= seg$X2 & occ$X2 >= seg$X1)) {
                      overlap <- TRUE
                      break
                    }
                  }
                }
                if (!overlap) {
                  for (j in seq_len(nrow(pep))) {
                    seg <- pep[j,]
                    rowid <- as.character(seg$Y1)
                    lanes[[k]][[rowid]] <- if (is.null(lanes[[k]][[rowid]])) {
                      data.frame(X1 = seg$X1,
                                 X2 = seg$X2)
                    } else {
                      rbind(lanes[[k]][[rowid]],
                            data.frame(
                              X1 = seg$X1,
                              X2 = seg$X2))
                    }
                  }
                  pepinfo$offset[i] <- k - 1L
                  assigned <- TRUE
                  break
                }
              }
              ## Create a new offset if necessary
              if (!assigned) {
                k <- length(lanes) + 1L
                lanes[[k]] <- list()
                for (j in seq_len(nrow(pep))) {
                  seg <- pep[j,]
                  lanes[[k]][[as.character(seg$Y1)]] <- data.frame(X1 = seg$X1,
                                                                   X2 = seg$X2)
                }
                pepinfo$offset[i] <- k - 1L
              }
            }
            sub$offset <- pepinfo$offset[match(sub$ID, pepinfo$ID)]
            return(sub)
          })
          matches <- do.call(rbind, matches_list)
          #
          # ---- Code chunk developed with istaGPT and chatGPT - END ---- #
          #
          interPepDist <- 1/(max(matches$offset)*1.5)
          matches$Y <- matches$Y - matches$offset*interPepDist - 0.05
          Yextent <- 0.4*interPepDist
          if (Yextent < 0.01) {
            # (heuristic: a solution would be to change the aspect ratio of the plot... maybe for a future edit?)
            warning("There are so many peptides that it is becoming hard to plot all! Setting peptide segment height to 0.001!")
            Yextent <- 0.01
          }
          matches$X1 <- matches$X1 - Xextent
          matches$X2 <- matches$X2 + Xextent
          matches$Y1 <- matches$Y + Yextent
          matches$Y2 <- matches$Y - Yextent
        }
        #
        # Restore correct I/L sequence
        wAA <- which(align.temp$AA %in% AA)
        align.temp$AA[wAA] <- align.temp$Amino_acid[wAA]
        align.temp$Amino_acid <- NULL
        #
        if (Mode == "Align2") {
          yxRat <- yxRat2 <- 3
          covPlot <- ggplot2::ggplot(align.temp) + ggplot2::coord_fixed(yxRat) +
            ggplot2::geom_text(ggplot2::aes(label = AA, x = X, y = Y, fontface = Face, hjust = Hjust),
                               colour = "white", cex = size, vjust = 0) + 
            ggplot2::geom_rect(data = matches,
                               ggplot2::aes(xmin = X1, xmax = X2, ymin = Y2, ymax = Y1, fill = Intensity, peptide = peptide),
                               show.legend = (length(unique(matches$Intensity)) > 1L)) +
            ggplot2::scale_fill_gradient(low = colscale[1L], high = colscale[2L], na.value = na,
                                         breaks = brks, labels = brkstxt, limits = c(0, maxInt)) +
            ggplot2::scale_x_discrete(breaks = NULL) + 
            ggplot2::scale_y_discrete(breaks = NULL#, labels = YlabZ
            ) +
            ggplot2::labs(x = "", y = "Position") +
            ggplot2::ggtitle(title, subtitle = covStr) +
            ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                           panel.grid.minor = ggplot2::element_blank(),
                           panel.background = ggplot2::element_rect(fill = bgcol),
                           plot.background = ggplot2::element_rect(fill = bgcol),
                           legend.title = ggplot2::element_blank(),
                           axis.title.y = ggplot2::element_text(colour = "white"),
                           axis.text.y = ggplot2::element_blank(),
                           axis.ticks.y = ggplot2::element_blank(),
                           legend.background = ggplot2::element_rect(fill = bgcol),
                           legend.text = ggplot2::element_text(colour = "white"),
                           title = ggplot2::element_text(colour = "white"),
                           plot.margin = ggplot2::margin(1L, 4L, 1L, 1L, "cm"))
          plotCreated <- TRUE
          #poplot(covPlot, new.window = new.window)
        }
        wMods <- grep("\\(", matches$Match.seq)
        if (length(wMods)) {
          matches$mods <- as.list(rep(NA, nrow(matches)))
          tmp <- data.frame(seq = paste0("_", matches$Match.seq[wMods], "_"),
                            X = matches$X1[wMods])
          if (Mode == "Align2") {
            tmp$X <- tmp$X+Xextent
          }
          matches$mods[wMods] <- lapply(1L:nrow(tmp), \(x) {#x <- 1L
            x <- tmp[x,]
            loc <- x[[2L]]
            x <- annot_to_tabl(x[[1L]])[[1L]]
            l <- nrow(x)-2L
            wh <- which(x$Annotations != "")
            if (length(wh)) {
              x <- x[wh, , drop = FALSE]
              x$Which <- wh-1L
              x$Type <- vapply(x$Which, \(y) { c("Nterm", "Internal", "Cterm")[unlist(which(c(y == 0L,
                                                                                              (y > 0L)&(y < l+1L),
                                                                                              y == l+1L)))] }, "")
              x$Which <- as.character(x$Which)
              x <- apply(x[, c("Annotations", "Which", "Type")], 1L, paste, collapse = "___")
              x <- paste0(x, "___", loc)
            } else { x <- NA_real_ }
            return(x)
          })
          mods <- listMelt(matches$mods[wMods], wMods, c("value", "Match"))
          mods$Y <- matches$Y[mods$Match]
          mods$Intensity <- matches$Intensity[mods$Match]
          mods[, c("Modification", "Xb", "Type", "Xa")] <- Isapply(strsplit(mods$value, "___"), unlist)
          mods$X <- as.numeric(mods$Xa) + as.numeric(mods$Xb) - 1L
          mods$X <- mods$X + c(0.5, 0, -0.5)[match(mods$Type, c("Nterm", "Internal", "Cterm"))]
          if (Mode == "Align2") {
            mods$Y <- as.numeric(mods$Y) + interPepDist*0.2
          }
          mods$Modification <- gsub(" \\(Internal\\)$", "", apply(mods[, c("Modification", "Type")], 1L, \(x) {
            paste0(x[[1L]], " (", x[[2L]], ")")
          }))
          #tmpshp <- unique(mods$Modification)
          #myShapes <- setNames(substr(tmpshp, 1L, 1L), tmpshp)
          #tstshp <- aggregate(myShapes, list(myShapes), length)
          #w <- which(tstshp$x > 1L)
          #if (length(w)) {
          #  w <- which(myShapes %in% tstshp$Group.1[w])
          #  myShapes[w] <- 1L:length(w)
          #}
          if (Mode == "Align2") {
            covPlot <- covPlot +
              ggplot2::geom_point(data = mods,
                                  ggplot2::aes(x = X, y = Y,
                                               colour = Modification, shape = Modification,
                                               PTM = Modification),
                                  size = size*0.33) +
              ggplot2::theme(legend.background = ggplot2::element_blank()) # +
            #ggplot2::geom_text(data = mods,
            #                   ggplot2::aes(x = X, y = Y, label = Modification),
            #                   size = size*0.8)
            #poplot(covPlot, new.window = new.window)
          }
        }
        temp <- setNames(lapply(nmz, \(n)  {
          wm <- which(matches$Protein == n)
          x <- if (length(wm)) { aggregate(matches$ID[wm], list(matches$ID[wm]), length) } else { NA }
          return(x)
        }), nmz)
        whinna <- which(!is.na(temp))
        if (length(whinna)) {
          dotstest <- plyr::rbind.fill(temp[whinna])
          dw <- which(dotstest$x > 1L)
          if (length(dw)) {
            dotstest <- dotstest[dw, , drop = FALSE]
            dotstest[, c("X", "Y")] <- Isapply(dotstest$Group.1, \(x) {
              matches[which(matches$ID == x), c("X1", "Y")]
            })
            dotstest[, c("Left", "Right")] <- as.data.frame(t(apply(dotstest[, c("X", "Y")],
                                                                    1L,
                                                                    \(x) { # x <- dotstest[1L, c("X", "Y")]
                                                                      w <- which(unlist(x[[1L]]) == 1L-Xextent)
                                                                      w <- w[which(w >= 2L)]
                                                                      return(c(unlist(x[[2L]])[w], unlist(x[[2L]])[w-1L]))
                                                                    })))
            Left <- unlist(dotstest$Left); ll <- length(Left)
            Right <- unlist(dotstest$Right); lr <- length(Right)
            dots <- list()
            if (ll) { dots$"LL" <- data.frame(Y = Left,
                                              X = 1-Xextent-0.5,
                                              Label = "...") }
            if (lr) { dots$"RR" <- data.frame(Y = Right,
                                              X = scale+Xextent+0.5,
                                              Label = "...") }
            dots <- plyr::rbind.fill(dots)
            if ((ll+lr) && (Mode == "Align2")) {
              covPlot <- covPlot +
                ggplot2::geom_text(data = dots,
                                   ggplot2::aes(x = X, y = Y, label = Label),
                                   vjust = 0, colour = "white", cex = size)
              #poplot(covPlot, new.window = new.window)
            }
          }
        }
        if (Mode == "Align2") {
          Xleft <- -3
          Xright <- ceiling(max(align.temp$X, na.rm = TRUE)) + 4L
          Ytop <- max(align.temp$Y, na.rm = TRUE) + 3L
          Ybottom <- min(align.temp$Y, na.rm = TRUE) - 3L
          covPlot <- covPlot + ggplot2::xlim(Xleft, Xright) + ggplot2::ylim(Ybottom, Ytop)
        }
        if (Mode == "Heat") {
          align.temp2 <- align.temp
          align.temp2$Intensity <- NA_real_
          matches2 <- matches
          matches2$Match <- as.integer(matches2$Match)
          matches2$Length <- as.integer(matches2$Length)
          w <- which(!is.na(matches2$Intensity))
          if (length(w)) {
            for (i in w) { #i <- 1L
              rg <- matches2$Match[i]-1L+(matches2$Start[i]:matches2$End[i])
              rg <- which((align.temp2$N %in% rg)&(align.temp2$Protein == matches2$Protein[i]))
              stopifnot(length(rg) == matches2$Length[i])
              #stopifnot(paste(align.temp2$AA[rg], collapse = "") == matches2$stripOrig[i])
              if (I_eq_L) {
                stopifnot(gsub("I", "L", paste(align.temp2$AA[rg], collapse = "")) == gsub("I", "L", matches2$stripOrig[i]))
              } else {
                stopifnot(paste(align.temp2$AA[rg], collapse = "") == matches2$stripOrig[i])
              }
              #align.temp2$AA[rg]
              tmp <- data.frame(A = align.temp2$Intensity[rg],
                                B = matches2$Intensity[i])
              tmp <- rowSums(tmp, na.rm = TRUE)
              align.temp2$Intensity[rg] <- tmp
            }
          }
          # We need to calculate maxInt here for this type of plot!
          maxInt <- if (misFun(maxInt)) { NA } else { suppressWarnings(as.numeric(maxInt)) }
          maxInt <- max(c(maxInt, align.temp2$Intensity), na.rm = TRUE)
          maxInt <- ceiling(maxInt)
          l <- nchar(as.character(maxInt))
          n <- min(c(2, l))
          maxInt2 <- ceiling(as.numeric(substr(as.character(maxInt), 1L, n))/(10L^(n-1L)))*(10L^(l-1L))
          brks <- sort(unique(c(c(0L:4L)*(maxInt2/4), maxInt)))
          #brks <- brks[which(brks <= maxInt)]
          brkstxt <- as.character(brks)
          w <- which(brkstxt == as.character(maxInt))
          brkstxt[w] <- paste0(" ---> max. int.: ", brkstxt[w])
          #
          nSubRows <- 3L # One each for the sequence, coverage boxes and some spacing...
          Xleft <- -1L
          if (length(wMods)) { # ... and add to this one per PTM
            uMods <- unique(mods$Modification)
            nMods <- length(uMods)
            #unique(align.temp2$Y)
            nSubRows <- nSubRows + nMods
            Xleft <- -max(nchar(uMods), na.rm = TRUE)/3-2
          }
          align.temp2$Y <- align.temp2$Y*(nSubRows)
          #unique(align.temp2$Y)
          wAA <- which(align.temp2$AA %in% AA)
          yxRat <- 1
          yxRat2 <- 3
          wFnd <- wAA[which((!is.na(align.temp2$Intensity[wAA]))
                            &(align.temp2$Legend[wAA] == "found"))]
          covPlot <- ggplot2::ggplot(align.temp2) +
            ggplot2::coord_fixed(yxRat) +
            ggplot2::geom_text(ggplot2::aes(label = AA, x = X, y = Y, fontface = Face, hjust = Hjust),
                               colour = "white", cex = size, vjust = 0) +
            # geom_rtile(data = align.temp2[wFnd,],
            #            ggplot2::aes(x = X, y = Y-1, fill = Intensity),
            #            width = 0.9, height = 0.9,
            #            radius = ggplot2::unit(3, "pt"))
            ggplot2::geom_tile(data = align.temp2[wFnd,],
                               ggplot2::aes(x = X, y = Y-1L, fill = Intensity),
                               width = 0.9, height = 0.9)
          #poplot(covPlot, 12, 20, new.window)
          if (length(wMods)) {
            mods2 <- mods
            mods2$Y <- mods2$Y*(nSubRows) - match(mods2$Modification, uMods)
            uModsDf <- plyr::rbind.fill(lapply(unique(align.temp2$Y[wAA]), \(y) {
              data.frame(Mod = uMods, Y = y-1L-1L:nMods, X = -1L)
            }))
            covPlot <- covPlot +
              ggplot2::geom_text(data = uModsDf,
                                 ggplot2::aes(label = Mod, x = X, y = Y),
                                 hjust = 1, vjust = 0.5, fontface = "italic",
                                 colour = "white", cex = size*1.1) +
              # geom_rtile(data = mods2,
              #            ggplot2::aes(x = X, y = Y-1, fill = Intensity, colour = Modification),
              #            width = 0.9, height = 0.9,
              #            radius = ggplot2::unit(3, "pt"))
              ggplot2::geom_tile(data = mods2,
                                 ggplot2::aes(x = X, y = Y-1L, fill = Intensity, colour = Modification),
                                 width = 0.9, height = 0.9#, cex = 0.5
              )
            #poplot(covPlot, 12, 20, new.window)
            Xright <- ceiling(max(c(align.temp2$X, mods2$X, uModsDf$X), na.rm = TRUE)) + 2L
            Ybottom <- floor(min(c(align.temp2$Y, mods2$Y, uModsDf$Y, 0L), na.rm = TRUE)) - 3L
            Ytop <- ceiling(max(c(align.temp2$Y, mods2$Y, uModsDf$Y), na.rm = TRUE)) + 1L
          } else {
            Xright <- ceiling(max(align.temp2$X, na.rm = TRUE)) + 2L
            Ybottom <- floor(min(c(0L, align.temp2$Y), na.rm = TRUE)) - 1L
            Ytop <- ceiling(max(align.temp2$Y, na.rm = TRUE)) + 1L
          }
          Ytop <- max(align.temp2$Y, na.rm = TRUE) + 1L
          covPlot <- covPlot +
            ggplot2::scale_fill_gradient(low = colscale[1L], high = colscale[2L], na.value = na,
                                         breaks = brks, labels = brkstxt,
                                         limits = c(0L, max(as.numeric(brks)))) +
            ggplot2::scale_x_continuous(limits = c(Xleft, Xright),
                                        breaks = NULL) +
            ggplot2::scale_y_continuous(limits = c(Ybottom, Ytop), breaks = NULL#, labels = YlabZ
            ) +
            ggplot2::labs(x = "", y = "") +
            ggplot2::ggtitle(title, subtitle = covStr) +
            ggplot2::guides(colour = "none") +
            ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                           panel.grid.minor = ggplot2::element_blank(),
                           panel.background = ggplot2::element_rect(fill = bgcol),
                           plot.background = ggplot2::element_rect(fill = bgcol),
                           legend.title = ggplot2::element_blank(),
                           axis.title.y = ggplot2::element_text(colour = "white"),
                           legend.background = ggplot2::element_rect(fill = bgcol),
                           legend.text = ggplot2::element_text(colour = "white"),
                           title = ggplot2::element_text(colour = "white"),
                           plot.margin = ggplot2::margin(1L, 4L, 1L, 1L, "cm"))
          plotCreated <- TRUE
          #poplot(covPlot, 12, 20, new.window)
        }
      } else {
        warning(paste0("Not a single peptide matches the sequence(s)!"))
      }
    }
    if (plotCreated) {
      wdth <- 20L
      xRg <- Xright-Xleft
      yRg <- Ytop-Ybottom
      hght <- max(c(wdth*yRg/(xRg*yxRat2), 10))
      wdth <- max(c(wdth, hght*0.5))
      if (display) {
        poplot(covPlot, height = hght, width = wdth, new.window = new.window)
      }
      if (!save.grph) { return() }
      for (svpth in unlist(tbl$Path)) { #svpth <- unlist(tbl$Path)[1L]
        ext <- rev(unlist(strsplit(svpth, "\\.")))[1L]
        suppressMessages({
          if (ext == "html") {
            covPlotLy <- plotly::ggplotly(covPlot, tooltip = c("peptide", "PTM"))
            if (!P %in% names(plotLyst)) { plotLyst[[P]] <- list() }
            plotLyst[[P]][[Ttl]] <- covPlotLy
            htmlwidgets::saveWidget(covPlotLy, svpth, selfcontained = TRUE)
          } else {
            ggplot2::ggsave(svpth, covPlot,
                            dpi = 300L,
                            width = wdth,
                            height = hght,
                            units = "in",
                            limitsize = FALSE)
          }
        })
        #openwd(dirname(svpth))
      }
    }
    if (Mode == "XML") {
      tst <- grep("^Protein name: ", align.temp$AA)
      XML_Cov <- lapply(1L:length(tst), \(x) { #x <- 1L
        x <- align.temp[(tst[x]+2L):(c(tst, nrow(align.temp)+1L)[x+1L]-1L),
                        c("Amino_acid", "Legend")]
        x <- x[which(x$Amino_acid != ""),]
        u <- unique(x$Legend)
        w <- c(0L, which(vapply(1L:(nrow(x)-1L), \(y) {
          x$Legend[y] != x$Legend[y+1L] # These are the ends of each stretch
        }, TRUE)), length(x$Legend))
        rs <- as.data.frame(t(sapply(1L:(length(w)-1L), \(y) {
          rg <- (w[y]+1L):w[y+1L]
          c(paste0(x$Amino_acid[rg], collapse = ""), x$Legend[rg[1L]])
        })))
        rs$V2 <- c(1L, 2L)[match(rs$V2, c("not found", "found"))]
        rs <- apply(rs[, c("V1", "V2")], 1L, \(y) {
          #y <- rs[1, c("V1", "V2")]
          z <- as.numeric(y[[2L]])
          list(openxlsx2::fmt_txt(y[[1L]], bold = as.logical(z-1L),
                                  color = openxlsx2::wb_color(hex = c("grey", colour)[as.numeric(z)])))
        })
        txt <- paste0("rs <- ", paste0("rs[[", as.character(1L:length(rs)), "]][[1L]]", collapse = " + "))
        eval(parse(text = txt))
        return(rs)
      })
      if (length(XML_Cov) > 1L) {
        rs2 <- openxlsx2::fmt_txt(";")
        txt <- paste0("XML_Cov <- ", paste0("XML_Cov[[", as.character(1L:length(XML_Cov)), "]][[1L]]", collapse = " + rs2 + "))
        eval(parse(text = txt))
      }
      xmList[[Ttl]] <- XML_Cov
    }
  }
  if (Mode == "XML") {
    return(XML_Cov)
  }
  if ((Mode == "Align2") && length(plotLyst)) {
    return(plotLyst)
  }
}
