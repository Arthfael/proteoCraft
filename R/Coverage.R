#' Coverage
#'
#' @description
#' This function can be used in 5 different modes:
#' - Coverage: calculates the proportion (from 1 to 0) of parent protein sequence(s) covered by the peptide sequence(s) provided.
#' - Align: plot the protein(s), highlighting amino acids covered by the peptides sequence(s) provided.
#' - Align2: as "Align", but peptides are displayed individually as bars; accepts PTM-modified sequences and will display PTM positions on the map; in addition, peptide colour can be mapped to the "intensities" argument.
#' - Heat: Amino acid color reflects the sum of intensities for all peptides observed covering said position. There is in addition one row for each PTM - also with color mapped to intensity.
#' - XML: formats covered sequence as XML to display coverage in Excel
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
#' @param colscale Only used if Mode = "Align2". Which colour scale do you want to use for peptide intensities? Default = 1 ("darkblue"->"orange"). Other values are 2: "darkblue"->"red", 3: "darkblue"->"green", 4: "black"->"blue", 5: "black"->"orange", 6: "black"->"red", 7: "black"->"green" and 8: "red"->"green"
#' @param size Only used if Mode = "Align", "Align2" or "Heat". The character size. Default = 2.5, you may want to reduce it for large proteins.
#' @param save Only used if Mode = "Align", "Align2" or "Heat". Preferred alternative to "save.path", though less flexible: vector of ggsave compatible file extensions to add to the "title" argument to save the graph to. Allows saving as multiple formats through a single function call. Default = "jpeg"; set to FALSE to not save it.
#' @param save.path Only used if Mode = "Align", "Align2" or "Heat". Prefer the "save" argument. If set, must be a character vector with the name of the file(s) in which to save each plot, so must have the same length as \"title\". The name must end in a ggsave-compatible extension. Default = FALSE
#' @param intensities Only used if Mode = "Align2" or "Heat". Values to map to the "colscale" argument so that peptides can be printed with different colours, e.g. mapped to abundance.
#' @param na Only used if Mode = "Align2" or "Heat". Colour for NA values.
#' @param maxInt If provided, the maximum of the intensity scale is not detected from the values in intensities but provided externally. Useful if drawing maps from several intensity vectors and wanting to apply a single scale to all.
#' @param bgcol Only used if Mode = "Align2" or "Heat". Which colour should the background be? Default = "black".
#' @param I_eq_L Should we consider I and L identical? Currently, by default, TRUE for both DIA and DDA: see https://github.com/vdemichev/DiaNN/discussions/1631
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
#' # Notes:
#' # * Never use Mode "Align" with sapply at it will create a graph for each protein!
#' #   Instead, this function works with a vector of several proteins sequences as input.
#' #   For the Align Mode, the display argument (ignored in Coverage Mode) should be a logical determining which proteins should be displayed.
#' # * In general, since at present the Align Mode makes some assumptions about the length of the protein sequence(s) submitted, I would avoid applying it to:
#' #   - too many or too large proteins;
#' #   - too many peptides.
#' @export

Coverage <- function(proteins,
                     peptides,
                     Mode = "Coverage",
                     new.window = TRUE, 
                     display = TRUE,
                     scale = 60,
                     title = "Coverage",
                     colour = "red",
                     colscale = 1,
                     size = 2.5,
                     save = "jpeg",
                     save.path,
                     intensities = NULL,
                     na = "red",
                     maxInt = NULL,
                     bgcol = "black",
                     I_eq_L) {
  TESTING <- FALSE
  #proteoCraft::DefArg(proteoCraft::Coverage);TESTING = TRUE
  #proteins = P; peptides = tmp$"Modified sequence"; Mode = "Align2"; title = ttl; save = c("jpeg", "pdf"); intensities = tmp$`log10(Intensity)`
  #proteins = x[[1]]; peptides = x[[2]]; Mode = "XML"; colour = "green"
  #proteins = sq; peptides = tmpSq$"Modified sequence_verbose"; intensities = tmpSq$Intensity; Mode = "Align2"; save = FALSE
  #proteins = sq; peptides = tmpSq$"Modified sequence_verbose"; intensities = tmpSq$Intensity; Mode = "Heat"; save = FALSE
  #proteins = seq; peptides = p1$Sequence; Mode = "Align2"; title = paste0("Coverage map - ", nm); save = c("jpeg", "pdf"); intensities = p1$Intensity; display = FALSE
  if (TESTING) {
    # Note:
    # This is not a perfect alternative to missing but will work in most cases, unless x matches a function imported by a package 
    misFun <- function(x) { return(!exists(deparse(substitute(x)))) }
  } else { misFun <- missing }
  #
  if ((misFun(I_eq_L))||(!is.logical(I_eq_L))||(is.na(I_eq_L))) {
    if ((exists("isDIA"))&&(is.logical(isDIA))&&(!is.na(isDIA))) {
      I_eq_L <- !idDIA
    } else {
      I_eq_L <- TRUE
    }
  }
  if (!Mode %in% c("Coverage", "Align", "Align2", "XML", "Heat")) {
    stop("Accepted value for \"Mode\" argument:\n - \"Coverage\": only return percentage of sequence coverage\n - \"Align\": create simple peptide coverage map\n - \"Align2\": create peptide coverage map with individual peptide intensities\n - \"XML\": writes xml formated coverage text for writing into Excel tables.")
  }
  if ((Mode == "Align")&&(colour == "black")) {
    warning("Any R-compatible colour but \"black\" is allowed - though I would advise something shiny...\nDefaulting to \"red\".\nYou're welcome!")
    colour <- "red"
  }
  # if (Mode == "Heat") {
  #   # Good idea in principle, but doesn't work for saving
  #   if (!require(statebins)) { install.packages("statebins") }
  #   require(statebins)
  #   # Thanks to https://stackoverflow.com/questions/64355877/round-corners-in-ggplots-geom-tile-possible("statebins") !!!
  #   `%||%` <- function(a, b) { if (is.null(a)) { b } else { a } }
  #   GeomRtile <- ggplot2::ggproto("GeomRtile", 
  #                        statebins:::GeomRrect, # 1) only change compared to ggplot2:::GeomTile
  #                        extra_params = c("na.rm"),
  #                        setup_data = function(data, params) {
  #                          data$width <- data$width %||% params$width %||% ggplot2::resolution(data$x, FALSE)
  #                          data$height <- data$height %||% params$height %||% ggplot2::resolution(data$y, FALSE)
  #                          
  #                          transform(data,
  #                                    xmin = x - width / 2,  xmax = x + width / 2,  width = NULL,
  #                                    ymin = y - height / 2, ymax = y + height / 2, height = NULL
  #                          )
  #                        },
  #                        default_aes = ggplot2::aes(
  #                          fill = "grey20", colour = NA, size = 0.1, linetype = 1,
  #                          alpha = NA, width = NA, height = NA
  #                        ),
  #                        required_aes = c("x", "y"),
  #                        
  #                        # These aes columns are created by setup_data(). They need to be listed here so
  #                        # that GeomRect$handle_na() properly removes any bars that fall outside the defined
  #                        # limits, not just those for which x and y are outside the limits
  #                        non_missing_aes = c("xmin", "xmax", "ymin", "ymax"),
  #                        draw_key = ggplot2::draw_key_polygon
  #   )
  #   geom_rtile <- function(mapping = NULL,
  #                          data = NULL,
  #                          stat = "identity",
  #                          position = "identity",
  #                          radius = grid::unit(6, "pt"), # 2) add radius argument
  #                          ...,
  #                          na.rm = FALSE,
  #                          show.legend = NA,
  #                          inherit.aes = TRUE) {
  #     ggplot2::layer(
  #       data = data,
  #       mapping = mapping,
  #       stat = stat,
  #       geom = GeomRtile, # 3) use ggproto object here
  #       position = position,
  #       show.legend = show.legend,
  #       inherit.aes = inherit.aes,
  #       params = rlang::list2(
  #         radius = radius,
  #         na.rm = na.rm,
  #         ...
  #       )
  #     )
  #   }
  # }
  if (is.null(names(proteins))) { names(proteins) <- 1:length(proteins) }
  namez <- names(proteins)
  if (length(peptides)) {
    if (Mode != "Coverage") {
      Table <- data.frame(Title = title)
      if (length(title) == 1) { Table$"Protein(s)" <- list(namez) } else {
        if (length(title) != length(proteins)) {
          if (length(title) == 1) { title <- paste0(title, " - ", namez) } else {
            warning("\"title\"'s length must be either one of the same length as the protein vector, creating a default title!")
            title <- paste0("Coverage - ", namez)
            Table <- data.frame(Title = gsub(":|\\*|\\?|<|>|\\||/", "-", title))
          }
        }
        Table$"Protein(s)" <- as.list(namez)
      }
      if ((length(save) > 1)||(as.character(save) != "FALSE")) {
        if ((misFun(save.path))||(length(save.path) != length(title))) {
          if ((!misFun(save.path))&&(length(save.path) != length(title))) {
            warning("Argument \"save.path\" will be ignored as its length is not compatible with that of the \"title\" argument!")
          }
          if (as.character("save") == "FALSE") { save.grph <- FALSE } else {
            Table$Path <- lapply(Table$Title, function(x) { paste0(x, ".", save) })
            save.grph <- TRUE
          }
        } else {
          if (as.character("save") == "FALSE") {
            Table$Path <- save.path
            save.grph <- TRUE
          } else {
            #warning("Both the \"save.path\" and \"save\" arguments were provided. I will use \"save.path\" for the file name but \"save\" for the extension(s).")
            Table$Path <- lapply(gsub("\\.[^A-Z,a-z,0-9]+$", "", save.path), function(x) { paste0(x, ".", save) })
            save.grph <- TRUE
          }
        }
      } else { save.grph <- FALSE }
    }
    if ((misFun(intensities))||(is.null(intensities))) { intensities <- rep(1, length(peptides)) } else {
      intensities <- as.numeric(unlist(intensities))
      if (length(intensities) != length(peptides)) {
        stop("The \"intensities\" vector should be the same length as the \"peptides\" vector!")
      }
      w2 <- which(!proteoCraft::is.all.good(intensities, 2))
      intensities[w2] <- NA
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
      if (misFun(maxInt)) { maxInt <- max(intensities, na.rm = TRUE) }
      if ((is.na("maxInt"))||(!is.finite(maxInt))) { maxInt <- max(intensities, na.rm = TRUE) }
      maxInt <- ceiling(as.numeric(maxInt))
      l <- nchar(as.character(maxInt))
      maxInt2 <- round(maxInt, 1-l)
      while (maxInt2 < maxInt) { maxInt2 <- maxInt2 + 10^(l-1) }
      brks <- sort(unique(c(c(0:4)*(maxInt2/4), maxInt)))
      #brks <- brks[which(brks <= maxInt)]
      brkstxt <- as.character(brks)
      w <- which(brkstxt == as.character(maxInt))
      brkstxt[w] <- paste0(" ---> max. int.: ", brkstxt[w])
    }
    peptides <- data.frame(Original.seq = peptidesOrig,
                           Match.seq = peptides,
                           Sequence = gsub("\\([^\\)]+\\)", "", peptides),
                           id = 1:length(peptides),
                           stripped.Original.seq = gsub("\\([^\\)]+\\)", "", peptidesOrig))
    if (!is.null(intensities)) { peptides$Intensity <- intensities }
    peptides$Length <- nchar(peptides$Sequence)
    seq <- unique(peptides$Sequence)
    coverage <- c()
    pos <- list()
    pos2 <- data.frame(Seq = seq)
    for (P in namez) { #P <- namez[1]
      protein <- paste0("_", proteins[P], "_")
      proteinOrig <- paste0("_", proteinsOrig[P], "_")
      pos2[[P]] <- lapply(seq, function(sq) { #sq <- seq[1]
        temp <- unlist(strsplit(protein, sq))
        l <- length(temp)
        if (l > 1) { res <- cumsum(nchar(temp[1:l-1])) + cumsum(c(0, rep(nchar(sq), l-2))) } else { res <- NA }
        return(res)
      })
      temp <- setNames(lapply(1:length(seq), function(x) {
        sq <- seq[x]
        x <- data.frame(start = pos2[[match(sq, pos2$Seq), P]])
        w <- which(!is.na(x$start))
        x <- x[w, , drop = FALSE]
        x$end <- x$start + nchar(sq) - 1
        x$spread <- apply(x[,c("start", "end")], 1, function(x) { list(x[[1]]:x[[2]])})
        return(x)
      }), seq)
      spread <- unique(unlist(sapply(seq, function(sq) { temp[[sq]]$spread })))
      res2 <- rep(0, nchar(protein)-2)
      res2[spread] <- 1
      #setNames(res2, unlist(strsplit(substr(protein, 2, nchar(protein)-1), "")))
      pos[[P]] <- res2
      coverage[P] <- mean(pos[[P]])
    }
    if (Mode != c("Coverage")) {
      if (Mode %in% c("Align2", "Heat")) {
        colscales <- rbind(c("darkblue", "orange"),
                           c("darkblue", "red"),
                           c("darkblue", "green"),
                           c("black", "blue"),
                           c("black", "orange"),
                           c("black", "red"),
                           c("black", "green"),
                           c("red", "green"))
        if (!suppressWarnings(as.numeric(colscale[1])) %in% 1:nrow(colscales)) {
          warning(paste0("Argument \"colscale\" must be a (single) number between 1 and ", nrow(colscales), ", defaulting to 1!"))
          colscale <- 1
        }
        colscale <- colscales[colscale,]
      }
      for (Ttl in unique(title)) { #Ttl <- unique(title)[1]
        plotCreated <- FALSE
        tbl <- Table[which(Table$Title == Ttl),]
        nmz <- unlist(tbl$Protein)
        prots <- proteins[nmz]
        protsOrig <- proteinsOrig[nmz]
        poZ2 <- pos2[,c("Seq", nmz)]
        align.temp <- list()
        YlabZ <- list()
        for (P in nmz) { #P <- nmz[1]
          protein <- prots[P]
          proteinOrig <- protsOrig[P]
          tmpAlgn <- data.frame(AA = unlist(strsplit(protein, "")))
          tmpAlgn$N <- 1:nrow(tmpAlgn)
          tmpAlgn$Amino_acid <- unlist(strsplit(proteinOrig, ""))
          tmpAlgn$Protein <- P
          tmpAlgn$Legend <- c("not found", "found")[pos[[P]] + 1]#poZ[[P]] + 1]
          tmpAlgn$X <- ((1:length(pos[[P]])) - 1) %% scale + 1#poZ[[P]])) - 1) %% scale + 1
          tmpAlgn$Y <- -sapply(1:length(pos[[P]]) - 1, function(x) {x - (x %% scale)})/scale#Z[[P]]) - 1, function(x) {x - (x %% scale)})/scale
          tmpAlgn$Y <- -min(tmpAlgn$Y) + 1 + tmpAlgn$Y
          tmpAlgn$Face <- sapply(pos[[P]], function(x) {c("plain", "bold")[x+1]})#poZ[[P]], function(x) {c("plain", "bold")[x+1]})
          tmpAlgn$Hjust <- 0.5
          YlabZ[[P]] <- c(rev(c(c(0:(max(tmpAlgn$Y, na.rm = TRUE) - 1)) * scale + 1)), "")
          covStr <- paste0("Coverage: ", round(100*coverage[P], 1), "%")
          align.temp[[P]] <- rbind(tmpAlgn,
                                   data.frame(AA = paste0("Protein name: ",  P),
                                              N = NA,
                                              Amino_acid = "",
                                              Protein = "",
                                              Legend = "not found",
                                              X = scale/2,
                                              Y = max(tmpAlgn$Y + 1.2),
                                              Face = "bold.italic",
                                              Hjust = 0.5),
                                   data.frame(AA = covStr,
                                              N = NA,
                                              Amino_acid = "",
                                              Protein = "",
                                              Legend = "not found",
                                              X = scale/2,
                                              Y = max(tmpAlgn$Y + 0.6),
                                              Face = "italic",
                                              Hjust = 0.5))
        }
        if (length(prots) > 1) {
          offset <- sapply(prots[2:length(prots)], function(x) { (nchar(x)/scale) }) + 2
          if (length(offset) > 1) {
            offset <- sapply(c(1:length(offset)), function(x) { sum(offset[x:length(offset)]) })
          }
          for (i in 1:length(offset)) { align.temp[[i]]$Y <- align.temp[[i]]$Y + offset[i] }
          for (i in 1:length(offset)) {
            align.temp[[1]] <- rbind(align.temp[[1]], align.temp[[i + 1]])
            YlabZ[[length(offset) + 1]] <- c(YlabZ[[length(offset) + 1]], "", YlabZ[[length(offset) + 1 - i]])
          }
          YlabZ <- c(YlabZ[[length(offset) + 1]], "")
        } else { YlabZ <- c(YlabZ[[1]], "") }
        align.temp <- align.temp[[1]]
        align.temp <- align.temp[order(align.temp$Y, decreasing = TRUE),]
        if (Mode == "Align") {
          myColors <- setNames(c("black", colour), c("not found", "found"))
          colScale <- ggplot2::scale_colour_manual(name = "colour", values = myColors)
          Xleft <- -3
          Xright <- ceiling(max(align.temp$X, na.rm = TRUE)) + 4
          Ytop <- max(align.temp$Y, na.rm = TRUE) + 3
          Ybottom <- min(align.temp$Y, na.rm = TRUE) - 3
          yxRat <- yxRat2 <- 2
          covplot <- ggplot2::ggplot(align.temp) + ggplot2::coord_fixed(yxRat) +
            ggplot2::geom_text(ggplot2::aes(label = AA,  x = X, y = Y, colour = Legend, fontface = Face,
                                            hjust = Hjust), cex = size) +
            colScale + ggplot2::scale_x_discrete(breaks = NULL) +
            ggplot2::scale_y_discrete(breaks = NULL, #labels = YlabZ
                                      ) +
            ggplot2::labs(x = "", y = "Position") +
            ggplot2::xlim(Xleft, Xright) + ggplot2::ylim(Ybottom, Ytop) +
            ggplot2::ggtitle(title, subtitle = covStr) +
            ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                           panel.grid.minor = ggplot2::element_blank(),
                           panel.background = ggplot2::element_blank(),
                           legend.title = ggplot2::element_blank(),
                  #aspect.ratio = 1.5
            )
          #proteoCraft::poplot(covplot)
          plotCreated <- TRUE
        }
        if (Mode %in% c("Align2", "Heat")) {
          matches <- lapply(nmz, function(n) { #n <- nmz[1]
            mtchs <- list(Outcome = FALSE)
            peptides[[paste0("Matches_", n)]] <<- lapply(peptides$Sequence, function(x) { #x <- peptides$Sequence[1]
              x <- poZ2[[match(x, poZ2$Seq), n]]
              return(x[which(!is.na(x))])
            })
            tst <- which(sapply(peptides[[paste0("Matches_", n)]], length) > 0)
            if (length(tst)) {
                mtchs <- peptides[tst, c(paste0("Matches_", n), "Match.seq", "Sequence", "id", "stripped.Original.seq")]
                colnames(mtchs) <- c("Match", "Match.seq", "Sequence", "pepid", "stripOrig")
                mtchs$Protein <- n
                mtchs <- list(Outcome = TRUE, Matches = mtchs)
            }
            return(mtchs)
          })
          matches <- matches[which(sapply(matches, function(x) { x$Outcome }))]
          if (length(matches)) {
            matches <- lapply(matches, function(x) { x$Matches })
            matches <- plyr::rbind.fill(matches)
            matches$temp <- do.call(paste, c(matches[, c( "Match.seq", "Sequence", "pepid", "Protein", "stripOrig")],
                                             sep = "_____"))
            matches <- proteoCraft::listMelt(matches$Match, matches$temp)
            colnames(matches)[1] <- "Match"
            matches[, c( "Match.seq", "Sequence", "pepid", "Protein", "stripOrig")] <- proteoCraft::Isapply(strsplit(matches$L1, "_____"), unlist)
            matches$L1 <- NULL
            matches$Length <- nchar(matches$Sequence)
            matches$pepid <- as.numeric(matches$pepid)
            matches[, c("X1", "Y1")] <- as.data.frame(t(apply(matches[,c("Protein", "Match")], 1, function(x) {
              unlist(align.temp[which(align.temp$Protein == as.character(x[[1]]))[as.numeric(x[[2]])], c("X", "Y")])
            })))
            matches[, c("X2", "Y2")] <- as.data.frame(t(apply(matches[,c("Length", "X1", "Y1")], 1, function(x) {
              x2 <- x[[2]]+x[[1]]-1
              y2 <- x[[3]] - floor(x2/scale)+1*(x2 %% scale == 0)
              x2 <- (x2 %% scale) + scale*(x2 %% scale == 0)
              return(c(x2, y2))
            })))
            if (!is.null(intensities)) { matches$Intensity <- peptides$Intensity[match(matches$pepid, peptides$id)] } else {
              matches$Intensity <- 1
            }
            # Deal with peptides spanning several rows
            matches$Start <- 1 # Not the match in the protein sequence! Used only for peptides spanning several rows!
            matches$End <- matches$Length # As above
            matches$Is_start <- TRUE
            matches$Is_end <- TRUE
            w <- which(matches$Y1 != matches$Y2) # Peptides spanning several rows
            if (length(w)) {
              for (i in w) { #i <- w[5]
                n <- matches$Y1[i]-matches$Y2[i] #(we almost never expect this to be more than 1: this would mean a super long peptide!!!)
                nr <- n+1
                temp <- matches[rep(i, nr),]
                temp$Y1 <- temp$Y2 <- max(temp$Y1)-(0:n)
                temp$X1[2:nr] <- 1
                temp$X2[1:(nr-1)] <- scale
                temp$Length <- temp$X2 - temp$X1 + 1
                temp2 <- as.data.frame(proteoCraft::annot_to_tabl(unique(temp$Match.seq),
                                                            Nterm = FALSE, Cterm = FALSE))
                colnames(temp2) <- c("A", "B")
                temp2$B[which(temp2$B != "")] <- paste0("(", temp2$B[which(temp2$B != "")], ")")
                temp2$A[which(temp2$A == "_")] <- ""
                rrange <- cumsum(temp$Length)+1
                temp$Match.seq <- sapply(1:nr, function(x) {
                  paste(apply(temp2[c(1, rrange+1)[x]:rrange[x], , drop = FALSE], 1, paste, collapse = ""), collapse = "")
                })
                temp$Sequence <- gsub("\\([^\\)]+\\)", "", temp$Match.seq)
                temp$End <- cumsum(temp$Length)
                temp$Start <- c(0, temp$End[1:(nr-1)]) + 1
                temp$stripOrig <- substr(temp$stripOrig, temp$Start, temp$End)
                temp$Is_start[2:nr] <- FALSE
                temp$Is_end[1:(nr-1)] <- FALSE
                matches[i,] <- temp[1,]
                matches <- rbind(matches, temp[2:(n+1),])
              }
            }
            matches <- matches[order(matches$Protein, matches$Match, decreasing = FALSE),]
            matches$Y <- matches$Y1
            matches$pepid <- as.character(matches$pepid)
            matches$Match <- as.character(matches$Match)
            matches$ID <- do.call(paste, c(matches[, c("pepid", "Match", "Protein")], sep = "_"))
            matches$Match <- as.integer(matches$Match)
            #olr <- 0
            olr <- 0.5
            # Let's calculate vertical offsets for overlapping peptides
            uniq <- setNames(lapply(1:nrow(matches), function(x) { #x <- 1
              x <- matches[x, c("X1", "X2", "Y")]
              paste0(as.character((x[[1]]-olr):(x[[2]]+olr)), "_", as.character(x[[3]]))
            }), matches$ID)
            uniq <- proteoCraft::listMelt(uniq)
            uniq <- aggregate(uniq$L1, list(uniq$value), list)
            uniq$Length <- sapply(uniq$x, length)
            offset0 <- 1.2
            matches$offset <- 0
            extents <- c()
            for (n in nmz) { #n <- nmz[1]#
              mt <- unique(matches$ID[which(matches$Protein == n)])
              for (i in mt) { #i <- mt[1] #i <- mt[2] #i <- mt[3]
                #w <- which((matches$ID == i)&(matches$Protein == n))
                w <- which(matches$ID == i)
                m <- matches[w, , drop = FALSE]
                rrange <- unique(unlist(apply(m[,c("X1", "X2", "Y")], 1, function(x) { paste0((x[[1]]-olr):(x[[2]]+olr), "_", x[[3]]) })))
                ol <- unique(unlist(uniq$x[which(uniq$Group.1 %in% rrange)]))
                ol <- ol[which(match(ol, mt) < match(i, mt))]
                l <- length(ol)
                if (l) {
                  #w2 <- which((matches$ID %in% ol)&(matches$Protein == n))
                  w2 <- which(matches$ID %in% ol)
                  exts <- lapply(0:max(matches$offset[w2]), function(x) {
                    unlist(apply(m[, c("X1", "X2", "Y"), drop = FALSE], 1, function(y) {
                      paste0((y[[1]]-olr):(y[[2]]+olr), "_", y[[3]], "_", x)
                    }))
                  })
                  tstext <- sapply(exts, function(x) { sum(x %in% extents) })
                  w3 <- which(tstext == 0)
                  if (length(w3)) { l <- w3[1]-1 }
                }
                matches$offset[w] <- l
                extents <- unique(c(extents,
                                    unlist(apply(matches[w, c("X1", "X2", "Y", "offset"), drop = FALSE], 1, function(x) { 
                                      paste0((x[[1]]-olr):(x[[2]]+olr), "_", x[[3]], "_", x[[4]])
                                    }))))
              }
            }
            matches$offset <-  matches$offset + offset0
            scale2 <- 1/(max(matches$offset)+4)
            matches$Y <- matches$Y-(matches$offset)*scale2
            Xextension <- 0.4
            Yextension <- 0.4*scale2
            matches$X1 <- matches$X1 - Xextension
            matches$X2 <- matches$X2 + Xextension
            matches$Y1 <- matches$Y + Yextension
            matches$Y2 <- matches$Y - Yextension
            #
            # Restore correct I/L sequence
            align.temp$AA <- NULL
            colnames(align.temp)[match("Amino_acid", colnames(align.temp))] <- "AA"
            if (Mode == "Align2") {
              yxRat <- yxRat2 <- 3
              covplot <- ggplot2::ggplot(align.temp) + ggplot2::coord_fixed(yxRat) +
                ggplot2::geom_text(ggplot2::aes(label = AA, x = X, y = Y, fontface = Face, hjust = Hjust),
                                   colour = "white", cex = size, vjust = 0) + 
                ggplot2::geom_rect(data = matches,
                                   ggplot2::aes(xmin = X1, xmax = X2, ymin = Y2, ymax = Y1, fill = Intensity),
                                   show.legend = (length(unique(matches$Intensity)) > 1)) +
                ggplot2::scale_fill_gradient(low = colscale[1], high = colscale[2], na.value = na,
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
                               plot.margin = ggplot2::margin(1, 4, 1, 1, "cm"))
              plotCreated <- TRUE
              #proteoCraft::poplot(covplot, new.window = new.window)
            }
            wMods <- grep("\\(", matches$Match.seq)
            if (length(wMods)) {
              matches$mods <- as.list(rep(NA, nrow(matches)))
              tmp <- data.frame(seq = paste0("_", matches$Match.seq[wMods], "_"),
                                X = matches$X1[wMods]+Xextension)
              matches$mods[wMods] <- lapply(1:nrow(tmp), function(x) {#x <- 1
                x <- tmp[x,]
                loc <- x[[2]]
                x <- proteoCraft::annot_to_tabl(x[[1]])[[1]]
                l <- nrow(x)-2
                wh <- which(x$Annotations != "")
                if (length(wh)) {
                  x <- x[wh, , drop = FALSE]
                  x$Which <- wh-1
                  x$Type <- sapply(wh-1, function(y) { c("Nterm", "Internal", "Cterm")[unlist(which(c(y == 0, (y > 0)&(y < l+1), y == l+1)))] })
                  x$Which <- as.character(x$Which)
                  x <- apply(x[,c("Annotations", "Which", "Type")], 1, paste, collapse = "___")
                  x <- paste0(x, "___", loc)
                } else { x <- NA }
                return(x)
              })
              mods <- proteoCraft::listMelt(matches$mods[wMods], wMods, c("value", "Match"))
              mods$Y <- matches$Y[mods$Match]
              mods$Intensity <- matches$Intensity[mods$Match]
              mods[, c("Modification", "Xb", "Type", "Xa")] <- proteoCraft::Isapply(strsplit(mods$value, "___"), unlist)
              mods$X <- as.numeric(mods$Xa) + as.numeric(mods$Xb) - 1
              mods$X <- mods$X + c(0.5, 0, -0.5)[match(mods$Type, c("Nterm", "Internal", "Cterm"))]
              mods$Y <- as.numeric(mods$Y) + scale2*0.2
              mods$Modification <- gsub(" \\(Internal\\)$", "", apply(mods[,c("Modification", "Type")], 1, function(x) {
                paste0(x[[1]], " (", x[[2]], ")")
              }))
              #tmpshp <- unique(mods$Modification)
              #myShapes <- setNames(substr(tmpshp, 1, 1), tmpshp)
              #tstshp <- aggregate(myShapes, list(myShapes), length)
              #w <- which(tstshp$x > 1)
              #if (length(w)) {
              #  w <- which(myShapes %in% tstshp$Group.1[w])
              #  myShapes[w] <- 1:length(w)
              #}
              if (Mode == "Align2") {
                covplot <- covplot +
                  ggplot2::geom_point(data = mods,
                                      ggplot2::aes(x = X, y = Y,
                                                   colour = Modification, shape = Modification),
                                      size = size*0.33) +
                  ggplot2::theme(legend.background = ggplot2::element_blank()) # +
                #ggplot2::geom_text(data = mods,
                #                   ggplot2::aes(x = X, y = Y, label = Modification),
                #                   size = size*0.8)
                #proteoCraft::poplot(covplot, new.window = new.window)
              }
            }
            temp <- setNames(lapply(nmz, function(n)  {
              wm <- which(matches$Protein == n)
              if (length(wm)) {
                x <- aggregate(matches$ID[wm], list(matches$ID[wm]), length)
              } else { x <- NA }
              return(x)
            }), nmz)
            whinna <- which(!is.na(temp))
            if (length(whinna)) {
              dotstest <- plyr::rbind.fill(temp[whinna])
              dw <- which(dotstest$x > 1)
              if (length(dw)) {
                dotstest <- dotstest[dw, , drop = FALSE]
                dotstest[, c("X", "Y")] <- proteoCraft::Isapply(dotstest$Group.1, function(x) {
                  matches[which(matches$ID == x), c("X1", "Y")]
                })
                dotstest[, c("Left", "Right")] <- as.data.frame(t(apply(dotstest[, c("X", "Y")], 1, function(x) { # x <- dotstest[1, c("X", "Y")]
                  w <- which(unlist(x[[1]]) == 1-Xextension)
                  w <- w[which(w >= 2)]
                  return(c(unlist(x[[2]])[w], unlist(x[[2]])[w-1]))
                })))
                Left <- unlist(dotstest$Left); ll <- length(Left)
                Right <- unlist(dotstest$Right); lr <- length(Right)
                dots <- list()
                if (ll) { dots$"LL" <- data.frame(Y = Left,
                                                  X = 1-Xextension-0.5,
                                                  Label = "...") }
                if (lr) { dots$"RR" <- data.frame(Y = Right,
                                                  X = scale+Xextension+0.5,
                                                  Label = "...") }
                dots <- plyr::rbind.fill(dots)
                if ((ll+lr)&&(Mode == "Align2")) {
                  covplot <- covplot +
                    ggplot2::geom_text(data = dots,
                                       ggplot2::aes(x = X, y = Y, label = Label),
                                       vjust = 0, colour = "white", cex = size)
                  #proteoCraft::poplot(covplot, new.window = new.window)
                }
              }
            }
            if (Mode == "Align2") {
              Xleft <- -3
              Xright <- ceiling(max(align.temp$X, na.rm = TRUE))+4
              Ytop <- max(align.temp$Y, na.rm = TRUE) + 3
              Ybottom <- min(align.temp$Y, na.rm = TRUE) - 3
              covplot <- covplot + ggplot2::xlim(Xleft, Xright) + ggplot2::ylim(Ybottom, Ytop)
            }
            if (Mode == "Heat") {
              align.temp2 <- align.temp
              align.temp2$Intensity <- NA
              matches2 <- matches
              matches2$Match <- as.integer(matches2$Match)
              matches2$Length <- as.integer(matches2$Length)
              w <- which(!is.na(matches2$Intensity))
              if (length(w)) {
                for (i in w) { #i <- 1
                  rg <- matches2$Match[i]-1+(matches2$Start[i]:matches2$End[i])
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
              if (misFun(maxInt)) {
                maxInt <- ceiling(max(align.temp2$Intensity, na.rm = TRUE))
              }
              maxInt <- as.numeric(maxInt)
              if ((is.na("maxInt"))||(!is.finite(maxInt))) {
                maxInt <- ceiling(max(align.temp2$Intensity, na.rm = TRUE))
              }
              l <- nchar(as.character(maxInt))
              n <- min(c(2, l))
              maxInt2 <- ceiling(as.numeric(substr(as.character(maxInt), 1, n))/(10^(n-1)))*(10^(l-1))
              brks <- sort(unique(c(c(0:4)*(maxInt2/4), maxInt)))
              #brks <- brks[which(brks <= maxInt)]
              brkstxt <- as.character(brks)
              w <- which(brkstxt == as.character(maxInt))
              brkstxt[w] <- paste0(" ---> max. int.: ", brkstxt[w])
              #
              offset2 <- 3
              Xleft <- -1
              if (length(wMods)) {
                mods2 <- mods
                uMods <- unique(mods2$Modification)
                nMods <- length(uMods)
                #unique(align.temp2$Y)
                offset2 <- offset2 + nMods
                Xleft <- -max(nchar(uMods), na.rm = TRUE)/3-2
              }
              align.temp2$Y <- align.temp2$Y*(offset2)
              #unique(align.temp2$Y)
              wAA <- which(align.temp2$AA %in% proteoCraft::AA)
              yxRat <- 1
              yxRat2 <- 3
              wFnd <- wAA[which((!is.na(align.temp2$Intensity[wAA]))
                                &(align.temp2$Legend[wAA] == "found"))]
              covplot <- ggplot2::ggplot(align.temp2) +
                ggplot2::coord_fixed(yxRat) +
                ggplot2::geom_text(ggplot2::aes(label = AA, x = X, y = Y, fontface = Face, hjust = Hjust),
                                   colour = "white", cex = size, vjust = 0) +
                # geom_rtile(data = align.temp2[wFnd,],
                #            ggplot2::aes(x = X, y = Y-1, fill = Intensity),
                #            width = 0.9, height = 0.9,
                #            radius = ggplot2::unit(3, "pt"))
                ggplot2::geom_tile(data = align.temp2[wFnd,],
                                   ggplot2::aes(x = X, y = Y-1, fill = Intensity),
                                   width = 0.9, height = 0.9)
              #proteoCraft::poplot(covplot, 12, 20, new.window)
              if (length(wMods)) {
                mods2$Y <- ceiling(mods2$Y)*(offset2) - match(mods2$Modification, uMods)
                uModsDf <- plyr::rbind.fill(lapply(unique(align.temp2$Y[wAA]), function(y) {
                  data.frame(Mod = uMods, Y = y-1-1:nMods, X = -1)
                }))
                covplot <- covplot +
                  ggplot2::geom_text(data = uModsDf,
                                     ggplot2::aes(label = Mod, x = X, y = Y),
                                     hjust = 1, vjust = 0.5, fontface = "italic",
                                     colour = "white", cex = size*1.1) +
                  # geom_rtile(data = mods2,
                  #            ggplot2::aes(x = X, y = Y-1, fill = Intensity, colour = Modification),
                  #            width = 0.9, height = 0.9,
                  #            radius = ggplot2::unit(3, "pt"))
                  ggplot2::geom_tile(data = mods2,
                                     ggplot2::aes(x = X, y = Y-1, fill = Intensity, colour = Modification),
                                     width = 0.9, height = 0.9#, cex = 0.5
                                     )
                #proteoCraft::poplot(covplot, 12, 20, new.window)
                Xright <- ceiling(max(c(align.temp2$X, mods2$X, uModsDf$X), na.rm = TRUE))+2
                Ybottom <- floor(min(c(align.temp2$Y, mods2$Y, uModsDf$Y, 0), na.rm = TRUE)) - 3
                Ytop <- ceiling(max(c(align.temp2$Y, mods2$Y, uModsDf$Y), na.rm = TRUE)) + 1
              } else {
                Xright <- ceiling(max(align.temp2$X, na.rm = TRUE))+2
                Ybottom <- floor(min(c(0, align.temp2$Y), na.rm = TRUE)) - 1
                Ytop <- ceiling(max(align.temp2$Y, na.rm = TRUE)) + 1
              }
              Ytop <- max(align.temp2$Y, na.rm = TRUE) + 1
              covplot <- covplot +
                ggplot2::scale_fill_gradient(low = colscale[1], high = colscale[2], na.value = na,
                                             breaks = brks, labels = brkstxt,
                                             limits = c(0, max(as.numeric(brks)))) +
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
                               plot.margin = ggplot2::margin(1, 4, 1, 1, "cm"))
              plotCreated <- TRUE
              #proteoCraft::poplot(covplot, 12, 20, new.window)
            }
          } else {
            warning(paste0("Not a single peptide matches the sequence(s)!"))
          }
        }
        if (plotCreated) {
          wdth <- 20
          xRg <- Xright-Xleft
          yRg <- Ytop-Ybottom
          hght <- max(c(wdth*yRg/(xRg*yxRat2), 10))
          wdth <- max(c(wdth, hght*0.5))
          if (display) {
            proteoCraft::poplot(covplot, height = hght, width = wdth, new.window = new.window)
          }
          if (save.grph) {
            for (svpth in unlist(tbl$Path)) {
              #ext <- rev(unlist(strsplit(svpth, "\\.")))[1]
              # if (ext == "pdf") {
              #   ggplot2::ggsave(svpth, covplot)
              # } else {
                ggplot2::ggsave(svpth, covplot, dpi = 300,
                                width = wdth,
                                height = hght,
                                units = "in",
                                limitsize = FALSE)
              #}
            }
          }
        }
        if (Mode == "XML") {
          tst <- grep("^Protein name: ", align.temp$AA)
          XML_Cov <- lapply(1:length(tst), function(x) {
            #x <- 1
            x <- align.temp[(tst[x]+2):(c(tst, nrow(align.temp)+1)[x+1]-1),
                            c("Amino_acid", "Legend")]
            x <- x[which(x$Amino_acid != ""),]
            u <- unique(x$Legend)
            w <- c(0, which(sapply(1:(nrow(x)-1), function(y) {
              x$Legend[y] != x$Legend[y+1] # These are the ends of each stretch
            })), length(x$Legend))
            rs <- as.data.frame(t(sapply(1:(length(w)-1), function(y) {
              rg <- (w[y]+1):w[y+1]
              c(paste0(x$Amino_acid[rg], collapse = ""), x$Legend[rg[1]])
            })))
            rs$V2 <- c(1, 2)[match(rs$V2, c("not found", "found"))]
            rs <- apply(rs[, c("V1", "V2")], 1, function(y) {
              #y <- rs[1, c("V1", "V2")]
              z <- as.numeric(y[[2]])
              list(openxlsx2::fmt_txt(y[[1]], bold = as.logical(z-1),
                                      color = openxlsx2::wb_color(hex = c("grey", colour)[as.numeric(z)])))
            })
            txt <- paste0("rs <- ", paste(sapply(1:length(rs), function(x) {
              paste0("rs[[", x, "]][[1]]")
            }), collapse = " + "))
            eval(parse(text = txt))
            return(rs)
          })
          if (length(XML_Cov) > 1) {
            rs2 <- openxlsx2::fmt_txt(";")
            txt <- paste0("XML_Cov <- ", paste(sapply(1:length(XML_Cov), function(x) {
              paste0("XML_Cov[[", x, "]][[1]]")
            }), collapse = " + rs2 + "))
            eval(parse(text = txt))
          }
          return(XML_Cov)
        }
      }
    } else { return(coverage) }
  } else {
    if (Mode == "Coverage") { cat("No peptides provided -> coverage = 0%\n") }
    if (Mode %in% c("Align", "Align2")) { warning("No peptides to plot!") }
    if (Mode == "XML") { return(fmt_txt(paste(proteins, collapse = ";"), color = wb_color(hex = c("grey")))) }
  }
}
