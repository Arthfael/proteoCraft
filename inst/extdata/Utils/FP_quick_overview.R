require(proteoCraft)
wd <- "D:/groups_temp/QCs/FP_20260310"

RPath %<o% as.data.frame(library()$results)
RPath <- normalizePath(RPath$LibPath[match("proteoCraft", RPath$Package)], winslash = "/")
libPath %<o% paste0(RPath, "/proteoCraft")
if (!exists("N.clust")) { N.clust <- max(c(round(parallel::detectCores()*0.95)-1, 1)) }
parSrc %<o% paste0(libPath, "/extdata/Sources/make_check_Cluster.R")
source(parSrc)


rDir <- "H:/aRmel_package/proteoCraft/R"
source(paste0(rDir, "/bind_worker.R"))
source(paste0(rDir, "/FP2MQ_modSeqWrkr1.R"))
source(paste0(rDir, "/FP2MQ_modSeqWrkr2.R"))
source(paste0(rDir, "/FP2MQ_modSeqWrkr3.R"))
source(paste0(rDir, "/FP2MQ_modSeqWrkr4.R"))
source(paste0(rDir, "/FP2MQ_modSeqWrkr5.R"))
source(paste0(rDir, "/FP_to_MQ.R"))
fls <- list.files(wd, full.names = TRUE)
wrkflw <- grep("\\.workflow$", fls, value = TRUE)
mnfst <- grep("\\.fp-manifest$", fls, value = TRUE)
tmpFP <- FP_to_MQ(wrkflw, mnfst, cl = parClust)

tmp <- tmpFP$Evidence
#
require(ggplot2)
tmp <- aggregate(tmp$id, list(tmp$Experiment), length)
colnames(tmp) <- c("Run", "PSMs")
plot <- ggplot(tmp) +
  geom_bar(aes(Run, PSMs, fill = Run), stat = "identity") +
  theme_bw()
poplot(plot)
