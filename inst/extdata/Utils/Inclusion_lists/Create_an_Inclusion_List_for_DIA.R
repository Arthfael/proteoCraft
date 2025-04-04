# Create an Inclusion List for DIA
# Header
options(stringsAsFactors = FALSE)
library(utils)
library(ggplot2)
library(reshape2)
library(proteoCraft)
library(svDialogs)

# Choose destination directory:
destd <- choose.dir("...Search_Folder/")

N.windows <- as.numeric(dlgInput("How many DIA windows do you want?", 40)$res)
Lower <- as.numeric(dlgInput("Enter the lower M/Z boundary:", 400)$res)
Upper <- as.numeric(dlgInput("Enter the upper M/Z boundary:", 1200)$res)
Overlap <- as.numeric(dlgInput("Enter the overlap between adjacent windows:", 1)$res)
#NB: I thought that OpenSWATH does not support overlap, but it does!!!
## Simply, while the acquisition windows list can have overlaps, the analysis windows list should not!
RT_start <- as.numeric(dlgInput("Enter the time (minutes) the method starts:", 0)$res)
RT_end <- as.numeric(dlgInput("Enter the time (minutes) the method ends:", 220)$res)


Width <- (Upper-Lower+(N.windows-1)*Overlap)/N.windows
Width <- round(Width, 0) # No need to be too precise (and more than 1 decimals is not accepted anyway)
SWATH.windows <- data.frame(Left_boundary = Lower+((1:N.windows)-1)*(Width-Overlap))
SWATH.windows$Right_boundary <- SWATH.windows$Left_boundary+Width
SWATH.windows$"Mass [m/z]" <- (SWATH.windows$Left_boundary+SWATH.windows$Right_boundary)/2
SWATH.windows$Width <- Width
SWATH.windows$"Start [min]" <- RT_start
SWATH.windows$"End [min]" <- RT_end
SWATH.windows2 <- SWATH.windows[,c("Mass [m/z]", "Start [min]", "End [min]")]
SWATH.windows2$"Polarity" <- "Positive"
for (i in c("Formula [M]", "Formula type", "Species", "CS [z]", "(N)CE", "(N)CE type", "MSX ID", "Comment")) {
  SWATH.windows2[[i]] <- ""
}
SWATH.windows2 <- SWATH.windows2[,c("Mass [m/z]", "Formula [M]", "Formula type", "Species", "CS [z]", "Polarity",
                                    "Start [min]", "End [min]", "(N)CE", "(N)CE type", "MSX ID", "Comment")]
write.csv(SWATH.windows2, file = paste0(destd, "\\DIA_Inclusion_list_-_centre.txt"), row.names = FALSE, quote = FALSE)
# Also write .tsv files for OpenSwath analysis
## Acquisition tsv
SWATH.windows3 <- SWATH.windows[,c("Left_boundary", "Right_boundary")]
#colnames(SWATH.windows3) <- c("start", "end")
colnames(SWATH.windows3) <- c("lower_offset", "upper_offset")
write.table(SWATH.windows3, file = paste0(destd, "\\DIA_Inclusion_list_-_acquisition.tsv"), row.names = FALSE,
            #col.names = FALSE,
            quote = FALSE, sep = "\t")
## Analysis tsv
SWATH.windows4 <- SWATH.windows3
for (i in 1:(nrow(SWATH.windows4)-1)) {
  m <- mean(c(SWATH.windows4[i, 2],SWATH.windows4[i+1, 1]))
  SWATH.windows4[i, 2] <- SWATH.windows4[i+1, 1] <- m
}
write.table(SWATH.windows4, file = paste0(destd, "\\DIA_Inclusion_list_-_analysis.tsv"), row.names = FALSE,
            #row.names = FALSE,
            quote = FALSE, sep = "\t")

shell.exec(destd)

