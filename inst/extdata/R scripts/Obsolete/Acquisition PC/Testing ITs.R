require(proteoCraft)

wd <- "...QC_File_Archive/20211130 - test IT/combined/txt"
wd <- "...My_Dataset/combined/txt_6Fr"
wd <- "...My_Dataset/combined/txt_12Fr"

setwd(wd)
MQ.load()

Mods <- cor_mod_seq(ev)
ev <- Mods$Peptides
Mods <- Mods$PTMs
aggregate(ev$Sequence, list(ev$`Raw file`), length)

tst <- MQ.summary(wd, ev, PG, mods = setNames(Mods$Mark, Mods$`Full name`))

#openwd()
