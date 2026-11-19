if (!exists("bckpsDir")) { bckpsDir %<o% paste0(wd, "/_backups") }
if (!dir.exists(bckpsDir)) { dir.create(bckpsDir, recursive = TRUE) }
backupFls %<o% setNames(paste0(bckpsDir, "/",
                               c("Backup", "QC", "evmatch", "AnalysisParam", "PG_assembly", "Contrasts", "pep_intens_norm",
                                 "PG_quant", "PG_quant_reNorm", "PG_quant_imput", "covPlots", "ratioPlots", "HeatMaps",
                                 "HeatMaps_data", "profilePlots", "quantPlots", "Venn_plotly", "xlDat", "volcPlotly",
                                 "DimRedPlots", "GSEA_plotly", "GO_plot_ly", "GO_mappings", "GO_terms", "Parsed_annotations",
                                 "CustPG", "limmaFits", "MSqRob", "ROTS"),
                               ".RDS"),
                        c("Main", "QC", "evmatch", "Param", "PG_assembly", "Contrasts", "pep_intens_norm",
                          "PG_quant", "PG_quant_reNorm", "PG_quant_imput", "Coverage", "Ratios", "HeatMaps",
                          "HeatMaps_data", "Profile", "Ranked", "Venn", "Excel report data", "Vulcan",
                          "dimRed", "GSEA", "GO plots", "GO mappings", "GO terms", "Annotations",
                          "custom PGs", "limma", "MSqRob", "ROTS"))
