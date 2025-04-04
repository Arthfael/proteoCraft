#####################################################
#                                                   #
#                    proteoCraft                    #
#                                                   #
#####################################################

This package is meant for analysis of proteomics data, with the following types of experiments in mind:
 - Replicates available (2 workflows: protein groups or modified peptides only):
	- "Regulation analysis", aka DEP (Differential Expression of Proteins)
	- Pull down analysis (co-IPs, BioID, TurboID, APEX...), in which case also saintEXPRESS analysis is performed
	- Regulation analysis of post-translationally modified peptides
	- Analysis of changes to subcellular localisation
	- Time course analysis
 - No replicates available (1 workflow):
	- Samples characterisation
	- Gel band ID
	- Basic, no-statistics involved comparisons between samples

All workflows could previously be run through a single master script, which ran code chunks sequentially and allowed partial processing.
I have had some issues with it and it hasn't been tested/maintained in a while, so I would recommend not to use it for now.

Though where possible the scripts have been written as generally as possible, they have only been tested on data from:
 - ThermoScientific Q-Exactive HF (a few datasets on Exploris with FAIMS or Fusion Tribrid)
 - Bruker timsTOF HT
Input data can be the output from MaxQuant, diaNN, FragPipe (using "reasonable" search parameters... FragPipe is very flexible) or Proteome Discoverer (largelly deprecated but could be easily "restored") searches.
We would like to add Spectronaut/Peaks/AlphaPept support too.
Depending on search parameters you may encounter unexpected situations, as we have not tested every possible situation but only ones occurring in our facility.

This is meant to run through RStudio in Windows and will probably not work without significant edits on Linux/Unix.
This is because the scripts will sometimes:
 - use R code which is Windows-specific, e.g. makes assumptions about the type of UI/GUI available.
 - build and run windows command line or powershell-specific commands.

The package has been somewhat optimised for speed, but probably many steps could be performed much faster.
Current bottlenecks include:
 - Extracting MS methods using ScanHeadsMan,
 - Checking peptide-to-protein assignments provided by the search software (recommended),
 - Extracting TIC/Base peak information using the rawrr package,
 - Cytoscape-dependent steps.

In addition to R and RStudio, for a flawless experience you should also install:
 - ScanHeadsMan (run by MatMet_LCMS, extracts methods from Thermo raw files)
 - Cytoscape (for ClueGO GO-terms enrichment analysis and network visualisation)
 - Word (for checking the Materials and Methods text)
 - Excel (to visualize reports)
 - Optional: a web browser for opening interactive plots (.html files)
Normally, the scripts provided should automatically install any R packages required.
In case this fails just manually install the missing ones.


Troubleshooting:
################
 - In our hands RStudio can become very slow after some time.
   To solve the issue, go to Tools > Global Options..., click on Code then Diagnostics and untick "Show diagnostics for R"
 - Up to v_6.4.0.0, some shiny apps in these workflows (specifically filling down dropdown menu selections) did not work if shiny_1.8.0 or later was used. This should be fixed now, and with the way the code was rewritten should also be faster and still work with older version. Still, kt keep an eye on it...

