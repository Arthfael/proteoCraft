# proteoCraft

<b><i>A set of scripts and functions for R-based MS-Proteomics analysis.</i></b>

The following types of experiments are covered:
 1) Replicates available (2 workflows: protein groups or modified peptides only):
	- "Regulation analysis", aka DEP (Differential Expression of Proteins)
	- Pull down analysis (co-IPs, BioID, TurboID, APEX...), in which case also saintEXPRESS analysis is performed
	- Regulation analysis of post-translationally modified peptides
	- Analysis of changes to subcellular localisation
	- Time course analysis
 2) No replicates available (1 workflow):
	- Samples characterisation
	- Gel band ID
	- Basic, no-statistics involved comparisons between samples
 3) Histones modified peptides analysis.

<i>(NB: All workflows could previously be run through a single master script, which ran code chunks sequentially and allowed partial processing. It is still included, but hasn't been tested/maintained in a while, so we would recommend not to use it for now.)</i>


### Input data
Though where possible the scripts have been written as generally as possible, they have only been tested on data from:
 - ThermoScientific Q-Exactive HF (a few datasets on Exploris with FAIMS or Fusion Tribrid)
 - Bruker timsTOF HT
Input data can be the output from:
 - MaxQuant
 - diaNN
 - FragPipe (using "reasonable" search parameters... FragPipe is very flexible)
There is burgeoning support for:
 - Proteome Discoverer: A conversion function exists but integration with the analysis scripts is currently deprecated.
 - Skyline: A basic conversion function exists and is usable with the Histones script. It probably needs fuller development and integration with the scripts but for now allows to deal with .tsv/.csv exports assuming a set of core columns were exported.
 - alphaDIA: A basic conversion function exists but is not integrated with the analysis scripts - for this parsing of the alphaDIA log would need to be introduced.

Depending on search parameters you may encounter unexpected situations, as we have not tested every possible situation but only ones occurring in our facility. It is likely that some options in some software will break the conversion functions, especially for the more versatile engines such as FragPipe. If that happens, get in touch and we will add support asap.
We would like to add Spectronaut/Peaks/alphaDIA/SAGE support too eventually...


### Environment and requisites
This is meant to run through RStudio in Windows and will probably not work without significant edits on Linux/Unix... though here and there some work has been done towards achieving that goal.
This limitation is because the scripts will sometimes:
 - use R code which is Windows-specific, e.g. makes assumptions about the type of (G)UI available.
 - build and run windows command line or powershell-specific commands.
Run this ideally on a multi-core PC optimized for data analysis with a lot of memory. Many steps are parallelized. We run this on a Windows Server PC with 56 vCPUs and 128 GB RAM.

In addition to R and RStudio, for a smoother experience you should also install:
 - Python: some python scripts are run here and there.
 - ScanHeadsMan (run by MatMet_LCMS for extracting methods from Thermo raw files)
 - Cytoscape (for ClueGO GO-terms enrichment analysis and STRINGdb networks visualisation)
 - Word or any software which can open .docx files (to check the Materials and Methods template text)
 - Excel or any software which can open .xlsx files (to check the final Report file)
 - Optional: a web browser for opening interactive plots (.html files)
Normally, the scripts provided should automatically install any R packages required.
In case this fails just manually install the missing ones.


### Troubleshooting:
 - In our hands RStudio can become very slow after some time. To solve the issue, go to Tools > Global Options..., click on Code then Diagnostics and untick "Show diagnostics for R"
 - We create a parallel cluster then reuse it. It is known that after some time the cluster gets corrupted. We tend to delete/recreate it regularly but this has a cost in time. If at some point the outcome of a parallel function has a length not matching its input, your should run the following lines then restart where the workflow stopped:
```
> stopCluster(parClust)
> source(parSrc)
# (also rerun any parallel::clusterCall(...) or parallel::clusterExport(...) calls before where this failed)
```

#### <i>Note on speed:</i>
<i>This package has been somewhat optimised for speed, but many steps could still probably be faster. The current main bottlenecks include:<i/>
 - <i>Extracting MS methods using ScanHeadsMan,<i/>
 - <i>Checking peptide-to-protein assignments provided by the search software (recommended),<i/>
 - <i>Extracting TIC/Base peak information from Thermo raw files using the rawrr package,<i/>
 - <i>Cytoscape-dependent steps.<i/>
