% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Prot.Quant.R
\name{Prot.Quant}
\alias{Prot.Quant}
\title{Prot.Quant}
\usage{
Prot.Quant(
  Prot,
  Peptide.IDs = "Peptide IDs",
  Mode = "Classic",
  Pep,
  id = "id",
  Summary.method = "median",
  Summary.weights,
  Intensity.weights = FALSE,
  Priority = "Intensities",
  Skip.ratios = FALSE,
  experiments.map,
  ref.groups,
  ratio.groups,
  sample.groups,
  param,
  aggregate.map = Aggregate.map,
  aggregate.list = Aggregate.list,
  aggregates = Aggregates,
  Pep.Intens.root,
  Pep.Ratios.root,
  log.Pep.Intens = FALSE,
  log.Pep.Ratios = 2,
  Prot.LFQ.to.log = TRUE,
  Prot.Ratios.to.log = TRUE,
  Mods,
  Mods.to.Exclude,
  Mod.Nms = "Modified sequence",
  Discard.unmod = TRUE,
  Prot.Seq = "Sequence (1st accession)",
  Min.N = 2,
  Max.N = 50,
  Ratios.SD = FALSE,
  Ratios.Pvalue = FALSE,
  N.clust,
  N.reserved = 1,
  Unique.peptide.IDs = "Unique peptide IDs",
  Min.Unique = 3,
  Refs_Mode = "2",
  Refs_AllGroups = NULL,
  cl
)
}
\arguments{
\item{Prot}{Protein/Protein groups table. A data.frame.}

\item{Peptide.IDs}{Name of the Protein/Protein groups table's peptide IDs column. Default = "Peptide IDs"}

\item{Mode}{In "Classic" mode (current default), the function just uses the "Peptide.IDs" argument. In "PreferUnique" mode, the "Peptide.IDs" column should contain razor or shared peptide IDs, and a second column name of unique peptides should be provided using the "Unique.peptide.IDs" argument. If at least "Min.Unique" unique peptides are present, then only those will be used for quantitation, otherwise as many razor/shared peptides as necessary will be added (sorted by decreasing average intensities).}

\item{Pep}{Peptides table. A data.frame.}

\item{id}{The name of the Peptides table's IDs column. Default = "id"}

\item{Summary.method}{The summary method used for ratios (the Levenberg-Marquardt algorithm is used for intensities). One of "mean", "median", or "weighted.mean".}

\item{Summary.weights}{If a "weighted.mean" summary method is chosen, then a vector of weights must be provided. This should be the name of a column of the peptides table}

\item{Intensity.weights}{default = FALSE; if TRUE, will take into account individual peptide intensities when calculating average profile if "Summary.method" is either mean "mean" or "weighted.mean" (thus, it will actually be a weighted mean regardless).}

\item{Priority}{One of "Ratios" or "Intensities" (default). Some flexibility in spelling is allowed. You want to prioritize ratios for SILAC because in this case MaxQuant measures peptides ratios more precisely than intensities. Otherwise, you want to prioritize Intensities and re-calculate ratios from them.}

\item{Skip.ratios}{Default = FALSE. If TRUE, ratios will not be calculated.}

\item{experiments.map}{The experiments map.}

\item{ref.groups}{Defines which samples are paired to which references. May alternatively (preferred solution) be provided indirectly through the param argument.}

\item{ratio.groups}{Defines groups within which ratios are calculated. May alternatively (preferred solution) be provided indirectly through the param argument.}

\item{sample.groups}{Defines groups of samples which are replicates of the same category. Used to calculate within group reference ratios.}

\item{param}{The experiment's parameters object. If provided, the ref.groups argument is not required.}

\item{aggregate.map}{The aggregate map. Default = Aggregate.map}

\item{aggregate.list}{The named list of aggregates. Default = Aggregate.list}

\item{aggregates}{The aggregates themselves. Default = Aggregates}

\item{Pep.Intens.root}{Root of the peptides intensity column(s) names}

\item{Pep.Ratios.root}{Root of the peptides ratios column(s) names}

\item{log.Pep.Intens}{Set to 0 or FALSE if input peptide intensities are linear scale. If the data is already log scale, set to the relevant scale's base. Default = FALSE}

\item{log.Pep.Ratios}{Set to 0 or FALSE if input peptide ratios are linear scale. If the data is already log scale, set to the relevant scale's base. Default = 2}

\item{Prot.LFQ.to.log}{Should the output protein LFQ values be log-scale or not? Can ve set to the log base desired. If set to TRUE, this will be the same base as the input intensities log scale, or 10 by default.}

\item{Prot.Ratios.to.log}{Should the output protein ratios be log-scale or not? Can ve set to the log base desired. If set to TRUE, this will be the same base as the input ratios log scale, or its default, 2}

\item{Mods}{Which modifications (2 lowercase letters PTM code) should be included? If set to FALSE, will not filter any modifications.}

\item{Mods.to.Exclude}{Alternative/complementary to Mods. Which modifications should be excluded? (use argument "Discard.unmod" to discard unmodified counterpart peptides.) A data.frame with columns "Mark" (code) and "Where", see argument "Discard.unmod".}

\item{Mod.Nms}{Default = "Modified sequence". The name of the column containing the modified sequence in the peptides table. Can be set to a non-modified sequence if "Mods" = FALSE and "Mods.to.Exclude" is empty.}

\item{Discard.unmod}{Default = TRUE. Should we discard those unmodified peptides whose primary sequence is the same as that of some modified peptides we will not use? If set to 2, will use the "Where" column in "Mods.to.Exclude" to identify (and exclude) peptides which could be modified even if the modified form was not identified. Requires knowledge of protein sequence ("Prot.Seq" argument)! Be careful! This will likely massively reduce the number of peptides available for quantitation!}

\item{Prot.Seq}{Default = "Sequence (1st accession)", used if "Discard.unmod" is set to 2}

\item{Min.N}{How many peptides should at least be present? Should be at the very least 1 (default = 2).}

\item{Max.N}{How many peptides can we use at most for the Levenberg-Marquardt procedure (used only for LFQs)? Using too many peptides can be an issue, e.g. with huge proteins like Titin. Default = 50. The most intense peptides will be selected.}

\item{Ratios.SD}{Deprecated (Should the output ratios also include peptide level standard deviations? Default = FALSE)}

\item{Ratios.Pvalue}{Deprecated (Should the Output also include peptide level P-values? Default = FALSE)}

\item{N.clust}{A limit on the number of vCPUs to use. If left as NULL (default), uses the number of available clusters - 1, to a minimum of 1.}

\item{N.reserved}{Default = 1. Number of reserved vCPUs the function is not to use. Note that for obvious reasons it will always use at least one.}

\item{Unique.peptide.IDs}{Used only if mode = "PreferUnique". Name of the Protein/Protein groups table's unique peptide IDs column. Default = "Unique peptide IDs"}

\item{Min.Unique}{Used only if mode = "PreferUnique". Minimum number of unique peptides to consider before adding shared peptides. Default = 3}

\item{Refs_Mode}{How are reference ratios calculated?\cr
- If set to "1", only references are considered (i.e. it compares individual references either to each other, or if available to the average reference for the group).\cr
- If set to "2" (default), for each ratios group, reference ratios are based on comparing every possible pair of samples within the group.\cr}

\item{Refs_AllGroups}{Deprecated, use Refs_Mode instead.}

\item{cl}{Already have a cluster handy? Why spend time making a new one, which on top of that may invalidate the old one. Just pass it along!}
}
\description{
A function to calculate estimated protein group Expression (Absolute Quant) and Ratios (optional) from individual peptide values.
The input is assumed to be normalized.

A note on ref.groups and ratio.groups:
These parameters are the method I have used to define if replicates are paired or not.
For an experiment with one master Experiment, several Conditions (incl. one reference) and Replicates:
- in a paired (= "nested") setup, you would set ref.groups to "ExpRep" and ratio.groups to "Exp"
- in an unpaired setup, you would set ref.groups to "Exp" to and ratio.groups to "Exp"
Essentially, ref.groups are the groups within which ratios are calculated to all available references,
while ratio.groups are groups within which reference-to-reference values are calculated:
- between individual references in the paired setup
- from individual references to the average reference in the unpaired setup
}
\examples{
temp <- Prot.Quant(Prot = PG, Peptide.IDs = Pep4Quant, Pep = pep, id = "New Peptide ID",
                   experiments.map = Exp.map,
                   ref.groups = Ratios.Ref.Groups,
                   Pep.Intens.root = pep.ref, Pep.Ratios.root = pep.ratios.ref,
                   log.Pep.Intens = FALSE, log.Pep.Ratios = 2,
                   Prot.LFQ.to.log = TRUE, Prot.Ratios.to.log = TRUE,
                   Mods = Mod4Quant, Mod.Nms = "Modified sequence",
                   Min.N = 2, Ratios.SD = FALSE, Ratios.Pvalue = FALSE)

}
