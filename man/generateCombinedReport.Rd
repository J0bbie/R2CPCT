% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysisWGS_generateCombinedReport.R
\name{generateCombinedReport}
\alias{generateCombinedReport}
\title{Generate a gene-level overview of all overlapping mutations (incl. SV and CNA) and dN/dS + GISTIC2 results.}
\usage{
generateCombinedReport(
  data.Cohort,
  dNdS,
  GISTIC2,
  nThreads = 40,
  mutantsOnly = TRUE
)
}
\arguments{
\item{data.Cohort}{(list): Cohort-wide data obtained from \link[R2CPCT]{importWGSOfCohort}.}

\item{dNdS}{(list): Output of the dN/dS analysis obtained from \link[R2CPCT]{rundNdS}.}

\item{GISTIC2}{(list): Output of the GISTIC2 analysis obtained from \link[R2CPCT]{importGISTIC2}.}

\item{nThreads}{(integer): Number of cores over which to parallelize (when needed).}

\item{mutantsOnly}{(logical): Only output records with mutations (CNA, LOH or Muts; TRUE) or return all records (FALSE).}
}
\value{
(tibble) Returns a per-sample gene-level overview of all reported somatic aberrations on protein-coding genes.
}
\description{
The gene-level overview will be performed on a per-sample basis using only protein-coding genes.
The combined report will only contain protein-coding genes which have at least a single somatic aberration or were detected in the dN/dS or GISTIC2 analysis.

Use the 'selectedGenes' parameter to retrieve your genes-of-interests, even if no aberrations were detected.
}
\examples{
\donttest{

	# Generate a combined gene-level report of an imported WGS cohort (which could be a subset).
	generateCombinedReport(data.Cohort, results.Cohort$dNdS, results.Cohort$GISTIC2)

}

}
\author{
Job van Riet \email{j.vanriet@erasmusmc.nl}
}
