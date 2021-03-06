% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add_limma.R
\name{add_limma_to_fdata}
\alias{add_limma_to_fdata}
\title{Run "limma" analysis and add to SummarizedExperiment}
\usage{
add_limma_to_fdata(
  object,
  contrasts = autonomics.find::default_contrasts(object),
  design = autonomics.import::create_design_matrix(object)
)
}
\arguments{
\item{object}{SummarizedExperiment}

\item{contrasts}{a character vector of contrasts (of subgroup levels)}

\item{design}{design matrix}

\item{overwrite}{whether to overwrite existing results in fdata (logical)}
}
\value{
updated SummarizedExperiment with limma results in fdata
        added columns: rank.xxx value.xxx p.xxx fdr.xxx bonf.xxx \cr
}
\description{
Performs a limma analysis for the specified contrasts.
}
\note{
Features with a single observation can get a significant p value
(though not a significant q value) when the exprs value for that feature is exceptionally large.
This is a consequence of the shrinkage estimation of the limma method, as
\href{https://support.bioconductor.org/p/4932/}{explained by Gordon Smyth on Bioconductor Support}.
}
\examples{
library(magrittr)

# STEM CELL COMPARISON
if (require(autonomics.data)){
   object <- autonomics.data::stemcomp.proteinratios
   object \%>\% autonomics.find::add_limma() \%>\% 
              autonomics.import::limma()   \%>\% 
              names()
}

if (require(subramanian.2016)){
   object <- subramanian.2016::metabolon
   object \%>\% autonomics.find::add_limma() \%>\% autonomics.import::fvars()
   contrasts <- subramanian.2016::contrasts.metabolon
   object \%<>\% autonomics.find::add_limma_to_fdata(contrasts)
}
if (require(billing.differentiation.data)){
   contrasts <- billing.differentiation.data::contrasts
   billing.differentiation.data::rna.counts \%>\% autonomics.import::logcpm() \%>\%
      autonomics.find::add_limma_to_fdata(contrasts) \%>\%  autonomics.import::fdata() \%>\%
      extract2('bonf.EM0.8_0') \%>\% is_less_than(0.05) \%>\% sum(na.rm=TRUE)
   billing.differentiation.data::rna.voomcounts \%>\%
      autonomics.find::add_limma_to_fdata(contrasts) \%>\%  autonomics.import::fdata() \%>\%
      extract2('bonf.EM0.8_0') \%>\% is_less_than(0.05) \%>\% sum(na.rm=TRUE)
}
}
\author{
Aditya Bhagwat
}
