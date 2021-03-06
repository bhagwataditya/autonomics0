% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_feature_sets.R
\name{extract_ora_query}
\alias{extract_ora_query}
\title{Extract ora query}
\usage{
extract_ora_query(object, contrast_name, topdef, direction)
}
\arguments{
\item{object}{SummarizedExperiment}

\item{contrast_name}{character}

\item{topdef}{character}

\item{direction}{'pos' or 'neg'}
}
\value{
chracter vector
}
\description{
Extract ora query
}
\examples{
require(magrittr)
if (require(autonomics.data)){
   object <- autonomics.data::billing2016
   object \%>\% autonomics.ora::extract_ora_query('BM_E', 'fdr<0.05', 'pos')          \%>\% head(2)
}
}
