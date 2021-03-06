% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter_contrast_features.R
\name{arrange_features_by_rank}
\alias{arrange_features_by_rank}
\title{Order features on rank}
\usage{
arrange_features_by_rank(object, contrast_name)
}
\arguments{
\item{object}{SummarizedExperiment}

\item{contrast_name}{contrast name}
}
\value{
arranged eSet
}
\description{
Order features on rank
}
\examples{
require(magrittr)
if (require(autonomics.data)){
   (object <- autonomics.data::stemdiff.proteinratios)
   contrast_name <-  names(autonomics.import::contrastdefs(object))[1]
   (object \%<>\% autonomics.find::arrange_features_by_rank(contrast_name))
   object \%>\% autonomics.import::limma() \%>\% magrittr::extract(1:2, contrast_name, 1:2)
}
}
