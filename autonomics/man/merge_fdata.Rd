% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/merge.R
\name{merge_fdata}
\alias{merge_fdata}
\title{Merge fdata into SummarizedExperiment}
\usage{
merge_fdata(object, newdata, by = "feature_id", newdata_first = FALSE)
}
\arguments{
\item{object}{SummarizedExperiment}

\item{newdata}{dataframe}

\item{by}{variable by which to merge}

\item{newdata_first}{logical: whether to have newdata first in sdata}
}
\description{
Merge fdata into SummarizedExperiment
}
\details{
Common variables (other than by) are removed from newdata prior to merge.
}
\examples{
require(magrittr)
if (require(autonomics.data)){
   object <- system.file('extdata/stemcomp/soma/stemcomp.adat',
                          package = 'autonomics.data') \%>\%
             autonomics::read_somascan()
   newdata <- autonomics.import::fdata(object)           \%>\%
              magrittr::extract(, 'feature_id', drop = FALSE) \%>\%
              cbind(Letter = sample(letters, nrow(.), TRUE))
   object \%>\% autonomics.import::fdata() \%>\% head()
   object \%>\% autonomics::merge_fdata(newdata) \%>\%
              autonomics.import::fdata() \%>\% head()
   object \%>\% autonomics.import::merge_fdata(newdata, newdata_first = TRUE) \%>\%
              autonomics.import::fdata() \%>\% head()
}
}
