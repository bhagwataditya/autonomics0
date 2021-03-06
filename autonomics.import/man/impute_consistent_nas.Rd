% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/impute_consistent_nas.R
\name{impute_common_nas}
\alias{impute_common_nas}
\alias{impute_consistent_nas}
\title{Impute consistent NA values}
\usage{
impute_common_nas(
  object,
  imputefun = function(x) impute_around_zero(x),
  verbose = FALSE
)

impute_consistent_nas(
  object,
  imputefun = function(x) impute_around_zero(x),
  svar = "subgroup",
  verbose = FALSE
)
}
\arguments{
\item{object}{SummarizedExperiment}

\item{imputefun}{imputation function}

\item{verbose}{TRUE or FALSE}

\item{svar}{string}
}
\value{
SummarizedExperiment with updated exprs
}
\description{
Impute values missing in all (subgroup) samples.
}
\examples{
require(magrittr)
if (require(autonomics.data)){

   # Read object
   object <- 'extdata/glutaminase/glutaminase.xlsx'     \%>\%
              system.file(package = 'autonomics.data')  \%>\%
              read_metabolon()

   # Common NA values - missing in all samples
   object \%>\% split_by_svar() \%>\%
              magrittr::extract2(2) \%>\%
              impute_common_nas(verbose = TRUE)

   # Consistent NA values - missing in all subgroup samples
   object \%>\% impute_consistent_nas(verbose = TRUE)
   # Use a different imputation method
   object \%>\% impute_consistent_nas(
      verbose = TRUE,
      imputefun = function(x) imputeLCMD::impute.QRILC(x)[[1]])
}
}
