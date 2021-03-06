% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_design_matrix.R
\name{create_design_matrix}
\alias{create_design_matrix}
\title{Create design matrix for statistical analysis}
\usage{
create_design_matrix(
  object,
  intercept = length(unique(sdata(object)$subgroup)) == 1,
  confounders = character(0)
)
}
\arguments{
\item{object}{SummarizedExperiment}

\item{intercept}{TRUE or FALSE: include an intercept in the design?}

\item{confounders}{confounder svars (character)}
}
\value{
design matrix
}
\description{
Create design matrix for statistical analysis
}
\examples{
if (require(autonomics.data)){
   require(magrittr)

   # STEM CELL COMPARISON
   object <- autonomics.data::stemcomp.proteinratios
   object \%>\% create_design_matrix()

   # GLUTAMINASE
   object <- autonomics.data::glutaminase
   object \%>\% create_design_matrix() \%>\% extract(1:10, 1:10)
}
}
