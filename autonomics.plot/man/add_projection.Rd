% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pca.R
\name{add_projection}
\alias{add_projection}
\alias{add_pca}
\alias{add_sma}
\alias{add_lda}
\alias{add_pls}
\title{Add PCA/LDA/PLS/SMA results to object}
\usage{
add_projection(
  object,
  method = c("pca", "lda", "pls")[1],
  implementation = NULL,
  na.impute = FALSE,
  ndim = 2
)

add_pca(object, ...)

add_sma(object, ...)

add_lda(object, ...)

add_pls(object, ...)
}
\arguments{
\item{object}{SummarizedExperiment}

\item{method}{'pca' or 'lda'}

\item{implementation}{string}

\item{na.impute}{TRUE or FALSE}

\item{ndim}{number of dimensions to include}

\item{...}{passed to add_projection}
}
\description{
Add PCA/LDA/PLS/SMA results to object
}
\examples{
if (require(autonomics.data)){

   # STEM CELL COMPARISON
   require(magrittr)
   object <- autonomics.data::stemcomp.proteinratios
   object \%>\%  add_pca(ndim = 3)
   object \%<>\% add_pca()
   object \%>\%  autonomics.import::svars()
   object \%>\%  autonomics.import::fvars()
}

}
