% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getters_setters.R
\name{fnames}
\alias{fnames}
\alias{fnames,SummarizedExperiment-method}
\alias{fnames,eSet-method}
\alias{fnames,EList-method}
\alias{fnames<-}
\alias{fnames<-,SummarizedExperiment,character-method}
\alias{fnames<-,eSet,character-method}
\alias{fnames<-,EList,character-method}
\title{Get/Set fnames}
\usage{
fnames(object)

\S4method{fnames}{SummarizedExperiment}(object)

\S4method{fnames}{eSet}(object)

\S4method{fnames}{EList}(object)

fnames(object) <- value

\S4method{fnames}{SummarizedExperiment,character}(object) <- value

\S4method{fnames}{eSet,character}(object) <- value

\S4method{fnames}{EList,character}(object) <- value
}
\arguments{
\item{object}{SummarizedExperiment, eSet, or EList}

\item{value}{character vector with feature names}
}
\value{
feature name vector (get) or updated object (set)
}
\description{
Get/Set feature names
}
\examples{
if (require(autonomics.data)){
   require(magrittr)
   object <- autonomics.data::stemcomp.proteinratios
   fnames(object) \%<>\% paste0('PG', .)
   object
}
}
