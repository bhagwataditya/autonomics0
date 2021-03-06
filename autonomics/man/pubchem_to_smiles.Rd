% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotate.R
\name{pubchem_to_smiles}
\alias{pubchem_to_smiles}
\alias{pubchem_to_smiles.character}
\alias{pubchem_to_smiles.factor}
\alias{pubchem_to_smiles.SummarizedExperiment}
\title{Map a vector of pubchemids to (canonical) smiles}
\usage{
pubchem_to_smiles(x, ...)

\method{pubchem_to_smiles}{character}(x, ...)

\method{pubchem_to_smiles}{factor}(x, ...)

\method{pubchem_to_smiles}{SummarizedExperiment}(x, pubchem_var = "PUBCHEM", smiles_var = "SMILES", ...)
}
\arguments{
\item{x}{charactervector, factorvector, or SummarizedExperiment}

\item{...}{provided to enable S3 dispatch}

\item{pubchem_var}{pubchem fvar}

\item{smiles_var}{smiles fvar}
}
\value{
charactervector or SummarizedExperiment
}
\description{
Map a vector of pubchemids to (canonical) smiles
}
\examples{
require(magrittr)
# charactervector
   x <- c(NA_character_, "11988421", "10236635", "5283147", "91477", NA_character_)
   x \%>\% autonomics::pubchem_to_smiles()
   
# SummarizedExperiment
if (require(autonomics.data)){
   x <- autonomics.data::glutaminase
   x \%>\% autonomics::pubchem_to_smiles()
}
}
\references{
Backend is \href{https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest-tutorial}{PUBCHEM's REST interface}
}
