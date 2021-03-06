% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_contrast_matrix.R
\name{create_contrast_matrix}
\alias{create_contrast_matrix}
\title{Create contrast matrix}
\usage{
create_contrast_matrix(design, contrast_defs = colnames(design))
}
\arguments{
\item{design}{design matrix}

\item{contrast_defs}{vector of contrast defs}
}
\value{
contrast matrix
}
\description{
Create contrast matrix
}
\examples{

# STEM CELL COMPARISON
if (require(autonomics.data)){
   require(magrittr)
   autonomics.data::stemcomp.proteinratios   \%>\% 
   autonomics.import::create_design_matrix()   \%>\% 
   autonomics.find::create_contrast_matrix()
}

if (require(billing.differentiation.data)){
   billing.differentiation.data::rna.voomcounts \%>\% 
      autonomics.import::create_design_matrix() \%>\% 
      autonomics.find::create_contrast_matrix(
         c(EM0.8_0 = 'EM0.8 - EM0.0'))
}
if (require(autonomics.data)){
    object <- autonomics.data::glutaminase
    autonomics.import::create_design_matrix(object)
    autonomics.find::create_contrast_matrix(object)
}
}
