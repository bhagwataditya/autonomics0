% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/targetscan.R
\name{mir_to_gsymbol}
\alias{mir_to_gsymbol}
\title{Map vector of microRNAs to target gsymbols}
\usage{
mir_to_gsymbol(mir)
}
\arguments{
\item{mir}{vector of microRNAs}
}
\value{
named list which maps each mir to its target gene symbols
}
\description{
Map vector of microRNAs to target gsymbols
}
\examples{
require(magrittr)
mir <- c('hsa-mir-7-5p', 'hsa-mir-29a-3p', 'hsa-mir-199a-5p ')
mir \%>\% mir_to_gsymbol()
}
