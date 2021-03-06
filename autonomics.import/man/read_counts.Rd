% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_omics.R
\name{read_counts}
\alias{read_counts}
\title{Read rnaseq counts}
\usage{
read_counts(file, fid_var, fname_var = character(0))
}
\arguments{
\item{file}{string: path to rnaseq counts file}

\item{fid_var}{string or number: feature id variable}

\item{fname_var}{string or number: feature name variable}
}
\description{
Read tsv file with rnaseq counts into SummarizedExperiment
}
\details{
File format: header row
             feature annotations in first few columns
             feature counts      in next columns
}
\examples{
if (require(autonomics.data)){
   require(magrittr)
   file <- 'extdata/stemcomp/rnaseq/gene_counts.txt' \%>\%
            system.file(package = 'autonomics.data')
   file \%>\% read_counts(fid_var = 'gene_id', fname_var = 'gene_name')
}
}
\seealso{
merge_sdata, merge_fdata
}
