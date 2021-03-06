% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_omics.R
\name{download_gtf}
\alias{download_gtf}
\title{Download GTF file}
\usage{
download_gtf(
  organism,
  release = 95,
  gtffile = sprintf("~/.autonomics/gtf/\%s", basename(make_gtf_url(organism, release)
    \%>\% substr(1, nchar(.) - 3)))
)
}
\arguments{
\item{organism}{'Homo sapiens', 'Mus musculus' or 'Rattus norvegicus'}

\item{release}{GTF release. By default release 95 selected}

\item{gtffile}{string: path to local GTF file}
}
\description{
Download GTF file with feature annotations
}
\examples{
\dontrun{ # requires internet and does not always work: https://stackoverflow.com/questions/55532102
   download_gtf(organism = 'Homo sapiens')
   download_gtf(organism = 'Mus musculus')
   download_gtf(organism = 'Rattus norvegicus')
}
}
