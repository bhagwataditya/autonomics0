% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/uniprot.R
\name{clean_uniprot_reviewed_values}
\alias{clean_uniprot_reviewed_values}
\title{Clean uniprot "reviewed" values}
\usage{
clean_uniprot_reviewed_values(values)
}
\arguments{
\item{values}{uniprot "reviewed" values (character vector)}
}
\value{
cleaned values (character vector)
}
\description{
Clean uniprot "reviewed" values
}
\examples{
\dontrun{
require(magrittr)
values <- c("A0A024R4M0", "M0R210", "G3HSF3")
up <- values \%>\% autonomics.annotate::connect_to_uniprot()
annotations <- values \%>\% autonomics.annotate::fetch_uniprot_annotations(up)
annotations$REVIEWED \%T>\% print() \%>\% autonomics.annotate::clean_uniprot_reviewed_values()
}
}
