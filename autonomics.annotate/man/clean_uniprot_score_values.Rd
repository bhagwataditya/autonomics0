% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/uniprot.R
\name{clean_uniprot_score_values}
\alias{clean_uniprot_score_values}
\title{Clean uniprot "score" values}
\usage{
clean_uniprot_score_values(values)
}
\arguments{
\item{values}{character vector}
}
\value{
character vectors
}
\description{
Clean uniprot "score" values
}
\examples{
\dontrun{
require(magrittr)
values <- c("A0A024R4M0", "M0R210", "G3HSF3")
up <- values \%>\% autonomics.annotate::connect_to_uniprot()
annotations <- values \%>\% autonomics.annotate::fetch_uniprot_annotations(up)
annotations$SCORE \%T>\% print() \%>\% autonomics.annotate::clean_uniprot_score_values()
}
}
