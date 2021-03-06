% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/map_ids.R
\name{entrezg_set_to_uniprot_set}
\alias{entrezg_set_to_uniprot_set}
\title{Map entrezg ids to uniprot accessions}
\usage{
entrezg_set_to_uniprot_set(x, annotation_map)
}
\arguments{
\item{x}{vector with entrezg ids}

\item{annotation_map}{org.xx.eg.db object}
}
\value{
vector with uniprot accessions
}
\description{
Map entrezg ids to uniprot accessions
}
\examples{
\dontrun{
   entrezg_set_to_uniprot_set(c('1', '2'), org.Hs.eg.db::org.Hs.eg.db)
}
}
