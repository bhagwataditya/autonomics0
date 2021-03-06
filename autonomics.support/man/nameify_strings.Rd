% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nameify_strings.R
\name{nameify_strings}
\alias{nameify_strings}
\title{Nameify strings}
\usage{
nameify_strings(x, verbose = TRUE)
}
\arguments{
\item{x}{character vector}

\item{verbose}{TRUE or FALSE}
}
\value{
character vector
}
\description{
Update character vector values so that they can serve as R variable names
}
\examples{
x <- c('266S', '270A', 'DM')
nameify_strings(x)
}
