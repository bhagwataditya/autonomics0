% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract.R
\name{extract_first_from_collapsed}
\alias{extract_first_from_collapsed}
\alias{extract_first_from_collapsed.character}
\alias{extract_first_from_collapsed.factor}
\title{Extract first from collapsed values}
\usage{
extract_first_from_collapsed(x, ...)

\method{extract_first_from_collapsed}{character}(x, sep = guess_sep(x), ...)

\method{extract_first_from_collapsed}{factor}(x, sep = guess_sep(x), ...)
}
\arguments{
\item{x}{charactervector or factorvector}

\item{...}{to allow for S3 method dispatch}

\item{sep}{collapsed string separator, e.g. ';'}
}
\value{
Updated x
}
\description{
Extract first from collapsed values
}
\examples{
require(magrittr)
x <- c('a;b;c', '1;2;3', 'alpha;beta;gamma')
x \%>\% extract_first_from_collapsed(sep = ';')
}
