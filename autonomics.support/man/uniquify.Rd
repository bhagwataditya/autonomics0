% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/factors.R
\name{uniquify}
\alias{uniquify}
\title{Make vector components unique by appending spaces}
\usage{
uniquify(x, method = "make.unique.spaces", verbose = TRUE)
}
\arguments{
\item{x}{character or factor vector}

\item{method}{'make.unique' or 'make.unique.spaces'}

\item{verbose}{TRUE (default) or FALSE}
}
\value{
character vector
}
\description{
Make vector components unique by appending spaces
}
\examples{
x <- c('A', 'B', 'C', 'A', 'D')
uniquify(x, 'make.unique')
uniquify(x, 'make.unique.spaces')
}
\seealso{
\code{\link[base]{make.unique}}
}
