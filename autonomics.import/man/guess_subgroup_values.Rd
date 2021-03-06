% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/guess.R
\name{guess_subgroup_values}
\alias{guess_subgroup_values}
\alias{guess_subgroup_values.character}
\alias{guess_subgroup_values.SummarizedExperiment}
\title{Guess subgroup values}
\usage{
guess_subgroup_values(x, ...)

\method{guess_subgroup_values}{character}(
  x,
  sep = x \%>\% guess_sep(),
  invert = FALSE,
  verbose = FALSE,
  ...
)

\method{guess_subgroup_values}{SummarizedExperiment}(
  x,
  sep = x \%>\% guess_sep(),
  invert = FALSE,
  verbose = FALSE,
  ...
)
}
\arguments{
\item{x}{charactervector, SummarizedExperiment}

\item{...}{used for proper S3 method dispatch}

\item{sep}{character(1)}

\item{invert}{FALSE (default) or TRUE: whether to guess "non-subgroup" component}

\item{verbose}{logical(1)}
}
\value{
character(n)
}
\description{
Guess subgroup values
}
\examples{
require(magrittr)

# charactervector
   # No sep: subgroup = x
      x <- c("EM00", "EM01", "EM02")
      x \%>\% guess_subgroup_values()

   # Sep: subgroup = head components of x
      x <- c("UT_10h_R1", "UT_10h_R2", "UT_10h_R3")
      x \%>\% guess_subgroup_values()
      x \%>\% guess_subgroup_values(invert = TRUE)

      x <- c("EM00_STD.R1", "EM01_STD.R1", "EM01_EM00.R1")
      x \%>\% guess_subgroup_values()
      x \%>\% guess_subgroup_values(invert = TRUE)

}
