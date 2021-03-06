% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_features.R
\name{gg_xvalues}
\alias{gg_xvalues}
\title{Get x values}
\usage{
gg_xvalues(b)
}
\arguments{
\item{b}{built ggplot (result of ggplot_build)}
}
\value{
x values
}
\description{
Get x values from faceted ggplot object
}
\examples{
if (require(autonomics.data)){
   require(magrittr)
   object <- autonomics.data::stemcomp.proteinratios \%>\% magrittr::extract(1:9, )
   b <- object \%>\% plot_features('point') \%>\% ggplot2::ggplot_build()
   b \%>\% gg_xvalues()
}
}
