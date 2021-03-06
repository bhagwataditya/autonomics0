% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocess_eset2.R
\name{preprocess_smart}
\alias{preprocess_smart}
\title{Sample normalize object automatically}
\usage{
preprocess_smart(
  object,
  plot = TRUE,
  result_dir = NULL,
  color_var = autonomics.plot::default_color_var(object),
  shape_var = autonomics.plot::default_shape_var(object),
  txt_var = autonomics.plot::default_txt_var(object)
)
}
\arguments{
\item{object}{eset}

\item{plot}{whether to plot sample distributions for different prepro steps}

\item{result_dir}{if non-NULL, directory to which plots will be printed}

\item{color_var}{svar mapped to color}

\item{shape_var}{svar mapped to shape}

\item{txt_var}{svar mapped to txt}
}
\value{
object after sample normalization
}
\description{
Sample normalize object automatically
}
\examples{
# max quant data
# --------------
library(magrittr)
if (require(autonomics.data)){
   autonomics.data::billing2016 \%>\% autonomics::preprocess_smart()
   autonomics.data::ALL \%>\% autonomics::preprocess_smart()
}
if (require(billing.differentiation.data)){
   billing.differentiation.data::protein.ratios \%>\%
      autonomics::preprocess_smart()
   billing.differentiation.data::phospho.occupancies \%>\%
      autonomics::preprocess_smart()
}
if (require(atkin.2014)){
   atkin.2014::soma  \%>\% 
      autonomics::preprocess_smart()
}
}
