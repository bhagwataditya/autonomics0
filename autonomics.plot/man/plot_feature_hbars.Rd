% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_feature_hbars.R
\name{plot_feature_hbars}
\alias{plot_feature_hbars}
\title{Plot feature bars}
\usage{
plot_feature_hbars(
  object,
  fvars = default_fvars(object),
  color_var = default_color_var(object),
  color_values = default_color_values(object, color_var),
  facet_def = "~ sample",
  alpha_var = NULL,
  file = NULL,
  width = NULL,
  height = NULL,
  xlab = NULL,
  title = "",
  legend.position = "right",
  x = NULL,
  shape_var = NULL
)
}
\arguments{
\item{object}{\code{SummarizedExperiment}}

\item{fvars}{fvars used to annotate in plot}

\item{color_var}{svar mapped to color}

\item{color_values}{color value vector (names = subgroups, contents = colours)}

\item{facet_def}{facet definition string}

\item{alpha_var}{svar mapped to transparancy}

\item{file}{name of file to write to}

\item{width}{width (inches)}

\item{height}{height (inches)}

\item{xlab}{xlab annotation of plot}

\item{title}{plot title}

\item{legend.position}{legend position}

\item{x}{NULL. Generifies interface of plot_xxx functions.}

\item{shape_var}{NULL. Generifies interface of plot_xxx functions.}
}
\description{
Create bar plot of exprs(object). Annotate it with fvars. Color per 'color_var'. Print to 'file'
}
\examples{
require(magrittr)
if (require(autonomics.data)){
   small_eset <- autonomics.data::stemcomp.proteinratios[1:10, ]
   plot_feature_hbars(small_eset)
   plot_feature_hbars(small_eset, fvars = 'Gene names')
}
}
