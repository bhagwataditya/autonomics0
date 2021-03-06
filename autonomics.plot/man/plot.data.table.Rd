% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_sample_distributions.R
\name{plot.data.table}
\alias{plot.data.table}
\title{Plot}
\usage{
\method{plot}{data.table}(
  plotdt,
  geom,
  stat,
  ...,
  coord_flip = FALSE,
  xlab = NULL,
  ylab = NULL,
  title = NULL
)
}
\arguments{
\item{geom}{string: ggplot2 geom}

\item{stat}{string: ggplot2 stat}

\item{...}{ggplot2 aesthetic specifications: color = subgroup, facet = condition, ...}

\item{coord_flip}{TRUE or FALSE: whether to flip coordinates}

\item{xlab}{if non-NULL, x axis label}

\item{ylab}{if non-NULL, y axis label}

\item{title}{if non-NULL, title text}

\item{object}{data.table}

\item{facet}{string: plotdt variable mapped to facet}
}
\description{
Plot
}
\examples{
if (require(autonomics.data)){
   require(magrittr)
   object <- autonomics.data::glutaminase \%>\% prepare_plot_dt()
   plot(object,
        geom  = 'line',
        stat  = 'density',
        x     = value,
        color = subgroup,
        group = sample_id,
        facet = subgroup)
}
}
