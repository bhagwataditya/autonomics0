% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/impute_consistent_nas.R
\name{impute_around_zero}
\alias{impute_around_zero}
\title{Impute around zero}
\usage{
impute_around_zero(x)
}
\arguments{
\item{x}{exprs matrix}
}
\value{
exprs matrix
}
\description{
Impute around zero
}
\examples{

# Example data
#-------------
x <- matrix(c(20, 21, 22, 23,
              30, 31, 32, 33,
              NA, NA, NA, NA,
              40, 41, NA, NA,
              NA, NA, NA, NA,
              50, 51, NA, NA), ncol=4, byrow=TRUE)
object <- SummarizedExperiment::SummarizedExperiment(assays = list(exprs=x))
object$subgroup <- c('EV', 'EV', 'MS', 'MS')

# Impute
#-------
x
impute_around_zero(x)

exprs(object)
exprs(impute_consistent_nas(object, imputefun = impute_around_zero))
exprs(impute_common_nas(object, imputefun = impute_around_zero))
}
