# NOTE: roxygenize('.') gives errors, but document('') works!!!
# Best to use devtools::document() I think!


#' Annotated and autonomics-ready version of ALL data
#' 
#' This is a clean and annotated version of the well-known \code{\link[ALL]{ALL}}
#' dataset from the Bioconductor package 'ALL' (which in its recent version 
#' doesn't contain feature annotations).
#' 
#' Modifications in the pData:
#' \enumerate{
#'   \item rm samples with an NA value for pvar 'sex'
#'   \item adding pvar 'subgroup' as a combination of B/T and sex
#' }
#' 
#' Modifications in the fData:
#' \enumerate{
#'   \item adding feature annotations using the hguav2.db package
#'   \item rm features with an NA gene symbol
#'   \item rm duplicate gene symbol features
#' }
#' @seealso \code{\link[ALL]{ALL}}
"ALL"