
#' Impute
#'
#' Data from proteomics and metabolomics LCMS experiments are missing not at random (MNAR).
#' They are well imputed with the QRILC (quantile regression imputation of left censored data).
#'
#' See \href{https://cran.r-project.org/package=imputeLCMD}{imputeLCMD package on CRAN}
#' See \href{https://www.biorxiv.org/content/early/2017/08/17/171967}{Wei et al. 2017}
#'
#' @param object SummarizedExperiment, eSet, or EList
#' @return updated object
#' @importFrom magrittr %<>% %>%
#' @export
impute <- function(object){
   autonomics.import::exprs(object) %<>% imputeLCMD::impute.QRILC() %>%
                                       magrittr::extract2(1)
   object
}


