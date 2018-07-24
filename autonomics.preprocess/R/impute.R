
#' Impute
#'
#' Data from proteomics and metabolomics LCMS experiments are missing not at random (MNAR).
#' They are well imputed with the QRILC (quantile regression imputation of left censored data).
#'
#' See \href{https://cran.r-project.org/package=imputeLCMD}{imputeLCMD package on CRAN}
#' See \href{https://www.biorxiv.org/content/early/2017/08/17/171967}{Wei et al. 2017}
#'
#' @param object SummarizedExperiment, eSet, or EList
#' @param random_seed Integer used with \code{\link[R.utils]{withSeed}} for reproducible
#' imputation (the used \code{\link[imputeLCMD]{impute.QRILC}} implys random draws from a distribution)
#' @return updated object
#' @importFrom magrittr %<>% %>%
#' @export
impute <- function(object, random_seed = NULL){
   if(is.null(random_seed))
   {
     autonomics.import::exprs(object) %<>%
       imputeLCMD::impute.QRILC() %>%
       magrittr::extract2(1)
   } else {
     assertive.types::assert_is_an_integer(random_seed)
     autonomics.import::exprs(object) <- R.utils::withSeed(
       { 
         autonomics.import::exprs(object) %>%
           imputeLCMD::impute.QRILC() %>%
           magrittr::extract2(1)
       },
       seed = random_seed)
   }
   object
}


