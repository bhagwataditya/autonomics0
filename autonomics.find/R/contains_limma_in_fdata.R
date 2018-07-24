#' Does eset contain limma results
#' @param object eset
#' @return logical
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    object <- autonomics.data::stemcomp.proteinratios
#'    object %>% autonomics.find::contains_limma_in_fdata()
#'    object %>% autonomics.find::add_limma_to_fdata() %>% 
#'               autonomics.find::contains_limma_in_fdata()
#' }
#' if (require(billing.differentiation.data)){
#'    billing.differentiation.data::rna.voomcounts %>% 
#'       contains_limma_in_fdata()
#' }
#' @importFrom magrittr  %>% 
#' @export
contains_limma_in_fdata <- function(object){
   autonomics.import::fvars(object) %>% 
      stringi::stri_detect_fixed('coef.') %>% 
      any()
}
