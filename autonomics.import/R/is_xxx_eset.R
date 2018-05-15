#' Is max quant eSet?
#' @param object eSet
#' @return logical
#' @examples
#' library(magrittr)
#' if (require(autonomics.data)){
#'    object <- autonomics.data::billing2016
#'    object %>% is_maxquant_eset()
#' }
#' @export
is_maxquant_eset <- function(object){
   autonomics.import::assert_is_valid_eset(object)
   if (!autonomics.import::contains_prepro(object))     return(FALSE)
   autonomics.import::prepro(object)$software == 'maxquant'
}

#' Is RNA seq eset
#' @param object  eset
#' @examples
#' if (require(billing.differentiation.data)){
#'    require(magrittr)
#'    object <- rna.voomcounts
#'    object %>% is_rnaseq_eset()
#' }
#' @export
is_rnaseq_eset <- function(object){
   autonomics.import::assert_is_valid_eset(object)
   if (!autonomics.import::contains_prepro(object))     return(FALSE)
   autonomics.import::prepro(object)$assay == 'rnaseq'
}

#' Is object a soma eset?
#' @param object eset
#' @return logical
#' @export
is_soma_eset <- function(object){
   autonomics.import::assert_is_valid_eset(object)
   if (!autonomics.import::contains_prepro(object))     return(FALSE)
   autonomics.import::prepro(object)$assay == 'somascan'
}

#' Is object a metabolon eset?
#' @param object eset
#' @return logical
#' @examples
#' require(magrittr)
#' if (require(subramanian.2016)){
#'    object <- subramanian.2016::metabolon
#'    object %>% autonomics.import::is_metabolon_eset()
#' }
#' @export
is_metabolon_eset <- function(object){
   autonomics.import::assert_is_valid_eset(object)
   if (!autonomics.import::contains_prepro(object))     return(FALSE)
   autonomics.import::prepro(object)$software == 'metabolon'
}

#' Is object an exiqon eset?
#' @param object   eset
#' @return logical
#' @examples
#' require(magrittr)
#' if (require(subramanian.2016)){
#'    object <- subramanian.2016::exiqon
#'    object %>% autonomics.import::is_exiqon_eset()
#' }
#' @export
is_exiqon_eset <- function(object){
   autonomics.import::assert_is_valid_eset(object)
   if (!autonomics.import::contains_prepro(object))     return(FALSE)
   autonomics.import::prepro(object)$assay == 'exiqon'
}
