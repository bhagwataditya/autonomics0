#' Is max quant eSet?
#' @param object eSet
#' @return logical
#' @examples
#' if (require(autonomics.data)){
#'    library(magrittr)
#'    object <- autonomics.data::stemcomp.proteinratios
#'    object %>% is_maxquant_eset()
#' }
#' @export
is_maxquant_eset <- function(object){
   assert_is_valid_object(object)
   all(c('Contaminant', 'Reverse') %in% fvars(object))
}

#' Is RNA seq eset
#' @param object  eset
#' @export
is_rnaseq_eset <- function(object){
   assert_is_valid_object(object)
   if (!contains_prepro(object))     return(FALSE)
   prepro(object)$assay == 'rnaseq'
}

#' Is object a soma eset?
#' @param object eset
#' @return logical
#' @examples
#' if (require(autonomics.data)){
#' require(magrittr)
#'    object <- autonomics.data::stemcomp.soma
#'    object %>% is_soma_eset()
#' }
#' @export
is_soma_eset <- function(object){
   assert_is_valid_object(object)
   ('SeqId' %in% fvars(object))
}

#' Is object a metabolon eset?
#' @param object eset
#' @return logical
#' @examples
#' if (require(autonomics.data)){
#'    require(magrittr)
#'    object <- autonomics.data::glutaminase
#'    object %>% is_metabolon_eset()
#' }
#' @export
is_metabolon_eset <- function(object){
   assert_is_valid_object(object)
   all(c('COMP_ID', 'BIOCHEMICAL') %in% fvars(object))
}

#' Is object an exiqon eset?
#' @param object   eset
#' @return logical
#' @export
is_exiqon_eset <- function(object){
   assert_is_valid_object(object)
   if (!contains_prepro(object))     return(FALSE)
   prepro(object)$assay == 'exiqon'
}
