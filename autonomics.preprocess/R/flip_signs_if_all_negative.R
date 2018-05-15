#' Flip sign if all expr values are negative
#' @param object SummarizedExperiment, eSet, or EList
#' @return updated object
#' @importFrom magrittr %>% %<>%
#' @export
flip_sign_if_all_exprs_are_negative <- function(object){
   idx <- !is.na(autonomics.import::exprs(object))
   if (all(sign(autonomics.import::exprs(object)[idx])==-1)){
      autonomics.support::cmessage('\t\tAll values negative: flip signs to prevent singularities.')
      autonomics.import::exprs(object) %<>% magrittr::multiply_by(-1)
   }
   object
}
