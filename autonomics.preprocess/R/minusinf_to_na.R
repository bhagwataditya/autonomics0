#' Replace -inf by NA
#'
#' Minus inf can be created as log(0). This function replaces such minus inf by NA
#'
#' @param object SummarizedExperiment, eSet, or EList
#' @return updated object
#' @export
minusinf_to_na <- function(object){
   idx <- autonomics.import::exprs(object)==-Inf
   if (sum(idx, na.rm=TRUE)>0){
      autonomics.support::cmessage('\t\tReplace -Inf by NA')
      autonomics.import::exprs(object)[idx] <- NA
   }
   object
}

#' Convert NAs to zeroes
#' @param object    SummarizedExperiment
#' @return Updated SummarizedExperiment
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    object <- autonomics.data::glutaminase
#'    object %>% na_to_zero()
#'    object %>% na_to_zero(condition)
#' }
#' @importFrom magrittr %>%
#' @export
na_to_zero <- function(object){
   selector <- autonomics.import::exprs(object) %>% is.na()
   autonomics.support::cmessage('\t\tReplace NA -> 0 for %d/%d features per sample (avg)',
                                selector %>% colSums(na.rm=TRUE) %>% mean() %>% round(),
                                nrow(object))
   autonomics.import::exprs(object)[selector] <- 0
   object

}

#' Convert zeroes to NAs
#' @param object    SummarizedExperiment
#' @return Updated SummarizedExperiment
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    object <- autonomics.data::glutaminase %>% na_to_zero()
#'    object %>% autonomics.preprocess::zero_to_na()
#' }
#' @importFrom magrittr %>%
#' @export
zero_to_na <- function(object){
   selector <- autonomics.import::exprs(object) == 0
   autonomics.support::cmessage('\t\tReplace 0 -> NA for %d/%d features per sample (avg)',
                                selector %>% colSums(na.rm=TRUE) %>% mean() %>% round(),
                                nrow(object))
   autonomics.import::exprs(object)[selector] <- NA_real_
   object
}
