#' Split SummarizedExperiment by svar
#' @param object SummarizedExperiment
#' @param svar 'subgroup'
#' @return list
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    autonomics.data::stemdiff.proteinratios %>%
#'    split_by_svar()
#' }
#' @export
split_by_svar <- function(object, svar = 'subgroup'){
   subgroup <- NULL
   Map(function(sg) object %>% autonomics.import::filter_samples(subgroup %in% sg),
       autonomics.import::slevels(object, svar))
}
