#' Split SummarizedExperiment by svar
#' @param object SummarizedExperiment
#' @param svar 'subgroup'
#' @return list
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    autonomics.data::stemdiff.proteinratios %>%
#'    autonomics.import::split_by_svar()
#' }
#' @export
split_by_svar <- function(object, svar = 'subgroup'){
   Map(function(sg) object %>% autonomics.import::filter_samples_(sprintf("%s %%in%% '%s'", svar, sg)),
       autonomics.import::slevels(object, svar))
}
