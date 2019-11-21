#' Split SummarizedExperiment by svar
#' @param object SummarizedExperiment
#' @param svar 'subgroup'
#' @return list of SummarizedExperiments
#' @examples
#' if (require(autonomics.data)){
#'    require(magrittr)
#'    object <- autonomics.data::stemcomp.proteinratios
#'    object %>% split_by_svar()
#' }
#' @importFrom magrittr %>%
#' @export
split_by_svar <- function(object, svar = 'subgroup'){

   # Return object if null svar
   if (is.null(svar)) return(list(object))

   # Split
   Map(function(sg) object %>% filter_samples_(sprintf("%s %%in%% '%s'", svar, sg), record = FALSE),
       slevels(object, svar))
}

