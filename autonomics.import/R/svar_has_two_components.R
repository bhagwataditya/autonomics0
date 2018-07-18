#' Does svar have two components?
#'
#' @param object SummarizedExperiment, eSet, or EList
#' @param svar sample var
#' @return logical
#' @examples
#' require(magrittr)
#' if (require(subramanian.2016)){
#'    object <- subramanian.2016::metabolon
#'    object %>% autonomics.import::svar_has_two_components('subgroup')
#'    object %>% autonomics.import::subgroup_has_two_components()
#' }
#' if (require(autonomics.data)){
#'    file <- system.file('extdata/glutaminase/glutaminase.xlsx',
#'                         package = 'autonomics.data')
#'    object <- autonomics.import::load_metabolon(file)
#'    object %>% autonomics.import::subgroup_has_two_components()
#' }
#' if (require(billing.differentiation.data)){
#'    billing.differentiation.data::protein.ratios %>%
#'       autonomics.import::subgroup_has_two_components()
#' }
#' @export
svar_has_two_components <- function(object, svar){
   object %>%
   autonomics.import::scomponents(svar) %>%
   ncol() %>%
   magrittr::equals(2)
}

#' @rdname svar_has_two_components
#' @importFrom magrittr  %>%
#' @export
subgroup_has_two_components <- function(object){
   object %>% svar_has_two_components('subgroup')
}
