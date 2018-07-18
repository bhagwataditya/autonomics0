

#' Get subgroup components
#'
#' Get data.table in which subgroup components have been decomposed
#'
#' @param object  SummarizedExperiment, eSet, or EList
#' @param svar    sample var
#' @return data.table with decomposed subgroup components
#' @examples
#' require(magrittr)
#' if (require(billing.differentiation.data)){
#'    billing.differentiation.data::protein.ratios %>% scomponents('subgroup')
#'    billing.differentiation.data::protein.ratios %>% subgroup_components()
#' }
#' if (require(autonomics.data)){
#'    file <- system.file('extdata/glutaminase/glutaminase.xlsx',
#'                         package = 'autonomics.data')
#'    object <- autonomics.import::load_metabolon(file)
#'    object  %>% scomponents('subgroup')
#'    object  %>% subgroup_components()
#' }
#' if (require(subramanian.2016)){
#'    subramanian.2016::metabolon    %>% scomponents('subgroup')
#'    subramanian.2016::metabolon    %>% subgroup_components()
#' }
#' if (require(graumann.lfq)){
#'    graumann.lfq::lfq.intensities  %>% scomponents('subgroup')
#'    graumann.lfq::lfq.intensities  %>% subgroup_components()
#' }
#' @importFrom magrittr %>%
#' @export
scomponents <- function(object, svar){

   # Extract
   svar_levels <- object %>% autonomics.import::slevels(svar)
   sep         <- object %>% autonomics.import::ssep(svar)

   # Single component
   if (is.null(sep)) return(data.table::data.table(V1 = svar_levels))

   # Multiple components
   y <- svar_levels %>% stringi::stri_split_fixed(sep)
   n.component <- length(y[[1]])
   1:n.component %>% lapply(function(z) vapply(y, magrittr::extract, character(1), z)) %>%
      data.table::as.data.table()
}

#' @rdname scomponents
#' @importFrom magrittr %>%
#' @export
subgroup_components <- function(object){
   object %>% scomponents('subgroup')
}
