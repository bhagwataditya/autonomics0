

#' Split (subgroup) components
#'
#' Get data.table in which subgroup components have been decomposed
#'
#' @param values character vector with subgroup values
#' @param sep character(1): separator
#' @return data.table with decomposed subgroup components
#' @examples
#' require(magrittr)
#' if (require(billing.differentiation.data)){
#'    object <- billing.differentiation.data::protein.ratios
#'    object %>% subgroup_values() %>% split_components()
#' }
#' if (require(autonomics.data)){
#'    autonomics.data::glutaminase %>% subgroup_values()  %>% split_components()
#' }
#' if (require(subramanian.2016)){
#'    subramanian.2016::metabolon %>% subgroup_values()   %>% split_components()
#' }
#' if (require(graumann.lfq)){
#'    graumann.lfq::lfq.intensities %>% subgroup_values() %>% split_components()
#' }
#' @importFrom magrittr %>%
#' @export
split_components <- function(
   values,
   sep = autonomics.import::infer_design_sep(values, verbose = FALSE)
){

   # Single component
   if (is.null(sep)) return(data.table::data.table(V1 = values))

   # Multiple components
   y <- values %>% stringi::stri_split_fixed(sep)
   n.component <- length(y[[1]])
   1:n.component %>% lapply(function(z) vapply(y, magrittr::extract, character(1), z)) %>%
                     data.table::as.data.table()
}

#' Split (subgroup) components
#'
#' Get data.table in which subgroup components have been decomposed
#'
#' @param object SummarizedExperiment
#' @param svar character
#' @return data.table with decomposed subgroup components
#' @examples
#' require(magrittr)
#' if (require(billing.differentiation.data)){
#'    object <- billing.differentiation.data::protein.ratios
#'    object %>% subgroup_components()
#'    object %>% scomponents('subgroup')
#' }
#' if (require(autonomics.data)){
#'    object <- autonomics.data::glutaminase
#'    object %>% subgroup_components()
#'    object %>% scomponents('subgroup')
#' }
#' if (require(subramanian.2016)){
#'    object <- subramanian.2016::metabolon
#'    object %>% subgroup_components()
#'    object %>% scomponents('subgroup')
#' }
#' if (require(graumann.lfq)){
#'    object <- graumann.lfq::lfq.intensities
#'    object %>% subgroup_components()
#'    object %>% scomponents('subgroup')
#' }
#' @importFrom magrittr %>%
#' @export
subgroup_components <- function(object){
   object %>% autonomics.import::subgroup_values() %>% autonomics.import::split_components()
}

#' @rdname subgroup_components
#' @importFrom magrittr %>%
#' @export
scomponents <- function(object, svar){
   object %>% autonomics.import::svalues(svar) %>% autonomics.import::split_components()
}
