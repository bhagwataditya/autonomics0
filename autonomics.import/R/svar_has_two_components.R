#' Do values consist of two components?
#' @param x character vector
#' @return logical(1)
#' @export
has_two_components <- function(x){

   dt <- x %>% autonomics.import::split_components()

   # FALSE if less than two columns
   if (ncol(dt) < 1) return(FALSE)

   # FALSE if only one V2 component per V1 component
   if (nrow(dt[, .SD[.N>1], by = 'V1'])==0) return(FALSE)

   # TRUE if dt has two components
   return(ncol(dt)==2)

}

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
#'
#' # GLUTAMINASE (METABOLON)
#' if (require(autonomics.data)){
#'    file <- system.file('extdata/glutaminase/glutaminase.xlsx',
#'                         package = 'autonomics.data')
#'    object <- autonomics.import::load_metabolon(file)
#'    object %>% autonomics.import::subgroup_has_two_components()
#' }
#'
#' # STEM CELL DIFFERENTIATION (MAXQUANT)
#' if (require(autonomics.data)){
#'    autonomics.data::stemdiff.proteinratios %>%
#'    autonomics.import::subgroup_has_two_components()
#' }
#' @export
svar_has_two_components <- function(object, svar){
   object %>%
   autonomics.import::slevels(svar) %>%
   autonomics.import::has_two_components()
}

#' @rdname svar_has_two_components
#' @importFrom magrittr  %>%
#' @export
subgroup_has_two_components <- function(object){
   object %>%
   svar_has_two_components('subgroup')
}
