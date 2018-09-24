#' Get svar separator
#' @param object SummarizedExperiment, eSet, or EList
#' @param svar sample variable
#' @return string or NULL (if no separator found)
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    file <- system.file('extdata/glutaminase/glutaminase.xlsx',
#'                         package = 'autonomics.data')
#'    object <- autonomics.import::load_metabolon(file)
#'    object  %>% ssep()
#' }
#' if (require(autonomics.data)){
#'    autonomics.data::stemcomp.proteinratios %>% autonomics.import::ssep()
#' }
#' if (require(subramanian.2016)){
#'    subramanian.2016::metabolon   %>% autonomics.import::ssep()
#' }
#' if (require(graumann.lfq)){
#'    graumann.lfq::lfq.intensities %>% autonomics.import::ssep()
#' }
#' @importFrom magrittr %>%
#' @export
ssep <- function(object, svar = 'subgroup'){
   possible_separators <- if (autonomics.import::contains_ratios(object)) c('.', ' ') else c('.', '_', ' ')
   object %>% autonomics.import::slevels(svar) %>%
              autonomics.import::infer_design_sep(possible_separators = possible_separators, verbose = FALSE)
}

#' @rdname ssep
#' @importFrom  magrittr %>%
#' @export
subgroup_sep <- function(object){
   object %>% ssep('subgroup')
}

