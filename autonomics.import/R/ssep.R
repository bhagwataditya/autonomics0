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
#' if (require(subramanian.2016)){
#'    subramanian.2016::metabolon   %>% autonomics.import::ssep()
#' }
#' if (require(graumann.lfq)){
#'    graumann.lfq::lfq.intensities %>% autonomics.import::ssep()
#' }
#' @importFrom magrittr %>%
#' @export
ssep <- function(object, svar = 'subgroup'){
   object %>%
   autonomics.import::slevels(svar) %>%
   autonomics.import::infer_design_sep()
}

#' @rdname ssep
#' @importFrom  magrittr %>%
#' @export
subgroup_sep <- function(object){
   object %>% ssep('subgroup')
}

