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
   svar_levels <- object %>% autonomics.import::slevels(svar)
   n_components <- 1
   # Note: ' ' will not occur: subgroups have been undergone make.names() in autonomics.import
   possible_seps <- if(autonomics.import::is_maxquant_eset(object)) '.' else  c('.', '_') # in LCMS proteomics, this indicates ratios
   for (sep in possible_seps){
      subgroup_components <- svar_levels %>% stringi::stri_split_fixed(sep)
      n_components <- subgroup_components %>%
         vapply(length, numeric(1)) %>%
         (function(x) if (all(x[1]==x))   x[1]   else  1 )
      if (n_components > 1)   return(sep)
   }
   return(NULL)
}

#' @rdname ssep
#' @importFrom  magrittr %>%
#' @export
subgroup_sep <- function(object){
   object %>% ssep('subgroup')
}

