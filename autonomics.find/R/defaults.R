

#' Validify contrast names
#' @param x vector with contrast names
#' @return vector with valid names
#' @importFrom magrittr %>%
#' @export
validify_contrast_names <- function(x){
   x %>% gsub(' ', '',  ., fixed = TRUE) %>%
         gsub('-', '_', ., fixed = TRUE) %>%
         make.names()
}


#' Default contrasts
#' @param object eset
#' @return contrast definitions
#' @examples
#' require(magrittr)
#' 
#' # STEM CELL COMPARISON
#' if (require(autonomics.data)){
#'    object <- autonomics.data::stemcomp.proteinratios 
#'    object %>% autonomics.find::default_contrasts()
#' }
#' 
#' # STEM CELL DIFFERENTIATION
#' if (require(billing.differentiation.data)){
#'    billing.differentiation.data::rna.voomcounts %>%
#'       autonomics.find::default_contrasts() %>% data.frame()
#' }
#' 
#' # GLUTAMINASE
#' if (require(autonomics.data)){
#'    object <- autonomics.data::glutaminase
#'    object %>% autonomics.import::filter_samples(TREATMENT != 0)  %>%
#'               autonomics.find::default_contrasts()               %>% 
#'               data.frame()
#' }
#' 
#' if (require(subramanian.2016)){
#'    subramanian.2016::metabolon %>%
#'       autonomics.find::default_contrasts() %>%
#'       data.frame()
#' }
#' # if (require(atkin.2014)){
#' #    atkin.2014::metabolon %>%
#' #      autonomics.find::default_contrasts() %>%
#' #      data.frame()
#' # }
#' if (require(graumann.lfq)){
#'    graumann.lfq::lfq.intensities %>%
#'       autonomics.find::default_contrasts() %>%
#'       data.frame()
#' }
#' @importFrom magrittr   %>%
#' @export
default_contrasts <- function(object){
   
   # If contrasts present in object, use these
   if (!is.null(autonomics.import::contrastdefs(object))){
      autonomics.support::cmessage('\tUse contrasts in autonomics.import::contrastdefs(object)')
      return(autonomics.import::contrastdefs(object))
   }

   # Extract subgroup levels
   subgroup_values <- autonomics.import::sdata(object)$subgroup
   if (is.character(subgroup_values))  subgroup_values %<>% autonomics.support::factorify()

   # Ratios for ratio esets
   if (autonomics.import::contains_ratios(object)){
      autonomics.support::cmessage('\tGenerate contrasts from ratios')
      contrasts <- subgroup_values %>% levels() %>% magrittr::set_names(validify_contrast_names(.))

   # Reference contrasts for intensity esets
   } else {
      autonomics.support::cmessage('\tGenerate contrasts to reference level')
      contrasts <- diff_contrasts(object)
   }

   # Return
   contrasts
}

#' Default nplot
#' @param object   eset
#' @return number of top features to be plotted
#' @export
#' @examples
#' require(magrittr)
#' if (require(billing.differentiation.data)){
#'    billing.differentiation.data::rna.voomcounts %>% 
#'       autonomics.find::default_nplot()
#' }
default_nplot <- function(object){
   16
   # 500
}

#' default topdef
#' @param object eset
#' @return default value of topdef
#' @export
#' @examples
#' require(magrittr)
#' if (require(billing.differentiation.data)){
#'    billing.differentiation.data::rna.voomcounts %>% 
#'       autonomics.find::default_topdef()
#' }
default_topdef <- function(object){

   subgroups             <- !is.null(autonomics.import::sdata(object)$subgroup)
   replicates            <- any(duplicated(autonomics.import::sdata(object)$subgroup)) # anyDuplicated sometimes returns a number!

   if (subgroups & replicates)   'p < 0.05'
   else                          '(effect < 1 | effect > 1)'
}

