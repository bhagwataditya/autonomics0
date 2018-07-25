

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
#' if (require(autonomics.data)){
#'    autonomics.data::billing2016 %>% autonomics.find::default_contrasts()
#' }
#' if (require(billing.differentiation.data)){
#'    billing.differentiation.data::rna.voomcounts %>%
#'       autonomics.find::default_contrasts() %>% data.frame()
#' }
#' if (require(halama.2016)){
#'    halama.2016::cell.metabolites %>%
#'       autonomics.import::filter_samples(TREATMENT != 0) %>%
#'       autonomics.find::default_contrasts() %>% data.frame()
#' }
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
      autonomics.support::cmessage('\t\tUse contrasts in autonomics.import::contrastdefs(object)')
      return(autonomics.import::contrastdefs(object))
   }

   # Extract subgroup levels
   subgroup_values <- autonomics.import::sdata(object)$subgroup
   if (is.character(subgroup_values))  subgroup_values %<>% autonomics.support::factorify()

   # Ratios for ratio esets
   if (autonomics.import::contains_ratios(object)){
      autonomics.support::cmessage('\t\tGenerate contrasts from ratios')
      contrasts <- subgroup_values %>% levels() %>% magrittr::set_names(validify_contrast_names(.))

   # Reference contrasts for intensity esets
   } else {
      autonomics.support::cmessage('\t\tGenerate contrasts to reference level')
      contrasts <- object %>% autonomics.find::make_ref_contrasts()
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

#' default top_definition
#' @param object eset
#' @return default value of top_definition
#' @export
#' @examples
#' require(magrittr)
#' if (require(billing.differentiation.data)){
#'    billing.differentiation.data::rna.voomcounts %>% 
#'       autonomics.find::default_top_definition()
#' }
default_top_definition <- function(object){

   subgroups             <- !is.null(autonomics.import::sdata(object)$subgroup)
   replicates            <- anyDuplicated(autonomics.import::sdata(object)$subgroup)
   sufficient_features   <- nrow(object) > 50

   if (subgroups & replicates & sufficient_features)   'p < 0.05'
   else                                                '(quantile < 0.05 | quantile > 0.95)'
}

