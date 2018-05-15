#' Infer contrast names
#' @param object eset
#' @return character vector
#' @examples 
#' require(magrittr)
#' if (require(subramanian.2016)){
#'    subramanian.2016::exiqon %>%  
#'       autonomics.find::infer_contrast_names()
#' }
#' @importFrom magrittr  %>% 
#' @export
infer_contrast_names <- function(object){
   
   if (!autonomics.find::contains_limma_in_fdata(object))   return(character(0))
   
   autonomics.import::fvars(object) %>% 
      magrittr::extract(stringi::stri_detect_fixed(., 'rank.')) %>% 
      stringi::stri_replace_first_fixed('rank.', '')
   
}