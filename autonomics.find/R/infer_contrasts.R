#' Infer contrast names
#' @param object eset
#' @return character vector
#' @examples 
#' require(magrittr)
#' if (require(subramanian.2016)){
#'    subramanian.2016::exiqon %>%  
#'    autonomics.find::infer_contrast_names()
#' }
#' @importFrom magrittr  %>% 
#' @export
infer_contrast_names <- function(object){

   limma <- autonomics.import::limma(object)
      
   if (is.null(limma))   return(character(0))
   
   colnames(limma)

}