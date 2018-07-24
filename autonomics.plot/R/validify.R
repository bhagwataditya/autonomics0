#' Validify shape_var values
#' @param object SummarizedExperiment
#' @param shape_var svar mapped to shape
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    autonomics.data::stemcomp.proteinratios %>%
#'    autonomics.plot::validify_shape_values('replicate')
#' }
#' @importFrom magrittr %<>%
#' @export
validify_shape_values <- function(object, shape_var){
   if (is.null(shape_var)) return(object)
   if (is.numeric(autonomics.import::svalues(object, shape_var))){
      autonomics.import::svalues(object, shape_var) %<>% as.character()
   }
   return(object)
}
