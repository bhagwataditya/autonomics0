#' Extract features
#' @param object SummarizedExperiment
#' @param extractor logical/numeric vector
#' @examples 
#' require(magrittr)
#' if (require(autonomics.data)){
#'    object <- autonomics.data::stemdiff.proteinratios
#'    (object %<>% autonomics.import::extract_features(c(1,2)))
#'    (object %<>% autonomics.import::extract_features(1))
#'    object %>% autonomics.import::limma()
#' }
#' @importFrom magrittr %<>%
#' @export
extract_features <- function(object, extractor){
   object %<>% magrittr::extract(extractor, )
   if (!is.null(autonomics.import::limma(object)))   autonomics.import::limma(object) %<>% magrittr::extract(extractor, , , drop = FALSE)
   object
}