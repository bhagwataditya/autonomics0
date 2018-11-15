#' Create arrow diagram
#' @param object SummarizedExperiment
#' @param topdef character(1): top definition
#' @return character(1)
#' @examples
#' if (require(autonomics.data)){
#'    object <- autonomics.data::stemdiff.proteinratios
#'    topdef <- 'fdr < 0.05'
#'    cat(create_arrow_diagram(object, topdef))
#' }
#' @importFrom magrittr %>%
#' @export
create_arrow_diagram <- function(object, topdef){
   posdef <- paste0(topdef, '& effect > 0')
   negdef <- paste0(topdef, '& effect < 0')
   s <- '$$'
   for (contrast_name in names(autonomics.import::contrastdefs(object))){
      s <- sprintf('%s\n\\xrightarrow[%d]{%d} %s', s,
                   autonomics.find::n_contrast_features(object, contrast_name, negdef),
                   autonomics.find::n_contrast_features(object, contrast_name, posdef),
                   contrast_name)
   }
   s <- sprintf('%s\n$$', s)
   s
}
