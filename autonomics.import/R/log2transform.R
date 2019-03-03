
#' Log2 transform
#' @param object SummarizedExperiment
#' @param verbose logical(1)
#' @return Updated SummarizedExperiment
#' @importFrom magrittr %<>%
#' @export
log2transform <- function(object, verbose = FALSE){
   if (verbose) message('\t\tLog2 transform')
   autonomics.import::exprs(object) %<>% log2()
   return(object)
}
