#' Replace -inf by NA
#'
#' Minus inf can be created as log(0). This function replaces such minus inf by NA
#'
#' @param object SummarizedExperiment, eSet, or EList
#' @return updated object
#' @export
minusinf_to_na <- function(object){
   idx <- autonomics.import::exprs(object)==-Inf
   if (sum(idx, na.rm=TRUE)>0){
      autonomics.support::cmessage('\t\tReplace -Inf by NA')
      autonomics.import::exprs(object)[idx] <- NA
   }
   object
}
