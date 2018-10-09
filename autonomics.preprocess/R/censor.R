#' Censor exprs
#' @param object          SummarizedExperiment
#' @param threshold       numeric scalar: threshold beyond which to censor values
#' @param censored_value  scalar
#' @return SummarizedExperiment
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    object <- autonomics.data::glutaminase
#'    object %>% autonomics.preprocess::left_censor(13)
#'    object %>% autonomics.preprocess::right_censor(30)
#' }
#' @importFrom magrittr %>%
#' @export
left_censor <- function(
   object,
   threshold,
   censored_value = threshold
){
   idx <- autonomics.import::exprs(object) %>% autonomics.support::is_less_than(threshold)
   autonomics.import::exprs(object)[idx] <- censored_value
   autonomics.support::cmessage('\t\tReplace exprs < %s with %s in %d/%d features of %d/%d samples',
                                as.character(threshold),
                                as.character(censored_value),
                                sum(matrixStats::rowAnys(idx)), nrow(idx),
                                sum(matrixStats::colAnys(idx)), ncol(idx))
   object
}

#' @rdname left_censor
#' @importFrom magrittr %>%
#' @export
right_censor <- function(
   object,
   threshold,
   censored_value = threshold
){
   idx <- autonomics.import::exprs(object) %>% autonomics.support::is_greater_than(threshold)
   autonomics.import::exprs(object)[idx] <- censored_value
   autonomics.support::cmessage('\t\tReplace exprs exprs > %s with %s in %d/%d features of %d/%d samples',
                                as.character(threshold),
                                as.character(censored_value),
                                sum(matrixStats::rowAnys(idx)), nrow(idx),
                                sum(matrixStats::colAnys(idx)), ncol(idx))
   object
}
