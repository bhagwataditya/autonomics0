#' Censor exprs
#' @param object     SummarizedExperiment
#' @param threshold  numeric scalar: threshold beyond which to censor values
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
left_censor <- function(object, threshold){
   idx <- autonomics.import::exprs(object) %>% autonomics.support::is_less_than(threshold)
   autonomics.import::exprs(object)[idx] <- threshold
   autonomics.support::cmessage('\t\tLeft censor exprs < %s in %d/%d features of %d/%d samples',
                                as.character(threshold),
                                sum(matrixStats::rowAnys(idx)), nrow(idx),
                                sum(matrixStats::colAnys(idx)), ncol(idx))
   object
}

#' @rdname left_censor
#' @importFrom magrittr %>%
#' @export
right_censor <- function(object, threshold){
   idx <- autonomics.import::exprs(object) %>% autonomics.support::is_greater_than(threshold)
   autonomics.import::exprs(object)[idx] <- threshold
   autonomics.support::cmessage('\t\tRight censor exprs > %s in %d/%d features of %d/%d samples',
                                as.character(threshold),
                                sum(matrixStats::rowAnys(idx)), nrow(idx),
                                sum(matrixStats::colAnys(idx)), ncol(idx))
   object
}
