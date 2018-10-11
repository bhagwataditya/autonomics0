#' NA exprs (weakly) beyond threshold
#' @param object SummarizedExperiment
#' @param threshold numeric(1)
#' @return Updated SummarizedExperiment
#' @export
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    object <- autonomics.data::glutaminase
#'    object %>% autonomics.preprocess::na_exprs_weakly_gt(30)
#'    object %>% autonomics.preprocess::na_exprs_weakly_lt(13)
#' }
na_exprs_weakly_gt <- function(object, threshold){
   idx <- autonomics.import::exprs(object) >= threshold
   autonomics.support::cmessage('\t\tNA exprs >= %s. Done for %d/%d features in %d/%d samples',
                                as.character(threshold),
                                sum(matrixStats::rowAnys(idx, na.rm=TRUE)), nrow(idx),
                                sum(matrixStats::colAnys(idx, na.rm=TRUE)), ncol(idx))
   autonomics.import::exprs(object)[idx] <- NA
   object
}

#' @rdname na_exprs_weakly_gt
#' @export
na_exprs_weakly_lt <- function(object, threshold){
   idx <- autonomics.import::exprs(object) <= threshold
   autonomics.support::cmessage('NA exprs <= %s. Done for %d/%d features in %d/%d samples',
                                as.character(threshold),
                                sum(matrixStats::rowAnys(idx, na.rm=TRUE)), nrow(idx),
                                sum(matrixStats::colAnys(idx, na.rm=TRUE)), ncol(idx))
   autonomics.import::exprs(object)[idx] <- NA
   object
}

