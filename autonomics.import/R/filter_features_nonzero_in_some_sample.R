#' Keep features that are non-zero (and non-NA) for some sample
#' @param   object  eSet
#' @return  filtered eSet
#' @importFrom  magrittr %>%
#' @export
filter_features_nonzero_in_some_sample <- function(object){
   # . != 0 needed due to stupid behaviour of rowAnys
   # https://github.com/HenrikBengtsson/matrixStats/issues/89
   selector <- matrixStats::rowAnys(autonomics.import::exprs(object) != 0, na.rm = TRUE)
   autonomics.support::cmessage('\t\tRetain %d/%d features that are non-zero (and non-NA) for some sample',
                                sum(selector), length(selector))
   object %>% magrittr::extract(selector, )
}
