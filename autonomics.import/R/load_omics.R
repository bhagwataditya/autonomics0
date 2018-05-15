
#' Validify sample/feature names
#' @param x character vector with sample ids
#' @return character vector
#' @importFrom magrittr %<>%
#' @export
validify_sample_ids <- function(x){
   . <- NULL
   selector <- substr(x,1,1) %in% 0:9
   x[selector] %<>% paste0('s', .)
   x
}


