#' NA-accepting is_weakly_greater_then
#' @param a numeric scalar, vector, or matrix
#' @param b numeric scalar, vector, or matrix
#' @return logical scalar, vector, or matrix
#' @importFrom magrittr %>%
#' @export
is_weakly_greater_than <- function(a,b){
   (a >= b) %>%
   (function(x){ x[is.na(x)] <- FALSE; x})
}
