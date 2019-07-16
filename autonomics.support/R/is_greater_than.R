#' NA-accepting versions of comparator operators
#' @param a numeric scalar, vector, or matrix
#' @param b numeric scalar, vector, or matrix
#' @return logical scalar, vector, or matrix
#' @importFrom magrittr %>%
#' @export
is_greater_than <- function(a,b){
  (a > b) %>%
    (function(x){ x[is.na(x)] <- FALSE; x})
}

