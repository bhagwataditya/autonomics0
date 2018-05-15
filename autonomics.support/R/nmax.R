#' Return nth max (min) value in vector
#' 
#' Orders a vector and returns n'th ordered value.
#' When vector length is smaller than n, returns last value.
#' 
#' @return value
#' @param x numeric vector
#' @param n integer
#' @importFrom magrittr %>% 
#' @export
nmax <- function(x, n){
  x      %>% 
  sort(decreasing = TRUE) %>% 
  magrittr::extract(min(length(.), n))
}

#' @rdname nmax
#' @importFrom magrittr %>% 
#' @export
nmin <- function(x, n){
  x      %>% 
  sort() %>% 
  magrittr::extract(min(length(.), n))
}