#' Return nth max (min) value in vector
#' 
#' Orders a vector and returns n'th ordered value.
#' When vector length is smaller than n, returns last value.
#' 
#' @param x numeric vector
#' @param n integer
#' @return value
#' @examples 
#' nmax(c(1,2,3,4,5), 2)
#' nmin(c(1,2,3,4,5), 2)
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