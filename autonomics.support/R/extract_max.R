#' Convenient equals operator
#' 
#' Performs x == y, but returns FALSE rather than NA for NA elements of x.
#' @param x numeric vector or scalar
#' @param y numeric scalar
#' @examples
#' require(magrittr)
#' x <- c(A=1,B=3,C=2,D=3, E=NA)
#' y <- 3
#' x %>% equals(y)
#' @export
na_aware_equals <- function(x,y){
  result <- rep(FALSE, length(x)) %>% 
    magrittr::set_names(names(x))
  if (is.na(y)){
    result[ is.na(x)] <- TRUE
    result[!is.na(x)] <- FALSE
  } else {
    result[ is.na(x)] <- FALSE
    result[!is.na(x)] <- x[!is.na(x)] == y
  }
  result
}

#' Is maximal
#' @param x numeric vector
#' @examples 
#' require(magrittr)
#' x <- c(A=1,B=3,C=2,D=3, E=NA)
#' x %>% is_max()
#' @export
is_max <- function(x){
  na_aware_equals(x, max(x, na.rm = TRUE))
}

