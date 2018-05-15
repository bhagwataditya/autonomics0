#' Is maximal
#' @param x numeric vector
#' @examples 
#' require(magrittr)
#' x <- c(A=1,B=3,C=2,D=3, E=NA)
#' x %>% is_max()
#' @export
is_max <- function(x){
  autonomics.support::equals(x, max(x, na.rm = TRUE))
}

#' Extract portion of vector that contains maximum values
#' @param x numeric vector
#' @examples
#' require(magrittr) 
#' x <- c(A=1,B=3,C=2,D=3, E=NA)
#' x %>% autonomics.support::extract_max()
#' @importFrom magrittr %>% 
#' @export
extract_max <- function(x){
  assertive.types::assert_is_numeric(x)
  x %>% magrittr::extract(autonomics.support::is_max(.))
}