#' All elements of vector are identical
#' @param x vector
#' @return TRUE or FALSE
#' @examples
#' x <- c(2,2,1,2)
#' has_identical_values(x)
#' @export
has_identical_values <- function(x){
  length(unique(x))==1
}

#' All have value
#' @param x vector
#' @param y scalar
#' @return TRUE or FALSE
#' @examples 
#' require(magrittr)
#' c(1,1,1,1) %>% all_have_value(1)
#' @export
all_have_value <- function(x, y){
  all(x==y)
}