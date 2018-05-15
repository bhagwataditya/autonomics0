#' All elements of vector are identical
#' @param x vector
#' @return logical(1)
#' @examples
#' require(magrittr)
#' x <- c(2,2,1,2)
#' x %>% has_identical_values()
#' @export
has_identical_values <- function(x){
  length(unique(x))==1
}

#' All have value
#' @param x vector
#' @param y scalar
#' @examples 
#' require(magrittr)
#' c(1,1,1,1) %>% all_have_value(1)
#' @return logical(1)
#' @export
all_have_value <- function(x, y){
  all(x==y)
}