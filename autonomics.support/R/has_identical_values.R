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
