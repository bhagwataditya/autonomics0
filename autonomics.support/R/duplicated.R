#' Convenient (two way) duplicated
#' @param x vector
#' @return logical vector
#' @examples 
#' require(magrittr)
#' c(1,2,3,4,5,2) %>% cduplicated()
#' @export
cduplicated <- function(x){
  duplicated(x) | duplicated(x, fromLast = TRUE)
}