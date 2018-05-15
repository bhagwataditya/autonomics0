#' Is even/odd?
#' @param x integer
#' @return logical
#' @examples
#' is_even(13)
#' is_even(12)
#' is_odd(13)
#' is_odd(12)
#' @importFrom magrittr %>% 
#' @export
is_even <- function(x){
  x %% 2 %>% magrittr::equals(0)
}

#' @rdname is_even
#' @importFrom magrittr %>% 
#' @export
is_odd <- function(x){
  x %>% is_even() %>% magrittr::not()
}

#' Has even/odd length?
#' @param x vector
#' @return logical
#' @examples
#' has_even_length(1:2)
#' has_odd_length(1:2)
#' has_even_length(1:3)
#' has_even_length(1:3)
#' @importFrom magrittr %>% 
#' @export
has_even_length <- function(x){
  length(x) %>% is_even()
}

#' @rdname has_even_length
#' @importFrom magrittr %>% 
#' @export
has_odd_length <- function(x){
  length(x) %>% is_odd()
}

#' Evenify upwards
#' @param x integer
#' @export
evenify_upwards <- function(x){
  if (is_odd(x)) x+1 else x
}
