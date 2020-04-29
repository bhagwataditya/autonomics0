#' Is even/odd?
#' @param x integer
#' @return TRUE or FALSE
#' @examples
#' is_even(13)
#' is_even(12)
#' is_odd(13)
#' is_odd(12)
#' @importFrom magrittr %>% 
#' @export
is_even <- function(x)   (x %% 2) == 0


#' @rdname is_even
#' @export
is_odd <- function(x)    !is_even(x)


#' Has even/odd length?
#' @param x vector
#' @return logical
#' @examples
#' has_even_length(1:2)
#' has_odd_length(1:2)
#' has_even_length(1:3)
#' has_odd_length(1:3)
#' @export
has_even_length <- function(x)   is_even(length(x))


#' @rdname has_even_length
#' @export
has_odd_length <- function(x)   is_odd(length(x))


#' Evenify upwards
#' @param x integer
#' @return integer
#' @examples 
#' evenify_upwards(3)
#' @export
evenify_upwards <- function(x)   if (is_odd(x)) x+1 else x
