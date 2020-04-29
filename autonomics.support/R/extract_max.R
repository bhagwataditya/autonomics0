#' Convenient equals operator
#' 
#' Performs x == y, but returns FALSE rather than NA for NA elements of x.
#' @param x numeric vector or scalar
#' @param y numeric scalar
#' @return logical vector
#' @examples
#' x <- c(A=1,B=3,C=2,D=3, E=NA)
#' y <- 3
#' equals(x, y)
#' @export
cequals <- function(x,y){
    result <- rep(FALSE, length(x)) %>% magrittr::set_names(names(x))
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
#' @return logical vector
#' @examples 
#' x <- c(A=1,B=3,C=2,D=3, E=NA)
#' is_max(x)
#' @export
is_max <- function(x){
    cequals(x, max(x, na.rm = TRUE))
}


#' Extract portion of vector that contains maximum values
#' @param x numeric vector
#' @return subvector 
#' @examples
#' x <- c(A=1,B=3,C=2,D=3, E=NA)
#' extract_max(x)
#' @importFrom magrittr %>% 
#' @export
extract_max <- function(x){
    assertive.types::assert_is_numeric(x)
    x %>% magrittr::extract(is_max(.))
}