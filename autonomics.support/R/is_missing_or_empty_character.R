
#' @rdname all_are_missing_or_empty_character
#' @export
is_missing_or_empty_character <- function(x){
    x[is.na(x)] <- ''
    x == ''
}


#' Are all missing or empty character
#' @param x character vector
#' @return TRUE or FALSE
#' @examples
#' x <- c('a', 'b', NA_character_, 'd', '')
#' all_are_missing_or_empty_character(x)
#' @importFrom magrittr %>%
#' @export
all_are_missing_or_empty_character <- function(x){
    x %>% is_missing_or_empty_character() %>% all()
}

