
#' Is missing or empty character
#' @param x character vector
#' @return logical vector
#' @export
is_missing_or_empty_character <- function(x){
  x[is.na(x)] <- ''
  x == ''
}


#' Are all missing or empty character
#' @param x character vector
#' @return logical(1)
#' @examples
#' require(magrittr)
#' x <- c('a', 'b', NA_character_, 'd', '')
#' x %>% autonomics.support::all_are_missing_or_empty_character()
#' @importFrom magrittr %>%
#' @export
all_are_missing_or_empty_character <- function(x){
  x %>% is_missing_or_empty_character() %>% all()
}



#' Are all svalues available?
#' @param object SummarizedExperiment
#' @param svar   character(1)
#' @return logical(1)
#' @importFrom magrittr %>%
#' @export
all_svalues_available <- function(object, svar){
   svalues1 <- object %>% autonomics.import::svalues(svar)
   if (is.null(svalues1))                            return(FALSE)
   if (all_are_missing_or_empty_character(svalues1)) return(FALSE)
   TRUE
}
