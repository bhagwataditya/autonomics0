#' Characterify a dataframe
#' 
#' Transform numeric variables into character variables in scientific notation
#' Limit length of character variables
#' 
#' @param x             dataframe
#' @param char_length   allowed length of strings
#' @return compactified dataframe
#' @importFrom magrittr          %>%
#' @export
characterify <- function(x, char_length = 40){
  
  assertive.types::assert_is_data.frame(x)
  assertive.types::assert_is_numeric(char_length)
  
  # limit length of character variables
  selector <- x %>% vapply(is.character, logical(1))
  x[selector] %<>% apply(2, substring, 1, char_length)
  
  # numeric -> characters in scientific notation
  selector <- x %>% vapply(is.numeric, logical(1))
  x[selector] %<>% signif(2) %>% apply(2, function(x){sprintf('%1.2e', x)})
  x
}
