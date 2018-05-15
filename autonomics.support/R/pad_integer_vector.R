#' Pad integer vector
#' @param x integer vector
#' @param pad padding character
#' @return character vector
#' @importFrom magrittr %>%
#' @export
pad_integer_vector <- function(x, pad = '0'){
  assertive.base::assert_all_are_not_na(x)
  # Avoid breakdown when x has class factor (but still contains integers)
  x %<>% as.character() %>% as.integer()
  n_digits <- max(x) %>% log10 %>% ceiling()
  formatC(x, width = n_digits, format = "d", flag = pad) 
}