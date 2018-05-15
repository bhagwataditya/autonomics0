#' Strsplit and extract first
#' @param x collapsed character vector
#' @param sep separator
#' @return uncollapsed character vector 
#' @importFrom magrittr %>% 
#' @export
strsplit_extract_first <- function(x, sep){
  
  x %>% as.character() %>% # factor requires conversion to character
        stringi::stri_split_fixed(sep) %>% 
        vapply(magrittr::extract, character(1), 1)
  
}