#' @rdname commonify_strings
#' @importFrom magrittr %>% 
#' @export
commonify_two_strings <- function(x,y){
  ta <- drop(attr(utils::adist(x,y, counts=TRUE), "trafos"))
  s <- ''
  if (stringi::stri_detect_regex(ta, 'M+')){
     s <- x %>% stringi::stri_sub(stringi::stri_locate_all_regex(ta, 'M+')[[1]])
  }
  if (stringi::stri_detect_regex(ta, '[ISD]+')){
     s <- sprintf('%s%s | %s', 
                   s,
                   x %>% stringi::stri_sub(stringi::stri_locate_all_regex(ta, '[ISD]+')[[1]]),
                   y %>% stringi::stri_sub(stringi::stri_locate_all_regex(ta, '[ISD]+')[[1]]))
  }
  return(s)
}

#' Commonify strings from their heads
#' @param x string
#' @param y string
#' @examples
#' require(magrittr)
#' # Example 1: tail slightly different
#' x <- "Golgi apparatus membrane protein TVP23 homolog B"
#' y <- "Golgi apparatus membrane protein TVP23 homolog C"
#' z <- "Golgi apparatus membrane protein TVP23 homolog D"
#' commonify_two_strings(x,y)
#' commonify_strings(c(x,y,z))
#' 
#' # Example 2: tail different
#' x <- 'Histone H2B type 1-K'
#' y <- 'Histone H2B type 1-C/E/F/G/I'
#' commonify_two_strings(x,y)
#' 
#' # Example 3: all common
#' x <- 'Retrotransposon Gag-like protein 8B'
#' y <- 'Retrotransposon Gag-like protein 8B'
#' commonify_two_strings(x,y)
#' 
#' @references https://stackoverflow.com/questions/28261825/longest-common-substring-in-r-finding-non-contiguous-matches-between-the-two-str
#' @importFrom magrittr %>% 
#' @export
commonify_strings <- function(x){
  Reduce(commonify_two_strings, x)
}
