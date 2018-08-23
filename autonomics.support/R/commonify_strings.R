x <- "heartspecific HLA class I histocompatibility antigen, A-68 alpha chain"
y <- "HLA class I/II histocompatibility antigen, A-2 alpha chain (isoform 3)"

#' @rdname commonify_strings
#' @importFrom magrittr %>% 
#' @export
commonify_two_strings <- function(x,y){
  ta <- drop(attr(utils::adist(x,y, counts=TRUE), "trafos"))
  tt <- ta %>% strsplit('') %>% unlist()
  
  x
  tx <- ta %>% stringi::stri_replace_all_fixed('I', '') %>% strsplit('') %>% unlist()
  
  y
  ty <-  ta %>% stringi::stri_replace_all_fixed('D', '') %>% strsplit('') %>% unlist()
  
  x2 <- x %>% strsplit('') %>% unlist() %>%  (function(tmp){ tmp[tx=='D'] <- '.'; tmp}) %>% paste0(collapse='')
  y2 <- y %>% strsplit('') %>% unlist() %>%  (function(tmp){ tmp[ty=='I'] <- '.'; tmp}) %>% paste0(collapse='')
  
  t2 <- rep('', length(tt))
  t2[tt=='D'] <- unlist(strsplit(x, ''))[tx=='D']                      # add "()"
  t2[tt=='I'] <- unlist(strsplit(y, ''))[ty=='I']     # add "()" before|after chunk
  t2[tt=='S'] <- sprintf('(%s|%s)', unlist(strsplit(x, ''))[tx=='S'],
                                    unlist(strsplit(y, ''))[ty=='S'])
  t2[tt=='M'] <- unlist(strsplit(y, ''))[ty=='M']
  t2 %>% paste0(collapse = '')
  
  
  
  
  tt <- ta %>% strsplit('') %>% unlist()
  xx <- x  %>% strsplit('') %>% unlist()
  yy <- y  %>% strsplit('') %>% unlist()
  
  ss <- rep('', length(tt))
  ss[tt=='D'] <- 
  
  
  s <- ''
  if (stringi::stri_detect_regex(ta, 'M+')){
       s <- x %>% stringi::stri_sub(stringi::stri_locate_all_regex(ta, 'M+')[[1]])
  }
  if (stringi::stri_detect_regex(ta, '[ISD]+')){
     s <- sprintf('%s%s|%s', 
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
