
#' Extract common substring
#' @param a first string
#' @param b second string
#' @return character(1)
#' @importFrom magrittr %>% 
#' @export
#' @examples 
#' # Sequences
#'   a <- "heart-specific Fatty acid binding protein"
#'   b <- "Fatty acid binding protein, isoform 3"
#'   extract_common_substr(a, b)
#' @references https://stackoverflow.com/questions/28261825/longest-common-substring-in-r-finding-non-contiguous-matches-between-the-two-str
extract_common_substr <- function(a,b){
  
  tt <- drop(attr(utils::adist(a,b, counts=TRUE), "trafos"))
  aa <- tt %>% stringi::stri_sub(stringi::stri_locate_all_regex(tt, '[DM]+')[[1]])
  bb <- tt %>% stringi::stri_sub(stringi::stri_locate_all_regex(tt, '[IM]+')[[1]])
  
  a %>% stringi::stri_sub(stringi::stri_locate_all_regex(aa, 'M+')[[1]])
  # different  = c(a %>% stringi::stri_sub(stringi::stri_locate_all_regex(aa, 'D+')[[1]])  %>% trimws(),
  # b %>% stringi::stri_sub(stringi::stri_locate_all_regex(bb, 'I+')[[1]])) %>% trimws())
  
}

#' Commonify strings
#' @param x string
#' @param y string
#' @examples
#' require(magrittr)
#' 
#' # NO DIFFERENCES
#' x <- c('Retrotransposon Gag-like protein 8B',
#'        'Retrotransposon Gag-like protein 8B')
#'x %>% commonify_strings()
#'
#' # TAILS DIFFER
#' x <- c('Histone H2B type 1-K',
#'        'Histone H2B type 1-C/E/F/G/I')
#' x %>% commonify_strings()
#' 
#' # MORE COMPLEX DIFFERENCES
#' x <- c('Fatty acid binding protein, isoform 3',
#'        'Fatty acid binding protein', 
#'        'heart-specific Fatty acid binding protein',
#'        'heart-specific Fatty acid binding protein, isoform 3')
#' x %>% commonify_strings()
#' 
#' @importFrom magrittr %>% 
#' @export
commonify_strings <- function(x){
  
  common <- Reduce(extract_common_substr, x)
  
  alternate <- x %>% stringi::stri_replace_first_fixed(common, '')  %>% 
                     stringi::stri_replace_first_fixed(', ', '')    %>% 
                     trimws()
  
  if (all(alternate == '')) return(common)
  
  alternate %>% magrittr::extract(.!='') %>% 
                paste0(collapse=' | ') %>% 
                            #paste0('(', ., ')')    %>% 
                            paste0(common, ' : ', .)
}
