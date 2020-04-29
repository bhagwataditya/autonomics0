
#' Extract common substring
#' @param a first string
#' @param b second string
#' @return  string
#' @importFrom magrittr %>% 
#' @export
#' @examples 
#' # Sequences
#'   a <- "heart-specific Fatty acid binding protein"
#'   b <- "Fatty acid binding protein, isoform 3"
#'   extract_common_substr(a, b)
#'   
#'   a <- "Small nuclear ribonucleoprotein-associated proteins B and B'" 
#'   b <- "Small nuclear ribonucleoprotein-associated protein N"
#'   extract_common_substr(a, b)
#' @references https://stackoverflow.com/questions/28261825/longest-common-substring-in-r-finding-non-contiguous-matches-between-the-two-str
extract_common_substr <- function(a,b){
    
    tt <- drop(attr(utils::adist(a,b, counts=TRUE), "trafos"))
    
    # Nothing in common
    if (!stringi::stri_detect_regex(tt, 'M+')) return('')
    
    # Something in common
    aa <- tt %>% 
        stringi::stri_sub(stringi::stri_locate_all_regex(tt, '[DM]+')[[1]]) %>% 
        paste0(collapse = '') %>% trimws() 
            # paste is required because multiple substrings can be found
    #bb <- tt %>% 
    # stringi::stri_sub(stringi::stri_locate_all_regex(tt, '[IM]+')[[1]]) %>% 
    # paste0(collapse = '')
    
    a %>% 
    stringi::stri_sub(stringi::stri_locate_all_regex(aa, 'M+')[[1]]) %>% 
    paste0(collapse = '') %>% 
    trimws()
    
    # different  = c(a %>% 
    #    stringi::stri_sub(stringi::stri_locate_all_regex(aa, 'D+')[[1]]) %>% 
    #    trimws(),
    # b %>% 
    #    stringi::stri_sub(stringi::stri_locate_all_regex(bb, 'I+')[[1]])) %>% 
    #    trimws())
    
}

#' Commonify strings
#' @param x character vector
#' @examples
#' # NO DIFFERENCES
#' require(magrittr)
#'    x <- c( 'Retrotransposon Gag-like protein 8B',
#'            'Retrotransposon Gag-like protein 8B')
#'    x %>% commonify_strings()
#' # TAILS DIFFER
#'    x <- c( 'Histone H2B type 1-K',
#'            'Histone H2B type 1-C/E/F/G/I')
#'    x %>% commonify_strings()
#'    x <- c("Small nuclear ribonucleoprotein-associated proteins B and B'",
#'           "Small nuclear ribonucleoprotein-associated protein N")
#'    x %>% commonify_strings()
#' # MORE COMPLEX DIFFERENCES
#'    x <- c( 'Fatty acid binding protein, isoform 3',
#'            'Fatty acid binding protein', 
#'            'heart-specific Fatty acid binding protein',
#'            'heart-specific Fatty acid binding protein, isoform 3')
#'    x %>% commonify_strings()
#' # NOTHING IN COMMON
#'    x <- c('ABC1', 'DEF2')
#'    x %>% commonify_strings()
#' @importFrom magrittr %>% 
#' @export
commonify_strings <- function(x){
    
    common <- Reduce(extract_common_substr, x)
    alternate <- if (common==''){ x
    } else {         x %>% stringi::stri_replace_first_fixed(common, '')  %>% 
            stringi::stri_replace_first_fixed(', ', '')    %>% 
            trimws()
    }
    
    if (all(alternate == '')) return(common)
    
    alternate %>% unique() %>% 
        (function(s){s[s==''] <- '.'; s}) %>%
        sort() %>% 
        #magrittr::extract(.!='') %>% 
        paste0(collapse=' | ') %>% 
        paste0('( ', ., ' )')    %>% 
        paste0(common, ' ', .)
}
