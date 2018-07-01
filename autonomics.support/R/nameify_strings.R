#' Nameify strings
#' 
#' Update character vector values so that they can serve as R variable names
#' @param x character vector
#' @param verbose character
#' @return validified subgroups
#' @examples
#' require(magrittr)
#' x <- c('266S', '270A', 'DM')
#' x %>% autonomics.support::nameify_strings()
#' require(magrittr)
#' @importFrom magrittr   %>%    %<>%
#' @export
nameify_strings <- function(x, verbose = TRUE){
  
  # Satisfy CHECK
  . <- NULL
  
  old_values <- x
  new_values <- old_values %>% make.names()
  
  old_values %<>% setdiff(new_values)
  new_values %<>% setdiff(old_values)
  
  if (length(old_values) > 0){
    for (i in seq_along(old_values)){
      x   %<>% stringi::stri_replace_first_fixed(old_values[i], new_values[i])
      if (verbose) autonomics.support::cmessage('\t\t%s -> %s', old_values[i], new_values[i])
    }
  }
  return(x)
  
}

