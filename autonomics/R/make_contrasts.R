utils::globalVariables('.')

#' Generate "A <-> B" contrasts
#' @examples 
#' make_contrasts_A.min.B(c('A', 'B', 'C'))
#' @importFrom magrittr   %<>%   %>%
#' @noRd
make_contrasts_A.min.B <- function(subgroup_levels){
   utils::combn(subgroup_levels, 2)                           %>%
   apply(2, paste, collapse = ' - ')                   %>%
   magrittr::set_names(stringr::str_replace(., ' - ', '__'))
   
}

#' Generate contrasts of subgroup levels
#' 
#' @param object ProtSet object
#' @param contrast_types types of contrasts to be generated
#' @importFrom magrittr %<>%
#' @export
make_contrasts <- function(object, contrast_types = c('A-B')){
  
   contrast_types %<>% stringr::str_replace(' ', '')
   subgroup_levels <- as.character(unique(object$subgroup))
   nlevels <- length(subgroup_levels)
   contrasts <- c()
   
   if (length(contrast_types)>0){
     if (nlevels > 1 & 'A-B' %in% contrast_types){
       contrasts %<>% c(make_contrasts_A.min.B(subgroup_levels))
     }
   }
   return(contrasts)
}
