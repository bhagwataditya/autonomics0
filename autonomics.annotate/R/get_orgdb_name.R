#' Get orgdb name
#' @param organism any value in \code{names(ANNOTATED_ORGANISMS)}
#' @examples
#' get_orgdb_name('Homo sapiens')
#' @importFrom magrittr %>%
#' @export
get_orgdb_name <- function(organism){
   assertive.sets::assert_is_subset(organism, names(ANNOTATED_ORGANISMS))
   organism <- ANNOTATED_ORGANISMS %>% magrittr::extract2(organism)
   orgdb_name <- sprintf('org.%s.%s.db', organism['abrev'], organism['source'])
   return(orgdb_name)
}
