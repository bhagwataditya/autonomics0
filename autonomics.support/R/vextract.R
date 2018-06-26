#' Vectorized extract
#' 
#' Extract vector from list_of_vectors
#' 
#' @param  list_of_vectors  list of n vectors
#' @param  extractor_list   list of n integer or character vector
#' @param  simplify         logical: whether to simplify list output
#' @return list of n vectors
#' @examples 
#' require(magrittr)
#' list_of_vectors <- list(c(L="STD", M="EM15", H="EM30"), 
#'                         c(L="STD", M="EM15", H="EM30"), 
#'                         c(L="STD", M="EM15", H="EM30" ))
#' extractor_list <- list(c('H', 'L'), c('H', 'L'), c('H', 'L'))
#' list_of_vectors %>% autonomics.support::vextract(extractor_list)
#' list_of_vectors %>% autonomics.support::vextract(extractor_list, simplify = TRUE)
#' @export
vextract <- function(
  list_of_vectors, 
  extractor_list, 
  simplify = FALSE
){
  extractor_list %<>% purrr::transpose() %>% lapply(unlist)
  result <- lapply(extractor_list, function(extractor_vector){
                           mapply(function(x, extractor) x[[extractor]], 
                                  list_of_vectors, extractor_vector)
                         }
            )
  result %<>% purrr::transpose() %>% lapply(unlist)
  if (simplify)  result %<>% do.call(rbind, .) %>% as.data.frame()
  result
}