#' Get unique tails
#' @param x character vector
#' @return  character vector
#' @examples
#' x <- c("E_1", "E_2", "E_3")
#' get_unique_tails(x)
#' @importFrom magrittr %>%
#' @export
get_unique_tails <- function(x){
   if (is.factor(x)) x %<>% as.character()
   for (k in 0:max(nchar(x))){
      xtail <- x %>% substr(max(1, nchar(.)-k), nchar(.))
      if (!any(duplicated(xtail))) return(xtail)
   }
}
