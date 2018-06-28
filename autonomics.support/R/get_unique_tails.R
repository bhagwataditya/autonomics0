#' Get unique unique tails
#' @param x character vector
#' @examples
#' x <- c("E_1", "E_2", "E_3")
#' x %>% get_unique_tail()
#' @importFrom magrittr %>%
#' @export
get_unique_tails <- function(x){
  if (is.factor(x)) x %<>% as.character()
   for (k in 0:max(nchar(x))){
      xtail <- x %>% substr(max(1, nchar(.)-k), nchar(.))
      if (!any(duplicated(xtail))) return(xtail)
   }
}
