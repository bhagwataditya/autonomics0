#' Get unique unique tails
#' @param x character vector
#' @examples
#' x <- c("E_1", "E_2", "E_3")
#' x %>% get_unique_tail()
#' @importFrom magrittr %>%
#' @export
get_unique_tails <- function(x){
   for (k in 0:max(nchar(x))){
      xtail <- sample_ids %>% substr(max(1, nchar(.)-k), nchar(.))
      if (!any(duplicated(xtail))) return(xtail)
   }
}
