#' Uniquify collapsed strings
#' @param x vector of collapsed strings
#' @param sep string separator
#' @examples 
#' require(magrittr)
#' a <- paste0("GO:0000786; GO:0002227; GO:0003677; GO:0005615; ",
#'             "GO:0005634; GO:0005654; GO:0005829; GO:0006334; ", 
#'             "GO:0016567; GO:0019731; GO:0031640; GO:0046982; ",
#'             "GO:0050829; GO:0050830; GO:0061844")
#' b <- paste0("GO:0000786; GO:0002227; GO:0003677; GO:0005615; ",
#'             "GO:0005634; GO:0005654; GO:0005829; GO:0006334; ",
#'             "GO:0016567; GO:0019731; GO:0042802; GO:0046982; ",
#'             "GO:0050830; GO:0061844; GO:0070062")
#' x <- c(a, b)
#' x %>% uniquify_collapsed_strings(';')
#' @importFrom magrittr %>% 
#' @export
uniquify_collapsed_strings <- function(x, sep){
  x %>% strsplit(sep) %>% unlist() %>% unique() %>% paste0(collapse = sep) 
}