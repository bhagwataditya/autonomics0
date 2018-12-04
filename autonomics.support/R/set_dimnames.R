#' Set dimnames of a matrix
#' @param x matrix
#' @param dim_names list
#' @examples
#' require(magrittr)
#' x <- matrix(1:4, nrow=2)
#' dim_names <- list(feature = c('F1', 'F2'), sample = c('S1', 'S2'))
#' x %>% autonomics.support::set_dimnames(dim_names)
#' @return matrix with set dimnames
#' @export
set_dimnames <- function(x, dim_names){
   dimnames(x) <- dim_names
   x
}
