#' Convenient grid.draw
#' @param x grid object
#' @importFrom magrittr %>%
#' @export
cdraw <- function(x){
  grid::grid.newpage()
  x %>% grid::grid.draw()
}