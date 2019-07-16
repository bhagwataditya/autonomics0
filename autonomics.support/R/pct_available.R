

#' Percentage of TRUE value
#' @param x logical vector
#' @export
#' @return rounded percentage 
#' @examples
#' x <- c(TRUE, FALSE, TRUE, TRUE)
#' pct_true(x)
#' @importFrom magrittr %>%
pct_true <- function(x){
  assertive.types::assert_is_logical(x)
  pct <- 100*sum(x)/length(x)
  floor(pct)
}