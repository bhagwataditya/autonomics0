#' Percentage of available values
#' @param x vector
#' @return rounded percentage
#' @examples
#' x <- c(NA, 2, 3, NA, 5, 6)
#' pct_available(x)
#' @importFrom magrittr %>%
#' @export
pct_available <- function(x){
  selector <- !is.na(x)
  pct <- 100*sum(selector)/length(selector)
  pct %>% floor()
}

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
  pct %>% floor()
}