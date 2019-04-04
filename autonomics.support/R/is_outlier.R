#' Identify outliers in vector
#' @param x numeric
#' @return logical
#' @examples 
#' x <- c(1:25,  200)
#' autonomics.support::is_outlier(x)
#' @importFrom magrittr %>% 
#' @export
is_outlier <- function(x){
   med <- x %>% stats::median(na.rm=TRUE)
   iqr <- x %>% stats::IQR(na.rm = TRUE)
   lower <- med - 1.5*iqr
   upper <- med + 1.5*iqr
   (x < lower) | (x > upper)
}
