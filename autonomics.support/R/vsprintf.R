#' Vectorized sprintf
#' @param fmt format string
#' @param ... passed to sprintf
#' @param first_slowest logical
#' @return character vector
#' @examples
#' vsprintf(fmt = '%s.%s', letters[1:3], 1:3)
#' vsprintf(fmt = '%s.%s', letters[1:3], 1:3, first_slowest = FALSE)
#' vsprintf('%s [%s]%s', 
#'          c("E(L).EM(M).BM(H).R1", "BM(L).E(M).EM(H).R2", "EM(L).BM(M).E(H).R3"), 
#'          c("H/L", "H/M", "M/L"), 
#'          '')
#' @importFrom magrittr %>%
#' @export
vsprintf <- function(fmt, ..., first_slowest = TRUE){
   arguments <- list(...)
   n <- arguments %>% vapply(length, integer(1)) %>% Reduce(max, .)
   argdf <- arguments %>%
           (function(x) if (first_slowest) rev(x) else x) %>%
            expand.grid() %>%
           (function(x) if (first_slowest) x %>% magrittr::extract(, ncol(.):1) else x)
   fun <- function(...) sprintf(fmt, ...)
   do.call(fun, argdf)
}
