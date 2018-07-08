#' Vectorized sprintf
#' @param fmt format string
#' @param ... arguments
#' @return character vector
#' @examples
#' vsprintf(fmt = '%s.%s', letters[1:3], 1:3)
#' vsprintf(fmt = '%s.%s', letters[1:3], 1:3, first_slowest = FALSE)
#' @importFrom magrittr %>%
#' @export
vsprintf <- function(fmt, ..., first_slowest = TRUE){
   arguments <- list(...)
   n <- arguments %>% vapply(length, integer(1)) %>% Reduce(max, .)
   arguments %<>% lapply(function(x) if (length(x)==1) rep(x,n) else x)
   assertive.base::assert_all_are_true(vapply(arguments, length, integer(1)) == length(arguments[[1]]))
   argdf <- arguments %>%
           (function(x) if (first_slowest) rev(x) else x) %>%
            expand.grid() %>%
           (function(x) if (first_slowest) x %>% magrittr::extract(, ncol(.):1) else x)
   fun <- function(...) sprintf(fmt, ...)
   do.call(fun, argdf)
}
