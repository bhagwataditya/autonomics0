#' @importFrom magrittr %<>%
rm_trailing_slashes <- function(x){
  if (substr(x, nchar(x), nchar(x)) == '/'){
    x %<>% substr(1, nchar(x)-1)
  }
  x
}

#' convenient interface to file.path
#' 
#' removes redundant slashes before running file.path
#' @param ... sent to file.path
#' @examples 
#' cfile.path('dir/', 'subdir')
#' @importFrom magrittr %>%
#' @export
cfile.path <- function(...){
  list(...) %>%
  lapply(rm_trailing_slashes) %>%
  (function(x){do.call(file.path, x)})
}