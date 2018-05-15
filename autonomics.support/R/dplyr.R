#' Left join while preserving rownames
#' 
#' Convenient left_join function that preserves row names
#' and avoids warning messages when by values are of different type.
#' 
#' @param x first object
#' @param y second object
#' @param by by
#' @param ... passed to dplyr::left_join
#' @return merged dataframe
#' @seealso \code{\link[dplyr]{left_join}}
#' @importFrom magrittr %>% 
#' @export
left_join_keeping_rownames <- function(x, y, by, ...){
  
  # Avoid that by values are of different type
  # This throws a warning in dplyr::left_join
  x[[by]] %<>% as.character()
  y[[by]] %<>% as.character()
  
  # Left join  
  dplyr::left_join(x,y,by,...) %>% 
    magrittr::set_rownames(rownames(x))
}

#' Mutate while preserving rownames
#' @param .data dataframe
#' @param ... passed to dplyr::mutate
#' @return data.frame
#' @importFrom magrittr  %>% 
#' @export
mutate_keeping_rownames <- function(.data, ...){
  dplyr::mutate(.data, ...) %>% 
    magrittr::set_rownames(rownames(.data))
}