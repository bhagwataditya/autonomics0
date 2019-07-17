


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
#' @seealso dplyr::left_join
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


#' Pull columns in a dataframe to the front
#' @param df dataframe
#' @param first_cols columns that need to be pulled to the front
#' @param verbose logical
#' @return dataframe with re-ordered columns
#' @examples 
#' require(magrittr)
#' df <- data.frame(
#'    symbol = c('A1BG', 'A2M'), 
#'    id     = c('1',    '2'),
#'    name   = c('alpha-1-B glycoprotein', 'alpha-2-macroglobulin'), 
#'    type   = c('proteincoding', 'proteincoding'))
#' first_cols <- c('id', 'symbol', 'location', 'uniprot') 
#' df %>% autonomics.support::pull_columns(first_cols)
#' @importFrom magrittr %>% 
#' @export
pull_columns <- function(df, first_cols, verbose = TRUE){
  
  assertive.types::assert_is_data.frame(df)
  assertive.types::assert_is_character(first_cols)
  
  idx <- first_cols %in% names(df)
  if (any(!idx)){
    if (verbose) autonomics.support::cmessage(
                   'pull_columns: ignore absent columns %s', 
                    first_cols[!idx] %>% sprintf("'%s'", .) %>% paste0(collapse = ', '))
    first_cols %<>% magrittr::extract(idx)
  }
  
  df %>% magrittr::extract(, c(first_cols, setdiff(names(df), first_cols)), drop = FALSE)
}



