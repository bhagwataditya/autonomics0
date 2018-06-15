#' Pull columns in a dataframe to the front
#' @param df dataframe
#' @param first_cols columns that need to be pulled to the front
#' @return dataframe with re-ordered columns
#' @examples 
#' require(magrittr)
#' df <- data.frame(
#'    symbol = c('A1BG', 'A2M'), 
#'    id     = c('1',    '2'),
#'    name   = c('alpha-1-B glycoprotein', 'alpha-2-macroglobulin'), 
#'    type   = c('proteincoding', 'proteincoding'))
#' first_cols <- c('id', 'symbol') 
#' df %>% autonomics.support::pull_columns(first_cols)
#' @importFrom magrittr %>% 
#' @export
pull_columns <- function(df, first_cols){
  
  assertive.types::assert_is_data.frame(df)
  assertive.types::assert_is_character(first_cols)
  
  df %>% 
  magrittr::extract(, c(first_cols, setdiff(names(df), first_cols)), drop = FALSE)
}