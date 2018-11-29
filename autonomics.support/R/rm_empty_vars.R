
#' rm empty vars from list
#' @param x list, dataframe, etc.
#' @return list without empty vars
#' @importFrom magrittr %>%
#' @export
rm_empty_vars <- function(x){
  assertive.types::assert_is_list(x)

  for (curvar in names(x)){
    isempty <- x[[curvar]] %>% as.character() %>% assertive.strings::is_missing_or_empty_character() %>% all()
    if (isempty){
      x[[curvar]] <- NULL
    }
  }
  x
}
