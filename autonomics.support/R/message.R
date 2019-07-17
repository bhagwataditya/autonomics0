#' Conveniently message dataframe
#' 
#' Conveniently message dataframe using sprintf syntax.
#' Use place holder '%s' for data.frame.
#' 
#' @param format_string sprintf style format string
#' @param x data.frame
#' @examples
#' require(magrittr)
#' x <- data.frame(feature_id = c('F001', 'F002'), symbol = c('FEAT1', 'FEAT2'))
#' autonomics.support::cmessage_df('\t%s', x)
#' 
#' x <- c(rep('PASS', 25), rep('FAIL', 25))
#' autonomics.support::cmessage_df(format_string = '%s', table(x))
#' @importFrom magrittr %>% 
#' @export
cmessage_df <- function(format_string, x){
  format_string %>% 
    sprintf(utils::capture.output(print(x))) %>% 
    paste0(collapse = '\n') %>% 
    enc2utf8() %>% 
    message()
}


#' Conveniently message 
#' 
#' Print message to screen with sprintf syntax
#' 
#' @param format_string sprintf format string
#' @param ... additional arguments  passed to sprintf
#' @examples 
#' cmessage('\t%s\t%s', 'Hi', 'there')
#' @importFrom magrittr %>%
#' @export
cmessage <- function(format_string, ...){
  format_string %>% 
    sprintf(...)%>% 
    message()
}

