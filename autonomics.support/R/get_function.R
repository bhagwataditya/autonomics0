get_function <- function(x){
  x %>% 
  stringi::stri_split_fixed('::') %>% 
  unlist() %>% 
  (function(y) utils::getFromNamespace(y[2],y[1]))
}