#' DIRECTIONS
#'
#' Allowed values for diresction argument
#' @export
DIRECTIONS <- c('neg', 'pos', 'both')

#' Convert direction into sign
#' @param direction any value in autonomics.find::DIRECTIONS
#' @return sign (character)
#' @export
direction_to_sign <- function(direction){
   switch(direction, 
          neg = '<', 
          pos = '>', 
          both = '<=>')
}