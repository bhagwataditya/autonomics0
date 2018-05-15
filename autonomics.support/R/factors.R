# note: earlier name was 'order_pres_factor'
#' Create factor with levels in order of appearance
#' 
#' Creates a factor from a vector, where the levels are in (possibly reverse)
#' order of appearance in the vector, rather than being alphabetically sorted.
#' @param avector An atomic vector.
#' @param reverse A logical value. Should factor levels be in reverse order of
#' appearance.
#' @return A factor.
#' @export
factorify <- function(avector, reverse = FALSE){
  myLevels <- unique(avector)
  if(reverse){myLevels <- rev(myLevels)}
  factor(avector, myLevels)
}

#' Make vector components unique by appending spaces
#' @param avector character or factor vector
#' @seealso \code{\link[base]{make.unique}}
#' @export
uniquify <- function(avector){
  avector <- as.character(avector)
  while(any(duplicated(avector))){
    idxDuplicates <- which(duplicated(avector))
    avector[idxDuplicates] <- paste0(avector[idxDuplicates], ' ')
  }
  return(avector)
}
