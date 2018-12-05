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


#' @importFrom magrittr %<>% 
make.unique.spaces <- function(x, verbose=TRUE){
  x %<>% as.character()
  while(any(duplicated(x))){
    idx_duplicates <- which(duplicated(x))
    x[idx_duplicates] <- paste0(x[idx_duplicates], ' ')
  }
  return(x)
}


#' Make vector components unique by appending spaces
#' @param x character or factor vector
#' @param method 'uniquify'
#' @param verbose logical
#' @seealso \code{\link[base]{make.unique}}
#' @examples 
#' require(magrittr)
#' x <- c('A', 'B', 'C', 'A', 'D')
#' x %>% autonomics.support::uniquify('make.unique')
#' x %>% autonomics.support::uniquify('make.unique.spaces')
#' @export
uniquify <- function(x, method = 'make.unique.spaces', verbose = TRUE){
  idx <- autonomics.support::cduplicated(x)
  if (any(idx)){
    uniquefun <- get(method)
    uniquex <- uniquefun(x)
    if (verbose){
      autonomics.support::cmessage(   '\t\tUniquify ( %s -> %s ) duplicates of', x[idx][1], uniquefun(x[idx][c(1, 1)])[2])
      autonomics.support::cmessage_df("\t\t\t  %s", table(x[idx]) %>% as.list() %>% unlist)                    # unlist(as.list(.)) is to prevent empty line
    }
  } else {
    uniquex <- x
  }
  uniquex
}
