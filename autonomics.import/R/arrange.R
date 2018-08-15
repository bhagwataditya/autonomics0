#' @rdname arrange_features
#' @importFrom magrittr   %>%    %<>%
#' @export
arrange_features_ <- function(x, fvars){
   idx <- do.call(order, autonomics.import::fdata(x)[, fvars, drop = FALSE])
   x %<>% magrittr::extract(idx, )
   return(x)
}

#' Arrange features of eset
#'
#' Arrange eset on specified fvars
#' @param x       eset
#' @param fvars   fvars (SE)
#' @param ...     fvars (NSE)
#' @return arranged eset
#' @examples
#' if(require(autonomics.data)){
#'    head(autonomics.import::fdata(ALL)[, -2])
#'    head(autonomics.import::fdata(arrange_features(ALL, gene_symbols))[, -2])
#'    head(autonomics.import::fdata(arrange_features_(ALL, "gene_symbols"))[, -2])
#' }
#' @importFrom magrittr   %>%
#' @export
arrange_features <- function(x, ...){
   lazy_list <- lazyeval::lazy_dots(...)
   fvars  <- lapply(lazy_list, getElement, 'expr') %>% as.character()
   arrange_features_(x, fvars)
}

#' @rdname arrange_samples
#' @importFrom magrittr   %<>%
#' @export
arrange_samples_ <- function(x, svars){
   idx <- do.call(order, autonomics.import::sdata(x)[, svars, drop = FALSE])
   x %<>% magrittr::extract(, idx)
   return(x)
}

#' Arrange samples of eset
#'
#' Arrange eset on specified svars
#' @param x       eset
#' @param svars   svars (SE)
#' @param ...     svars (NSE)
#' @return arranged eset
#' @examples
#' if(require(autonomics.data)){
#'    head(autonomics.import::sdata(ALL)[, -2])
#'    head(autonomics.import::sdata(arrange_samples(ALL, age)))
#'    head(autonomics.import::sdata(arrange_samples_(ALL, "age")))
#' }
#' @importFrom magrittr   %>%
#' @export
arrange_samples <- function(x, ...){
   lazy_list <- lazyeval::lazy_dots(...)
   svars  <- lapply(lazy_list, getElement, 'expr') %>% as.character()
   arrange_samples_(x, svars)
}

#' Order features in eSet on magnitude of given contrast
#' @param object eset
#' @param contrast_name contrast_name
#' @importFrom magrittr extract %>%
#' @export
arrange_on_magnitude <- function(object, contrast_name){
   idx <- autonomics.import::fdata(object) %>%
      magrittr::extract2(sprintf('value.%s', contrast_name)) %>%
      abs(.) %>% order(decreasing = TRUE)
   object %>% magrittr::extract(idx, )
}
