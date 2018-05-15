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
      magrittr::extract2(sprintf('coef.%s', contrast_name)) %>%
      abs(.) %>% order(decreasing = TRUE)
   object %>% magrittr::extract(idx, )
}




#' Filter features on condition
#' @param x eSet
#' @param condition filter condition
#' @return filtered eSet
#' @examples
#' if (require(autonomics.data)){
#'    filter_features( ALL,  gene_symbols %in% c('LIG4', 'MAPK12', 'MAPK1') )
#'    filter_features_(ALL, "gene_symbols %in% c('LIG4', 'MAPK12', 'MAPK1')")
#' }
#' @export
filter_features_ <- function(x, condition){
  idx <- lazyeval::lazy_eval(condition, autonomics.import::fdata(x))
  idx <- idx & !is.na(idx)
  # sprintf('\t%d/%d features pass filter', sum(idx), length(idx)) %>% message()
  x[idx, ]
}

#' @rdname filter_features_
#' @export
filter_features <- function(x, condition){
  filter_features_(x, lazyeval::lazy(condition))
}

#' Identify reference features
#' @param object   eset
#' @param fvar       fvar
#' @param split      split
#' @param fvalues    fvalues
#' @examples
#' if (require(billing.differentiation.data)){
#'    require(magrittr)
#'    fvalues <- c("A0MZ66", "A0JLT2", "Q8NG66", "Q9NPI5")
#'    object <- billing.differentiation.data::protein.ratios
#'    object %>% is_fvalue_feature(fvar = 'Uniprot accessions',
#'                                   split = ';', fvalues = fvalues) %>%
#'                 sum()
#' }
#' @importFrom magrittr  %>%
#' @export
is_fvalue_feature <- function(object, fvar, split, fvalues){
  autonomics.import::fdata(object) %>%
    magrittr::extract2(fvar) %>%
    as.character() %>%
    strsplit(split) %>%
    vapply(function(x) any(x %in% fvalues), logical(1))
}


#' Filter subet for reference features
#' @param  object            eSet
#' @param  fvar                fvar
#' @param  split               split
#' @param  fvalues             fvalues
#' @examples
#' if (require(billing.differentiation.data)){
#'    require(magrittr)
#'    fvalues <- c("A0MZ66", "A0JLT2", "Q8NG66", "Q9NPI5")
#'    object <- billing.differentiation.data::protein.ratios
#'    object %>% filter_features_on_fvalues(fvar    = 'Uniprot accessions',
#'                                            split   = ';',
#'                                            fvalues = fvalues)
#' }
#' @return subset of object with reference features
#' @importFrom magrittr       %>%     %<>%
#' @export
filter_features_on_fvalues <- function(object, fvar, split, fvalues){
  fvalues %<>% unique()
  idx <- object %>% is_fvalue_feature(fvar, split, fvalues)
  object %>% magrittr::extract(idx, )
}


utils::globalVariables('.')
#'Select portion of eSet with only up (or down) features
#' @param object eSet
#' @param contrast contrast (named string)
#' @param direction 'neg' or 'pos'
#' @importFrom magrittr %<>%
#' @export
filter_features_on_direction <- function(object, contrast, direction){
  if (direction == 'pos'){
    object %<>% magrittr::extract(which(autonomics.import::fdata(.)[[paste0('coef.', names(contrast))]] > 0))
  } else if (direction == 'neg'){
    object %<>% magrittr::extract(which(autonomics.import::fdata(.)[[paste0('coef.', names(contrast))]] < 0))
  }
  return(object)
}

#' Filter features with min expr for all samples
#' @param object  exprs object
#' @param min_expr  minimum expression
#' @importFrom magrittr %>%
#' @export
filter_features_min_expr <- function(object, min_expr){
   selector <- matrixStats::rowAnys(autonomics.import::exprs(object) > min_expr)
   autonomics.support::cmessage('\t\tRetain %d/%d features: expr > %0.2f for all samples', sum(selector), length(selector), min_expr)
   object %>% magrittr::extract(selector, )
}


#' Filter samples on condition
#' @param x eSet
#' @param condition filter condition
#' @return filtered eSet
#' @examples
#' if (require(autonomics.data)){
#'    filter_samples( ALL,  sex == 'M')
#'    filter_samples_(ALL, "sex == 'M'")
#' }
#' @importFrom magrittr %<>%
#' @export
filter_samples_ <- function(x, condition){
  idx <- lazyeval::lazy_eval(condition, autonomics.import::sdata(x))
  idx <- idx & !is.na(idx)
  # sprintf('\t%d/%d samples pass filter', sum(idx), length(idx)) %>% message()
  x %<>% magrittr::extract(, idx)
  autonomics.import::sdata(x) %<>% droplevels()
  x
}

#' @rdname filter_samples_
#' @export
filter_samples <- function(x, condition){
  filter_samples_(x, lazyeval::lazy(condition))
}

