
#' @rdname filter_features
#' @export
filter_features_ <- function(object, condition, verbose = FALSE){
   if (is.null(condition)) return(object)
   idx <- lazyeval::lazy_eval(condition, autonomics.import::fdata(object))
   idx <- idx & !is.na(idx)
   if (verbose) message('\t\tRetain ', sum(idx), '/', length(idx), ' features: ', if (class(condition)=='lazy') deparse(condition$expr) else condition)
   object %>% autonomics.import::extract_features(idx)
}

#' Filter features on condition
#' @param object SummarizedExperiment
#' @param condition filter condition
#' @param verbose logical
#' @return filtered eSet
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    autonomics.data::ALL %>%
#'    filter_features(gene_symbols %in% c('LIG4', 'MAPK12', 'MAPK1') )
#'
#'    autonomics.data::ALL %>%
#'    filter_features(gene_symbols %in% c('LIG4', 'MAPK12', 'MAPK1'), verbose = TRUE)
#'
#'    autonomics.data::ALL %>%
#'    filter_features_("gene_symbols %in% c('LIG4', 'MAPK12', 'MAPK1')")
#'
#'    autonomics.data::ALL %>%
#'    filter_features_("gene_symbols %in% c('LIG4', 'MAPK12', 'MAPK1')", verbose = TRUE)
#' }
#' @export
filter_features <- function(object, condition, verbose = FALSE){
   filter_features_(object, lazyeval::lazy(condition), verbose = verbose)
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
      object %<>% magrittr::extract(which(autonomics.import::fdata(.)[[paste0('value.', names(contrast))]] > 0))
   } else if (direction == 'neg'){
      object %<>% magrittr::extract(which(autonomics.import::fdata(.)[[paste0('value.', names(contrast))]] < 0))
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


#' @rdname filter_samples
#' @importFrom magrittr %<>%
#' @export
filter_samples_ <- function(object, condition, verbose = FALSE){
   if (is.null(condition)) return(object)
   idx <- lazyeval::lazy_eval(condition, autonomics.import::sdata(object))
   idx <- idx & !is.na(idx)
   if (verbose) if (verbose) message('\t\tRetain ', sum(idx), '/', length(idx), ' features: ', if (class(condition)=='lazy') deparse(condition$expr) else condition)
   object %<>% magrittr::extract(, idx)
   autonomics.import::sdata(object) %<>% droplevels()
   object
}

#' Filter samples on condition
#' @param object SummarizedExperiment
#' @param condition filter condition
#' @param verbose logical
#' @return filtered eSet
#' @examples
#' if (require(autonomics.data)){
#'    filter_samples( ALL,  sex == 'M')
#'    filter_samples_(ALL, "sex == 'M'")
#' }
#' @export
filter_samples <- function(object, condition, verbose = FALSE){
   filter_samples_(object, lazyeval::lazy(condition), verbose = verbose)
}

