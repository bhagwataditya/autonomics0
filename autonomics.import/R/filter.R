#=================
# FILTER EXPRS
#=================

utils::globalVariables('.')


#' Filter features with replicated expression in some subgroup
#' @param object      SummarizedExperiment
#' @param comparator  '>' or '!='
#' @param lod         numeric(1)
#' @return Filtered SummarizedExperiment
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    object <- autonomics.data::stemcomp.proteinratios
#'    object %>% autonomics.import::filter_exprs_replicated_in_some_subgroup()
#'    object <- autonomics.data::glutaminase
#'    object %>% autonomics.import::filter_exprs_replicated_in_some_subgroup()
#' }
#' @importFrom magrittr %>%
#' @export
filter_exprs_replicated_in_some_subgroup <- function(
   object,
   comparator = if (autonomics.import::contains_ratios(object)) '!=' else '>',
   lod = 0
){

   # Datatablify
   fidvar <- object %>% autonomics.import::fid_var()
   dt     <- object %>% autonomics.import::sumexp_to_long_dt(svars = 'subgroup')

   # Is expr replicated in its subgroup?
   if (comparator ==  '>') dt %>% magrittr::extract(, replicated_in_its_subgroup := sum(value  > lod, na.rm=TRUE) > 1,   by = c(fidvar, 'subgroup'))
   if (comparator == '!=') dt %>% magrittr::extract(, replicated_in_its_subgroup := sum(value != lod, na.rm=TRUE) > 1,   by = c(fidvar, 'subgroup'))

   # Is it replicated in any subgroup
   dt %>% magrittr::extract(, replicated_in_any_subgroup := any(replicated_in_its_subgroup),    by =  fidvar)

   # Keep only replicated features
   replicated_features <- dt %>% magrittr::extract(replicated_in_any_subgroup == TRUE, get(fidvar))
   idx <- autonomics.import::fid_values(object) %in% replicated_features
   autonomics.support::cmessage('\t\tFilter %d/%d features: expr %s %s, for at least two samples in some subgroup', sum(idx), length(idx), comparator, as.character(lod))
   object %>% magrittr::extract(idx, )
}




#=================
# FILTER FEATURES
#=================

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
#' \dontrun{
#' if (require(autonomics.data)){
#'    autonomics.data::ALL %>%
#'    autonomics.import::filter_features(gene_symbols %in% c('LIG4', 'MAPK12', 'MAPK1') )
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
#' if (require(autonomics.data)){
#'    require(magrittr)
#'    fvalues <- c("A0MZ66", "A0JLT2", "Q8NG66", "Q9NPI5")
#'    object <- autonomics.data::stemdiff.proteinratios
#'    object %>% is_fvalue_feature(fvar = 'Uniprot accessions',
#'                                 split = ';', fvalues = fvalues) %>%
#'               sum()
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
#' if (require(autonomics.data)){
#'    require(magrittr)
#'    fvalues <- c("A0MZ66", "A0JLT2", "Q8NG66", "Q9NPI5")
#'    object <- autonomics.data::stemdiff.proteinratios
#'    object %>% filter_features_on_fvalues(fvar    = 'Uniprot accessions',
#'                                          split   = ';',
#'                                          fvalues = fvalues)
#' }
#' @return subset of object with reference features
#' @importFrom magrittr       %>%     %<>%
#' @export
filter_features_on_fvalues <- function(object, fvar, split, fvalues){
   fvalues %<>% unique()
   idx <- object %>% is_fvalue_feature(fvar, split, fvalues)
   object %>% magrittr::extract(idx, )
}




#=========================
# FILTER SAMPLES
#=========================

#' @rdname filter_samples
#' @importFrom magrittr %<>%
#' @export
filter_samples_ <- function(object, condition, verbose = FALSE){
   if (is.null(condition)) return(object)
   idx <- lazyeval::lazy_eval(condition, autonomics.import::sdata(object))
   idx <- idx & !is.na(idx)
   if (verbose) if (verbose) message('\t\tRetain ', sum(idx), '/', length(idx), ' samples: ', if (class(condition)=='lazy') deparse(condition$expr) else condition)
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

