#' @rdname are_contrast_features
#' @export
are_top_features <- function(object, topdef, contrast_name){
   .Deprecated('are_contrast_features')
   are_contrast_features(object, contrast_name = contrast_name, topdef = topdef)
}

#' Are contrast features?
#'
#' @param object         SummarizedExperiment
#' @param contrast_name  character(1): contrast name 
#' @param topdef         character(1): top definition
#' @param ...            for backward compatibility only
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    object <- autonomics.data::stemdiff.proteinratios
#'    topdef <- 'p < 0.05'
#'    contrast_name <- object %>% autonomics.import::contrastdefs() %>% magrittr::extract(1) %>% names()
#'    object %>% autonomics.find::are_contrast_features('p<0.05',            contrast_name) %>% sum()
#'    object %>% autonomics.find::are_contrast_features('p<0.05 & effect>0', contrast_name) %>% sum()
#' }
#' @return integer vector with indices of top features
#' @importFrom magrittr   %<>%
#' @export
are_contrast_features <- function(
   object, 
   contrast_name,
   topdef  = autonomics.find::default_topdef(object)
){
   autonomics.import::assert_is_valid_eset(object)
   assertive.types::assert_is_a_string(contrast_name)
   
   limma_df <- object %>% 
               autonomics.import::limma() %>% 
               magrittr::extract(, contrast_name, ) %>% 
               data.frame(stringsAsFactors = FALSE)

   selector <- lazyeval::lazy_eval(topdef, limma_df)
   selector[is.na(selector)] <- FALSE
   
   selector
}


#' Get number of contrast features
#'
#' @param object         SummarizedExperiment
#' @param contrast_name  character(1): contrast name 
#' @param topdef         character(1): top definition
#' @param ...            for backward compatibility only
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    object <- autonomics.data::stemdiff.proteinratios
#'    topdef <- 'p < 0.05'
#'    contrast_name <- object %>% autonomics.import::contrastdefs() %>% magrittr::extract(1) %>% names()
#'    object %>% autonomics.find::n_contrast_features('p<0.05',            contrast_name) %>% sum()
#'    object %>% autonomics.find::n_contrast_features('p<0.05 & effect>0', contrast_name) %>% sum()
#' }
#' @return integer vector with indices of top features
#' @importFrom magrittr   %<>%
#' @export
n_contrast_features <- function(object, contrast_name, topdef = autonomics.find::default_topdef(object)){
   object %>% 
   autonomics.find::are_contrast_features(contrast_name = contrast_name, topdef = topdef) %>% 
   sum()
}


#' @rdname are_contrast_features
#' @export
get_top_features <- function(object, topdef, contrast_name){
   .Deprecated('get_contrast_features')
   get_contrast_features(object, contrast_name = contrast_name, topdef = topdef)
}

#' @rdname are_contrast_features
#' @importFrom magrittr   %>% 
#' @export
get_contrast_features <- function(
   object,
   contrast_name,
   topdef = autonomics.find::default_topdef(object)
){
   idx <- object %>% are_contrast_features(contrast_name = contrast_name, topdef = topdef)
   autonomics.import::fnames(object)[idx]
}

#' @rdname filter_contrast_features
#' @export
filter_top_features <- function(...){
   .Deprecated('filter_top_features')
   filter_contrast_features(...)
}

#' Filter (and arrange) contrast features
#' @param object           SummarizedExperiment
#' @param contrast_name    character(1): contrast name
#' @param topdef           character(1): top definition.
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    (object <- autonomics.data::stemdiff.proteinratios)
#'    contrast_name <-  names(autonomics.import::contrastdefs(object))[1]
#'    topdef <- 'bonf < 0.05'
#'    (object %<>% autonomics.find::filter_contrast_features(contrast_name, topdef))
#'    object %>% autonomics.import::limma()
#' }
#' @importFrom magrittr %>%
#' @export
filter_contrast_features <- function(object, contrast_name, topdef){
   idx <- are_contrast_features(object, contrast_name = contrast_name, topdef = topdef)
   object %>% autonomics.import::extract_features(idx)
}

#' Order features on rank
#' @param object SummarizedExperiment
#' @param contrast_name contrast name
#' @return arranged eSet
#' @importFrom magrittr %<>%
#' @export
#' @examples 
#' require(magrittr)
#' if (require(autonomics.data)){
#'    (object <- autonomics.data::stemdiff.proteinratios)
#'    contrast_name <-  names(autonomics.import::contrastdefs(object))[1]
#'    (object %<>% autonomics.find::arrange_features_by_rank(contrast_name))
#'    object %>% autonomics.import::limma() %>% magrittr::extract(1:2, contrast_name, 1:2)
#' }
arrange_features_by_rank <- function(object, contrast_name){
   idx <- object %>% autonomics.import::limma() %>% 
                     magrittr::extract(, contrast_name, 'rank') %>%
                     order()
   object %>% autonomics.import::extract_features(idx)
}

#' @rdname filter_n_arrange_contrast_features
#' @export
filter_n_arrange_top_features <- function(...){
   .Deprecated('filter_n_arrange_contrast_features')
   filter_n_arrange_contrast_features(...)
}
   
#' Filter and arrange contrast features
#' @param object         SummarizedExperiment
#' @param contrast_name  contrast name
#' @param topdef         top definition
#' @param nmax           max number of features
#' @param ...            for backward compatibility only
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    (object <- autonomics.data::stemdiff.proteinratios)
#'    contrast_name <- names(autonomics.import::contrastdefs(object))[1]
#'    topdef <- 'bonf < 0.05 & effect<0'
#'    object %<>% filter_n_arrange_contrast_features(contrast_name, topdef)
#'    object %>% autonomics.import::limma()
#' }
#' @importFrom magrittr  %>%
#' @export
filter_n_arrange_contrast_features <- function(object, contrast_name, topdef, nmax = Inf){
   object                                                      %>%
   autonomics.find::filter_contrast_features(contrast_name, topdef) %>%
   autonomics.find::arrange_features_by_rank(contrast_name)    %>%
   autonomics.import::extract_features(seq_len(min(nmax, nrow(.))))
}

#' Filter significant features
#' 
#' Filter features which are significant for at least one contrast
#' 
#' @param object eset
#' @param n number of top features to select
#' @examples
#' library(magrittr)
#' if (require(autonomics.data)){
#'    autonomics.data::stemcomp.proteinratios %>% autonomics.find::filter_significant_features()
#' }
#' if (require(subramanian.2016)){
#'    subramanian.2016::metabolon %>% autonomics.find::filter_significant_features()
#'    subramanian.2016::exiqon    %>% autonomics.find::filter_significant_features()
#'    subramanian.2016::rnaseq    %>% autonomics.find::filter_significant_features() 
#' }
#' @return filtered eset
#' @importFrom data.table   data.table   :=   .N
#' @importFrom magrittr     %>%
#' @export
filter_significant_features <- function(object){

   # Return if no limma results in object
   if (is.null(autonomics.import::limma(object))) return(object)
   
   # Return if no p values in object
   contains_pvalues <- autonomics.import::limma(object) %>% dimnames() %>% magrittr::extract2(3) %>% magrittr::is_in('p', .)
   if (!contains_pvalues) return(object)
   
   # Filter significant
   selector <- object                                 %>% 
               autonomics.import::limma()             %>% 
               magrittr::extract(, , 'p')             %>%
               magrittr::is_less_than(0.05)           %>% 
              (function (x){x[is.na(x)] <- FALSE; x}) %>% 
               matrixStats::rowAnys()
   autonomics.support::cmessage('\t\tFilter %d/%d features significant for at least one contrast (p < 0.05)', sum(selector), length(selector))
   object %<>% autonomics.import::extract_features(selector)

   # Return
   object
   
}

#' @rdname filter_significant_features
#' @export
filter_top_ranked_features <- function(object, n = 1000){
   
   # Return eset if no limma in fdata
   if (is.null(autonomics.import::limma(object))){
      autonomics.support::cmessage('\t\tNo limma in object - abort')
      return(object) 
   } 
   
   # Order on rank
   selector <- autonomics.import::limma(object) %>% 
               magrittr::extract(,, 'rank')     %>%
               matrixStats::rowMins()           %>% 
               order()
   object %<>% autonomics.import::extract_features(selector)
   
   # Filter top ranked
   if (n < nrow(object)){
      autonomics.support::cmessage('\t\tRestrict to %d/%d top ranked features for performance reasons', n, nrow(object))
      object %<>% autonomics.import::extract_features(1:n)
   }
}
