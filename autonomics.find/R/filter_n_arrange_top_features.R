
#' Which are the top features?
#'
#' Which features fulfill top definition and direction for specified contrast?
#' @param object       eset
#' @param top_definition top definition
#' @param contrast_name  contrast name (string)
#' @param direction      'both', 'pos', or 'neg'
#' @examples
#' require(magrittr)
#' if (require(billing.differentiation.data)){
#'    object <- billing.differentiation.data::rna.voomcounts
#'    contrast <-  c(EM0.8_0 = 'EM0.8 - EM0.0')
#'    topdef <- 'bonf < 0.05'
#'    autonomics.find::are_top_features(
#'       object, topdef, names(contrast), direction = 'neg')  %>% 
#'       sum()
#'    autonomics.find::are_top_features(
#'       object, topdef, names(contrast), direction = 'pos')  %>% 
#'       sum()
#'    autonomics.find::are_top_features(
#'       object, topdef, names(contrast), direction = 'both') %>% 
#'       sum()
#' }
#' @return integer vector with indices of top features
#' @importFrom magrittr   %<>%
#' @export
are_top_features <- function(
   object, 
   top_definition  = autonomics.find::default_top_definition(object), 
   contrast_name, 
   direction
){
   autonomics.import::assert_is_valid_eset(object)
   assertive.types::assert_is_a_string(contrast_name)
   assertive.types::assert_is_a_string(direction)
   assertive.sets::assert_is_subset(direction, c('both', 'pos', 'neg'))
   
   top_definition %<>% autonomics.find::complete_top_definition(contrast_name, direction)
   selector <- lazyeval::lazy_eval(top_definition, autonomics.import::fdata(object))
   selector[is.na(selector)] <- FALSE
   
   selector
}

#' @rdname are_top_features
#' @importFrom magrittr   %>% 
#' @export
get_top_features <- function(
   object,
   top_definition = autonomics.find::default_top_definition(object), 
   contrast_name, 
   direction
){
   idx <- object %>% are_top_features(top_definition = top_definition, contrast_name = contrast_name, direction = direction)
   autonomics.import::fnames(object)[idx]
}

#' Filter (and arrange) top features
#' @param object           eset
#' @param contrast_name    contrast name
#' @param top_definition   definition of 'top features'.
#' @param direction        'both', 'neg', or 'pos'
#' @examples
#' require(magrittr)
#' if (require(billing.differentiation.data)){
#'    object <- billing.differentiation.data::rna.voomcounts
#'    contrast <-  c(EM0.8_0 = 'EM0.8 - EM0.0')
#'    topdef <- 'bonf < 0.05'
#'    autonomics.find::filter_top_features(object, names(contrast), topdef, 'neg') %>% dim()
#'    autonomics.find::filter_top_features(object, names(contrast), topdef, 'both') %>% dim()
#' }
#' @importFrom magrittr %>%
#' @export
filter_top_features <- function(object, contrast_name, top_definition, direction){
   idx <- are_top_features(object, top_definition, contrast_name, direction)
   object %>% magrittr::extract(idx, )
}

#' Order features on rank
#' @param object eSet
#' @param contrast_name contrast name
#' @return arranged eSet
#' @importFrom magrittr %<>%
#' @export
#' @examples 
#' require(magrittr)
#' if (require(billing.differentiation.data)){
#'    object <- billing.differentiation.data::rna.voomcounts
#'    object %>% autonomics.import::fnames() %>% head()
#'    contrast <-  c(EM0.8_0 = 'EM0.8 - EM0.0')
#'    object %>% autonomics.find::arrange_features_by_rank(names(contrast)) %>% 
#'                 autonomics.import::fnames() %>% head()
#' }
arrange_features_by_rank <- function(object, contrast_name){
   rank_var <- sprintf('rank.%s', contrast_name)
   object %<>% autonomics.import::arrange_features_(rank_var)
   autonomics.import::fdata(object)[[rank_var]] <- seq_len(nrow(object))
   return(object)
}


#' Filter and arrange top features
#' @param object       eset
#' @param contrast_name  contrast name
#' @param top_definition top definition
#' @param direction      'neg' or 'pos'
#' @param nmax           max number of features
#' @examples
#' require(magrittr)
#' if (require(billing.differentiation.data)){
#'    object <- billing.differentiation.data::rna.voomcounts
#'    contrast <-  c(EM0.8_0 = 'EM0.8 - EM0.0')
#'    topdef <- 'bonf < 0.05'
#'    autonomics.find::filter_n_arrange_top_features(object, names(contrast), topdef, 'neg') %>% dim()
#'    autonomics.find::filter_n_arrange_top_features(object, names(contrast), topdef, 'both') %>% dim()
#' }
#' @importFrom magrittr  %>%
#' @export
filter_n_arrange_top_features <- function(object, contrast_name, top_definition, direction, nmax = Inf){
   object                                                                           %>%
   autonomics.find::filter_top_features(contrast_name, top_definition, direction)   %>%
   autonomics.find::arrange_features_by_rank(contrast_name)                         %>%
   magrittr::extract(seq_len(min(nmax, nrow(.))), )
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
#'    autonomics.data::billing2016               %>% 
#'       autonomics.find::add_limma_to_fdata()   %>% 
#'       autonomics.find::filter_significant_features()
#' }
#' if (require(subramanian.2016)){
#'    subramanian.2016::metabolon %>% 
#'       autonomics.find::add_limma_to_fdata(
#'          contrasts = subramanian.2016::contrasts.metabolon %>% extract(-length(.))) %>% 
#'       autonomics.find::filter_significant_features() %>% invisible()
#'    
#'    subramanian.2016::exiqon %>% 
#'       autonomics.find::add_limma_to_fdata(
#'          contrasts  = subramanian.2016::contrasts.exiqon %>% extract(-length(.))) %>% 
#'       autonomics.find::filter_significant_features() %>% invisible()
#'    
#'    subramanian.2016::rnaseq %>% 
#'       autonomics.find::add_limma_to_fdata(
#'          contrasts = subramanian.2016::contrasts.rnaseq %>% extract(-length(.))) %>%
#'       autonomics.find::filter_significant_features() %>% 
#'       autonomics.find::filter_top_ranked_features() %>% invisible()
#' }
#' @return filtered eset
#' @importFrom data.table   data.table   :=   .N
#' @importFrom magrittr     %>%
#' @export
filter_significant_features <- function(object){
   
   # Return eset if no limma in fdata
   if (!contains_limma_in_fdata(object)){
      autonomics.support::cmessage('\t\tNo limma in fdata - abort')
      return(object) 
   } 
   
   # Filter significant
   contains_p <- autonomics.import::fvars(object) %>% stringi::stri_detect_fixed('p.') %>% any()
   if (contains_p){
      selector <- autonomics.import::fdata(object) %>%
                  magrittr::extract(, stringi::stri_detect_fixed(names(.), 'p.')) %>%
                  # data.table::as.data.table() %>%
                  magrittr::is_less_than(0.05) %>% 
                 (function (x){x[is.na(x)] <- FALSE; x}) %>% 
                  matrixStats::rowAnys()
      autonomics.support::cmessage('\t\tFilter %d/%d features significant for at least one contrast (p < 0.05)', sum(selector), length(selector))
      object %<>% magrittr::extract(selector, )
   }
   
   # Return
   object
   
}

#' @rdname filter_significant_features
#' @export
filter_top_ranked_features <- function(object, n = 1000){
   
   # Return eset if no limma in fdata
   if (!contains_limma_in_fdata(object)){
      autonomics.support::cmessage('\t\tNo limma in fdata - abort')
      return(object) 
   } 
   
   # Order on rank
   selector <- autonomics.import::fdata(object) %>% 
               magrittr::extract(, stringi::stri_detect_fixed(names(.), 'rank.')) %>%
               data.matrix() %>% 
               matrixStats::rowMins() %>% 
               order()
   object %<>% magrittr::extract(selector, )
   
   # Filter top ranked
   if (n < nrow(object)){
      autonomics.support::cmessage('\t\tRestrict to %d/%d top ranked features for performance reasons', n, nrow(object))
      object %<>% magrittr::extract(1:n, )
   }
}
