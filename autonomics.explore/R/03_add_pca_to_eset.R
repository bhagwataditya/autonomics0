
#' Add PCA/LDA/PLS/SMA results to eset
#' @param object SummarizedExperiment, eSet, or EList
#' @param method 'pca' or 'lda'
#' @param implementation (character)
#' @param na.impute TRUE or FALSE
#' @param ndim number of dimensions to include
#' @param ... passed to add_projection_to_eset
#' @examples 
#' require(magrittr)
#' # STEM CELL COMPARISON
#' if (require(autonomics.data)){
#'    object <- autonomics.data::stemcomp.proteinratios
#'    object %>% autonomics.explore::add_pca_to_eset(ndim = 3)
#'    object %<>% autonomics.explore::add_pca_to_eset()
#'    object %>% autonomics.import::svars()
#'    object %>% autonomics.import::fvars()
#' }
#' 
#' # STEM CELL DIFFERENTIATION
#' if (require(billing.differentiation.data)){
#'    object <- billing.differentiation.data::rna.voomcounts
#'    object %>% autonomics.explore::add_pca_to_eset()
#' }
#' 
#' if (require(atkin.2014)){
#'    object <- atkin.2014::soma 
#'    object %>% autonomics.explore::add_pca_to_eset()
#' }
#' @importFrom magrittr  %>%
#' @export
add_projection_to_eset <- function(
   object, 
   method         = c('pca', 'lda', 'pls')[1],
   implementation = NULL, 
   na.impute      = FALSE, 
   ndim           = 2
){
   nsample <- ncol(object)
   npca    <- min(nsample-1, 2)
   if (npca == 0) return(object)
   
   # Wipe previous pca results, as they could be 
   #   - on a different data subset
   #   - for different ndim
   idx <- autonomics.import::fvars(object) %>% stringi::stri_detect_regex(paste0('^', method)) %>% which()
   if (length(idx)>0){
      autonomics.import::fdata(object) %<>% magrittr::extract(, -idx)
      # autonomics.support::cmessage('\t\tAbort %s - object already contains projection scores', method)
      # return(object)
   }
   
   projection <- object %>% autonomics.explore::project(method         = method, 
                                                        implementation = implementation, 
                                                        na.impute      = na.impute, 
                                                        ndim           = ndim)
   pca_features <- data.frame(feature_id = rownames(projection$features), projection$features)
   autonomics.import::fdata(object) %<>% autonomics.support::left_join_keeping_rownames(pca_features, by = 'feature_id')
   object
}

#' @rdname add_projection_to_eset
#' @export
add_pca_to_eset <- function(object, ...){
   add_projection_to_eset(object, method = 'pca', ...)
}

#' @rdname add_projection_to_eset
#' @export
add_sma_to_eset <- function(object, ...){
   add_projection_to_eset(object, method = 'sma', ...)
}

#' @rdname add_projection_to_eset
#' @export
add_lda_to_eset <- function(object, ...){
   add_projection_to_eset(object, method = 'lda', ...)
}

#' @rdname add_projection_to_eset
#' @export
add_pls_to_eset <- function(object, ...){
   add_projection_to_eset(object, method = 'pls', ...)
}
