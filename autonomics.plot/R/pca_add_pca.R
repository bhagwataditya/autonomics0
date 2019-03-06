
#' Add PCA/LDA/PLS/SMA results to object
#' @param object          SummarizedExperiment
#' @param method          'pca' or 'lda'
#' @param implementation  string
#' @param na.impute       TRUE or FALSE
#' @param ndim            number of dimensions to include
#' @param ...             passed to add_projection
#' @examples
#' if (require(autonomics.data)){
#'
#'    # STEM CELL COMPARISON
#'    require(magrittr)
#'    object <- autonomics.data::stemcomp.proteinratios
#'    object %>%  add_pca(ndim = 3)
#'    object %<>% add_pca()
#'    object %>%  autonomics.import::svars()
#'    object %>%  autonomics.import::fvars()
#' }
#'
#' @importFrom magrittr  %>%
#' @export
add_projection <- function(
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

   projection <- object %>% project(method         = method,
                                    implementation = implementation,
                                    na.impute      = na.impute,
                                    ndim           = ndim)
   pca_features <- data.frame(feature_id = rownames(projection$features), projection$features)
   autonomics.import::fdata(object) %<>% autonomics.support::left_join_keeping_rownames(pca_features, by = 'feature_id')
   object
}

#' @rdname add_projection
#' @export
add_pca <- function(object, ...){
   add_projection(object, method = 'pca', ...)
}

#' @rdname add_projection
#' @export
add_sma <- function(object, ...){
   add_projection(object, method = 'sma', ...)
}

#' @rdname add_projection
#' @export
add_lda <- function(object, ...){
   add_projection(object, method = 'lda', ...)
}

#' @rdname add_projection
#' @export
add_pls <- function(object, ...){
   add_projection(object, method = 'pls', ...)
}
