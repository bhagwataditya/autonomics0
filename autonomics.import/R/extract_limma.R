#' Extract limma matrix from object
#' @param object SummarizedExperiment
#' @param quantity 'effect', 'p', 'fdr', 'bonf'
#' @return matrix (nfeature x ncontrast)
#' @examples
#' if (require(autonomics.data)){
#'    require(magrittr)
#'    object <- autonomics.data::stemcomp.proteinratios
#'    object %>% extract_limma_matrix(quantity = 'effect') %>% extract(1:3, 1:3)
#'    object %>% effect() %>% extract(1:3, 1:3)
#'    object %>% p()      %>% extract(1:3, 1:3)
#'    object %>% fdr()    %>% extract(1:3, 1:3)
#' }
#' @importFrom magrittr %>%
#' @export
extract_limma_matrix <- function(object, quantity){

   # Extract limma res
   limma_res <- object %>% limma()

   # Assert
   assertive.sets::assert_is_subset(quantity, dimnames(limma_res)[[3]])

   # Return
   limma_res %>% magrittr::extract(fnames(object), , quantity)
}

#' @rdname extract_limma_matrix
#' @importFrom magrittr %>%
#' @export
effect <- function(object){
   object %>% extract_limma_matrix('effect')
}

#' @rdname extract_limma_matrix
#' @importFrom magrittr %>%
#' @export
p <- function(object){
   object %>% extract_limma_matrix('p')
}


#' @rdname extract_limma_matrix
#' @importFrom magrittr %>%
#' @export
fdr <- function(object){
   object %>% extract_limma_matrix('fdr')
}

#' Extract limma datatable
#' @param object SummarizedExperiment
#' @return melted data.table
#' @examples
#' if (require(autonomics.data)){
#'    require(magrittr)
#'    autonomics.data::stemcomp.proteinratios %>% extract_limma_dt()
#' }
#' @importFrom magrittr %>%
#' @export
extract_limma_dt <- function(object){

   # Satisfy CHECK
   . <- NULL

   c('effect', 'p', 'fdr', 'bonf') %>%
   lapply(function(quantity){
      data.table::data.table(
         fid   = object %>% fid_values(),
         fname = object %>% fname_values(),
         object %>% extract_limma_matrix(quantity)
      ) %>%
         data.table::melt.data.table(
            id.vars       = c('fid', 'fname'),
            variable.name = 'contrast',
            value.name    = quantity
         )
   }) %>%
   Reduce(function(x,y) merge(x,y, by = c('fid', 'fname', 'contrast')), .)

}
