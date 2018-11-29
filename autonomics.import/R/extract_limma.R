#' Extract limma matrix from object
#' @param object SummarizedExperiment
#' @param quantity 'value', p', 'fdr', 'bonf'
#' @return matrix (nfeature x ncontrast)
#' @examples
#' require(magrittr)
#' if (require(subramanian.2016)){
#'    object <- subramanian.2016::metabolon
#'    object %>% extract_limma_matrix(quantity = 'effect') %>% extract(1:3, 1:3)
#'    object %>% effect() %>% extract(1:3, 1:3)
#'    object %>% p()      %>% extract(1:3, 1:3)
#'    object %>% fdr()    %>% extract(1:3, 1:3)
#' }
#' @importFrom magrittr %>%
#' @export
extract_limma_matrix <- function(object, quantity){

   # Extract limma res
   limma_res <- object %>% autonomics.import::limma()

   # Assert
   assertive.sets::assert_is_subset(quantity, dimnames(limma_res)[[3]])

   # Return
   limma_res %>% magrittr::extract(autonomics.import::fnames(object), , quantity)
}

#' @rdname extract_limma_matrix
#' @importFrom magrittr %>%
#' @export
effect <- function(object){
   object %>% autonomics.import::extract_limma_matrix('effect')
}

#' @rdname extract_limma_matrix
#' @importFrom magrittr %>%
#' @export
p <- function(object){
   object %>% autonomics.import::extract_limma_matrix('p')
}


#' @rdname extract_limma_matrix
#' @importFrom magrittr %>%
#' @export
fdr <- function(object){
   object %>% autonomics.import::extract_limma_matrix('fdr')
}

#' Extract limma datatable
#' @param object SummarizedExperiment
#' @return melted data.table
#' @examples
#' \dontrun{
#' require(magrittr)
#' if (require(subramanian.2016)){
#'    subramanian.2016::metabolon %>% autonomics.import::extract_limma_dt()
#' }
#' }
#' @importFrom magrittr %>%
#' @export
extract_limma_dt <- function(object){

   # Satisfy CHECK
   . <- NULL

   c('effect', 'p', 'fdr', 'bonf') %>%
   lapply(function(quantity){
      data.table::data.table(
         fid   = object %>% autonomics.import::fid_values(),
         fname = object %>% autonomics.import::fname_values(),
         object %>% autonomics.import::extract_limma_matrix(quantity)
      ) %>%
         data.table::melt.data.table(
            id.vars       = c('fid', 'fname'),
            variable.name = 'contrast',
            value.name    = quantity
         )
   }) %>%
   Reduce(function(x,y) merge(x,y, by = c('fid', 'fname', 'contrast')), .)

}
