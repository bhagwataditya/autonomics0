#' Extract limma matrix from object
#' @param object SummarizedExperiment
#' @param quantity 'value', p', 'fdr', 'bonf'
#' @return matrix (nfeature x ncontrast)
#' @examples
#' require(magrittr)
#' if (require(subramanian.2016)){
#'    object <- subramanian.2016::metabolon
#'    object %>% autonomics.import::extract_limma_matrix(quantity = 'value')
#'    object %>% autonomics.import::value()
#'    object %>% autonomics.import::p()
#' }
#' @importFrom magrittr %>%
#' @export
extract_limma_matrix <- function(object, quantity){

   # Assert
   assertive.strings::assert_any_are_matching_fixed(autonomics.import::fvars(object), quantity)

   # Satisfy CHECK
   . <- NULL

   quantity %<>% paste0('.')
   autonomics.import::fdata(object) %>%
   dplyr::select(dplyr::contains(quantity)) %>%
   data.matrix() %>%
   magrittr::set_colnames(colnames(.) %>% stringi::stri_replace_first_fixed(quantity, '')) %>%
   magrittr::set_rownames(autonomics.import::fnames(object))
}

#' @rdname extract_limma_matrix
#' @importFrom magrittr %>%
#' @export
value <- function(object){
   object %>% autonomics.import::extract_limma_matrix('value')
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
#' require(magrittr)
#' if (require(subramanian.2016)){
#'    subramanian.2016::metabolon %>% autonomics.import::extract_limma_dt()
#' }
#' @importFrom magrittr %>%
#' @export
extract_limma_dt <- function(object){

   # Satisfy CHECK
   . <- NULL

   c('value', 'p', 'fdr', 'bonf') %>%
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
