#' @rdname get_median_exprs
#' @importFrom magrittr %>%
#' @export
get_geometric_median_exprs <- function(object){
   object %>%
autonomics.import::exprs() %>%
   t() %>%
   pcaPP::l1median()
}

#' @rdname get_median_exprs
#' @importFrom magrittr %>%
#' @export
get_arithmetic_median_exprs <- function(object){
   object %>%
autonomics.import::exprs() %>%
   matrixStats::rowMedians(na.rm = TRUE) %>%
   magrittr::set_names(rownames(object))
}

#' Get median (sample) exprs
#'
#' Returns exprs of median sample
#' The arithmetic median computes the median across samples separately for each feature.
#' The geometric median corresponds to the point in nfeature-dimensional space for which
#' the total distance to all other points is minimal (See details).
#'
#' The geometric median is robust against outliers and norm invariant.
#' It is also known as the spatial median or the L1 median.
#' We compute it using an implementation based on NLM, which
#' outperformed other implementations in accuracy and speed, and is
#' robust against NAS
#' \href{https://dx.doi.org/10.2139/ssrn.1690502}{Croux et al., 2010}
#' @param object SummarizedExperiment
#' @param method 'arithmetic' or 'geometric'
#' @examples
#' require(magrittr)
#' if (require(subramanian.2016)){
#'    subramanian.2016::exiqon %>% get_median_exprs()
#' }
#' @importFrom magrittr %>%
#' @export
get_median_exprs <- function(object, method = 'arithmetic'){

   # Assert
autonomics.import::assert_is_valid_eset(object)
   assertive.sets::assert_is_subset(method, c('arithmetic', 'geometric'))

   # Compute
   if      (method == 'arithmetic') object %>% autonomics.preprocess::get_arithmetic_median_exprs()
   else if (method == 'geometric')  object %>% autonomics.preprocess::get_geometric_median_exprs()

}

#' Get median reference exprs
#' @param object SummarizedExperiment
#' @param ref reference subgroup level
#' @param method 'arithmetic' or 'geometric'
#' @return vector with median reference exprs
#' @examples
#' require(magrittr)
#' if (require(subramanian.2016)){
#'    subramanian.2016::exiqon %>% get_median_ref_exprs()
#' }
#' @importFrom magrittr %>%
#' @export
get_median_ref_exprs <- function(
   object,
   ref = levels(autonomics.import::sdata(object)$subgroup)[1],
   method = 'arithmetic'
){
   # Satisfy CHECK
   subgroup <- NULL

   # Get median exprs
   object %>%
   autonomics.import::filter_samples(subgroup == ref) %>%
      autonomics.preprocess::get_median_exprs(method)
}

#' Subtract values from exprs
#' @param object SummarizedExperiment
#' @param values values
#' @importFrom magrittr %<>%
#' @export
subtract_exprs <- function(
   object,
   values
){
autonomics.import::exprs(object) %<>% magrittr::subtract(values)
   object
}



#' Subtract median ref exprs
#' @param object SummarizedExperiment
#' @param ref reference subgroup level
#' @param method 'arithmetic' or 'geometric'
#' @return SummarizedExperiment
#' @examples
#' require(magrittr)
#' if (require(subramanian.2016)){
#'    subramanian.2016::exiqon
#'    subramanian.2016::exiqon %>% autonomics.preprocess::subtract_median_ref_exprs()
#' }
#' @importFrom magrittr %>%
#' @export
#' @author Aditya Bhagwat
subtract_median_ref_exprs <- function(
   object,
   ref = levels(autonomics.import::sdata(object)$subgroup)[1],
   method = 'arithmetic'
){

   # Satisfy CHECK
   . <- NULL

   # Get median ref exprs
   median_ref_exprs <- object %>%
                       autonomics.preprocess::get_median_ref_exprs(ref = ref, method = method)

   # Subtract
   object %>% autonomics.preprocess::subtract_exprs(median_ref_exprs)
}


