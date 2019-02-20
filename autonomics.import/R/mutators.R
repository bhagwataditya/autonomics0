#==========================================================================
#' @title log2 transform assayData matrices
#' @description log2 transform assayData
#' @param x SummarizedExperiment
#' @rdname log2
#' @export
setMethod("log2", signature("SummarizedExperiment"), function(x){ # getGeneric('log2'): argument must be called 'x'!
   autonomics.import::exprs(x) <- log2(autonomics.import::exprs(x))
   x
})


#============================================================================
#' cbind two esets
#' @param x eSet
#' @param y eSet
#' @importFrom magrittr   %>%   %<>%
#' @importFrom methods cbind2
#' @export
setMethod("cbind2", signature("eSet", "eSet"), function(x, y){

   # Assert that samples are different
   assertive.sets::assert_are_disjoint_sets(autonomics.import::snames(x), autonomics.import::snames(y))

   # Assert that remaining portions are identical
   assertive.sets::assert_are_set_equal(autonomics.import::fnames(x), autonomics.import::fnames(y))
   assertive.sets::assert_are_set_equal(autonomics.import::fvars(x),  autonomics.import::fvars(y))
   assertive.sets::assert_are_set_equal(autonomics.import::svars(x),  autonomics.import::svars(y))
   assertive.sets::assert_are_set_equal(autonomics.import::prepro(x), autonomics.import::prepro(y))

   # cbind
   exprs1 <- cbind(autonomics.import::exprs(x), autonomics.import::exprs(y))
   eset1 <- Biobase::ExpressionSet(exprs1)
   autonomics.import::fdata(eset1)  <- autonomics.import::fdata(x)
   autonomics.import::sdata(eset1)  <- rbind(autonomics.import::sdata(x), autonomics.import::sdata(y))
   autonomics.import::prepro(eset1) <- autonomics.import::prepro(x)
   eset1
})


#===========================================================================================
#' Rename subgroup levels
#' @param object eSet
#' @param mapping named character vector (names = old, values = new)
#' @importFrom magrittr         %<>%
#' @export
rename_subgroups <- function(object, mapping){
   assertive.sets::assert_is_subset(names(mapping), unique(object$subgroup))

   for (i in seq_along(mapping)){
      old <- names(mapping[i])
      new <- unname(mapping[i])
      object$subgroup %<>% gsub(old, new, .)
      Biobase::sampleNames(object) %<>% gsub(old, new, .)
      object$sample_id %<>% gsub(old, new, .)
   }

   return(object)
}


#============================================================================================
#' Replace NAs with zeros
#' @param  object  eSet
#' @param  verbose   TRUE or FALSE
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    object <- autonomics.data::billing2016
#'    replace_nas_with_zeros(object)
#' }
#' if (require(billing.differentiation.data)){
#'    object <- billing.differentiation.data::rna.voomcounts
#'    object %>% autonomics.import::replace_nas_with_zeros()
#' }
#' @return eset
#' @importFrom magrittr %>%
#' @export
replace_nas_with_zeros <- function(object, verbose = TRUE){

   # Record
   na_features <- matrixStats::rowAnys(is.na(autonomics.import::exprs(object))) %>% sum()
   total_features <- nrow(object)

   # Replace
   selector <- is.na(autonomics.import::exprs(object))
   if (sum(selector)>0){
      autonomics.import::exprs(object)[selector] <- 0
   }

   # Report
   if (verbose){
      autonomics.support::cmessage('\t\tReplace NA with 0 in %s/%s features',
                                   na_features, total_features)
   }

   # Return
   object
}


#==========================================================================
#' Split on svar
#' @param object SummarizedExperiment
#' @param svar sample var
#' @return list of SummarizedExperiments
#' @examples
#' require(magrittr)
#' if (require(subramanian.2016)){
#'    object <- subramanian.2016::metabolon
#'    svar = 'condition'
#'    object %>% autonomics.import::split_on_svar(svar)
#' }
#' @importFrom magrittr %>%
#' @export
split_on_svar <- function(object, svar = NULL){

   # Return object if null svar
   if (is.null(svar)) return(list(object))

   # Split
   autonomics.import::slevels(object, svar) %>%
      Map(function(curlevel){
         object %>% autonomics.import::filter_samples_(sprintf("%s=='%s'", svar, curlevel))
      }, .)
}


