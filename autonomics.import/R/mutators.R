#==========================================================================
#' @title log2 transform assayData matrices
#' @description log2 transform assayData
#' @param x SummarizedExperiment
#' @rdname log2
#' @export
setMethod("log2", signature("SummarizedExperiment"), function(x){ # getGeneric('log2'): argument must be called 'x'!
   exprs(x) <- log2(exprs(x))
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
   assertive.sets::assert_are_disjoint_sets(snames(x), snames(y))

   # Assert that remaining portions are identical
   assertive.sets::assert_are_set_equal(fnames(x), fnames(y))
   assertive.sets::assert_are_set_equal(fvars(x),  fvars(y))
   assertive.sets::assert_are_set_equal(svars(x),  svars(y))
   assertive.sets::assert_are_set_equal(prepro(x), prepro(y))

   # cbind
   exprs1 <- cbind(exprs(x), exprs(y))
   eset1 <- Biobase::ExpressionSet(exprs1)
   fdata(eset1)  <- fdata(x)
   sdata(eset1)  <- rbind(sdata(x), sdata(y))
   prepro(eset1) <- prepro(x)
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
#' @param  object  SummarizedExperiment
#' @param  verbose TRUE or FALSE
#' @examples
#' if (require(autonomics.data)){
#'    require(magrittr)
#'    object <- autonomics.data::stemcomp.proteinratios
#'    replace_nas_with_zeros(object)
#' }
#' @return eset
#' @importFrom magrittr %>%
#' @export
replace_nas_with_zeros <- function(object, verbose = TRUE){

   # Record
   na_features <- matrixStats::rowAnys(is.na(exprs(object))) %>% sum()
   total_features <- nrow(object)

   # Replace
   selector <- is.na(exprs(object))
   if (sum(selector)>0){
      exprs(object)[selector] <- 0
   }

   # Report
   if (verbose){
      autonomics.support::cmessage('\t\tReplace NA with 0 in %s/%s features',
                                   na_features, total_features)
   }

   # Return
   object
}


