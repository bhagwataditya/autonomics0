

merge_limma_matrices <- function(limma_matrices){
   matrix_names  <- names(limma_matrices)
   comparisons   <- colnames(limma_matrices[[1]])
   n_matrices    <- length(matrix_names)
   n_comparisons <- length(comparisons)
   n_features    <- nrow(limma_matrices[[1]])
   fmat <- matrix(NA, nrow = n_features, ncol = n_comparisons*n_matrices)
   fcolnames <- rep(NA, n_comparisons*n_matrices)
   for (i in 1:n_matrices){
      selector <- seq(i, n_comparisons*n_matrices, n_matrices)
      fmat[, selector] <- limma_matrices[[i]]
      fcolnames[selector] <- paste0(matrix_names[[i]], '.', comparisons)
   }
   colnames(fmat) <- fcolnames
   fmat
}

#' Run "limma" analysis and add to SummarizedExperiment
#'
#' Performs a limma analysis for the specified contrasts.
#'
#' @param  object   SummarizedExperiment
#' @param  contrasts  a character vector of contrasts (of subgroup levels)
#' @param  design     design matrix
#' @param  overwrite  whether to overwrite existing results in fdata (logical)
#' @return updated SummarizedExperiment with limma results in fdata
#'         added columns: rank.xxx value.xxx p.xxx fdr.xxx bonf.xxx \cr
#' @author Aditya Bhagwat
#' @note Features with a single observation can get a significant p value
#' (though not a significant q value) when the exprs value for that feature is exceptionally large.
#' This is a consequence of the shrinkage estimation of the limma method, as
#' \href{https://support.bioconductor.org/p/4932/}{explained by Gordon Smyth on Bioconductor Support}.
#' @examples
#' library(magrittr)
#' 
#' # STEM CELL COMPARISON
#' if (require(autonomics.data)){
#'    object <- autonomics.data::stemcomp.proteinratios
#'    object %>% autonomics.find::add_limma_to_fdata() %>% 
#'               autonomics.import::fdata() %>%
#'               extract2('bonf.BM_EM') %>% 
#'               is_less_than(0.05) %>% 
#'               sum(na.rm = TRUE)
#' }
#' 
#' if (require(subramanian.2016)){
#'    object <- subramanian.2016::metabolon
#'    object %>% autonomics.find::add_limma_to_fdata() %>% autonomics.import::fvars()
#'    contrasts <- subramanian.2016::contrasts.metabolon
#'    object %<>% autonomics.find::add_limma_to_fdata(contrasts)
#' }
#' if (require(billing.differentiation.data)){
#'    contrasts <- billing.differentiation.data::contrasts
#'    billing.differentiation.data::rna.counts %>% autonomics.import::logcpm() %>%
#'       autonomics.find::add_limma_to_fdata(contrasts) %>%  autonomics.import::fdata() %>%
#'       extract2('bonf.EM0.8_0') %>% is_less_than(0.05) %>% sum(na.rm=TRUE)
#'    billing.differentiation.data::rna.voomcounts %>%
#'       autonomics.find::add_limma_to_fdata(contrasts) %>%  autonomics.import::fdata() %>%
#'       extract2('bonf.EM0.8_0') %>% is_less_than(0.05) %>% sum(na.rm=TRUE)
#' }
#' @importFrom magrittr   %>%   %<>%
#' @export
add_limma_to_fdata <- function(
   object,
   contrasts  = autonomics.find::default_contrasts(object),
   design     = autonomics.find::create_design_matrix(object),
   overwrite  = TRUE
){

  # Sumexps contains limma: overwrite or return
   autonomics.import::assert_is_valid_eset(object)
   if (autonomics.find::contains_limma_in_fdata(object)){
      if (overwrite){
         autonomics.import::fdata(object) %<>% magrittr::extract(, - grep('^(p|fdr|bonf|rank|quantile|value)[.]', names(.)), drop = FALSE)
      } else {
         message('\t\tSummarizedExperiment contains limma results already - returning')
         return(object)
      }
   }

  # Assert valid inputs
   assertive.types::assert_is_matrix(design)
   assertive.types::assert_is_numeric(design)
   assertive.base::assert_is_identical_to_true(unname(ncol(object)) == nrow(design))
   if (is.null(contrasts))    return(object)
   contrasts %<>% autonomics.find::validify_contrasts(design)

  # Return if no contrasts
  if (length(contrasts)==0) return(object)

  # Run limma
  limma_df <- object %>%
              # autonomics.preprocess::normalize_smart(plot = TRUE, result_dir = result_dir) %>%
              autonomics.find::run_limma(contrasts = contrasts, design = design) %>%
              merge_limma_matrices() %>%
              data.frame()

  # Update eset and return
  autonomics.import::fdata(object) %<>% cbind(limma_df)
  return(object)

}
