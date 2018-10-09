
#' Run limma and add results
#' 
#' Limma results can be easily accessed with limma(object).
#' This returns a list with components:
#' \itemize{
#'    \item {effect} matrix (ngene x ncontrast): effect sizes
#'    \item {rank}   matrix (ngene x ncontrast): effect ranks
#'    \item {t}      matrix (ngene x ncontrast): t    values (moderated t test)
#'    \item {se}     matrix (ngene x ncontrast): se   values (moderated t test)
#'    \item {p}      matrix (ngene x ncontrast): p    values (moderated t test)
#'    \item {fdr}    matrix (ngene x ncontrast): fdr  values (moderated t test)
#'    \item {bonf}   matrix (ngene x ncontrast): bonf values (moderated t test)
#'    \item {F}      vector (ngene)            : F    values (moderated F test)
#'    \item {F.p}    vector (ngene)            : p    values (moderated F test)
#' }
#' 
#' @param object        SummarizedExperiment
#' @param contrastdefs  contrastdef vector, preferably named (automatically generated names are not always intuitive)
#' @param design        design matrix
#' @return Updated SummarizedExperiment
#' @examples
#' require(magrittr)
#' #if (require(atkin.2014)){
#' #   results <- atkin.2014::soma  %>%  autonomics.find::run_limma()
#' #}
#' 
#' # STEM CELL DIFFERENTIATION
#' if (require(autonomics.data)){
#'    object <- autonomics.data::stemdiff.proteinratios
#'    contrastdefs <- autonomics.find::make_ref_contrasts(object)
#'    object %>% autonomics.find::add_limma(contrastdefs[1:2]) %>% 
#'               extract2('bonf')       %>%
#'               extract(, 'EM01_EM00') %>% 
#'               is_less_than(0.05)     %>% 
#'               sum(na.rm=TRUE)
#'
#' }
#' if (require(billing.differentiation.data)){
#'    contrastdefs <- billing.differentiation.data::contrasts
#'
#'    billing.differentiation.data::rna.counts %>%
#'       autonomics.import::logcpm() %>%
#'       autonomics.find::run_limma(contrastdefs) %>% extract2('bonf') %>%
#'       extract(, 'EM01_0') %>% is_less_than(0.05) %>% sum(na.rm=TRUE)
#'
#'    billing.differentiation.data::rna.counts %>%
#'       autonomics.import::voom_transform() %>%
#'       autonomics.find::run_limma(contrastdefs) %>% extract2('bonf') %>%
#'       extract(, 'EM01_0') %>% is_less_than(0.05) %>% sum(na.rm=TRUE)
#' }
#' if (require(subramanian.2016)){
#'    dir <- system.file('extdata/rnaseq', package = 'subramanian.2016')
#'    object <- autonomics.import::load_rnaseq_counts(dir, sample_file = NULL) %>% 
#'              autonomics.import::voom_transform()
#'    results <- object %>% autonomics.find::run_limma()
#'    results$p[1:3,1:3]
#' }
#' if (require(subramanian.2016)){
#'    object <- subramanian.2016::metabolon
#' }
#'
#' @importFrom  magrittr  %>%
#' @export
add_limma <- function(
   object,
   contrastdefs = autonomics.find::default_contrasts(object),
   design       = autonomics.import::create_design_matrix(object)
){
   # Assert
   assertive.types::assert_is_matrix(design)
   assertive.types::assert_is_numeric(design)
   assertive.base::assert_is_identical_to_true(unname(ncol(object)) == nrow(design))
   if (is.null(contrasts))    return(object)
   contrastdefs %<>% autonomics.find::validify_contrasts(design)
   
   # Return if no contrasts
   if (length(contrastdefs)==0) return(object)
   
   # Set block and correlation if required
   autonomics.support::cmessage('\t\tRun limma')
   block <- NULL
   correlation <- NULL
   my_sdata <- autonomics.import::sdata(object) # works for both eSet and EList
   if (autonomics.import::has_complete_block_values(object)){
      autonomics.support::cmessage("\t\tBlock on svar 'block'")
      block <- my_sdata$block
      correlation <- limma::duplicateCorrelation(autonomics.import::exprs(object),
                                                 design,
                                                 block = block) %>%
                     magrittr::extract2('consensus.correlation')
   }

   # Fit lm and compute contrasts
   contrast_matrix <- autonomics.find::create_contrast_matrix(design, contrastdefs)
   fit <- suppressWarnings( limma::lmFit(object      = autonomics.import::exprs(object),
                                         design      = design,
                                         block       = block,
                                         correlation = correlation, 
                                         weights     = autonomics.import::weights(object))) %>% 
                            limma::contrasts.fit(contrasts = contrast_matrix)
   
   limma_quantities <- if (all(fit$df.residual==0)){ c('effect', 'rank')
                       } else                      { c('effect', 'rank', 't', 'se', 'p', 'fdr', 'bonf') }
   S4Vectors::metadata(object)$limma <- array(dim      =    c(nrow(fit),               ncol(fit),             length(limma_quantities)), 
                                              dimnames = list(feature = rownames(fit), contrast = colnames(fit), quantity = limma_quantities))
   S4Vectors::metadata(object)$limma[, , 'effect']  <- fit$coefficients
   S4Vectors::metadata(object)$limma[, , 'rank'  ]  <- fit$coefficients %>% abs() %>% magrittr::multiply_by(-1) %>% apply(2, rank)
   
   # Perform moderated t test
   if (!all(fit$df.residual==0)){
      fit %<>% limma::eBayes()
      
      S4Vectors::metadata(object)$limma[, , 't'   ] <- fit$t
      S4Vectors::metadata(object)$limma[, , 'se'  ] <- sqrt(fit$s2.post) * fit$stdev.unscaled
      S4Vectors::metadata(object)$limma[, , 'p'   ] <- fit$p.value
      
      S4Vectors::metadata(object)$limma[, , 'fdr' ] <- fit$p.value %>% apply(2, stats::p.adjust, 'fdr')
      S4Vectors::metadata(object)$limma[, , 'bonf'] <- fit$p.value %>% apply(2, stats::p.adjust, 'bonferroni')
      
      autonomics.import::fdata(object)$F   <- fit$F
      autonomics.import::fdata(object)$F.p <- fit$F.p
   }

   # Return
   return(object)
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
#'    object %>% autonomics.find::add_limma() %>% 
#'               autonomics.import::limma()   %>% 
#'               names()
#' }
#' 
#' if (require(subramanian.2016)){
#'    object <- subramanian.2016::metabolon
#'    object %>% autonomics.find::add_limma() %>% autonomics.import::fvars()
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
   design     = autonomics.import::create_design_matrix(object)
){
   .Deprecated('add_limma')
   object %>% autonomics.find::add_limma(contrasts, design)
   return(object)
}