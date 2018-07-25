fit_contrasts <- function(fit, contrasts, design){

   contrast_matrix <- autonomics.find::create_contrast_matrix(design, contrasts)

   fit_single_contrast <- function(i.contrast){
      cur.contrast.matrix <- contrast_matrix[, i.contrast, drop = FALSE]
      cur.contrast.matrix %<>% magrittr::extract(cur.contrast.matrix!=0, , drop = FALSE)
      cur.subgroups <- rownames(cur.contrast.matrix)[cur.contrast.matrix != 0]
      cur.fit <- fit %>% magrittr::extract(, cur.subgroups)
      cur.fit %<>% limma::contrasts.fit(cur.contrast.matrix)
      if (!all(fit$df.residual == 0)){
         cur.fit %<>% limma::eBayes()
      }
      cur.fit
   }

   contrast_results <- lapply(seq(1,ncol(contrast_matrix)), fit_single_contrast)

   results <- list()
   results$coef     <- contrast_results %>% lapply(function(x)  x$coefficients) %>% do.call(cbind, .)
   results$rank     <- apply(-abs(results$coef), 2, rank)
   results %<>% magrittr::extract(c('rank', 'coef'))
   results$quantile <- results$coef %>% apply(2, function(x) stats::ecdf(x)(x))
   results$p        <- contrast_results %>% lapply(function(x){
                                                      if ('p.value' %in% names(x)){
                                                         x$p.value
                                                      } else {
                                                         matrix(NA, nrow = nrow(x$coef), ncol = ncol(x$coef), dimnames = dimnames(x$coef))
                                                      }
                                                   }) %>% do.call(cbind, .)
   if (all(is.na(results$p))){
      results$p <- NULL
   } else {
      results$fdr      <- results$p %>% apply(2, stats::p.adjust, 'fdr')
      results$bonf     <- results$p %>% apply(2, stats::p.adjust, 'bonferroni')
      results$rank     <- results$p %>% apply(2, rank)
   }
   results
}


#' Run limma
#'
#' Compare contrasts of object$subgroup levels
#'
#' @param object  eset
#' @param contrasts contrast vector, preferably named (automatically generated names are not always intuitive)
#' @param design    design matrix
#' @return list with four components: rank, fc, p, fdr, bonf
#' @examples
#' require(magrittr)
#' #if (require(atkin.2014)){
#' #   results <- atkin.2014::soma  %>%  autonomics.find::run_limma()
#' #}
#' 
#' # STEM CELL DIFFERENTIATION
#' if (require(autonomics.data)){
#'    object <- autonomics.data::stemdiff.proteinratios
#'    contrasts <- autonomics.find::make_ref_contrasts(object)
#'    object %>% autonomics.find::run_limma(contrasts[1:2]) %>% 
#'               extract2('bonf')       %>%
#'               extract(, 'EM01_EM00') %>% 
#'               is_less_than(0.05)     %>% 
#'               sum(na.rm=TRUE)
#'
#' }
#' if (require(billing.differentiation.data)){
#'    contrasts <- billing.differentiation.data::contrasts
#'
#'    billing.differentiation.data::rna.counts %>%
#'       autonomics.import::logcpm() %>%
#'       autonomics.find::run_limma(contrasts) %>% extract2('bonf') %>%
#'       extract(, 'EM01_0') %>% is_less_than(0.05) %>% sum(na.rm=TRUE)
#'
#'    billing.differentiation.data::rna.counts %>%
#'       autonomics.import::voom_transform() %>%
#'       autonomics.find::run_limma(contrasts) %>% extract2('bonf') %>%
#'       extract(, 'EM01_0') %>% is_less_than(0.05) %>% sum(na.rm=TRUE)
#' }
#' if (require(subramanian.2016)){
#'    dir <- system.file('extdata/rnaseq', package = 'subramanian.2016')
#'    object <- autonomics.import::load_rnaseq_counts(dir, sample_file = NULL) %>% 
#'              autonomics.import::voom_transform()
#'    results <- object %>% autonomics.find::run_limma()
#'    results$p[1:3,1:3]
#' }
#'
#' @importFrom  magrittr  %>%
#' @export
run_limma <- function(
   object,
   contrasts = autonomics.find::default_contrasts(object),
   design    = autonomics.find::create_design_matrix(object)
){

   # Set block and correlation if required
   autonomics.support::cmessage('\t\tRun limma')
   block <- NULL
   correlation <- NULL
   my_sdata <- autonomics.import::sdata(object) # works for both eSet and EList
   if (autonomics.import::contains_block(object)){
      autonomics.support::cmessage("\t\tBlock on svar 'block'")
      block <- my_sdata$block
      correlation <- limma::duplicateCorrelation(autonomics.import::exprs(object),
                                                 design,
                                                 block = block) %>%
                     magrittr::extract2('consensus.correlation')
   }

   # Run limma
   fit <- suppressWarnings(limma::lmFit(
      object      = autonomics.import::exprs(object),     # it will use the voom weights in an EList object
      design      = design,
      block       = block,
      correlation = correlation, 
      weights     = autonomics.import::weights(object)))

   # Analyze contrasts and significances
   results <- fit %>% fit_contrasts(contrasts, design)

   # Return
   return(results)
}
