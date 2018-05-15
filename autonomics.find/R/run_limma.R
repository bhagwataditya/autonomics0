
#'@importFrom magrittr  %>%
summarize_limma_analysis <- function(object, fit){
   .Deprecated('fit_contrasts')
  results <- list()
  results$rank     <- - abs(fit$coefficients) %>% apply(2, rank) # Note the - before the abs to ensure ranking in correct order!
  results$coef     <- fit$coefficients
  results$quantile <- fit$coefficients %>% apply(2, function(x) stats::ecdf(x)(x))
  if ('p.value' %in% names(fit)){
    results$p    <- fit$p.value
    results$fdr  <- fit$p.value %>% apply(2, stats::p.adjust, 'fdr')
    results$bonf <- fit$p.value %>% apply(2, stats::p.adjust, 'bonferroni')
    results$rank <- fit$p.value %>% apply(2, rank)
  }
  return(results)
}

run_limma_old <- function(
   object, 
   contrasts = autonomics.find::default_contrasts(object),
   design    = autonomics.find::create_design_matrix(object)
){
  .Deprecated('run_limma')
  # Set block and correlation if required
  autonomics.support::cmessage('\t\tRun limma')
  block <- NULL
  correlation <- NULL
  my_sdata <- autonomics.import::sdata(object) # works for both eSet and EList
  if ('block' %in% names(my_sdata)){
     autonomics.support::cmessage("\t\tBlock on svar 'block'")
     block <- my_sdata$block
     correlation <- limma::duplicateCorrelation(object, design, block = block)$consensus.correlation
  }
  
  # Run limma
  fit <- suppressWarnings(limma::lmFit(
            object      = object,     # it will use the voom weights in an EList object
            design      = design, 
            block       = block, 
            correlation = correlation))
  
  # Analyze contrasts and significances
  contrast_matrix <- autonomics.find::create_contrast_matrix(design, contrasts)
  #results <- fit %>% fit_contrasts_iteratively(contrast_matrix)
  fit %<>% limma::contrasts.fit(contrast_matrix)
  if (!all(fit$df.residual == 0)){
    fit %<>% limma::eBayes()
  }
  results <- summarize_limma_analysis(object, fit)
  
  # Return
  return(results)
}