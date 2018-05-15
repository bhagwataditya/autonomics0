#' @rdname project
#' @importFrom magrittr %>% 
#' @export
pca <- function(object, na.impute = FALSE, dims = 1:2, ...){
   
   # Rm all NA features
   object %<>% autonomics.preprocess::filter_features_nonzero_in_some_sample()
   
   # (Double) center and (global) normalize
   raw_x <- autonomics.import::exprs(object)
   row_means <- raw_x %>% rowMeans(na.rm=TRUE)
   col_means <- raw_x %>% matrixStats::colWeightedMeans(abs(row_means), na.rm = TRUE)
   global_mean <- mean(col_means)
   x <- raw_x %>% 
        apply(1, '-', col_means)                %>%  # Center columns
        apply(1, '-', row_means)                %>%  # Center rows
        magrittr::add(global_mean)              %>%  # Add doubly subtracted component
        magrittr::divide_by(stats::sd(., na.rm=TRUE))       # Normalize
   
   # SVD:       X = U*D*V' (U = L singular vectors, V = R singular vectors, D = singular values)
   #       Scores = U*D    (or U*sqrt(D)   in SMA)
   #     Loadings = V      (or   sqrt(D)*V in SMA)
   #                
   # NIPALS:    X = T*P'   (no further decomposition of scores, so SMA not possible)
   #       Scores = T = U*D
   #     Loadings = P = V
   pca_res  <- pcaMethods::pca(t(x), nPcs = max(dims, na.rm = TRUE), scale = 'none', center = FALSE, method = 'nipals')
   samples  <- pca_res@scores   %>% magrittr::set_colnames(sprintf('pca%d', 1:ncol(.)))
   features <- pca_res@loadings %>% magrittr::set_colnames(sprintf('pca%d', 1:ncol(.)))
   var <- round(100*pca_res@R2) %>% magrittr::set_names(   sprintf('pca%d', 1:length(.)))
   list(samples  = samples %>%
           magrittr::extract(,dims),
        features = features %>%
           magrittr::extract(,dims),
        var      = var %>%
           magrittr::extract(dims))
   
   #df=pca_res@scores %>% data.frame(sample = rownames(.), subgroup = rownames(.) %>% substr(1, nchar(.)-3))
   #ggplot2::ggplot(df) + 
   #   ggplot2::geom_point(ggplot2::aes(x=PC1, y=PC2, color = subgroup), shape = 15, size = 3) + 
   #   ggplot2::scale_color_manual(values = color_values)
   #color_values   = autonomics.plot::default_color_values(object, 'subgroup')
   
}
   

#' @rdname project
#' @importFrom magrittr %>% 
#' @export
sma <- function(object, na.impute = FALSE, dims = 1:2, ...){
   
   # Preprocess
   object %<>% autonomics.preprocess::minusinf_to_na()
   object %<>% autonomics.preprocess::flip_sign_if_all_exprs_are_negative() # else SVD singular 
   if (na.impute) object %<>% autonomics.preprocess::impute()
   object %<>% autonomics.preprocess::filter_features_available_in_all_samples()
   
   # Transform
   mpm_tmp <- data.frame(feature = rownames(object), autonomics.import::exprs(object)) %>% 
              mpm::mpm(logtrans = FALSE, closure = 'none',  center = 'double',
                       normal = 'global', row.weight = 'mean', col.weight = 'constant')
   ncomponents <- (100*mpm_tmp$contrib) %>% magrittr::is_greater_than(1) %>% sum() %>% autonomics.support::evenify_upwards()
   if(ncomponents < max(dims,na.rm = TRUE))
   {
      stop('\'dims\' includes \'', max(dims, na.rm = TRUE), '\', but only \'', ncomponents, 'can be provided.')
   }
   npairs      <-  ncomponents/2
   mpm_out <-  split(1:ncomponents, rep(1:npairs, each = 2)) %>% 
      lapply(function(x){
         y <- suppressPackageStartupMessages(mpm::plot.mpm(mpm_tmp, do.plot = FALSE, dim = x))
         list(features = y$Rows    %>% magrittr::extract(, c('X', 'Y')), 
              samples  = y$Columns %>% magrittr::extract(, c('X', 'Y')))
      })
   
   # Extract
   samples  <- mpm_out %>% lapply(magrittr::extract2, 'samples')  %>% do.call(cbind, .) %>% magrittr::set_colnames(sprintf('sma%d', 1:ncol(.)))
   features <- mpm_out %>% lapply(magrittr::extract2, 'features') %>% do.call(cbind, .) %>% magrittr::set_colnames(sprintf('sma%d', 1:ncol(.)))
   percentages <- round(100*mpm_tmp$contrib[1:ncomponents])                             %>% magrittr::set_names(   sprintf('sma%d', 1:length(.)))
   
   # Return
   list(samples  = samples %>%
           magrittr::extract(,dims),
        features = features %>%
           magrittr::extract(,dims),
        var      = percentages %>%
           magrittr::extract(dims))
}

#' @rdname project
#' @importFrom magrittr %>% 
#' @importFrom magrittr %<>% 
#' @export
lda <- function(
   object,
   na.impute = FALSE,
   dims      = 1:2,
   ...)
{
   # Assert
   assertive.numbers::assert_all_are_greater_than(
      object %>%
         autonomics.import::subgroup_levels() %>%
         length(),
      1)

   max_dim <- object %>%
      autonomics.import::subgroup_levels() %>%
      length() %>%
      magrittr::subtract(1)
   
   legal_dims_subsetter <- assertive.base::assert_all_are_true(
      dims %>%
         magrittr::is_weakly_less_than(max_dim))
   
   dims %<>%
      magrittr::extract(legal_dims_subsetter)

   # Preprocess
   object %<>% autonomics.preprocess::minusinf_to_na()
   object %<>% autonomics.preprocess::flip_sign_if_all_exprs_are_negative() # else SVD singular 
   if (na.impute) object %<>% autonomics.preprocess::impute()
   object %<>% autonomics.preprocess::filter_features_available_in_all_samples()
   
   # Transform
   exprs_t  <- autonomics.import::exprs(object) %>% t()
   lda_out  <- suppressWarnings(MASS::lda(exprs_t, grouping = autonomics.import::sdata(object)$subgroup))
   features <- lda_out$scaling 
   if (ncol(features)==1) features %<>% cbind(LD2 = 0)
   samples  <- scale(exprs_t, center = colMeans(lda_out$means), scale = FALSE) %*% features
   percentages <- round((lda_out$svd^2)/sum(lda_out$svd^2)*100)
   if (length(percentages)==1) percentages <- c(LD1 = percentages, LD2 = 0)
   
   # Rename
   features    %<>% (function(x){colnames(x) %<>% stringi::stri_replace_first_fixed('LD', 'lda'); x})
   samples     %<>% (function(x){colnames(x) %<>% stringi::stri_replace_first_fixed('LD', 'lda'); x})
   percentages %<>% (function(x){names(x) <- sprintf('lda%d', 1:ncol(samples)); x})

   # Subset
   features %<>%
      magrittr::extract(,dims)
   samples %<>%
      magrittr::extract(,dims)
   percentages %<>%
      magrittr::extract(dims)

   # Return
   list(samples = samples, features = features, var = percentages)
}

#' @rdname project
#' @importFrom magrittr %>% 
#' @export
pls <- function(
   object, 
   na.impute = FALSE, 
   implementation = NULL,
   dims = 1:2,
   ...
){
   
   # Assert
   assertive.base::assert_is_identical_to_true(length(autonomics.import::subgroup_levels(object)) > 1)
   
   x <- t(autonomics.import::exprs(object))
   y <- autonomics.import::subgroup_values(object)
   
   # mixOmics implementation
   if (is.null(implementation))  implementation <- 'mixOmics::plsda'   # much faster than the ropls implementation
   if (implementation %in% c('mixOmics::plsda', 'mixOmics::splsda')){
      # Project
      pls_out <- if (implementation=='mixOmics::plsda'){
         mixOmics::plsda( x, y, ncomp = max(dims, na.rm = TRUE)) 
      } else {
         mixOmics::splsda(x, y, ncomp = max(dims, na.rm = TRUE))
      }
      # Extract
      samples  <- pls_out$variates$X
      features <- pls_out$loadings$X
      var      <- round(100*pls_out$explained_variance$X)
      
      # Rename
      colnames(samples)  %<>% stringi::stri_replace_first_fixed('comp ', 'pls')
      colnames(features) %<>% stringi::stri_replace_first_fixed('comp ', 'pls')
      names(var)         %<>% stringi::stri_replace_first_fixed('comp ', 'pls')
      
   }
   
   # ropls implementation
   if (implementation == 'ropls::opls'){
      
      # Project
      pls_out <- ropls::opls(x, y, predI = max(dims, na.rm = TRUE), permI = 0, plotL = FALSE)
      
      # Extract
      samples  <- pls_out@scoreMN
      features <- pls_out@loadingMN
      var      <- round(pls_out@modelDF$R2X*100)
      
      # Rename
      colnames(samples)  %<>% stringi::stri_replace_first_fixed('p ', 'pls')
      colnames(features) %<>% stringi::stri_replace_first_fixed('p ', 'pls')
      names(var) <- c('pls1', 'pls2')
   }
   
   # Return
   return(
      list(
         samples  = samples %>%
            magrittr::extract(,dims),
         features = features %>%
            magrittr::extract(,dims),
         var      = var %>%
            magrittr::extract(dims)))
}

#' PCA, SMA, LDA, or PLS transform
#' 
#' Perform principal component analysis (PCA), spectral map analysis (SMA), 
#' linear discriminant analysis (LDA), or partial least squares analysis (PLS).
#' 
#' Compute SMA and LDA using SVD (utilizable for all methods, but restricted to non NA features).
#' Compute PCA and PLS using NIPALS (utilizable for PCA and PLS only, but allows features with NAS).
#' 
#' PCA and SMA are related: their SVD is identical, but score and loading computation differs 
#' (Wouters et al. 2003 - see references) \cr
#'    PCA Scores: U * D                              \cr
#'    SMA Scores: U * D^(1/2)                        \cr 
#'                                                   \cr
#'    PCA Loadings: V                                \cr
#'    SMA Loadings: D^(1/2)*V                        \cr
#'    
#' @param object          SummarizedExperiment
#' @param method         'pca', 'sma', 'lda', 'pls'
#' @param dims            \code{\link{numeric}} vector indicating what 'dimentions' (e.g. PCA: dimensions == components) to return
#' @param implementation  which implementation of the method to use
#' @param dims_conflict   how to react to a conflict between what \code{lda} can deliver (\code{dims = 1:2}) and what dims are requested
#' @param na.impute       TRUE or FALSE
#' @param ...             further parameters
#' @return list(samples, features, var)
#' @examples 
#' require(magrittr)
#' if (require(autonomics.data)){
#'    object <- autonomics.data::billing2016
#'    object %>% autonomics.explore::pca() %>% str()
#'    object %>% autonomics.explore::pca(dims = 1:6) %>% str()
#'    object %>% autonomics.explore::project(method = 'lda') %>% str()
#' \dontrun{
#'    # Fails, as max dim length(subgroup) - 1
#'    object %>%
#'      autonomics.explore::project(method = 'lda', dims = 1:3) %>%
#'      str()}
#' }
#' if (require(subramanian.2016)){
#'    object <- subramanian.2016::metabolon
#'    object %>% autonomics.explore::pca()  %>%  str()
#'    object %>% autonomics.explore::sma()  %>%  str()
#'    
#'    object %>% autonomics.explore::lda()  %>%  str()
#'    object %>% autonomics.import::filter_samples(
#'                  subgroup %in% autonomics.import::subgroup_levels(.)[1:2]) %>% 
#'               autonomics.explore::lda()
#'               
#'    object %>% autonomics.explore::pls()  %>%  str()
#'    object %>% autonomics.explore::pls(implementation = 'mixOmics::splsda') %>% str()
#'    \dontrun{object %>% autonomics.explore::pls(implementation = 'ropls::opls'     ) %>% str()}
#' }
#' if (require(halama.2016)){
#'    object <- halama.2016::cell.metabolites
#'    object %>% autonomics.explore::project('lda') %>% str()
#'    object %>% autonomics.explore::project('pca') %>% str()
#' }
#' @references
#' PCA 
#' SMA: Wouters et al (2003) Graphical exploration of gene expression data: a comparative study of 
#'      three multivariate methods
#' @export
project <- function(object, method = 'pca', na.impute = FALSE, implementation = NULL, dims = 1:2, dims_conflict = 'warn'){
   get(method)(object, na.impute = na.impute, implementation = implementation, dims = dims)
}
