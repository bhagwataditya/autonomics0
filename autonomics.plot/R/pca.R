#' Principal Component Analysis (PCA)
#'
#' Returns \code{PCA} sample scores, feature loadings, and component variances.
#'
#' \code{PCA} is performed using the (\code{NA}-friendly) \code{NIPALS} algorithm \cr
#' instead of the mathematically elegant, but \code{NA}-unfriendly \code{SVD} formulation. \cr
#' ```
#' SVD:       X = U*D*V' (U = L singular vectors, V = R singular vectors, D = singular values)
#'       Scores = U*D    (or U*sqrt(D)   in SMA)
#'     Loadings = V      (or   sqrt(D)*V in SMA)
#'
#' NIPALS:    X = T*P'   (no further decomposition of scores, so SMA not possible)
#'       Scores = T = U*D
#'     Loadings = P = V
#' ```
#' @param object SummarizedExperiment
#' @param ndim   Number of PCA components to include in results
#' @param ...    included only to keep \code{\link{project}} generic
#' @return
#' ```
#' list(samples  =   sample scores,       # matrix:  nsamples x ndim
#'      features =  feature loadings,     # matrix: nfeatures x ndim
#'      var      = variance percentages   # vector:         1 x ndim
#' ```
#' @examples
#' if (require(autonomics.data)){
#'
#'    # STEM CELL DIFFERENTIATION
#'       require(magrittr)
#'       object <- 'extdata/stemdiff/maxquant/proteinGroups.txt' %>%
#'                  system.file(package='autonomics.data')       %>%
#'                  autonomics.import::read_proteingroups()
#'       object %>% pca() %>% str()
#'
#'    # GLUTAMINASE
#'       object <- 'extdata/glutaminase/glutaminase.xlsx'    %>%
#'                  system.file(package = 'autonomics.data') %>%
#'                  autonomics.import::read_metabolon()
#'       object %>% pca() %>% str()
#' }
#' @seealso \code{\link[autonomics.plot]{sma}},
#'          \code{\link[autonomics.plot]{lda}},
#'          \code{\link[autonomics.plot]{pls}}
#' @author Aditya Bhagwat
#' @md
#' @importFrom magrittr %>%
#' @export
pca <- function(object, ndim = 2, ...){

   # Convert NaN into NA (otherwise pcaMethods::pca fails)
   autonomics.import::exprs(object) %<>% (function(x){x[is.nan(x)|is.infinite(x)] <- NA_real_; x})

   # Rm features which are NA in all samples
   object %<>% autonomics.preprocess::filter_features_nonzero_in_some_sample()

   # (Double) center and (global) normalize
   raw_x <- autonomics.import::exprs(object)
   row_means <- raw_x %>% rowMeans(na.rm=TRUE)
   col_means <- raw_x %>% matrixStats::colWeightedMeans(abs(row_means), na.rm = TRUE)
   global_mean <- mean(col_means)
   x <- raw_x %>%
        apply(1, '-', col_means)                %>%   # Center columns
        apply(1, '-', row_means)                %>%   # Center rows
        magrittr::add(global_mean)              %>%   # Add doubly subtracted component
        magrittr::divide_by(stats::sd(., na.rm=TRUE)) # Normalize

   # Perform PCA
   pca_res  <- pcaMethods::pca(t(x), nPcs = ndim, scale = 'none', center = FALSE, method = 'nipals')
   samples  <- pca_res@scores   %>% magrittr::set_colnames(sprintf('pca%d', 1:ncol(.)))
   features <- pca_res@loadings %>% magrittr::set_colnames(sprintf('pca%d', 1:ncol(.)))
   var <- round(100*pca_res@R2) %>% magrittr::set_names(   sprintf('pca%d', 1:length(.)))

   # Return
   list(samples  = samples  %>% magrittr::extract(,1:ndim, drop = FALSE), # drop = FALSE required when ndim=1!
        features = features %>% magrittr::extract(,1:ndim, drop = FALSE), # drop = FALSE required when ndim=1!
        var      = var      %>% magrittr::extract( 1:ndim))
}


#' Spectral Map Analysis (SMA)
#'
#' Returns \code{SMA} sample scores, feature loadings, and component variances.
#'
#' \code{SMA} is based on the (\code{NA}-unfriendly) \code{SVD}, so only \code{NA}-free features are used. \cr
#' \code{SMA} and \code{PCA} are related: their \code{SVD} is identical, but score and loading computation differ:
#' ```
#'    PCA Scores  : U * D
#'    SMA Scores  : U * D^(1/2)
#'
#'    PCA Loadings:         V
#'    SMA Loadings: D^(1/2)*V
#' ```
#' @param object     SummarizedExperiment
#' @param na.impute  whether to \code{\link[autonomics.preprocess]{impute}} prior to \code{SMA}
#' @param ndim       Number of SMA components to include in results
#' @param ...        used to keep \code{\link{project}} generic
#' @return
#' ```
#' list(samples  =   sample scores,       # matrix:  nsample x ndim
#'      features =  feature loadings,     # matrix: nfeature x ndim
#'      var      = variance percentages   # vector:        1 x ndim
#' ```
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'
#'    # STEM CELL DIFFERENTIATION
#'       require(magrittr)
#'       object <- 'extdata/stemdiff/maxquant/proteinGroups.txt' %>%
#'                  system.file(package='autonomics.data')       %>%
#'                  autonomics.import::read_proteingroups()
#'       object %>% sma() %>% str()
#'
#'    # GLUTAMINASE
#'       object <- 'extdata/glutaminase/glutaminase.xlsx'    %>%
#'                  system.file(package = 'autonomics.data') %>%
#'                  autonomics.import::read_metabolon()
#'       object %>% sma() %>% str()
#' }
#' @references
#' Wouters et al (2003) Graphical exploration of gene expression data: a comparative study of
#' three multivariate methods.
#' @seealso \code{\link[autonomics.plot]{pca}},
#'          \code{\link[autonomics.plot]{lda}},
#'          \code{\link[autonomics.plot]{pls}}
#' @author Aditya Bhagwat, Laure Cougnaud
#' @md
#' @importFrom magrittr %>%
#' @export
sma <- function(object, na.impute = FALSE, ndim = 2, ...){

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
   if(ncomponents < ndim){
      stop('\'ndim\' = \'', ndim, '\', but only \'', ncomponents, 'can be provided.')
   }
   npairs  <- ncomponents/2
   mpm_out <- split(1:ncomponents, rep(1:npairs, each = 2)) %>%
              lapply(function(x){
                        y <- suppressPackageStartupMessages(mpm::plot.mpm(mpm_tmp, do.plot = FALSE, dim = x))
                        list(features = y$Rows    %>% magrittr::extract(, c('X', 'Y')),
                             samples  = y$Columns %>% magrittr::extract(, c('X', 'Y')))})

   # Extract
   samples  <- mpm_out %>% lapply(magrittr::extract2, 'samples')  %>% do.call(cbind, .) %>% magrittr::set_colnames(sprintf('sma%d', 1:ncol(.)))
   features <- mpm_out %>% lapply(magrittr::extract2, 'features') %>% do.call(cbind, .) %>% magrittr::set_colnames(sprintf('sma%d', 1:ncol(.)))
   percentages <- round(100*mpm_tmp$contrib[1:ncomponents])                             %>% magrittr::set_names(   sprintf('sma%d', 1:length(.)))

   # Return
   list(samples  = samples     %>% magrittr::extract(,1:ndim, drop = FALSE), # drop = FALSE required when ndim=1!
        features = features    %>% magrittr::extract(,1:ndim, drop = FALSE), # drop = FALSE required when ndim=1!
        var      = percentages %>% magrittr::extract( 1:ndim))
}


#' Linear discriminant analysis (LDA)
#'
#' Returns \code{LDA} sample scores, feature loadings, and discriminant variances.
#'
#' \code{LDA} relies on the (\code{NA}-unfriendly) \code{SVD}, so only \code{NA}-free features are included in the analysis.
#'
#' @param object     SummarizedExperiment
#' @param na.impute  currently ignored. May be used in future.
#' @param ndim       Number of linear discriminants to include in results
#' @param ...        only included to keep \code{\link[autonomics.plot]{project}} generic
#' @return
#' ```
#' list(samples  =   sample scores,       # matrix:  nsample x ndim
#'      features =  feature loadings,     # matrix: nfeature x ndim
#'      var      = variance percentages   # vector:        1 x ndim
#' ```
#' @examples
#' if (require(autonomics.data)){
#'    require(magrittr)
#'
#'    # STEM CELL COMPARISON
#'       object <- 'extdata/stemcomp/maxquant/proteinGroups.txt' %>%
#'                  system.file(package = 'autonomics.data')     %>%
#'                  autonomics.import::read_proteingroups()
#'       object %>% autonomics.plot::lda() %>% str()
#'       \dontrun{ # Fails, as max dim length(subgroup) - 1
#'          object %>% autonomics.plot::lda(ndim = 3) %>% str()
#'        }
#'
#'    # STEM CELL DIFFERENTIATION
#'       require(magrittr)
#'       object <- 'extdata/stemdiff/maxquant/proteinGroups.txt' %>%
#'                  system.file(package='autonomics.data')       %>%
#'                  autonomics.import::read_proteingroups()
#'       object %>% lda() %>% str()
#'
#'    # GLUTAMINASE
#'       object <- 'extdata/glutaminase/glutaminase.xlsx'    %>%
#'                  system.file(package = 'autonomics.data') %>%
#'                  autonomics.import::read_metabolon()
#'       object %>% lda() %>% str()
#' }
#' @seealso \code{\link[autonomics.plot]{pca}},
#'          \code{\link[autonomics.plot]{sma}},
#'          \code{\link[autonomics.plot]{pls}}
#' @author Aditya Bhagwat, Laure Cougnaud
#' @md
#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%
#' @export
lda <- function(object, na.impute = FALSE, ndim = 2,  ...){

   # Assert
   nsubgroup <- object %>% autonomics.import::subgroup_levels() %>% length()
   assertive.numbers::assert_all_are_greater_than(nsubgroup, 1)
   if (ndim > (nsubgroup-1)){
      stop(sprintf('LDA requires ndim (%d) <= nsubgroup-1 (%d)', ndim, nsubgroup-1))
   }

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

   # Return
   list(samples  = samples     %>% magrittr::extract(,1:ndim, drop = FALSE), # drop = FALSE required when ndim=1!
        features = features    %>% magrittr::extract(,1:ndim, drop = FALSE), # drop = FALSE required when ndim=1!
        var      = percentages %>% magrittr::extract( 1:ndim))

}

#' Partial least squares analysis (PLS)
#'
#' Returns \code{PLS} sample scores, feature loadings, and latent variable variances.
#'
#' \code{PLS} is performed using the (\code{NA}-friendly) \code{NIPALS} algorithm.
#'
#' @param object SummarizedExperiment
#' @param ndim   Number of latent variables to include in results
#' @param implementation 'mixOmics::plsda', 'mixOmics::splsda', or 'ropls::opls'
#' @param ...    only inlcuded to keep \code{\link[autonomics.plot]{project}} generic
#' @return
#' ```
#' list(samples  =   sample scores,       # matrix:  nsample x ndim
#'      features =  feature loadings,     # matrix: nfeature x ndim
#'      var      = variance percentages   # vector:        1 x ndim
#' ```
#' @examples
#' if (require(autonomics.data)){
#'    require(magrittr)
#'
#'    # STEM CELL DIFFERENTIATION
#'    object <- 'extdata/stemdiff/maxquant/proteinGroups.txt' %>%
#'               system.file(package = 'autonomics.data') %>%
#'               autonomics.import::read_proteingroups()
#'    object %>% autonomics.plot::pls() %>% str()
#'
#'   # GLUTAMINASE
#'    object <- 'extdata/glutaminase/glutaminase.xlsx' %>%
#'               system.file(package='autonomics.data') %>%
#'               autonomics.import::read_metabolon()
#'    object %>% autonomics.plot::pls() %>% str()
#'    \dontrun{ # slow
#'       object %>% autonomics.plot::pls(implementation = 'ropls::opls') %>% str()
#'    }
#' }
#' @seealso \code{\link[autonomics.plot]{pca}},
#'          \code{\link[autonomics.plot]{sma}},
#'          \code{\link[autonomics.plot]{lda}}
#' @author Aditya Bhagwat
#' @md
#' @importFrom magrittr %>%
#' @export
pls <- function(object, implementation = NULL, ndim = 2, ...){

   # Assert
   assertive.base::assert_is_identical_to_true(length(autonomics.import::subgroup_levels(object)) > 1)

   x <- t(autonomics.import::exprs(object))
   y <- autonomics.import::subgroup_values(object)

   # mixOmics implementation
   if (is.null(implementation))  implementation <- 'mixOmics::plsda'   # much faster than the ropls implementation
   if (implementation %in% c('mixOmics::plsda', 'mixOmics::splsda')){
      # Project
      pls_out <- if (implementation=='mixOmics::plsda'){
         mixOmics::plsda( x, y, ncomp = ndim)
      } else {
         mixOmics::splsda(x, y, ncomp = ndim)
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
      pls_out <- ropls::opls(x, y, predI = ndim, permI = 0, plotL = FALSE)

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
   list(samples  = samples  %>% magrittr::extract(,1:ndim, drop = FALSE), # drop = FALSE required when ndim=1!
        features = features %>% magrittr::extract(,1:ndim, drop = FALSE), # drop = FALSE required when ndim=1!
        var      = var      %>% magrittr::extract( 1:ndim))
}

#' Project
#'
#' Return sample scores, feature loadings, and component variances.
#'
#' @param object          SummarizedExperiment
#' @param method         'pca', 'sma', 'lda', 'pls'
#' @param ndim            Number of 'dimensions' (e.g. PCA: dimensions == components) to return
#' @param implementation  which implementation of the method to use
#' @param na.impute       TRUE or FALSE
#' @return list(samples, features, var)
#' @examples
#' if (require(autonomics.data)){
#'    require(magrittr)
#'
#'    # STEM CELL DIFFERENTIATION
#'    object <- 'extdata/stemdiff/maxquant/proteinGroups.txt' %>%
#'               system.file(package = 'autonomics.data')     %>%
#'               autonomics.import::read_proteingroups()
#'    object %>% autonomics.plot::project() %>% str()
#'
#'    # GLUTAMINASE
#'    object <- 'extdata/glutaminase/glutaminase.xlsx'    %>%
#'               system.file(package = 'autonomics.data') %>%
#'               autonomics.import::read_metabolon()
#'    object %>% autonomics.plot::project() %>% str()
#' @seealso \code{\link[autonomics.plot]{pca}},
#'          \code{\link[autonomics.plot]{sma}},
#'          \code{\link[autonomics.plot]{lda}},
#'          \code{\link[autonomics.plot]{pls}}
#' @author Aditya Bhagwat
#' @export
project <- function(object, method = 'pca', na.impute = FALSE, implementation = NULL, ndim = 2){
   utils::getFromNamespace(method, 'autonomics.plot')(object,
                                                      na.impute      = na.impute,
                                                      implementation = implementation,
                                                      ndim           = ndim)
}
