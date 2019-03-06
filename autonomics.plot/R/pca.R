#=====================
# PROJECT
#=====================


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
#'    require(magrittr)
#'    object <- autonomics.data::glutaminase
#'    object %>% pca() %>% str()
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
#' @param na.impute  logical: \code{\link[autonomics.preprocess]{impute}} prior to \code{SMA}?
#' @param ndim       number: how many SMA components to include in results
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
#'    object <- autonomics.data::glutaminase
#'    object %>% sma() %>% str()
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
   object %<>% autonomics.import::minusinf_to_na()
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
#' @param na.impute  logical: \code{\link[autonomics.preprocess]{impute}} prior to \code{LDA}?
#' @param ndim       number: how many linear discriminants to include in results
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
#'       object <- autonomics.data::stemcomp.proteinratios
#'       object %>% lda() %>% str()
#'       \dontrun{ # Fails, as max dim length(subgroup) - 1
#'          object %>% lda(ndim = 3) %>% str()
#'        }
#'
#'    # GLUTAMINASE
#'       object <- autonomics.data::glutaminase
#'       object %>% lda() %>% str()
#' }
#' @seealso \code{\link[autonomics.plot]{pca}},
#'          \code{\link[autonomics.plot]{sma}},
#'          \code{\link[autonomics.plot]{pls}}
#' @author Aditya Bhagwat, Laure Cougnaud
#' @md
#' @importFrom magrittr %>% %<>%
#' @export
lda <- function(object, na.impute = FALSE, ndim = 2,  ...){

   # Assert
   nsubgroup <- object %>% autonomics.import::subgroup_levels() %>% length()
   assertive.numbers::assert_all_are_greater_than(nsubgroup, 1)
   if (ndim > (nsubgroup-1)){
      stop(sprintf('LDA requires ndim (%d) <= nsubgroup-1 (%d)', ndim, nsubgroup-1))
   }

   # Preprocess
   object %<>% autonomics.import::minusinf_to_na()
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
#' @param object          SummarizedExperiment
#' @param ndim            number: how many latent variables to include in results
#' @param implementation 'mixOmics::plsda', 'mixOmics::splsda', or 'ropls::opls'
#' @param ...             only inlcuded to keep \code{\link[autonomics.plot]{project}} generic
#' @return
#' ```
#' list(samples  =   sample scores,       # matrix:  nsample x ndim
#'      features =  feature loadings,     # matrix: nfeature x ndim
#'      var      = variance percentages   # vector:        1 x ndim
#' ```
#' @examples
#' if (require(autonomics.data)){
#'    require(magrittr)
#'    object <- autonomics.data::glutaminase
#'    object %>% pls() %>% str()
#'    \dontrun{ # slow
#'       object %>% pls(implementation = 'ropls::opls') %>% str()
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
#'    object <- autonomics.data::glutaminase
#'    object %>% project() %>% str()
#' }
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



#==================================================
# PLOT PROJECTED SAMPLES
#==================================================

#' Make datatable to plot PCA/PLS/LDA/SMA scores
#' @param projection  output from PCA/PLS/LDA/SMA
#' @param object      SummarizedExperiment
#' @param dims        PCA/LDA dimensions
#' @param color_var   svar mapped to color
#' @param shape_var   svar mapped to shape
#' @param size_var    svar mapped to size
#' @param txt_var     svar mapped to txt
#' @param split_var   svar on which object is splitted prior to lda analysis
#' @param facet_var   svar on which plot is faceted
#' @param group_var   svar mapped to group_var
#' @examples
#' if (require(autonomics.data)){
#'    require(magrittr)
#'    object <- autonomics.data::glutaminase
#'
#'    object %>% pca() %>%
#'               make_projected_samples_df(object,
#'                                         dims = c(1,2),
#'                                         color_var = 'subgroup',
#'                                         shape_var = NULL,
#'                                         size_var  = NULL,
#'                                         txt_var   = NULL,
#'                                         split_var = NULL,
#'                                         facet_var = NULL,
#'                                         group_var = NULL) %>%
#'               print()
#'
#'    object %>% pca(ndim=4) %>%
#'               make_projected_samples_df(object,
#'                                         dims      = 3:4,
#'                                         color_var = 'subgroup',
#'                                         shape_var = NULL,
#'                                         size_var  = NULL,
#'                                         txt_var   = NULL,
#'                                         split_var = NULL,
#'                                         facet_var = NULL,
#'                                         group_var = NULL) %>%
#'               print()
#' }
#' @importFrom data.table  data.table  :=
#' @importFrom magrittr    %<>%
#' @export
make_projected_samples_df <- function(
   projection,
   object,
   dims,
   color_var,
   shape_var,
   size_var,
   txt_var,
   split_var,
   facet_var,
   group_var
){
   # Satisfy CHECK
   color <- shape <- size <- label <- group <- facet <- NULL

   # Add MPA coordinates
   DT <- data.table::data.table(
            sample_id = autonomics.import::snames(object),
            x         = projection$samples[, dims[1]],
            y         = projection$samples[, dims[2]])

   # Add other aesthetics
   DT %<>% magrittr::extract(, color := if(is.null(color_var)) 'default'  else  autonomics.import::sdata(object)[[color_var]] )
   DT %<>% magrittr::extract(, shape := if(is.null(shape_var)) 'default'  else  autonomics.import::sdata(object)[[shape_var]] )
   DT %<>% magrittr::extract(, size  := if(is.null(size_var))  'default'  else  autonomics.import::sdata(object)[[size_var ]] )
   if (!is.null(txt_var))   DT %<>% magrittr::extract(, label := autonomics.import::sdata(object)[[txt_var]])
   if (!is.null(group_var)) DT %<>% magrittr::extract(, group := autonomics.import::sdata(object)[[group_var]])
   if (!is.null(facet_var)) DT %<>% magrittr::extract(, facet := autonomics.import::sdata(object)[[facet_var]])
   if (!is.null(split_var)) DT %<>% magrittr::extract(, facet := sprintf('%s %d%% %d%%', facet, round(projection$var[1]), round(projection$var[2])))
   DT
}




#' Plot PCA/LDA/PLS/SMA sample scores
#' @param object            SummarizedExperiment
#' @param ...               Additional SummarizedExperiment(s)
#' @param listed_objects    Additional SummarizedExperiment(s), wrapped in a \code{link{list}}
#' @param method            'pca', 'lda', 'pls', 'sma'
#' @param implementation    which implementation of pls to use (see \code{\link{project}})
#' @param dims              dimensions
#' @param color_var         svar mapped to color
#' @param color_values      color values vector
#' @param shape_var         svar mapped to shape
#' @param size_var          svar mapped to size
#' @param txt_var           svar mapped to txt
#' @param txt_outliers      logical(1): whether to annotate outliers
#' @param group_var         svar mapped to group
#' @param split_var         svar on which to split data prior to transformation
#' @param facet_var         svar on which to facet sample biplot
#' @param scales            'free' or 'fixed'
#' @param labs              named \code{\link{list}} handed on to \code{\link[ggplot2]{labs}}
#' @param nrow              integer
#' @param mention_method    Whether to use \code{method} in title and/or facet labels
#' @param title             title
#' @param na.impute         TRUE or FALSE
#' @param legend.position   character
#' @examples
#' if (require(autonomics.data)){
#'    require(magrittr)
#'
#'    # GLUTAMINASE
#'       object <- autonomics.data::glutaminase
#'       object %>% plot_pca_samples()
#'       object %>% plot_lda_samples()
#'       object %>% plot_pls_samples()
#'       object %>% plot_projected_samples(
#'                     method = c('pca', 'lda', 'sma', 'pls'),
#'                     facet_var = c('PCA', 'LDA', 'SMA', 'PLS'),
#'                     nrow = 2)
#'
#'    # STEM CELL COMPARISON
#'       object <- autonomics.data::stemcomp.proteinratios
#'       object %>% plot_pca_samples()
#'       object %>% plot_pca_samples(dims = c(3,4))
#'
#'       object %>% plot_projected_samples(
#'                     method = c('pca', 'lda', 'sma', 'pls'),
#'                     facet_var = c('OK PCA', 'Not so nice LDA', 'Nice SMA', 'Very nice PLS'),
#'                     nrow = 2)
#'
#' }
#'
#' @author Aditya Bhagwat, Johannes Graumann
#' @importFrom data.table   data.table   :=
#' @importFrom magrittr %>% %<>%
#' @export
plot_projected_samples <- function(
   object,
   ...,
   listed_objects    = NULL,
   method            = c('pca', 'lda', 'sma', 'pls')[1],
   implementation    = c('mixOmics::plsda', 'mixOmics::splsda', 'ropls::opls')[1],
   dims              = 1:2,
   color_var         = 'subgroup', # default_color_var(object) gives 'block' when present!
   color_values      = default_color_values(object, color_var),
   shape_var         = default_shape_var(object),
   size_var          = NULL,
   txt_var           = default_txt_var(object),
   txt_outliers      = FALSE,
   group_var         = NULL,
   split_var         = NULL,
   facet_var         = NULL,
   scales            = ifelse(is.null(facet_var), 'fixed' ,'free'),
   labs              = list(color = stringi::stri_trans_totitle(color_var), shape = stringi::stri_trans_totitle(shape_var)),
   nrow              = NULL,
   mention_method    = ifelse(is.null(facet_var), TRUE, FALSE),
   title             = NULL,
   na.impute         = FALSE,
   legend.position   = 'right'
){

# Check prerequisites -----------------------------------------------------
   # Objects
   ## Combine all handed in objects
   if(!is.null(listed_objects)) listed_objects %>% assertive.types::assert_is_list()
   obj_list <- list(object) %>% c(list(...), listed_objects)
   obj_list %<>% magrittr::extract(!sapply(obj_list, is.null))

   ## Assert all
   obj_list %>%  lapply(autonomics.import::assert_is_valid_object) %>%
                 lapply(function(x) x %>% autonomics.import::exprs() %>% assertive.properties::assert_is_non_empty())
   obj_list %<>% lapply(function(x) x %>% validify_shape_values(shape_var))
   dims %>% assertive.types::assert_is_numeric()              %>%
            assertive.properties::assert_is_of_length(2)      %>%
            assertive.numbers::assert_all_are_whole_numbers() %>%
            assertive.numbers::assert_all_are_greater_than(0)

   for(var in c(color_var, group_var, shape_var, size_var, split_var, txt_var)){
      obj_list %>% sapply(function(x) x %>% autonomics.import::svars() %>% assertive.sets::assert_is_superset(var))
   }

   if(length(facet_var) == 1){
      obj_list %>% sapply(function(x) x %>% autonomics.import::svars() %>% assertive.sets::assert_is_superset(facet_var))

   } else if(length(facet_var) > 1){
      if(length(obj_list) == 1){ obj_list <- obj_list[[1]] %>% list() %>% rep(times = length(facet_var))
      } else {                   facet_var %>% assertive.properties::assert_are_same_length(obj_list)
      }
   }

   # Must be left AFTER 'facet_var' - as obj_list may change depending on that
   method %<>% match.arg(choices = c('pca', 'lda', 'sma', 'pls'), several.ok = TRUE)
   if(length(method) > 1) method %>% assertive.properties::assert_are_same_length(obj_list)

   # color_values # how to check this?

   scales %<>% match.arg(choices = c('fixed', 'free_x', 'free_y', 'free'), several.ok = FALSE)

   if(!is.null(labs)){
      labs %>% assertive.types::assert_is_list() %>%
               assertive.properties::assert_has_names() # Leaky: 'all_have_names' does not exist
   }

   if(!is.null(nrow)){
      nrow %>% assertive.types::assert_is_a_number() %>%
               assertive.numbers::assert_all_are_whole_numbers() %>%
               assertive.numbers::assert_all_are_greater_than(0)
   }

   mention_method %>% assertive.types::assert_is_a_bool()

   if(!is.null(title)) title %>% assertive.types::assert_is_a_string()

   na.impute %>% assertive.types::assert_is_a_bool()

   if(length(legend.position) == 1){
     legend.position %<>% match.arg(choices    = c("none", "left", "right", "bottom", "top"), several.ok = FALSE)
   } else if(length(legend.position) == 2){
      legend.position %>% assertive.types::assert_is_numeric() %>%
                          assertive.numbers::assert_all_are_in_closed_range(c(0,1))
   } else {
      stop('\'legend.position\' must be of length [1,2].')
   }

   # Transform
   project_list <- seq_along(obj_list) %>%
                   lapply(function(x){ obj_list[[x]] %>%
                                       project(method         = ifelse(length(method) == 1, method, method[x]),
                                               implementation = implementation,
                                               ndim           = max(dims),
                                               na.impute      = na.impute)})

# Processing --------------------------------------------------------------
   # Are we dealing with an on-the-fly generated facetting variable or do we use one from the data?
   ## Augment the facet vars
   if(length(facet_var) > 1){
      facet_var <- seq_along(facet_var) %>%
                   sapply(function(x){ make_sample_scores_title2(
                                          object                = obj_list[[x]],
                                          project_result        = project_list[[x]],
                                          dims                  = dims,
                                          method                = ifelse(length(method) == 1, method, method[x]),
                                          manual_label          = facet_var[x],
                                          mention_method        = mention_method,
                                          mention_feature_count = TRUE,
                                          mention_variance      = TRUE,
                                          separator             = '\n')})
   }

   # Generate data.tables/frames
   plotDF <- seq_along(obj_list) %>%
             lapply(function(x){
                      if(length(facet_var) <= 1){
                         make_projected_samples_df(object     = obj_list[[x]],
                                                   projection = project_list[[x]],
                                                   dims       = dims,
                                                   color_var  = color_var,
                                                   shape_var  = shape_var,
                                                   size_var   = size_var,
                                                   txt_var    = txt_var,
                                                   facet_var  = facet_var,
                                                   split_var  = split_var,
                                                   group_var  = group_var)
                      } else {
                         tmp_df <- make_projected_samples_df(object = obj_list[[x]],
                                                             projection = project_list[[x]],
                                                             dims       = dims,
                                                             color_var  = color_var,
                                                             shape_var  = shape_var,
                                                             size_var   = size_var,
                                                             txt_var    = txt_var,
                                                             facet_var  = NULL,
                                                             split_var  = split_var,
                                                             group_var  = group_var)
                         tmp_df[['facet']] <- facet_var[x]
                         tmp_df  %<>% dplyr::mutate_(facet = ~facet %>% factor(levels = facet_var))
                         tmp_df
                      }}) %>%
             data.table::rbindlist()

   # Initialize plot
   p <- ggplot2::ggplot(plotDF) + ggplot2::theme_bw()

   # Points
   p <- p + ggplot2::geom_point(
            ggplot2::aes_string(
               x     = 'x',
               y     = 'y',
               color = 'color',
               shape = 'shape',
               size  = 'size'))
   if (is.null(shape_var))  p <- p + ggplot2::scale_shape_manual(values = c(default = 15), guide = FALSE)
   if (is.null(size_var))   p <- p + ggplot2::scale_size_manual( values = c(default = 3),  guide = FALSE)
   plot_guide <- if (is.null(color_var)) FALSE  else TRUE
   p <- p + ggplot2::scale_color_manual(values = color_values, guide = plot_guide)
   if (!is.null(color_var)) p <- p + ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(shape = 15, size = 2)))

   # Line
   if (!is.null(group_var)){
      p <- p + ggplot2::geom_path(ggplot2::aes_string(x = 'x', 'y', color = 'color', group = 'group'))
   }

   # txt
   if (!is.null(txt_var)){
      p <- p + ggrepel::geom_text_repel(ggplot2::aes_string(x='x',y='y', label='label', color = 'color'), show.legend = FALSE)
   }
   if (txt_outliers){
      is_an_outlier <- autonomics.support::is_outlier(plotDF$x) | autonomics.support::is_outlier(plotDF$y)
      p <- p + ggrepel::geom_text_repel(data    = plotDF[is_an_outlier, ],
                                        mapping = ggplot2::aes_string(x='x', y='y', color = 'color', label = 'sample_id'))
   }

   # Facet
   if (!is.null(facet_var)){
      p <- p + ggplot2::facet_wrap(
               stats::as.formula('~ facet'),
               scales = scales,
               nrow   = nrow)
      p <- p + ggplot2::theme(strip.text = ggplot2::element_text(hjust = 0.05))
   }

   # Axes
   p <- p + ggplot2::geom_vline(xintercept = 0, linetype = 'dashed') +
            ggplot2::geom_hline(yintercept = 0, linetype = 'dashed')

   # Axis labels
   xlab <- paste0('X', dims[1])
   ylab <- paste0('X', dims[2])
   if(length(facet_var) <= 1){
      if (is.null(split_var)){ var <- project_list[[1]]$var
                               xlab %<>% paste0(' (', var[dims[1]], '%)')
                               ylab %<>% paste0(' (', var[dims[2]], '%)')}
   }
   p <- p + ggplot2::xlab(xlab) + ggplot2::ylab(ylab)

   # Title
   if(length(facet_var) <= 1){
      default_title <- make_sample_scores_title2(
                          object                = obj_list[[1]],
                          method                = method,
                          dims                  = dims,
                          project_result        = project_list[[1]],
                          mention_method        = mention_method,
                          mention_feature_count = TRUE,
                          mention_variance      = FALSE,
                          separator             = '; ')
      title <- if (is.null(title)) default_title   else   sprintf('%s: %s', title, default_title)
      p <- p + ggplot2::ggtitle(title)
   }
   p <- p + ggplot2::theme(legend.position = legend.position)

   # Legends/labs
   if(!is.null(labs)){
      p <- p + do.call(what = ggplot2::labs, args = labs)
   }

   # Left align
   if(length(facet_var) > 1){
      gp <- ggplot2::ggplotGrob(p)

      # Edit the grob
      # Is throwing a warning - fix later - outcomment now
      # gp <- grid::editGrob(grid::grid.force(gp),
      #                      grid::gPath("GRID.stripGrob", "GRID.text"),
      #                      grep = TRUE,
      #                      global = TRUE,
      #                      just = "left",
      #                      x = grid::unit(0, "npc"))
      grid::grid.draw(gp)
      invisible(gp)
   } else {
      p
   }

}


make_sample_scores_title2 <- function(
   object,
   project_result,
   method,
   dims,
   manual_label          = NULL,
   mention_method        = TRUE,
   mention_feature_count = TRUE,
   mention_variance      = FALSE,
   separator             = ' on '
){
   # Manual label
   string1 <- ifelse(is.null(manual_label), NA_character_, manual_label)
   # Method used
   string2 <- ifelse(mention_method, toupper(method), NA_character_)
   # Dimensions used
   if(mention_feature_count){
      selector <- if (method %in% c('pca', 'pls')){
         selector <- matrixStats::rowAnys(!is.na(autonomics.import::exprs(object))) &
            matrixStats::rowAnys(!is.infinite(autonomics.import::exprs(object)))
      } else {
         selector <- matrixStats::rowAlls(!is.na(autonomics.import::exprs(object))) &
            matrixStats::rowAlls(!is.infinite(autonomics.import::exprs(object)))
      }
      string3 <- paste0(sum(selector), ' / ', length(selector), ' features')
   } else {
      string3 <- NA_character_
   }
   # Variances explained
   if(mention_variance){
      var <- project_result[['var']][dims]
      string4 <- 'Variance expl. ' %>%
                 paste0( var     %>%
                         names() %>%
                         stringi::stri_replace_all_regex('.*(\\d+)$', 'X$1') %>%
                         paste0(': ', var, '%', collapse = ', '))
   } else {
      string4 <- NA_character_
   }
   # Assemble
   all_strings <- c(string1, string2, string3, string4) %>%
                  stats::na.omit()
   paste0(all_strings, collapse = separator)
}


#' @rdname plot_projected_samples
#' @export
plot_pca_samples <- function(
   object,
   ...,
   listed_objects    = NULL,
   dims              = 1:2,
   color_var         = 'subgroup', # default_color_var(object) gives 'block' when present!
   color_values      = default_color_values(object, color_var),
   shape_var         = default_shape_var(object),
   size_var          = NULL,
   txt_var           = default_txt_var(object),
   txt_outliers      = FALSE,
   group_var         = NULL,
   split_var         = NULL,
   facet_var         = NULL,
   scales            = ifelse(is.null(facet_var), 'fixed' ,'free'),
   labs              = list(color = stringi::stri_trans_totitle(color_var), shape = stringi::stri_trans_totitle(shape_var)),
   nrow              = NULL,
   mention_method    = ifelse(is.null(facet_var), TRUE, FALSE),
   title             = NULL,
   na.impute         = FALSE,
   legend.position   = 'right'
){
   plot_projected_samples( object,
                           ...,
                           listed_objects    = listed_objects,
                           method            = 'pca',
                           implementation    = character(0),
                           dims              = dims,
                           color_var         = color_var,
                           color_values      = color_values,
                           shape_var         = shape_var,
                           size_var          = size_var,
                           txt_var           = txt_var,
                           txt_outliers      = txt_outliers,
                           group_var         = group_var,
                           split_var         = split_var,
                           facet_var         = facet_var,
                           scales            = scales,
                           labs              = labs,
                           nrow              = nrow,
                           mention_method    = mention_method,
                           title             = title,
                           na.impute         = na.impute,
                           legend.position   = legend.position)
}


#' @rdname plot_projected_samples
#' @export
plot_sma_samples <- function(
   object,
   ...,
   listed_objects    = NULL,
   dims              = 1:2,
   color_var         = 'subgroup', # default_color_var(object) gives 'block' when present!
   color_values      = default_color_values(object, color_var),
   shape_var         = default_shape_var(object),
   size_var          = NULL,
   txt_var           = default_txt_var(object),
   txt_outliers      = FALSE,
   group_var         = NULL,
   split_var         = NULL,
   facet_var         = NULL,
   scales            = ifelse(is.null(facet_var), 'fixed' ,'free'),
   labs              = list(color = stringi::stri_trans_totitle(color_var), shape = stringi::stri_trans_totitle(shape_var)),
   nrow              = NULL,
   mention_method    = ifelse(is.null(facet_var), TRUE, FALSE),
   title             = NULL,
   na.impute         = FALSE,
   legend.position   = 'right'
){
   plot_projected_samples( object,
                           ...,
                           listed_objects    = listed_objects,
                           method            = 'sma',
                           implementation    = character(0),
                           dims              = dims,
                           color_var         = color_var,
                           color_values      = color_values,
                           shape_var         = shape_var,
                           size_var          = size_var,
                           txt_var           = txt_var,
                           txt_outliers      = txt_outliers,
                           group_var         = group_var,
                           split_var         = split_var,
                           facet_var         = facet_var,
                           scales            = scales,
                           labs              = labs,
                           nrow              = nrow,
                           mention_method    = mention_method,
                           title             = title,
                           na.impute         = na.impute,
                           legend.position   = legend.position)
}


#' @rdname plot_projected_samples
#' @export
plot_lda_samples <- function(
   object,
   ...,
   listed_objects    = NULL,
   dims              = 1:2,
   color_var         = 'subgroup', # default_color_var(object) gives 'block' when present!
   color_values      = default_color_values(object, color_var),
   shape_var         = default_shape_var(object),
   size_var          = NULL,
   txt_var           = default_txt_var(object),
   txt_outliers      = FALSE,
   group_var         = NULL,
   split_var         = NULL,
   facet_var         = NULL,
   scales            = ifelse(is.null(facet_var), 'fixed' ,'free'),
   labs              = list(color = stringi::stri_trans_totitle(color_var), shape = stringi::stri_trans_totitle(shape_var)),
   nrow              = NULL,
   mention_method    = ifelse(is.null(facet_var), TRUE, FALSE),
   title             = NULL,
   na.impute         = FALSE,
   legend.position   = 'right'
){
   plot_projected_samples( object,
                           ...,
                           listed_objects    = listed_objects,
                           method            = 'lda',
                           implementation    = character(0),
                           dims              = dims,
                           color_var         = color_var,
                           color_values      = color_values,
                           shape_var         = shape_var,
                           size_var          = size_var,
                           txt_var           = txt_var,
                           txt_outliers      = txt_outliers,
                           group_var         = group_var,
                           split_var         = split_var,
                           facet_var         = facet_var,
                           scales            = scales,
                           labs              = labs,
                           nrow              = nrow,
                           mention_method    = mention_method,
                           title             = title,
                           na.impute         = na.impute,
                           legend.position   = legend.position)
}


#' @rdname plot_projected_samples
#' @export
plot_pls_samples <- function(
   object,
   ...,
   listed_objects    = NULL,
   implementation    = c('mixOmics::plsda', 'mixOmics::splsda', 'ropls::opls')[1],
   dims              = 1:2,
   color_var         = 'subgroup', # default_color_var(object) gives 'block' when present!
   color_values      = default_color_values(object, color_var),
   shape_var         = default_shape_var(object),
   size_var          = NULL,
   txt_var           = default_txt_var(object),
   txt_outliers      = FALSE,
   group_var         = NULL,
   split_var         = NULL,
   facet_var         = NULL,
   scales            = ifelse(is.null(facet_var), 'fixed' ,'free'),
   labs              = list(color = stringi::stri_trans_totitle(color_var), shape = stringi::stri_trans_totitle(shape_var)),
   nrow              = NULL,
   mention_method    = ifelse(is.null(facet_var), TRUE, FALSE),
   title             = NULL,
   na.impute         = FALSE,
   legend.position   = 'right'
){
   plot_projected_samples( object,
                           ...,
                           listed_objects    = listed_objects,
                           method            = 'pls',
                           implementation    = implementation,
                           dims              = dims,
                           color_var         = color_var,
                           color_values      = color_values,
                           shape_var         = shape_var,
                           size_var          = size_var,
                           txt_var           = txt_var,
                           txt_outliers      = txt_outliers,
                           group_var         = group_var,
                           split_var         = split_var,
                           facet_var         = facet_var,
                           scales            = scales,
                           labs              = labs,
                           nrow              = nrow,
                           mention_method    = mention_method,
                           title             = title,
                           na.impute         = na.impute,
                           legend.position   = legend.position)
}




#==================================================
# PLOT PROJECTION FEATURES
#==================================================


get_pca_var <- function(method, dim){
   sprintf('%s%d', method, dim)
}

#' Order on feature loadings
#' @param object    SummarizedExperiment
#' @param method   'pca', 'sma', 'lda', 'pls'
#' @param dim       dimension
#' @param na.impute whether to impute (if applicable)
#' @return re-ordered sumexp
#' @importFrom magrittr %>%
#' @export
order_on_feature_loadings <- function(object, method, dim, na.impute){
   pca_var <- get_pca_var(method, dim)
   loadings <- autonomics.import::fdata(object) %>% magrittr::extract2(pca_var)
   idx <- order(loadings, na.last = NA)
   object[idx, ]
}

#' Extract top and bottom of sumexp
#' @param object  SummarizedExperiment
#' @param n       number of top features to be selected
#' @return sumex with top n/2 and bottom n/2 rows
#' @importFrom magrittr   %>%
#' @export
extract_top_and_bottom <- function(object, n = 16){

   nfeature <- nrow(object)
   n1 <- ceiling(n/2)
   n2 <- n - n1

   top <- seq(1, n1)
   bottom <- seq(nfeature-n2+1, nfeature) # opposite direction

   object %>% magrittr::extract(c(top, bottom), )
}

#' Plot PCA/PLS/LDA/SMA features
#'
#' Plots top PCA/PLS/LDA/SMA features.
#' Uses factor loadings in object when available.
#'
#' @param object           SummarizedExperiment
#' @param method           'pca', 'lda', 'pls', 'sma'
#' @param implementation   'character' or NULL
#' @param geom             value in \code{\link[autonomics.plot]{FEATURE_PLOTS}}
#' @param fvars            fvars used for plot annotation
#' @param dim              principal component dimension
#' @param n                number of top features to plot
#' @param na.impute        TRUE or FALSE
#' @param title            title
#' @param file             file
#' @param ...              passed to \code{\link[autonomics.plot]{plot_features}}
#' @examples
#' require(magrittr)
#'
#' # STEM CELL COMPARISON
#' if (require(autonomics.data)){
#'    object <- autonomics.data::stemcomp.proteinratios
#'    object %>% plot_pca_features(n=9)
#'    object %>% plot_pca_features(geom = 'bar')
#' }
#'
#' # GLUTAMINASE
#' if (require(autonomics.data)){
#'    object <- autonomics.data::glutaminase
#'    object %>% plot_pca_features()
#' }
#'
#' @importFrom magrittr  %>%
#' @export
plot_projection_features <- function(
   object,
   method,
   implementation  = NULL,
   geom            = default_feature_plots(object)[1],
   fvars           = default_fvars(object),
   dim             = 1,
   n               = 9,
   na.impute       = FALSE,
   title           = sprintf('X%d', dim),
   file            = NULL,
   ...
){
   # Check input args
   autonomics.import::assert_is_valid_object(object)
   if (ncol(object) < 3){
      autonomics.support::cmessage('\tExit %s: only %d samples', method, ncol(object))
      return(invisible(NULL))
   }
   assertive.sets::assert_is_subset(geom, FEATURE_PLOTS)

   # Add projection to eset if required
   projection_dim_in_eset <- paste0(method, dim) %in% autonomics.import::fvars(object)
   if (!projection_dim_in_eset){
      # Wipe earlier accounts
      idx <- which(autonomics.import::fvars(object) %>% stringi::stri_detect_regex(paste0('^', method)))
      if (length(idx)>0) autonomics.import::fdata(object) %<>% magrittr::extract(, -idx)
      # Add projection
      object %<>% add_projection(method, na.impute = na.impute, ndim = dim)
   }

   # Order on projection and write to file
   object %<>% order_on_feature_loadings(method = method, dim = dim, na.impute = na.impute) %>%
               extract_top_and_bottom(n=n)

   # plot
   object %>% plot_features(geom = geom, fvars = fvars, title = title, ...)

}

#' @rdname plot_projection_features
#' @export
plot_pca_features <- function(
   object,
   implementation  = NULL,
   geom            = default_feature_plots(object)[1],
   fvars           = default_fvars(object),
   dim             = 1,
   n               = 9,
   na.impute       = FALSE,
   title           = sprintf('X%d', dim),
   file            = NULL,
   ...
){
   plot_projection_features(object,
                            method         = 'pca',
                            implementation = implementation,
                            geom           = geom,
                            fvars          = fvars,
                            dim            = dim,
                            n              = n,
                            na.impute      = na.impute,
                            title          = title,
                            file           = file,
                            ...)
}

#' @rdname plot_projection_features
#' @export
plot_sma_features <- function(object,
   implementation  = NULL,
   geom            = default_feature_plots(object)[1],
   fvars           = default_fvars(object),
   dim             = 1,
   n               = 9,
   na.impute       = FALSE,
   title           = sprintf('X%d', dim),
   file            = NULL,
   ...
){
   plot_projection_features(object,
                            method         = 'sma',
                            implementation = implementation,
                            geom           = geom,
                            fvars          = fvars,
                            dim            = dim,
                            n              = n,
                            na.impute      = na.impute,
                            title          = title,
                            file           = file,
                            ...)
}

#' @rdname plot_projection_features
#' @export
plot_lda_features <- function(object,
   implementation  = NULL,
   geom            = default_feature_plots(object)[1],
   fvars           = default_fvars(object),
   dim             = 1,
   n               = 9,
   na.impute       = FALSE,
   title           = sprintf('X%d', dim),
   file            = NULL,
   ...
){
   plot_projection_features(object,
                            method         = 'lda',
                            implementation = implementation,
                            geom           = geom,
                            fvars          = fvars,
                            dim            = dim,
                            n              = n,
                            na.impute      = na.impute,
                            title          = title,
                            file           = file,
                            ...)
}

#' @rdname plot_projection_features
#' @export
plot_pls_features <- function(object,
   implementation  = NULL,
   geom            = default_feature_plots(object)[1],
   fvars           = default_fvars(object),
   dim             = 1,
   n               = 9,
   na.impute       = FALSE,
   title           = sprintf('X%d', dim),
   file            = NULL,
   ...
){
   plot_projection_features(object,
                            method         = 'pls',
                            implementation = implementation,
                            geom           = geom,
                            fvars          = fvars,
                            dim            = dim,
                            n              = n,
                            na.impute      = na.impute,
                            title          = title,
                            file           = file,
                            ...)

}



#========================
# ADD
#========================


#' Add PCA/LDA/PLS/SMA results to object
#' @param object          SummarizedExperiment
#' @param method          'pca' or 'lda'
#' @param implementation  string
#' @param na.impute       TRUE or FALSE
#' @param ndim            number of dimensions to include
#' @param ...             passed to add_projection
#' @examples
#' if (require(autonomics.data)){
#'
#'    # STEM CELL COMPARISON
#'    require(magrittr)
#'    object <- autonomics.data::stemcomp.proteinratios
#'    object %>%  add_pca(ndim = 3)
#'    object %<>% add_pca()
#'    object %>%  autonomics.import::svars()
#'    object %>%  autonomics.import::fvars()
#' }
#'
#' @importFrom magrittr  %>%
#' @export
add_projection <- function(
   object,
   method         = c('pca', 'lda', 'pls')[1],
   implementation = NULL,
   na.impute      = FALSE,
   ndim           = 2
){
   nsample <- ncol(object)
   npca    <- min(nsample-1, 2)
   if (npca == 0) return(object)

   # Wipe previous pca results, as they could be
   #   - on a different data subset
   #   - for different ndim
   idx <- autonomics.import::fvars(object) %>% stringi::stri_detect_regex(paste0('^', method)) %>% which()
   if (length(idx)>0){
      autonomics.import::fdata(object) %<>% magrittr::extract(, -idx)
      # autonomics.support::cmessage('\t\tAbort %s - object already contains projection scores', method)
      # return(object)
   }

   projection <- object %>% project(method         = method,
                                    implementation = implementation,
                                    na.impute      = na.impute,
                                    ndim           = ndim)
   pca_features <- data.frame(feature_id = rownames(projection$features), projection$features)
   autonomics.import::fdata(object) %<>% autonomics.support::left_join_keeping_rownames(pca_features, by = 'feature_id')
   object
}

#' @rdname add_projection
#' @export
add_pca <- function(object, ...){
   add_projection(object, method = 'pca', ...)
}

#' @rdname add_projection
#' @export
add_sma <- function(object, ...){
   add_projection(object, method = 'sma', ...)
}

#' @rdname add_projection
#' @export
add_lda <- function(object, ...){
   add_projection(object, method = 'lda', ...)
}

#' @rdname add_projection
#' @export
add_pls <- function(object, ...){
   add_projection(object, method = 'pls', ...)
}


