
#' Transform vector to normal distribution
#' @param x numeric vector
#' @param mean  mean
#' @param sd    standard deviation
#' @return transformed vector
transform_to_normal <- function(x, mean, sd){
   selector <- !is.na(x) & !is.nan(x) & !is.infinite(x)
   pvals <- rank(x[selector]) / (length(x[selector]) + 1)
   y <- x
   y[selector] <- stats::qnorm(pvals, mean = mean, sd = sd)
   y
}

#' Transform vector to standard normal distribution
#' @param x numeric vector
#' @return transformed vector
#' @importFrom magrittr %>%
transform_to_standard_normal <- function(x){
   x %>% transform_to_normal(mean = 0, sd = 1)
}

#' Transform vector to fitting normal distribution
#' @param x numeric vector
#' @return transformed vector
#' @importFrom magrittr  %>%
transform_to_fitting_normal <- function(x){
   . <- NULL
   pars <- x %>% (function(x){
                     x %>%
                     magrittr::extract(!is.na(.) & !is.infinite(.)) %>%
                     MASS::fitdistr('normal') %>% magrittr::extract2('estimate')
                 })
   x %>% transform_to_normal(mean = pars[['mean']], sd = pars[['sd']])
}

#' Inverse normal transform samples
#' @param object eset
#' @return normalized eset
#' @importFrom magrittr  %<>%
#' @export
invnorm <- function(object){
   autonomics.import::exprs(object) %<>% apply(2, transform_to_fitting_normal)
   object
}

#' Quantile normalize samples
#' @param object eset
#' @return normalized eset
#' @importFrom  magrittr %<>%
#' @export
quantnorm <- function(object){
   autonomics.import::exprs(object) %<>% limma::normalizeBetweenArrays()
   object
}

#' Quantile normalize within subgroups
#' @param object eset
#' @return normalized eset
#' @importFrom magrittr   %>%
#' @export
quantnorm_within_subgroups <- function(object){
   for (cur_subgroup in object$subgroup){
      idx <- object$subgroup == cur_subgroup
      normalized_exprs <- autonomics.import::exprs(object)[, idx, drop = FALSE] %>%
                          limma::normalizeBetweenArrays()
      autonomics.import::exprs(object)[, idx] <- normalized_exprs
   }
   object
}

#' Normalize samples to standard normal
#' @param object eset
#' @importFrom  magrittr  %<>%
#' @export
normalize_samples_to_standard_normal <- function(object){
   autonomics.import::exprs(object) %<>% apply(2, transform_to_standard_normal)
   object
}

#' z score samples
#' @param object eset
#' @return z-scored eset
#' @importFrom magrittr %<>%
#' @export
z_score_samples <- function(object){
   autonomics.import::exprs(object) %<>% scale()
   object
}

#' Normalize samples on svar
#' @param object      SummarizedExperiment
#' @param normvar     svar on which to normalize
#' @param result_dir  directory where normalization plots should be printed
#' @return normalized SummarizedExperiment
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    object <- autonomics.data::glutaminase
#'    object %>% autonomics.preprocess::normalize_samples_on_svar('PROTEIN_CONTENT')
#' }
#' @importFrom magrittr  %<>%   %>%
#' @export
normalize_samples_on_svar <- function(object, normvar, result_dir = NULL){
   . <- NULL

   object[[normvar]] %<>% (function(x) if (is.factor(x))    as.character(x) else x)
   object[[normvar]] %<>% (function(x) if (is.character(x)) as.numeric(x)   else x)
   object[[normvar]] %<>% (function(x){  x[is.na(x)] <- 0; x})

   #   plotlist <- vector(mode = 'list', length = 4)
   #   plotlist[[1]] <- object %>% autonomics.plot::plot_sample_distributions(title = 'Original')
   #   plotlist[[2]] <- object %>% autonomics.plot::plot_pca_samples(title = 'Original')

   autonomics.import::exprs(object) %<>%  magrittr::subtract(object[[normvar]] %>%
                                                             replicate(nrow(object), .) %>%
                                                             t())
   #    title <-  paste0('Normalized on ', normvar)
   #    plotlist[[3]] <- object %>% autonomics.plot::plot_sample_distributions(title = title)
   #    plotlist[[4]] <- object %>% autonomics.plot::plot_pca_samples(title = title)
   #    if (is.null(result_dir)){
   #       autonomics.plot::multiplot(plotlist = plotlist, cols = 2)
   #    } else {
   #       file_name <- paste0(result_dir, '/sample_normalization.pdf')
   #       grDevices::pdf(file_name)
   #       autonomics.plot::multiplot(plotlist = plotlist, cols = 2, width = 20, height = 18)
   #       grDevices::dev.off()
   #    }

   object
}


# normalize_samples <- function(...){
#    .Deprecated('preprocess_eset')
#    preprocess_eset(...)
# }

