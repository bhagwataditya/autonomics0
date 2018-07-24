utils::globalVariables('.')

#' Write pca features
#' @param object SummarizedExperiment, eSet, or EList
#' @param result_dir result dir
#' @param method which character
#' @examples
#' require(magrittr)
#' if (require(subramanian.2016)){
#'    object <- subramanian.2016::metabolon
#'    result_dir <- tempdir() %T>% message()
#'    object %>% autonomics.explore::add_pca_to_eset() %>% 
#'               autonomics.explore::write_projected_features(result_dir = result_dir)
#' }
#' @importFrom magrittr   %>%
#' @export
write_projected_features <- function(object, result_dir, method = 'pca'){
   pattern <- sprintf('^(%s)([0-9]+)', method)
   projection_columns <- autonomics.import::fvars(object) %>% 
                         magrittr::extract(stringi::stri_detect_regex(., sprintf(pattern, method)))
   for (pca_col in projection_columns){
      which_features <- pca_col %>% stringi::stri_replace_first_regex(pattern, '$1')
      dim            <- pca_col %>% stringi::stri_replace_first_regex(pattern, '$2') %>% as.numeric()
      file_name <- sprintf('%s/%s/%s%s_features.txt', result_dir, method, method, dim, which_features)
      idx <- order(autonomics.import::fdata(object)[[pca_col]], na.last = NA)
      object %>% magrittr::extract(idx, ) %>% autonomics.import::write_fdata_to_file(file_name)
   }
}

#' Run pca and print results
#' @param object         eset
#' @param method        'pca', 'pls', or 'lda'
#' @param implementation character
#' @param result_dir     result dir
#' @param feature_plots  subset of \code{\link[autonomics.plot]{FEATURE_PLOTS}}
#' @param x              svar mapped to x in feature plots
#' @param color_var      svar mapped to color in feature plots
#' @param color_values        color vector (names = color_var levels, values = colors)
#' @param shape_var      svar mapped to shape in feature plots
#' @param group_var      svar mapped to group in feature plots
#' @param txt_var        svar mapped to txt   in feature plots
#' @param line           whether to connect points in feature plots with a line (logical)
#' @param na.impute      whether to impute missing values
#' @param ...            for backward compatibility
#' @return eset with pca results
#' @examples
#' require(magrittr)
#' result_dir <- tempdir() %T>% message()
#' 
#' # STEM CELL COMPARISON
#' if (require(autonomics.data)){
#'    object <- autonomics.data::stemcomp.proteinratios
#'    object %>% autonomics.explore::add_and_write_pca(result_dir = result_dir)
#' }
#' 
#' # STEM CELL DIFFERENTIATION
#' if (require(billing.differentiation.data)){
#'    require(magrittr)
#'    result_dir <- tempdir() %T>% message()
#'    billing.differentiation.data::rna.voomcounts %>% 
#'       autonomics.explore::add_and_write_pca(result_dir = result_dir)
#' }
#' 
#' # GLUTAMINASE
#' if (require(autonomics.data)){
#'    result_dir <- tempdir() %T>% message()
#'    autonomics.data::glutaminase %>% autonomics.explore::add_and_write_pca(
#'       result_dir   = result_dir,
#'       color_var    = 'GROUP_DESCRIPTION', 
#'       color_values = c(Control = 'orange', Vehicle = 'red', `Concentration 1` = 'green', 
#'                       `Concentration 2` = 'blue'))
#' }
#'                       
#' if (require(subramanian.2016)){
#'    object <- subramanian.2016::metabolon
#'    object %>% add_and_write_pca(result_dir = result_dir, color_var = 'condition')
#' }
#' @importFrom magrittr   %>%  %<>%
#' @export
add_and_write_projection <- function(
   object,
   method        = 'pca',
   implementation = NULL,
   result_dir,
   feature_plots = autonomics.plot::default_feature_plots(object),
   x             = autonomics.plot::default_x(object, feature_plots),
   color_var     = autonomics.plot::default_color_var(object),
   color_values  = autonomics.plot::default_color_values(object, color_var),
   shape_var     = autonomics.plot::default_shape_var(object),
   group_var     = autonomics.plot::default_group_var(object),
   txt_var       = autonomics.plot::default_txt_var(object),
   line          = autonomics.plot::default_line(object),
   na.impute     = FALSE
){
   # Check input args
   autonomics.import::assert_is_valid_eset(object)
   if (ncol(object) < 3){
      autonomics.support::cmessage('\tNo PCA (only %d samples)', ncol(object))
      return(object)
   }
   assertive.sets::assert_is_subset(feature_plots, autonomics.plot::FEATURE_PLOTS)
   
   # Run pca, add to fdata, write to file
   object %<>% autonomics.explore::add_projection_to_eset(
                  method         = method, 
                  implementation = implementation, 
                  na.impute      = na.impute)
   object %>% autonomics.explore::write_projected_features(result_dir, method = method)
   
   # Create sample + feature plots
   plot_args <- list(object     = object,
                     method       = method, 
                     implementation = implementation,
                     result_dir   = result_dir,
                     color_var    = color_var,
                     color_values = color_values,
                     shape_var    = shape_var,
                     group_var    = group_var,
                     txt_var      = txt_var,
                     line         = line, 
                     na.impute    = na.impute)
   for (iplot in seq_along(feature_plots)){
      cur_plot_args <- plot_args %>% c(list(feature_plot = feature_plots[iplot], x = x[iplot]))
      autonomics.explore::plot_projected_samples_and_features %>% do.call(cur_plot_args)
   }

   # Return
   object
}

#' @rdname add_and_write_projection
#' @export
add_and_write_pca <- function(object, ...){
   add_and_write_projection(object, method = 'pca', ...)
}

#' @rdname add_and_write_projection
#' @export
add_and_write_sma <- function(object, ...){
   add_and_write_projection(object, method = 'sma', ...)
}

#' @rdname add_and_write_projection
#' @export
add_and_write_pls <- function(object, ...){
   add_and_write_projection(object, method = 'pls', ...)
}

#' @rdname add_and_write_projection
#' @export
add_and_write_lda <- function(object, ...){
   add_and_write_projection(object, method = 'lda', ...)
}

#' @rdname add_and_write_projection
#' @export
add_and_print_pca <- function(...){
   .Deprecated('add_and_write_pca')
   add_and_write_pca(...)
}

