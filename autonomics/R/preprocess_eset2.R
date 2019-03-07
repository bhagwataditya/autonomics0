#' Preprocess
#' @param  object  eset
#' @param  normalize_on_svar                    svar on which to normalize (should have numeric values)
#' @param  invert_subgroups                     vector with subgroups that require inversion (relevant for ratio esets)
#' @param  channel_frac                         channel fractionator
#' @param  subgroup_frac                        subgroup fractionator
#' @param  mode_center                          whether to mode center each sample distribution (logical)
#' @param  displayed_features                   features to be displayed in the sample distributions (vector of numeric indexes or character \code{feature_id}s)
#' @param  feature_center                       features to center each sample distribution on (vector of numeric indexes or character \code{feature_id}s); median if \code{length(feature_center) > 1}
#' @param  normalize_within_samples             whether to normalize distribution within each sample (logical)
#' @param  normalize_between_subgroup_samples   whether to normalize distributions between samples of same subgroup (logical)
#' @param  normalize_between_all_samples        whether to normalize distributions of all samples (logical)
#' @param  filter_features_nonzero_for_two_samples_in_some_subgroup logical
#' @param  plot                                 whether to plot original and normalized sample distributions (logical)
#' @param  retain_plot_objects                  whether to return plot objects (grobs) in a \code{tmpplot} slot of the output object
#' @param  color_var                            svar mapped to color
#' @param  shape_var                            svar mapped to shape
#' @param  txt_var                              svar mapped to txt
#' @param  result_dir                           directory to which sample distributions (both before and after normalization) are printed
#' @return eset with updated exprs
#' @examples
#' library(magrittr)
#' result_dir <- tempdir() %T>% message()
#' if (require(autonomics.data)){
#'    object <- autonomics.data::billing2016
#'    object %>% autonomics::preprocess_eset(
#'                    normalize_between_subgroup_samples = TRUE, result_dir = result_dir)
#'    object %>% autonomics::preprocess_eset(
#'                    invert_subgroups                   = c('BM_E', 'BM_EM', 'EM_E'),
#'                    mode_center                        = TRUE,
#'                    normalize_within_samples           = TRUE,
#'                    normalize_between_subgroup_samples = TRUE,
#'                    normalize_between_all_samples      = TRUE,
#'                    result_dir               = result_dir)
#'    # Retain plotted object
#'    tmp_object <- object %>% autonomics::preprocess_eset(
#'                    normalize_between_subgroup_samples = TRUE, result_dir = result_dir,
#'                    retain_plot_objects = TRUE)
#'    grid::grid.draw(autonomics.import::tmpplot(tmp_object)[[1]])
#'    grid::grid.draw(autonomics.import::tmpplot(tmp_object)[[2]])
#' }
#' if (require(billing.differentiation.data)){
#'    object <- billing.differentiation.data::protein.ratios
#'    object %>% preprocess_eset(normalize_within_samples = TRUE, result_dir = result_dir)
#'    object %>% preprocess_eset(
#'                    mode_center                        = TRUE,
#'                    normalize_within_samples           = TRUE,
#'                    normalize_between_subgroup_samples = TRUE,
#'                    normalize_between_all_samples      = TRUE,
#'                    result_dir                         = result_dir)
#' }
#' if (require(atkin.2014)){
#'    object <- atkin.2014::soma
#'    object %>% preprocess_eset(normalize_within_samples           = TRUE, result_dir = result_dir)
#'    object %>% preprocess_eset(normalize_between_subgroup_samples = TRUE, result_dir = result_dir)
#'    object %>% preprocess_eset(normalize_between_all_samples      = TRUE, result_dir = result_dir)
#' }
#' if (require(halama.2016)){
#'    object <- halama.2016::cell.metabolites
#'    object %>% preprocess_eset(normalize_on_svar = 'ProteinContent')
#' }
#' if (require(alnoubi.2017)){
#'    object <- alnoubi.2017::proteins.replicates
#'    object %>% preprocess_eset(displayed_features = 'PG1396', feature_center = 'PG1396')
#' }
#' @importFrom magrittr   %>%   %<>%
#' @export
preprocess_eset <- function(
   object,
   normalize_on_svar        = character(0),
   invert_subgroups         = NULL,
   channel_frac             = '/',
   subgroup_frac            = '_',
   mode_center              = FALSE,
   displayed_features       = NULL,
   feature_center           = displayed_features,
   normalize_within_samples           = FALSE,
   normalize_between_subgroup_samples = FALSE,
   normalize_between_all_samples      = FALSE,
   filter_features_nonzero_for_two_samples_in_some_subgroup = FALSE,
   plot                     = TRUE,
   retain_plot_objects      = FALSE,
   color_var                = autonomics.plot::default_color_var(object),
   shape_var                = autonomics.plot::default_shape_var(object),
   txt_var                  = autonomics.plot::default_txt_var(object),
   result_dir               = NULL
){
   commonargs <- list(color_var = color_var, shape_var = shape_var, txt_var = txt_var)

   # Display individual features
   if(!is.null(displayed_features)){
      ## Check input validity
      displayed_features %<>% autonomics.import::assert_all_are_valid_features(object)
      ## Amend commonargs
      commonargs[['displayed_features']] <- displayed_features
   }

   # Original
   if (plot || retain_plot_objects){
      plotargs <- commonargs %>% c(list(object = object, descr = 'original'))
      plotout  <- do.call(autonomics.plot::plot_sample_distrs_n_pca, plotargs)
      distrs   <- plotout %>% magrittr::extract(1)
      pcas     <- plotout %>% magrittr::extract(2)
   }

   # Normalize on svar
   if (length(normalize_on_svar)>0){
      object %<>% autonomics.preprocess::normalize_samples_on_svar(normalize_on_svar)
      if (plot || retain_plot_objects){
         plotargs <- commonargs %>% c(list(object = object, descr = paste0('normalize.', normalize_on_svar)))
         plotout <- do.call(autonomics.plot::plot_sample_distrs_n_pca, plotargs)
         distrs %<>% c(plotout[1])
         pcas   %<>% c(plotout[2])
      }
   }

   # Invert subgroups
   if (!is.null(invert_subgroups)){
      object %<>% autonomics.import::invert(
                       invert_subgroups, channel_frac = channel_frac, subgroup_frac = subgroup_frac)
      if (plot || retain_plot_objects){
         plotargs <- commonargs %>% c(list(object = object, descr = 'invert.subgroups'))
         plotout <- do.call(autonomics.plot::plot_sample_distrs_n_pca, plotargs)
         distrs %<>% c(plotout[1])
         pcas   %<>% c(plotout[2])
      }
   }

   # Mode center
   if (mode_center){
      autonomics.support::cmessage('\t\tMode center')
      object %<>% autonomics.preprocess::mode_center()
      if (plot || retain_plot_objects){
         plotargs <- commonargs %>% c(list(object = object, descr = 'mode.center'))
         plotout <- do.call(autonomics.plot::plot_sample_distrs_n_pca, plotargs)
         distrs %<>% c(plotout[1])
         pcas   %<>% c(plotout[2])
      }
   }

   # Feature(s) center
   if(!is.null(feature_center)){
      ## Check input validity
      feature_center %<>% autonomics.import::assert_all_are_valid_features(object)
      ## Center on the feature(s)
      autonomics.support::cmessage('\t\tFeature center')
      object %<>% autonomics.preprocess::feature_center(feature_center)
      if (plot || retain_plot_objects){
         plotargs <- commonargs %>% c(list(object = object, descr = 'feature.center'))
         plotout <- do.call(autonomics.plot::plot_sample_distrs_n_pca, plotargs)
         distrs %<>% c(plotout[1])
         pcas   %<>% c(plotout[2])
      }
   }

   # Normalize within samples
   if (normalize_within_samples){
      autonomics.support::cmessage('\t\tNormalize within samples')
      object %<>% autonomics.preprocess::invnorm()
      if (plot || retain_plot_objects){
         plotargs <- commonargs %>% c(list(object = object, descr = 'normalize.within.each.sample'))
         plotout <- do.call(autonomics.plot::plot_sample_distrs_n_pca, plotargs)
         distrs %<>% c(plotout[1])
         pcas   %<>% c(plotout[2])
      }
   }

   # Normalize between subgroup samples
   if (normalize_between_subgroup_samples){
      autonomics.support::cmessage('\t\tNormalize between subgroup samples')
      object %<>% autonomics.preprocess::quantnorm_within_subgroups()
      if (plot || retain_plot_objects){
         plotargs <- commonargs %>% c(list(object = object, descr = 'normalize.between.subgroup.samples'))
         plotout <- do.call(autonomics.plot::plot_sample_distrs_n_pca, plotargs)
         distrs %<>% c(plotout[1])
         pcas   %<>% c(plotout[2])
      }
   }

   # Normalize between all samples
   if (normalize_between_all_samples){
      autonomics.support::cmessage('\t\tNormalize between all samples')
      object %<>% autonomics.preprocess::quantnorm()
      if (plot || retain_plot_objects){
         plotargs <- commonargs %>% c(list(object = object, descr = 'normalize.between.all.samples'))
         plotout <- do.call(autonomics.plot::plot_sample_distrs_n_pca, plotargs)
         distrs %<>% c(plotout[1])
         pcas   %<>% c(plotout[2])
      }
   }

   # Filter
   if (filter_features_nonzero_for_two_samples_in_some_subgroup){
      object %<>% autonomics.preprocess::filter_features_nonzero_for_two_samples_in_some_subgroup()
   }

   # Plot distrs
   if (plot){
      if (is.null(result_dir)){
         plotout <- c(distrs, pcas)
         autonomics.plot::multiplot(plotlist = plotout, layout = matrix(1:length(plotout), ncol = 2, byrow = FALSE))
      } else {
         grDevices::pdf(sprintf('%s/normalization_distrs.pdf', result_dir), width = 12, height = 9)
         length(distrs)
         ncols <- 2
         nrows <- ceiling(length(distrs) / ncols)
         autonomics.plot::multiplot(plotlist = distrs, layout = matrix(1:(ncols*nrows), ncol = ncols, byrow = FALSE))
         grDevices::dev.off()
         grDevices::pdf(sprintf('%s/normalization_pcas.pdf', result_dir), width = 12, height = 9)
         autonomics.plot::multiplot(plotlist = pcas, layout = matrix(1:(ncols*nrows), ncol = ncols, byrow = FALSE))
         grDevices::dev.off()
      }
   }
   if (retain_plot_objects){
      ncols <- 2
      nrows <- ceiling(length(distrs) / ncols)
      autonomics.import::tmpplot(object) <- list(
         grid::grid.grabExpr(
            autonomics.plot::multiplot(
               plotlist = distrs,
               layout   = matrix(
                  1:(ncols*nrows),
                  ncol  = ncols,
                  byrow = FALSE)),
            wrap = TRUE,
            warn = 1),
         grid::grid.grabExpr(
            autonomics.plot::multiplot(
               plotlist = pcas, 
               layout   = matrix(
                  1:(ncols*nrows),
                  ncol  = ncols,
                  byrow = FALSE)),
            wrap = TRUE,
            warn = 1))
   }

   # Return
   object
}

#' Preprocess max quant ratios
#'
#' Normalize max quant ratios and make them ready for statistical analysis
#' @param object   eset
#' @param plot       whether to plot effects on distribution and PCA (logical)
#' @param result_dir result dir
#' @param color_var  svar mapped to color
#' @param shape_var  svar mapped to shape
#' @param txt_var    svar mapped to txt
#' @return normalized eset
#' @export
preprocess_maxquant_ratios <- function(
   object,
   plot,
   result_dir = NULL,
   color_var = autonomics.plot::default_color_var(object),
   shape_var = autonomics.plot::default_shape_var(object),
   txt_var   = autonomics.plot::default_txt_var(object)
){

   # No normalization if # features too small
   # This happens e.g. for depleted samples
   within_samples <- if (nrow(object) > 500) TRUE else FALSE

   preprocess_eset(
      object                             = object,
      normalize_within_samples             = within_samples,
      normalize_between_subgroup_samples   = FALSE,
      normalize_between_all_samples        = FALSE,
      plot                                 = plot,
      result_dir                           = result_dir,
      color_var                            = color_var,
      shape_var                            = shape_var,
      txt_var                              = txt_var
   )
}

#' Normalize somascan data
#' @param object     eset
#' @param plot         whether to plot sample distributions and PCA (logical)
#' @param result_dir   result dir
#' @param color_var    svar mapped to color
#' @param shape_var    svar mapped to shape
#' @param txt_var      svar mapped to txt
#' @examples
#' if (require(atkin.2014)){
#'    require(magrittr)
#'    atkin.2014::soma %>% preprocess_somascan(plot = TRUE, color_var = 'time')
#' }
#' @return normalized eset
#' @export
preprocess_somascan <- function(
   object,
   plot,
   result_dir = NULL,
   color_var = autonomics.plot::default_color_var(object),
   shape_var = autonomics.plot::default_shape_var(object),
   txt_var   = autonomics.plot::default_txt_var(object)
){

   preprocess_eset(
      object,
      normalize_within_samples             = FALSE,
      normalize_between_subgroup_samples   = FALSE,
      normalize_between_all_samples        = TRUE,
      plot                                 = plot,
      result_dir                           = result_dir,
      color_var                            = color_var,
      shape_var                            = shape_var,
      txt_var                              = txt_var
   )
}

#' Normalize metabolon data
#' @param object   eset
#' @param plot       whether to plot sample distributions and PCA (logical)
#' @param result_dir result dir
#' @param color_var  svar mapped to color
#' @param shape_var  svar mapped to shape
#' @param txt_var    svar mapped to txt
#' @return normalized eset
#' @examples
#' if (require(subramanian.2016)){
#'    require(magrittr)
#'    object <- subramanian.2016::metabolon
#'    object %>% preprocess_metabolon(plot = TRUE, color_var = 'condition')
#' }
#' @export
preprocess_metabolon <- function(
   object,
   plot,
   result_dir = NULL,
   color_var  = autonomics.plot::default_color_var(object),
   shape_var  = autonomics.plot::default_shape_var(object),
   txt_var    = autonomics.plot::default_txt_var(object)
){

   preprocess_eset(
      object,
      normalize_within_samples           = FALSE,
      normalize_between_subgroup_samples = FALSE,
      normalize_between_all_samples      = FALSE,
      plot                               = plot,
      result_dir                         = result_dir,
      color_var                          = color_var,
      shape_var                          = shape_var,
      txt_var                            = txt_var
   )
}

#' Sample normalize object automatically
#' @param object          eset
#' @param plot              whether to plot sample distributions for different prepro steps
#' @param result_dir        if non-NULL, directory to which plots will be printed
#' @param color_var         svar mapped to color
#' @param shape_var         svar mapped to shape
#' @param txt_var           svar mapped to txt
#' @return object after sample normalization
#' @examples
#' # max quant data
#' # --------------
#' library(magrittr)
#' if (require(autonomics.data)){
#'    autonomics.data::billing2016 %>% autonomics::preprocess_smart()
#'    autonomics.data::ALL %>% autonomics::preprocess_smart()
#' }
#' if (require(billing.differentiation.data)){
#'    billing.differentiation.data::protein.ratios %>%
#'       autonomics::preprocess_smart()
#'    billing.differentiation.data::phospho.occupancies %>%
#'       autonomics::preprocess_smart()
#' }
#' if (require(atkin.2014)){
#'    atkin.2014::soma  %>% 
#'       autonomics::preprocess_smart()
#' }
#' @importFrom magrittr   %>%
#' @export
preprocess_smart <- function(
   object,
   plot       = TRUE,
   result_dir = NULL,
   color_var  = autonomics.plot::default_color_var(object),
   shape_var  = autonomics.plot::default_shape_var(object),
   txt_var    = autonomics.plot::default_txt_var(object)
){
   # maxquant eset
   if (autonomics.import::is_maxquant_eset(object) & autonomics.import::contains_ratios(object)){
      object %>% preprocess_maxquant_ratios(plot, result_dir, color_var = color_var, shape_var = shape_var, txt_var = txt_var)

   # metabolon eset
   } else if (autonomics.import::is_metabolon_eset(object)){
      object %>% preprocess_metabolon(      plot, result_dir, color_var = color_var, shape_var = shape_var, txt_var = txt_var)

   # soma eset
   } else if (autonomics.import::is_soma_eset(object)) {
      object %>% preprocess_somascan(       plot, result_dir, color_var = color_var, shape_var = shape_var, txt_var = txt_var)

   } else {
      object
   }
}
