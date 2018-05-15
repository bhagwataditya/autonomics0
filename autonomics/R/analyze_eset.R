####
#' @title Analyze eset with omics data
#' @description
#' \enumerate{
#'   \item Sample normalization
#'   \item Principal component analysis
#'   \item General linear model based contrast analysis (see details)
#'   \item Gene Set Over representation analysis
#' }
#' \cr
#'
#' @details
#' The data is analyzed with the model "~ 0 + subgroup", which fits a separate coefficient for each subgroup level.
#' Each element in the contrast vector can contain any linear combination of the model coefficients.
#' \itemize{
#'   \item the coefficients themselves:      contrasts = c('A', 'B', 'C')
#'   \item differences between coefficients: contrasts = c('A-B', 'A-C', 'B-C')
#'   \item the mean of all coefficients:     contrasts = c('(A + B + C)/3')
#'   \item or a combination of the previous types: contrasts = c('A', 'B', 'C', 'A-B', 'B-C', 'A-C', '(A+B+C)/3')
#' }
#' See the examples for more details
#'
#' @param object          eset
#' @param design            design matrix
#' @param contrasts         contrast definitions
#' @param top_definition    Definition of 'top features'.
#' @param result_dir        directory to which to write results
#' @param universe          'detectome', 'genome', or NULL (no ora)
#' @param feature_plots     subset of c('bars', 'profiles', 'distributions')
#' @param x                 svar mapped to x in feature plots
#' @param color_var         svar mapped to color in feature plots
#' @param color_values      color values vector
#' @param shape_var         svar mapped to shape in feature plots
#' @param group_var         svar mapped to group in feature plots
#' @param txt_var           svar mapped to txt in PCA plots (but not in feature plots!)
#' @param facet_var         svar used to facet feature plots
#' @param line              whether to join points in feature plot with line (logical)
#' @param fvars             fvars to use in plot annotation
#' @param do_pca            whether to perform PCA analysis (logical)
#' @param cluster_features  whether to cluster features (logical)
#' @param nplot             no of top features to be plotted
#' @param feature_plot_width in inches
#' @param feature_plot_height in inches
#' @param min_set_size      min set size (for ora)
#' @param max_set_size      max set size (for ora)
#' @param zip_results       zip analysis results (logical)?
#' @return eset with analysis results in fdata
#' @examples
#'
#' require(magrittr)
#' # Compare each subgroup level to zero
#' #------------------------------------
#' if (require(autonomics.data)){
#'    result_dir <- tempdir() %>% paste0('/analysis.results') %T>% message()
#'    dir.create(result_dir, showWarnings = FALSE)
#'    object <- autonomics.data::billing2016
#'    unique(object$subgroup)
#'    object %>% autonomics::analyze_eset(result_dir, universe = NULL)
#' }
#'
#' # Differences between subgroup levels
#' #------------------------------------
#' if (require(billing.differentiation.data)){
#'    result_dir <- tempdir() %>% paste0('/analysis.results') %T>% message()
#'    dir.create(result_dir, showWarnings = FALSE)
#'    autonomics::analyze_eset(
#'       object   = billing.differentiation.data::protein.ratios,
#'       result_dir = result_dir,
#'       # contrasts  = billing.differentiation.data::contrasts[10],
#'       universe   = NULL)
#'
#'    autonomics::analyze_contrasts(
#'       object   = billing.differentiation.data::rna.voomcounts,
#'       result_dir = result_dir,
#'       contrasts  = billing.differentiation.data::contrasts[1:2],
#'       universe   = NULL)
#' }
#'
#' # Commonalities across subgroup levels
#' #-------------------------------------
#' if (require(atkin.2014)){
#'    require(magrittr)
#'    result_dir <- tempdir() %>% paste0('/analysis.results') %T>% message()
#'    dir.create(result_dir, showWarnings = FALSE)
#'    object <- atkin.2014::soma
#'    contrasts <- atkin.2014::contrasts[1]
#'    object$subgroup %>% unique()
#'    object %>% autonomics::analyze_eset(
#'                    contrasts = contrasts, result_dir = result_dir, universe = NULL)
#' }
#'
#' # Metabolon data
#' #---------------
#' if (require(subramanian.2016)){
#'    require(magrittr)
#'    result_dir <- tempdir() %>% paste0('/analysis.results') %T>% message()
#'    dir.create(result_dir, showWarnings = FALSE)
#'    object <- subramanian.2016::metabolon
#'    contrasts <- subramanian.2016::contrasts.metabolon[1]
#'    object %>% autonomics::analyze_eset(
#'                    contrasts = contrasts, result_dir = result_dir,
#'                    color_var = 'condition', group_var = 'condition', line = TRUE)
#' }
#' @importFrom   magrittr   %>%   %<>%
#' @export
analyze_eset <- function(
  object,
  result_dir       = default_result_dir(object),
  design           = autonomics.find::create_design_matrix(object),
  contrasts        = autonomics.find::default_contrasts(object),
  top_definition   = autonomics.find::default_top_definition(object),
  universe         = autonomics.ora::default_universe(object),
  feature_plots    = autonomics.plot::default_feature_plots(object),
  x                = autonomics.plot::default_x(object, feature_plots),
  color_var        = autonomics.plot::default_color_var(object),
  color_values     = autonomics.plot::default_color_values(object),
  shape_var        = autonomics.plot::default_shape_var(object),
  group_var        = autonomics.plot::default_group_var(object),
  txt_var          = autonomics.plot::default_txt_var(object),
  facet_var        = autonomics.plot::default_facet_var(),
  line             = autonomics.plot::default_line(object),
  fvars            = autonomics.plot::default_fvars(object),
  do_pca           = TRUE,
  cluster_features = default_cluster_features(),
  nplot            = autonomics.find::default_nplot(object),
  feature_plot_width  = NULL, 
  feature_plot_height = NULL,
  min_set_size     = 5,
  max_set_size     = Inf,
  zip_results      = FALSE
){
   # Assert valid args
   autonomics.import::assert_is_valid_eset(object)
   assertive.types::assert_is_character(result_dir)
   assertive.types::assert_is_logical(zip_results)

   # Exit if no exprs
   if(assertive.properties::is_empty(autonomics.import::exprs(object))){
      message("\tNo expression data; skipping PCA and sample distribution plots, and differential expression.")
      return(object)
   }

   # Initialize
   ggplot2::theme_set(ggplot2::theme_bw())
   if(!file.exists(result_dir)){
     dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
   }

   # Sample distributions
   object %>% autonomics.plot::plot_sample_distributions(
                   color_var    = color_var,
                   color_values = color_values,
                   ylab         = 'sample distributions',
                   file         = paste0(result_dir, '/sample_distributions.pdf'))
   # # Preprocess
   # object %<>% autonomics.preprocess::normalize_smart(plot       = TRUE,
   #                                                      result_dir = result_dir,
   #                                                      color_var  = color_var,
   #                                                      shape_var  = shape_var,
   #                                                      txt_var    = txt_var)

   # PCA
   message('\tPerform PCA')
   if (do_pca){
      object %<>% autonomics.explore::add_and_write_pca(
                       result_dir    = result_dir,
                       feature_plots = feature_plots,
                       color_var     = color_var,
                       color_values  = color_values,
                       shape_var     = shape_var,
                       x             = x,
                       group_var     = group_var,
                       txt_var       = txt_var,
                       line          = line)
   }

   # Differential expression
   # replicates <- !is.null(object$subgroup) && anyDuplicated(object$subgroup)
   object %<>% autonomics::analyze_contrasts(
                    result_dir     = result_dir,
                    design         = design,
                    contrasts      = contrasts,
                    top_definition = top_definition,
                    universe       = universe,
                    feature_plots  = feature_plots,
                    x              = x,
                    color_var      = color_var,
                    color_values   = color_values,
                    shape_var      = shape_var,
                    group_var      = group_var,
                    facet_var      = facet_var,
                    line           = line,
                    fvars          = fvars,
                    nplot          = nplot,
                    feature_plot_height = feature_plot_height,
                    feature_plot_width  = feature_plot_width,
                    cluster_features = cluster_features,
                    min_set_size   = min_set_size,
                    max_set_size   = max_set_size)

   # Save eSet. Zip results. Return eSet.
   saveRDS(object, sprintf('%s/object.rds', result_dir))
   result_dir %>% autonomics.support::save_session_info()
   if (zip_results)   result_dir %>% autonomics.support::zip_dir()
   return(object)
}
