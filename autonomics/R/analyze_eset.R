####
#' @title Analyze SummarizeExperiment
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
#' @param object            SummarizedExperiment
#' @param design            design matrix
#' @param contrasts         contrast definitions
#' @param topdef            Definition of 'top features'.
#' @param result_dir        directory to which to write results
#' @param universe          'detectome', 'genome', or NULL (no ora)
#' @param do_pca            whether to perform PCA analysis (logical)
#' @param cluster_features  whether to cluster features (logical)
#' @param nplot             no of top features to be plotted
#' @param feature_plot_width in inches
#' @param feature_plot_height in inches
#' @param min_set_size      min set size (for ora)
#' @param max_set_size      max set size (for ora)
#' @param zip_results       zip analysis results (logical)?
#' @param ...               passed to autonomics.plot::plot_features(...)
#' @return Updated SummarizedExperiment
#' @examples
#'
#' require(magrittr)
#' # Compare each subgroup level to zero
#' #------------------------------------
#' if (require(autonomics.data)){
#'    object <- 'extdata/stemcomp/maxquant/proteinGroups.txt' %>% 
#'               system.file(package = 'autonomics.data')     %>% 
#'               autonomics::read_proteingroups()             %>%
#'               autonomics::prepare_proteingroups(invert_subgroups = c('E_BM', 'EM_BM', 'E_EM'))
#'               
#'    result_dir <- tempdir() %>% paste0('/analysis.results') %T>% message()
#'    dir.create(result_dir, showWarnings = FALSE)
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
  design           = autonomics.import::create_design_matrix(object),
  contrasts        = autonomics.find::default_contrasts(object),
  topdef           = autonomics.find::default_topdef(object),
  universe         = NULL, #autonomics.ora::default_universe(object),
  dodge_width      = 0,
  do_pca           = TRUE,
  cluster_features = default_cluster_features(),
  nplot            = autonomics.find::default_nplot(object),
  feature_plot_width  = NULL, 
  feature_plot_height = NULL,
  min_set_size        = 5,
  max_set_size        = Inf,
  zip_results         = FALSE, 
  ...
){
   # Assert valid args
   autonomics.import::assert_is_valid_object(object)
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
      object %<>% autonomics.plot::add_and_write_pca(result_dir = result_dir, ...)
   }

   # Analyze contrasts
   object %<>% autonomics::analyze_contrasts(result_dir          = result_dir,
                                             design              = design,
                                             contrasts           = contrasts,
                                             topdef              = topdef,
                                             universe            = universe,
                                             nplot               = nplot,
                                             feature_plot_height = feature_plot_height,
                                             feature_plot_width  = feature_plot_width,
                                             cluster_features    = cluster_features,
                                             min_set_size        = min_set_size,
                                             max_set_size        = max_set_size, 
                                             ...)

   # Save eSet. Zip results. Return eSet.
   saveRDS(object, sprintf('%s/object.rds', result_dir))
   result_dir %>% autonomics.support::save_session_info()
   if (zip_results)   result_dir %>% autonomics.support::zip_dir()
   return(object)
}
