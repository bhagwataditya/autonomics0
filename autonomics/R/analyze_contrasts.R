
#' @rdname analyze_eset
#' @importFrom magrittr %>%  %<>%
#' @export
analyze_contrasts <- function(
   object,
   result_dir       = default_result_dir(object),
   design           = autonomics.import::create_design_matrix(object),
   contrasts        = autonomics.find::default_contrasts(object),
   direction        = c('neg', 'pos'),
   topdef           = autonomics.find::default_topdef(object),
   universe         = NULL,#autonomics.ora::default_universe(object),
   cluster_features = autonomics::default_cluster_features(),
   nplot            = autonomics.find::default_nplot(object),
   feature_plot_width  = NULL,
   feature_plot_height = NULL,
   min_set_size     = 5,
   max_set_size     = Inf, 
   ...
){

  # Assert valid inputs
  autonomics.import::assert_is_valid_object(object)
  assertive.types::assert_is_character(result_dir)

  # Initialize
  current_theme <- ggplot2::theme_get()
  ggplot2::theme_set(ggplot2::theme_bw())
  on.exit(ggplot2::theme_set(current_theme))

  # Validify contrasts
  contrasts %<>% autonomics.find::validify_contrasts(design)
  if (assertive.properties::is_empty(contrasts))   return(object)

  # Run limma analysis
  message('\tRun significance analysis (limma)')
  object %<>% autonomics.find::add_limma(contrasts, design)

  # Create dirs
  subdirs <- autonomics.find::get_contrast_subdir(result_dir, names(contrasts))
  subdirs %>% plyr::l_ply(dir.create, showWarnings = FALSE)

  # Write to file
  message('\tWrite results to file')
  file <- paste0(result_dir, '/features.txt')
  autonomics.import::write_features(object, file)

  # Plot top features
  message('\tPlot features')
  object %>% autonomics.find::plot_per_contrast_features(
                 contrasts      = contrasts,
                 design         = design,
                 topdef         = topdef,
                 result_dir     = result_dir,
                 nplot          = nplot, 
                 width          = feature_plot_width, 
                 height         = feature_plot_height, 
                 ...)

  # Analyse over representation
  for (cur_universe in universe){
     object %>% autonomics.ora::run_ora_on_eset(
                     contrasts,
                     result_dir     = result_dir,
                     topdef         = topdef,
                     universe       = cur_universe,
                     min_set_size   = min_set_size,
                     max_set_size   = max_set_size)
  }


  # Run cluster analysis
  if (cluster_features){
     object %>% autonomics.find::cluster_top_features_on_subgroups(
                   contrasts    = contrasts,
                   result_dir   = result_dir,
                   ...)
  }

  # Return
  return(object)

}
