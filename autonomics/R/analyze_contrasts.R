
#' @rdname analyze_eset
#' @importFrom magrittr %>%  %<>%
#' @export
analyze_contrasts <- function(
   object,
   result_dir       = default_result_dir(object),
   design           = autonomics.find::create_design_matrix(object),
   contrasts        = autonomics.find::default_contrasts(object),
   top_definition   = autonomics.find::default_top_definition(object),
   universe         = autonomics.ora::default_universe(object),
   feature_plots    = autonomics.plot::default_feature_plots(object),
   x                = autonomics.plot::default_x(object, feature_plots[1]),
   color_var        = autonomics.plot::default_color_var(object),
   color_values     = autonomics.plot::default_color_values(object, color_var),
   shape_var        = autonomics.plot::default_shape_var(object),
   group_var        = autonomics.plot::default_group_var(object),
   facet_var        = autonomics.plot::default_facet_var(),
   line             = autonomics.plot::default_line(object),
   fvars            = autonomics.plot::default_fvars(object),
   cluster_features = autonomics::default_cluster_features(),
   nplot            = autonomics.find::default_nplot(object),
   feature_plot_width  = NULL,
   feature_plot_height = NULL,
   min_set_size     = 5,
   max_set_size     = Inf
){

  # Assert valid inputs
  autonomics.import::assert_is_valid_eset(object)
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
  object %<>% autonomics.find::add_limma_to_fdata(contrasts, design, overwrite = TRUE)

  # Create dirs
  subdirs <- autonomics.find::get_contrast_subdir(result_dir, names(contrasts))
  subdirs %>% plyr::l_ply(dir.create, showWarnings = FALSE)

  # Write to file
  message('\tWrite results to file')
  object %>% autonomics.find::write_features(design, contrasts, top_definition, result_dir)

  # Plot top features
  message('\tPlot top features')
  object %>% autonomics.find::plot_top_features_all_contrasts(
     design         = design,
     contrasts      = contrasts,
     result_dir     = result_dir,
     top_definition = top_definition,
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
     width  = feature_plot_width, 
     height = feature_plot_height)

  # Run sigb analysis
  # if(inherits(object, "ProtSet")){
  #    message('\tRun outlier analysis (sigb)')
  #    object %<>% autonomics.find::add_sigb_to_fdata(contrasts)
  # }

  # Analyse over representation
  for (cur_universe in universe){
     object %>% autonomics.ora::run_ora_on_eset(
                     contrasts,
                     result_dir     = result_dir,
                     top_definition = top_definition,
                     universe       = cur_universe,
                     min_set_size   = min_set_size,
                     max_set_size   = max_set_size)
  }


  # Run cluster analysis
  if (cluster_features){
     object %>% autonomics.find::cluster_top_features_on_subgroups(
        contrasts    = contrasts,
        result_dir   = result_dir,
        x            = x,
        color_var    = color_var,
        color_values = color_values,
        shape_var    = shape_var,
        group_var    = group_var,
        line         = line)
  }

  # Return
  return(object)

}
