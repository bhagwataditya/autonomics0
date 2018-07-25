#' @importFrom magrittr   %>%
create_feature_plot_title <- function(object, contrast, top_definition, direction){

  # Add no of top features
  n_top      <- sum(autonomics.find::are_top_features(object, top_definition, names(contrast), direction))
  down_or_up <- autonomics.find::direction_to_sign(direction)
  my_title <- sprintf("%s %s 0   |   %d %s", contrast, down_or_up, n_top, top_definition)

  # Add no of directional features
  coefvar        <- paste0('coef.',    names(contrast))
  coef_available <- coefvar %in% autonomics.import::fvars(object)
  if (coef_available){
     n_direction <- sum(autonomics.find::are_top_features(object, '1', names(contrast), direction))
     my_title <- sprintf('%s   |   %d %s', my_title, n_direction, autonomics.find::direction_to_sign(direction))
  }

  # Add no of significant features
  pvar        <- paste0('p.',    names(contrast))
  p_available <- pvar %in% autonomics.import::fvars(object)
  if (p_available){
     n_p    <- sum(autonomics.find::are_top_features(object, 'p   < 0.05',   names(contrast), direction))
     n_fdr  <- if (p_available)  sum(autonomics.find::are_top_features(object, 'fdr < 0.05',   names(contrast), direction)) else ''
     n_bonf <- if (p_available)  sum(autonomics.find::are_top_features(object, 'bonf < 0.05',  names(contrast), direction)) else ''
     my_title <- sprintf('%s  %d p   %d fdr   %d bonf (0.05)', my_title, n_p, n_fdr, n_bonf)
  }

  return(my_title)
}

infer_alpha <- function(object, design, contrast){
  contrast_mat <- create_contrast_matrix(design, contrast)
  are_relevant_samples(object  = object,
                       design    = design,
                       contrasts = contrast)
}

#' Collapse fvars
#' @param object   eset
#' @param fvars      fvars
#' @importFrom magrittr   %>%
#' @export
collapse_fvars <- function(object, fvars){

   # Assert
   autonomics.import::assert_is_valid_eset(object)
   assertive.sets::assert_is_subset(fvars, autonomics.import::fvars(object))

   # Extract
   fdf <- autonomics.import::fdata(object)     %>%
          magrittr::extract(fvars)             %>%
          lapply(as.character)                 %>%
          lapply(function(x){
                    x[is.na(x)] <- ''
                    x
                 }
          ) %>%
          data.frame(stringsAsFactors = FALSE)

   # Collapse
   do.call(paste, c(fdf, sep = '   ')) %>% trimws()
}

pvar_in_fdata <- function(object, contrast){
   get_pvar(object, contrast) %in% autonomics.import::fvars(object)
}

get_pvar <- function(object, contrast){
   sprintf('p.%s', names(contrast))
}

get_fdrvar <- function(object, contrast){
   sprintf('fdr.%s', names(contrast))
}

format_pvalues <- function(object, contrast){
   pvar   <- get_pvar(object, contrast)
   fdrvar <- get_fdrvar(object, contrast)
   if (pvar %in% autonomics.import::fvars(object)){
      pvalues   <- autonomics.import::fdata(object)[[pvar]]
      fdrvalues <- autonomics.import::fdata(object)[[fdrvar]]
      autonomics.import::fdata(object)[[pvar]] <-  sprintf('fdr = %2.1e  (p = %2.1e)', fdrvalues, pvalues)
   }
}

format_sigbvalues <- function(object, contrast){
   # Format sigb values
   sigbvar  <- sprintf('sigb.%s', names(contrast))
   if (sigbvar %in% autonomics.import::fvars(object)){
      autonomics.import::fdata(object)[[sigbvar]] %<>% sprintf('   sigb %2.1e', .)
   }
}



#' Plot top features
#'
#' This function visualizes the top features for the specified contrast.
#'
#' @param object \code{eSet}
#' @param design           design matrix
#' @param contrast         named contrast for which to select the top feature bars
#' @param top_definition   Definition of 'top features'.
#' @param direction        any value in autonomics.find::DIRECTIONS
#' @param xlab             passed to plotting
#' @param feature_plot     value in \code{\link[autonomics.plot]{FEATURE_PLOTS}}
#' @param x                svar mapped to x
#' @param color_var        svar mapped to color
#' @param color_values     color value vector (names = subgroups, contents = colors)
#' @param shape_var        svar mapped to shape
#' @param group_var        svar mapped to group
#' @param txt_var          svar mapped to txt
#' @param facet_var        svar used for faceting
#' @param line             whether to add line (logical)
#' @param fvars            fvars to use in plot
#' @param nplot            max no of top features to plot
#' @param file             file to which fesults are written
#' @param width            width in inches
#' @param height           height in inches
#' @examples
#' require(magrittr)
#' 
#' # STEMCELL COMPARISON
#'   file <- tempfile() %>% paste0('.pdf') %T>% message()
#'   if (require(autonomics.data)){
#'      object <- autonomics.data::stemcomp.proteinratios
#'      contrasts <- object %>% autonomics.find::default_contrasts()
#'      object %<>% autonomics.find::add_limma_to_fdata()
#'      object %>% plot_top_features(contrast = contrasts[1], feature_plot = 'hbars', file = file)
#'      object %>% plot_top_features(contrast = contrasts[1], feature_plot = 'bars', file = file)
#'      object %>% plot_top_features(contrast = contrasts[1], feature_plot = 'profiles', file = file)
#'      object %>% plot_top_features(contrast = contrasts[1], feature_plot = 'distributions', file = file)
#'      object %>% plot_top_features(contrast = contrasts[1], feature_plot = 'boxes', file = file)
#'      object %>% plot_top_features(direction  = 'neg', x = 'subgroup', file = file)
#'   }
#'
#'# GLUTAMINASE
#'   if (require(autonomics.data)){
#'      object <- autonomics.data::glutaminase
#'      object %>% autonomics.find::plot_top_features(
#'                    contrast       = halama.2016::contrasts[1],
#'                    top_definition = 'bonf < 0.05',
#'                    direction      = 'both',
#'         color_var      = 'GROUP_DESCRIPTION',
#'       # color_values   = c(Control = 'red', Vehicle = 'orange',
#'       #                   `Concentration 1` = 'green', `Concentration 2` = 'forestgreen'),
#'         x              = 'TIME_POINT',
#'         feature_plot   = 'boxes', 
#'         file = file)
#'   }
#'
#' # A somascan eset
#' if (require(atkin.2014)){
#'    object  <- atkin.2014::soma
#'    contrasts <- atkin.2014::contrasts
#'    object %>% autonomics.find::plot_top_features(
#'                    contrast = contrasts[1], direction = 'neg', nplot = 9,
#'                    feature_plot = 'profiles', x = 'time', color_var = 'condition',
#'                    facet_var = 'subject_id', fvars = 'TargetFullName', file = file)
#'    object %>% autonomics.find::plot_top_features(
#'                    contrast = contrasts[10], direction = 'neg', nplot = 9,
#'                    feature_plot = 'profiles', x = 'time', color_var = 'condition',
#'                    fvars = 'TargetFullName')
#'    object <- atkin.2014::soma.platelets %>%
#'                autonomics.import::filter_samples(time %in% c('t0', 't2'))
#'    object %>% autonomics.find::plot_top_features(
#'                    contrast  = atkin.2014::contrasts['t2_t0.D_C'],
#'                    color_var = 'condition',
#'                    fvars     = 'TargetFullName',
#'                    file      = file)
#' }
#'
#' # A metabolon eset
#' if (require(subramanian.2016)){
#'    contrasts <- subramanian.2016::contrasts.metabolon[1]
#'    object  <- subramanian.2016::metabolon %>%
#'                 autonomics.find::add_limma_to_fdata(contrasts = contrasts)
#'    object %>% autonomics.find::plot_top_features(
#'                    contrast      = contrasts,
#'                    direction     = 'neg',
#'                    feature_plot  = 'boxes',
#'                    color_var     = 'condition',
#'                    group_var     = 'condition',
#'                    line          = TRUE,
#'                    nplot         = 16, 
#'                    file          = file)
#' }
#'
#' if (require(billing.differentiation.data)){
#'    object <- billing.differentiation.data::rna.voomcounts
#'    contrast <- billing.differentiation.data::contrasts[1]
#'    object %>% autonomics.find::plot_top_features(contrast = contrast)
#' }
#'
#'
#' @author Aditya Bhagwat
#' @importFrom magrittr   %>%   %<>%
#' @export
plot_top_features <- function(
   object,
   design         = autonomics.find::create_design_matrix(object),
   contrast       = autonomics.find::default_contrasts(object)[1],
   top_definition = autonomics.find::default_top_definition(object),
   direction      = autonomics.find::DIRECTIONS[1],
#   result_dir     = NULL,
   xlab           = '',
   feature_plot   = autonomics.plot::default_feature_plots(object)[1],
   x              = autonomics.plot::default_x(object, feature_plot),
   color_var      = autonomics.plot::default_color_var(object),
   color_values   = autonomics.plot::default_color_values(object),
   shape_var      = autonomics.plot::default_shape_var(object),
   group_var      = autonomics.plot::default_group_var(object),
   txt_var        = autonomics.plot::default_txt_var(object),
   facet_var      = autonomics.plot::default_facet_var(),
   line           = autonomics.plot::default_line(object),
   fvars          = autonomics.plot::default_fvars(object),
   nplot          = autonomics.find::default_nplot(object),
   file           = NULL,
   width          = NULL,
   height         = NULL
){
  # Process and check args
  assertive.sets::assert_is_subset(feature_plot, autonomics.plot::FEATURE_PLOTS)
  assertive.strings::assert_is_a_non_empty_string(direction)
  assertive.sets::assert_is_subset(direction, autonomics.find::DIRECTIONS)
  assertive.base::assert_is_identical_to_true(autonomics.find::is_valid_contrast(contrast, design))

  # Limit eset to top features for chosen direction (abort if none)
  top <- object %>% autonomics.find::filter_n_arrange_top_features(names(contrast), top_definition, direction, nplot)
  if (nrow(top)==0){
     autonomics.support::cmessage('\t\t%s %s 0   no top features - abort',
                                  contrast, autonomics.find::direction_to_sign(direction))
     return(NULL)
  }

  # Prepare title & subdir
  my_title <- create_feature_plot_title(object, contrast, top_definition, direction)

  # Collapse fvars for hbars
  autonomics.import::sdata(top)[['plot.alpha']] <- top %>% infer_alpha(design, contrast)
  if (feature_plot == 'hbars'){
     autonomics.import::fdata(top)[['plot.annot']] <- top %>% collapse_fvars(rev(fvars))
     fvars <- 'plot.annot'
  }

  # Add p values
  if (pvar_in_fdata(object, contrast)){
     autonomics.import::fdata(top)[['plot.p'    ]] <- top %>% format_pvalues(contrast)
     fvars %<>% c('plot.p')
  }

  # Horizontal feature bars
  if (feature_plot == 'hbars'){
    facet_def <- if('replicate' %in% autonomics.import::svars(object)) '~ subgroup + replicate'  else '~ sample'
    top %>% autonomics.plot::plot_feature_hbars(fvars        = fvars,
                                                color_var    = color_var,
                                                color_values = color_values,
                                                facet_def    = facet_def,
                                                alpha_var    = 'plot.alpha',
                                                xlab         = 'log2 ratio',
                                                title        = my_title,
                                                file         = file,
                                                width        = width,
                                                height       = height)
  } else {
    top %>% autonomics.plot::plot_features(x            = x,
                                           color_var    = color_var,
                                           color_values = color_values,
                                           shape_var    = shape_var,
                                           group_var    = group_var,
                                           txt_var      = txt_var,
                                           facet_var    = facet_var,
                                           alpha_var    = 'plot.alpha',
                                           line         = line,
                                           title        = my_title,
                                           fvars        = fvars,
                                           file         = file,
                                           width        = width,
                                           height       = height,
                                           feature_plot = feature_plot)
  }

}

#' Plot top features
#'
#' Plot top features per contrast
#' @param object       eset, as returned by \link{add_limma_to_fdata}
#' @param design         design matrix
#' @param contrasts      named contrast vector
#' @param result_dir     directory where to store results
#' @param top_definition definition of top features
#' @param feature_plots  which feature plots to be created?
#' @param x              svar mapped to x in feature plots
#' @param color_var      svar mapped to color in feature plots
#' @param color_values   color vector (names = subgroups, values = colors)
#' @param shape_var      svar mapped to shape in feature plots
#' @param group_var      svar mapped to group in feature plots
#' @param txt_var        svar mapped to txt   in feature plots
#' @param facet_var      svar used to facet feature plots
#' @param line           whether to connect points in feature plots (logical)
#' @param fvars          fvars to use in plot
#' @param nplot          number of features to plot
#' @param width          figure width (inches)
#' @param height         figure height (inches)
#' @examples
#' require(magrittr)
#' if (require(atkin.2014)){
#'    object <- atkin.2014::soma
#'    contrasts <- atkin.2014::contrasts
#'    result_dir <- tempdir() %T>% message()
#'    object %>% autonomics.find::plot_top_features_all_contrasts(
#'                    contrasts = contrasts[1], result_dir = result_dir)
#'    object %>% autonomics.find::plot_top_features_all_contrasts(
#'                    contrasts = contrasts[1], result_dir = result_dir,
#'                    nplot = 9, feature_plot = 'profiles', x = 'time',
#'                    color_var = 'condition', facet_var = 'subject_id',
#'                    fvars = 'TargetFullName')
#' }
#'if (require(halama.2016)){
#'   object <- halama.2016::cell.metabolites
#'   result_dir <- tempdir() %T>% message()
#'   object %>% autonomics.find::plot_top_features_all_contrasts(
#'      contrast = halama.2016::contrasts[1],
#'      top_definition = 'fdr < 0.05',
#'      color_var = 'GROUP_DESCRIPTION',
#'      color_values = c(Control = 'red', Vehicle = 'orange',
#'                       `Concentration 1` = 'green', `Concentration 2` = 'forestgreen'),
#'      x = 'TIME_POINT',
#'      feature_plot = 'boxes',
#'      result_dir = result_dir)
#'}
#' @importFrom magrittr  %>%
#' @export
plot_top_features_all_contrasts <- function(
   object,
   design         = autonomics.find::create_design_matrix(object),
   contrasts      = autonomics.find::default_contrasts(object),
   result_dir,
   top_definition = autonomics.find::default_top_definition(object),
   feature_plots  = autonomics.plot::default_feature_plots(object),
   x              = autonomics.plot::default_x(object, feature_plots[1]),
   color_var      = autonomics.plot::default_color_var(object),
   color_values   = autonomics.plot::default_color_values(object),
   shape_var      = autonomics.plot::default_shape_var(object),
   group_var      = autonomics.plot::default_group_var(object),
   txt_var        = autonomics.plot::default_txt_var(object),
   facet_var      = autonomics.plot::default_facet_var(),
   line           = autonomics.plot::default_line(object),
   fvars          = autonomics.plot::default_fvars(object),
   nplot          = autonomics.find::default_nplot(object),
   width          = NULL,
   height         = NULL
){
  for (i in seq_along(contrasts)){
    contrast <- contrasts[i] %>% magrittr::set_names(names(contrasts)[i])
    subdir   <- autonomics.find::get_contrast_subdir(result_dir, names(contrast))
    dir.create(subdir, recursive = TRUE, showWarnings = FALSE)
    for (direction in c('neg', 'pos')){
       for (i_plot in seq_along(feature_plots)){
          cur_plot <- feature_plots[[i_plot]]
          if (length(x) > 1)   x %<>% magrittr::extract2(i_plot)
          my_file <- sprintf('%s/top_%s__%s__%s.pdf', subdir, cur_plot, names(contrast), direction)
          autonomics.support::cmessage('\t\t%s %s 0   %s',
                                        contrast,
                                        autonomics.find::direction_to_sign(direction), basename(my_file))
          object %>% autonomics.find::plot_top_features(  design         = design,
                                                            contrast       = contrast,
                                                            top_definition = top_definition,
                                                            direction      = direction,
                                                            feature_plot   = cur_plot,
                                                            x              = x,
                                                            color_var      = color_var,
                                                            color_values   = color_values,
                                                            shape_var      = shape_var,
                                                            group_var      = group_var,
                                                            txt_var        = txt_var,
                                                            facet_var      = facet_var,
                                                            line           = line,
                                                            fvars          = fvars,
                                                            nplot          = nplot,
                                                            file           = my_file,
                                                            width          = width,
                                                            height         = height)
       }
    }
  }
}
