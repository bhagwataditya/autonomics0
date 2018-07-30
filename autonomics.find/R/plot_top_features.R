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
#' @param geom             'point', 'boxplot', 'violin', 'bar', 'hbar'
#' @param fvars            fvars to use in plot
#' @param nplot            max no of top features to plot
#' @param ...              passed to autonomics.plot::plot_feature_xxx()
#' @examples
#' require(magrittr)
#' 
#' # STEMCELL COMPARISON
#'   file <- tempfile() %>% paste0('.pdf') %T>% message()
#'   if (require(autonomics.data)){
#'      object <- autonomics.data::stemcomp.proteinratios
#'      object %>% plot_top_features(geom = 'hbar')
#'      object %>% plot_top_features(geom = 'bar')
#'      object %>% plot_top_features(geom = 'point')
#'      object %>% plot_top_features(geom = 'violin')
#'      object %>% plot_top_features(geom = 'boxplot')
#'      object %>% plot_top_features(direction  = 'neg', x = 'subgroup')
#'   }
#'
#'# GLUTAMINASE
#'   if (require(autonomics.data)){
#'      object <- autonomics.data::glutaminase
#'      object %>% plot_top_features(top_definition = 'bonf < 0.05 & rank <= 4',
#'                                   direction      = 'both',
#'                                   geom           = 'boxplot')
#'   }
#'
#' # A somascan eset
#' if (require(atkin.2014)){
#'    object  <- atkin.2014::soma
#'    contrasts <- atkin.2014::contrasts
#'    object %>% plot_top_features(contrast = contrasts[1], direction = 'neg', nplot = 9,
#'                                 geom = 'point', x = 'time', color_var = 'condition',
#'                                 facet_var = 'subject_id', fvars = 'TargetFullName')
#' }
#'
#' # A metabolon eset
#' if (require(subramanian.2016)){
#'    contrasts <- subramanian.2016::contrasts.metabolon[1]
#'    object  <- subramanian.2016::metabolon %>% add_limma_to_fdata(contrasts = contrasts)
#'    object %>% plot_top_features(contrast  = contrasts,
#'                                 direction = 'neg',
#'                                 geom      = 'boxplot',
#'                                 color_var = 'condition',
#'                                 group_var = 'condition',
#'                                 line      = TRUE,
#'                                 nplot     = 4)
#' }
#'
#' if (require(billing.differentiation.data)){
#'    object <- billing.differentiation.data::rna.voomcounts
#'    contrast <- billing.differentiation.data::contrasts[1]
#'    object %>% autonomics.find::plot_top_features(contrast = contrast, 
#'                                                  nplot    = 4)
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
   geom           = autonomics.plot::default_feature_plots(object)[1],
   nplot          = autonomics.find::default_nplot(object), 
   fvars          = autonomics.plot::default_fvars(object),   
   ...
){
  # Process and check args
  assertive.sets::assert_is_subset(geom, autonomics.plot::FEATURE_PLOTS)
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
  if (geom == 'hbar'){
     autonomics.import::fdata(top)[['plot.annot']] <- top %>% collapse_fvars(rev(fvars))
     fvars <- 'plot.annot'
  }

  # Add p values
  if (pvar_in_fdata(object, contrast)){
     autonomics.import::fdata(top)[['plot.p'    ]] <- top %>% format_pvalues(contrast)
     fvars %<>% c('plot.p')
  }

  # Horizontal feature bars
  if (geom == 'hbar'){ facet_def <- if('replicate' %in% autonomics.import::svars(object)) '~ subgroup + replicate'  else '~ sample'
                                top %>% autonomics.plot::plot_feature_hbars(alpha_var = 'plot.alpha', title = my_title, fvars = fvars, xlab = 'log2 ratio', ...)
  } else {                      top %>% autonomics.plot::plot_features(     alpha_var = 'plot.alpha', title = my_title, fvars = fvars, geom = geom, ...)
  }

}

#' Plot top features
#'
#' Plot top features per contrast
#' @param object          SummarizedExperiment
#' @param design          design matrix
#' @param contrasts       named contrast vector
#' @param result_dir      directory where to store results
#' @param geoms           which feature plots to be created?
#' @param ...             passed to autonomics.plot::plot_features
#' @examples
#' require(magrittr)
#' if (require(atkin.2014)){
#'    object <- atkin.2014::soma
#'    contrasts <- atkin.2014::contrasts
#'    result_dir <- tempdir() %T>% message()
#'    object %>% plot_top_features_all_contrasts(
#'                    contrasts = contrasts[1], 
#'                     result_dir = result_dir)
#'    object %>% plot_top_features_all_contrasts(
#'                    contrasts  = contrasts[1], 
#'                    result_dir = result_dir,
#'                    nplot      = 9, 
#'                    geom       = 'point', 
#'                    x          = 'time',
#'                    color_var  = 'condition',
#'                    facet_var  = 'subject_id',
#'                    fvars      = 'TargetFullName')
#' }
#' # GLUTAMINASE
#' if (require(autonomics.data)){
#'   object <- autonomics.data::glutaminase
#'   result_dir <- tempdir() %T>% message()
#'   object %>% plot_top_features_all_contrasts(
#'                 contrasts      = autonomics.import::contrastdefs(.)[1:2],
#'                 top_definition = 'fdr < 0.05',
#'                 x              = 'TIME_POINT',
#'                 geom           = 'boxplot',
#'                 result_dir     = result_dir)
#'}
#' @importFrom magrittr  %>%
#' @export
plot_top_features_all_contrasts <- function(
   object,
   design         = autonomics.find::create_design_matrix(object),
   contrasts      = autonomics.find::default_contrasts(object),
   result_dir,
   geoms          = autonomics.plot::default_feature_plots(object),
   ...
){
  for (i in seq_along(contrasts)){
    contrast <- contrasts[i] %>% magrittr::set_names(names(contrasts)[i])
    subdir   <- autonomics.find::get_contrast_subdir(result_dir, names(contrast))
    dir.create(subdir, recursive = TRUE, showWarnings = FALSE)
    for (direction in c('neg', 'pos')){
       for (i_plot in seq_along(geoms)){
          cur_geom <- geoms[[i_plot]]
          #if (length(x) > 1)   x %<>% magrittr::extract2(i_plot)
          my_file <- sprintf('%s/top_%s__%s__%s.pdf', subdir, cur_geom, names(contrast), direction)
          autonomics.support::cmessage('\t\t%s %s 0   %s',
                                        contrast,
                                        autonomics.find::direction_to_sign(direction), basename(my_file))
          object %>% autonomics.find::plot_top_features( design    = design,
                                                         contrast  = contrast,
                                                         direction = direction,
                                                         geom      = cur_geom,
                                                         file      = my_file, 
                                                         ...)
       }
    }
  }
}
