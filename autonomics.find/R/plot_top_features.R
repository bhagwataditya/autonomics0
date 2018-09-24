#' @importFrom magrittr   %>%
create_feature_plot_title <- function(object, contrast, topdef, direction){
   
   # Note: I decided to use contrast names in the title, since some composite contrasts can be really huge
   
   limmamatrix <- autonomics.import::limma(object)[, names(contrast), ]
   effect <- limmamatrix  %>%  magrittr::extract(, 'effect')
   p      <- limmamatrix  %>%  magrittr::extract(, 'p'     )
   fdr    <- limmamatrix  %>%  magrittr::extract(, 'fdr'   )
   bonf   <- limmamatrix  %>%  magrittr::extract(, 'bonf'  )
   
   limmadf <- limmamatrix %>% data.table::data.table()
   
   n.total  <- nrow(object)
   n.top <- nrow(limmadf[eval(parse(text=topdef))])
   
   n.up      <- sum(effect > 0, na.rm = TRUE)
   n.up.p    <- sum(effect > 0  &  p    < 0.05, na.rm=TRUE)
   n.up.fdr  <- sum(effect > 0  &  fdr  < 0.05, na.rm=TRUE)
   n.up.bonf <- sum(effect > 0  &  bonf < 0.05, na.rm=TRUE)

   n.down      <- sum(effect < 0, na.rm = TRUE)
   n.down.p    <- sum(effect < 0  &  p    < 0.05, na.rm=TRUE)
   n.down.fdr  <- sum(effect < 0  &  fdr  < 0.05, na.rm=TRUE)
   n.down.bonf <- sum(effect < 0  &  bonf < 0.05, na.rm=TRUE)
   
   n.top      %<>% stringi::stri_pad(width = max(nchar(c(n.down,     n.up,   n.top)))    )
   n.down     %<>% stringi::stri_pad(width = max(nchar(c(n.down,     n.up,   n.top)))    )
   n.up       %<>% stringi::stri_pad(width = max(nchar(c(n.down,     n.up,   n.top)))    )
   n.down.p   %<>% stringi::stri_pad(width = max(nchar(c(n.down.p,   n.up.p)))  )
   n.up.p     %<>% stringi::stri_pad(width = max(nchar(c(n.down.p,   n.up.p)))  )
   n.down.fdr %<>% stringi::stri_pad(width = max(nchar(c(n.down.fdr, n.up.fdr))))
   n.up.fdr   %<>% stringi::stri_pad(width = max(nchar(c(n.down.fdr, n.up.fdr))))
   
   s0 <- sprintf('%s : %s %s', names(contrast), n.top, topdef)
   spaces <- paste0(rep(' ', nchar(names(contrast))), collapse = '')
   s1 <- sprintf('%s down > %s p > %s fdr > %s bonf', n.down, n.down.p, n.down.fdr, n.down.bonf)
   s2 <- sprintf('%s up   > %s p > %s fdr > %s bonf', n.up,   n.up.p,   n.up.fdr,   n.up.bonf)
   
   s <- sprintf('%s\n\n%s\n%s\n', s0, s1, s2) #%>% cat()

  return(s)
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


#' Plot top features
#'
#' This function visualizes the top features for the specified contrast.
#'
#' @param object \code{eSet}
#' @param design    design matrix
#' @param contrast  named contrast for which to select the top feature bars
#' @param topdef    definition of 'top features'.
#' @param geom     'point', 'boxplot', 'violin', 'bar', 'hbar'
#' @param fvars     fvars to use in plot
#' @param nplot     max no of top features to plot
#' @param ...       passed to autonomics.plot::plot_feature_xxx()
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
#'      object %>% plot_top_features(x = 'subgroup')
#'   }
#'
#'# GLUTAMINASE
#'   if (require(autonomics.data)){
#'      object <- autonomics.data::glutaminase
#'      object %>% plot_top_features(topdef = 'bonf < 0.05 & rank <= 4',
#'                                   geom           = 'boxplot')
#'   }
#'
#' # A somascan eset
#' if (require(atkin.2014)){
#'    object  <- atkin.2014::soma
#'    contrasts <- atkin.2014::contrasts
#'    object %>% plot_top_features(contrast = contrasts[1], nplot = 9,
#'                                 geom = 'point', x = 'time', color_var = 'condition',
#'                                 facet_var = 'subject_id', fvars = 'TargetFullName')
#' }
#'
#' # A metabolon eset
#' if (require(subramanian.2016)){
#'    contrasts <- subramanian.2016::contrasts.metabolon[1]
#'    object  <- subramanian.2016::metabolon %>% add_limma_to_fdata(contrasts = contrasts)
#'    object %>% plot_top_features(contrast  = contrasts,
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
   design         = autonomics.import::create_design_matrix(object),
   contrast       = autonomics.find::default_contrasts(object)[1],
   topdef         = autonomics.find::default_topdef(object),
   geom           = autonomics.plot::default_feature_plots(object)[1],
   nplot          = autonomics.find::default_nplot(object), 
   fvars          = autonomics.plot::default_fvars(object),
   color_var      = autonomics.plot::default_color_var(object),
   ...
){
  # Process and check args
  assertive.sets::assert_is_subset(geom, autonomics.plot::FEATURE_PLOTS)
  assertive.base::assert_is_identical_to_true(autonomics.find::is_valid_contrast(contrast, design))
  assertive.properties::is_not_null(rownames(autonomics.import::limma(object)))
  assertive.base::are_identical(autonomics.import::fnames(object), rownames(autonomics.import::limma(object)))
  

  # Prepare title & subdir
  my_title <- autonomics.find:::create_feature_plot_title(object, contrast, topdef)
  
  # Limit object to top features (abort if none)
  top <- object %>% autonomics.find::filter_n_arrange_top_features(contrast_name = names(contrast), topdef = topdef, nmax = nplot)
  if (nrow(top)==0){
     autonomics.support::cmessage('\t\tAbort, since no %s features: %s', names(contrast), topdef)
     return(NULL)
  }

  # Set alpha
  composite_colors_required <- top %>% autonomics.import::slevels(color_var) %>% autonomics.import::has_two_components()
  autonomics.import::sdata(top)[['plot.alpha']] <- if (composite_colors_required){ TRUE 
                                                   } else {                        top %>% autonomics.find::are_relevant_samples(design, contrast) }
  
  # Collapse fvars for hbars
  if (geom == 'hbar'){
     autonomics.import::fdata(top)[['plot.annot']] <- top %>% collapse_fvars(rev(fvars))
     fvars <- 'plot.annot'
  }

  # Add p & fdr values
  limma <- top %>% autonomics.import::limma()
  effect_values <- limma[, names(contrast), 'effect'] %>% round(digits = 1) %>% sprintf('%- .1f', .)
  fdr_values    <- limma[, names(contrast), 'fdr']    %>% sprintf('%.1e', .)
  p_values      <- limma[, names(contrast), 'p']      %>% sprintf('%.1e', .)
  autonomics.import::fdata(top)$plot.t <- sprintf('%s   %s   %s', effect_values, fdr_values, p_values)
  fvars %<>% c('plot.t')

  # Horizontal feature bars
  if (geom == 'hbar'){ 
     facet_def <- if('replicate' %in% autonomics.import::svars(object)) '~ subgroup + replicate'  else '~ sample'
     p <- top %>% autonomics.plot::plot_feature_hbars(alpha_var = 'plot.alpha', title = my_title, fvars = fvars, xlab = 'log2 ratio', color_var = color_var, ...)
  } else {
     p <- top %>% autonomics.plot::plot_features(     alpha_var = 'plot.alpha', title = my_title, fvars = fvars, geom = geom,         color_var = color_var, ...)
  }
  
  # Align title
  p <- p + ggplot2::theme(plot.title = ggplot2::element_text(family = 'mono', size = 11))
  p
}

#' Plot top features
#'
#' Plot top features per contrast
#' @param object          SummarizedExperiment
#' @param design          design matrix
#' @param contrasts       named contrast vector
#' @param result_dir      directory where to store results
#' @param geom           which feature plots to be created?
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
#'                 contrasts  = autonomics.import::contrastdefs(object)[1:2],
#'                 topdef     = 'fdr < 0.05',
#'                 x          = 'TIME_POINT',
#'                 geom       = 'boxplot',
#'                 result_dir = result_dir)
#'}
#' @importFrom magrittr  %>%
#' @export
plot_top_features_all_contrasts <- function(
   object,
   design         = autonomics.import::create_design_matrix(object),
   contrasts      = autonomics.find::default_contrasts(object),
   topdef         = autonomics.find::default_topdef(object),
   geom           = autonomics.plot::default_feature_plots(object),
   result_dir,
   ...
){
  for (i in seq_along(contrasts)){
    contrast <- contrasts[i] %>% magrittr::set_names(names(contrasts)[i])
    subdir   <- autonomics.find::get_contrast_subdir(result_dir, names(contrast))
    dir.create(subdir, recursive = TRUE, showWarnings = FALSE)
    for (direction in c('down', 'up')){
       for (i_plot in seq_along(geom)){
          cur_geom <- geom[[i_plot]]
          #if (length(x) > 1)   x %<>% magrittr::extract2(i_plot)
          my_file <- sprintf('%s/top_%ss__%s__%s.pdf', subdir, cur_geom, names(contrast), direction)
          comparator <- c(up='>', down='<')[[direction]]
          autonomics.support::cmessage('\t\t%s %s 0   %s', names(contrast), comparator, basename(my_file))
          object %>% autonomics.find::plot_top_features( design    = design,
                                                         contrast  = contrast,
                                                         topdef    = sprintf('effect %s 0  &  %s', comparator, topdef),
                                                         geom      = cur_geom,
                                                         file      = my_file, 
                                                         ...)
       }
    }
  }
}
