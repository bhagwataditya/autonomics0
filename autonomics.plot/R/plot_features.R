
#' Get no of facet rows/columns
#' @param b built ggplot (output of ggplot_build)
#' @return integer
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    object <- autonomics.data::stemcomp.proteinratios %>% magrittr::extract(1:9, )
#'    b <- object %>% autonomics.plot::plot_features('point') %>% ggplot2::ggplot_build()
#'    b %>% autonomics.plot::gg_nrow()
#'    b %>% autonomics.plot::gg_ncol()
#' }
#' @importFrom magrittr %>%
#' @export
gg_nrow <- function(b){
   assertive.sets::assert_are_set_equal(names(b), c('data', 'layout', 'plot'))
   b %>% magrittr::extract2('layout')       %>%
         magrittr::extract2('panel_layout') %>%
         magrittr::extract2('ROW')          %>%
         unique()                           %>%
         length()
}

#' @rdname gg_nrow
#' @importFrom magrittr %>%
#' @export
gg_ncol <- function(b){
   assertive.sets::assert_are_set_equal(names(b), c('data', 'layout', 'plot'))
   b %>% magrittr::extract2('layout')       %>%
         magrittr::extract2('panel_layout') %>%
         magrittr::extract2('COL')          %>%
         unique()                           %>%
         length()
}

#' Get x values
#'
#' Get x values from faceted ggplot object
#' @param b built ggplot (result of ggplot_build)
#' @return x values
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    object <- autonomics.data::stemcomp.proteinratios %>% magrittr::extract(1:9, )
#'    b <- object %>% autonomics.plot::plot_features('point') %>% ggplot2::ggplot_build()
#'    b %>% autonomics.plot::gg_xvalues()
#' }
#' @importFrom magrittr %>%
#' @export
gg_xvalues <- function(b){
   assertive.sets::assert_are_set_equal(names(b), c('data', 'layout', 'plot'))
   b %>% magrittr::extract2('data') %>%
         magrittr::extract2(which.max(vapply(., nrow, integer(1)))) %>%
         magrittr::extract2('x') %>%
         unique()
}

#' Get x labels
#'
#' Get x labels from ggplot object
#'
#' @param b built ggplot (result of ggplot_build)
#' @return x labels
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    object <- autonomics.data::stemcomp.proteinratios %>% magrittr::extract(1:9, )
#'    b <- object %>% autonomics.plot::plot_features('point') %>% ggplot2::ggplot_build()
#'    b %>% autonomics.plot::gg_xlabels()
#' }
#' @importFrom magrittr %>%
#' @export
gg_xlabels <- function(b){
   assertive.sets::assert_are_set_equal(names(b), c('data', 'layout', 'plot'))
   b %>% magrittr::extract2('layout') %>%
         magrittr::extract2('panel_ranges') %>%
         magrittr::extract2(1) %>%
         magrittr::extract2('x.labels')
         # ggplot2::ggplot_build(p) %>%
         # magrittr::extract2('data') %>%
         # magrittr::extract2(which.max(vapply(., nrow, integer(1)))) %>%
         # magrittr::extract2('x') %>% unique() %>%
         # length()
}

#' @title       Get column variables
#' @description Get column variables from faceted ggplot object
#' @param b built ggplot (result of ggplot_build)
#' @return vector with column variables
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    object <- autonomics.data::stemcomp.proteinratios %>% magrittr::extract(1:9, )
#'    b <- object %>% autonomics.plot::plot_features('point') %>% ggplot2::ggplot_build()
#'    b %>% autonomics.plot::gg_rowvars()
#'    b %>% autonomics.plot::gg_colvars()
#'    b %>% autonomics.plot::gg_nchar_rowvars()
#'    b %>% autonomics.plot::gg_panelvars()
#' }
#' if (require(atkin.2014)){
#'    object <- atkin.2014::soma %>% magrittr::extract(1:9, )
#'    b <- object %>%
#'         autonomics.plot::plot_features('point', x = 'time', facet_var = 'subject_id') %>%
#'         ggplot2::ggplot_build()
#'    b %>% autonomics.plot::gg_rowvars()
#'    b %>% autonomics.plot::gg_colvars()
#'    b %>% autonomics.plot::gg_nchar_rowvars()
#'    b %>% autonomics.plot::gg_panelvars()
#' }
#' @importFrom magrittr %>%
#' @export
gg_rowvars <- function(b){
   assertive.sets::assert_are_set_equal(names(b), c('data', 'layout', 'plot'))
   if ('rows' %in% names(b$plot$facet$params))   names(b$plot$facet$params$rows)
   else                                          return(character(0))
}

#' @rdname gg_rowvars
#' @importFrom magrittr %>%
#' @export
gg_colvars <- function(b){
   assertive.sets::assert_are_set_equal(names(b), c('data', 'layout', 'plot'))
   # Return 0 if absent
   if ('cols' %in% names(b$plot$facet$params))   names(b$plot$facet$params$cols)
   else                                          return(character(0))
}

#' @rdname gg_rowvars
#' @importFrom magrittr %>%
#' @export
gg_nchar_rowvars <- function(b){
   # Extract rowvars
   rowvars <- gg_rowvars(b)

   # Return nchar
   if (length(rowvars)==0)  return(0)
   b$plot$data %>% magrittr::extract(rowvars) %>%
                   apply(2, function(x) 2+ max(nchar(unique(x)))) %>% return()
}

#' @importFrom magrittr %>%
#' @export
#' @rdname gg_rowvars
gg_panelvars <- function(b){
   assertive.sets::assert_are_set_equal(names(b), c('data', 'layout', 'plot'))

   b %>% magrittr::extract2('layout') %>%
         magrittr::extract2('panel_layout') %>%
         names() %>%
         setdiff(c('PANEL', 'ROW', 'COL', 'SCALE_X', 'SCALE_Y')) %>%
         setdiff(gg_rowvars(b)) %>%
         setdiff(gg_colvars(b))
}

#' @importFrom magrittr %>%
gg_panelvar_width <- function(b){
   b$layout$panel_layout %>%
      magrittr::extract(, gg_panelvars(b), drop = FALSE) %>%
      unlist() %>% unname() %>% as.character() %>% nchar() %>%
      max()
}

#' Get name of variable mapped to a aesthetic
#' @param b built ggplot (result of ggplot_build)
#' @param aesthetic ggplot aesthetic (string)
#' @return character vector
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    object <- autonomics.data::stemcomp.proteinratios %>% magrittr::extract(1:9, )
#'    b <- object %>% autonomics.plot::plot_features('point') %>%
#'                    ggplot2::ggplot_build()
#'    b %>% autonomics.plot::gg_aesvar(aesthetic = 'color')
#'    b %>% autonomics.plot::gg_aesvar(aesthetic = 'x')
#'    b %>% autonomics.plot::gg_aesvar(aesthetic = 'shape')
#'    b <- object %>% autonomics.plot::plot_features('violin') %>%
#'                    ggplot2::ggplot_build()
#'    b %>% autonomics.plot::gg_aesvar(aesthetic = 'color')
#'    b %>% autonomics.plot::gg_aesvar(aesthetic = 'x')
#'    b %>% autonomics.plot::gg_aesvar(aesthetic = 'shape')
#' }
#' if (require(atkin.2014)){
#'    object <- atkin.2014::soma %>% magrittr::extract(1:10, )
#'    b <- object %>% autonomics.plot::plot_features('point',
#'                         x = 'time', color_var = 'condition', facet_var = 'subject_id',
#'                         line = TRUE, fvars = 'TargetFullName') %>%
#'                         ggplot2::ggplot_build()
#'    b %>% gg_aesvar(aesthetic = 'colour')
#' }
#' if (require(subramanian.2016)){
#'    b <- subramanian.2016::metabolon %>%
#'         magrittr::extract(1:10, ) %>%
#'         plot_features('boxplot',
#'                        color_var = 'condition',
#'                        group_var = 'condition',
#'                        line      = TRUE) %>%
#'         ggplot2::ggplot_build()
#'    aesthetic <- 'colour'
#'    gg_aesvar(b, aesthetic)
#' }
#' @importFrom magrittr %>%
#' @export
gg_aesvar <- function(b, aesthetic){
   layers <- b$plot$layers
   geom  <- layers %>% vapply(function(x)class(x$geom)[1], character(1))
   selected_layer <- which(geom %in% c('GeomBoxplot', 'GeomViolin', 'GeomPoint', 'GeomBar'))[1]
   selected_geom  <- geom[selected_layer]
   aesthetic %<>% stringi::stri_replace_first_fixed('color', 'colour')
   if (selected_geom %in% c('GeomBoxplot', 'GeomViolin', 'GeomBar')){
      aesthetic %<>% stringi::stri_replace_first_fixed('colour', 'fill')
   }
   b$plot$layers %>% magrittr::extract2(selected_layer) %>%
                     magrittr::extract2('mapping')      %>%
                     magrittr::extract2(aesthetic)      %>%
                     as.character()
}

#' Get levels of variable mapped to aesthetic
#' @param b built ggplot (result of ggplot_build)
#' @param aesthetic ggplot aesthetic (string)
#' @return character vector
#' @examples
#' require(magrittr)
#' require(Biobase)
#' if (require(autonomics.data)){
#'    object <- autonomics.data::stemcomp.proteinratios %>% magrittr::extract(1:9, )
#'    b <- object %>% autonomics.plot::plot_features('point') %>%
#'                    ggplot2::ggplot_build()
#'    b %>% autonomics.plot::gg_aeslevels(aesthetic = 'color')
#'    b %>% autonomics.plot::gg_aeslevels(aesthetic = 'shape')
#'    b %>% autonomics.plot::gg_aeslevels(aesthetic = 'x')
#' }
#' @importFrom magrittr %>%
#' @export
gg_aeslevels <- function(b, aesthetic){
   assertive.sets::assert_are_set_equal(names(b), c('data', 'layout', 'plot'))
   assertive.types::assert_is_a_string(aesthetic)
   aesthetic %<>% stringi::stri_replace_first_fixed('color', 'colour')
   aes_var <- autonomics.plot::gg_aesvar(b, aesthetic)
   b$plot %>% magrittr::extract2('data') %>%
              magrittr::extract2(aes_var) %>%
              unique() %>%
              as.character()
}

#' Compute appropriate feature plot dims
#' @param p ggplot object
#' @examples
#' require(magrittr)
#' require(SummarizedExperiment)
#' if (require(autonomics.data)){
#'    object <- autonomics.data::stemcomp.proteinratios %>% magrittr::extract(1:9, )
#'    object %>% plot_features('point')   %>% compute_feature_plot_dims()
#'    object %>% plot_features('violin')  %>% compute_feature_plot_dims()
#'    object %>% plot_features('boxplot') %>% compute_feature_plot_dims()
#' }
#' if (require(atkin.2014)){
#'    object <- atkin.2014::soma %>% magrittr::extract(1:10, )
#'    p <- object %>% plot_features('point', x = 'time', color_var = 'condition',
#'                         facet_var = 'subject_id', line = TRUE, fvars = 'TargetFullName')
#'    p %>% compute_feature_plot_dims()
#'
#' }
#' @export
compute_feature_plot_dims <- function(p){

   suppressWarnings(b <- p %>% ggplot2::ggplot_build())

   # n
   nrow <- b %>% gg_nrow()
   ncol <- b %>% gg_ncol()
   nx   <- b %>% gg_xlabels()        %>% length()

   rowvar_width    <- 0.01*(gg_nchar_rowvars(b) %>% magrittr::add(2) %>% sum())
   colvar_height   <- 0.1*length(gg_colvars(b))
   panelvar_height <- 0.1*(length(gg_panelvars(b)))
   panelvar_width  <- 0.055*gg_panelvar_width(b)

   # legends
   aesthetics <- b$plot$labels %>% (function(x)x[!is.null(x) & !x=='NULL']) %>%  names() %>% intersect(c('colour', 'fill', 'shape'))
   variables  <- aesthetics %>% lapply(function(cur_aes) gg_aesvar(   b, cur_aes)) %>% unlist()
   levels     <- aesthetics %>% lapply(function(cur_aes) gg_aeslevels(b, cur_aes)) %>% unlist()
   label_width  <- 0.4*max(c(nchar(variables), 3+nchar(levels))) # approx 3 characters for icon
   label_height <- 0.3*max(nchar(gg_aeslevels(b, 'x'))) #c(variables, levels) %>% length()          %>% magrittr::multiply_by(0.5)

   # delta
   dx <- 0.5
   aspect_ratio <- 1.2

   # width & height
   panel_width  <- max(dx*nx, panelvar_width)*aspect_ratio
   panel_height <- max(dx*nx, panelvar_width) + panelvar_height
   fig_width  <- panel_width * ncol + label_width  + rowvar_width
   fig_height <- panel_height* nrow + label_height + colvar_height  #max(dd*nx*nrow, vspace)

   # Resize very big figures to keep text readable
   if (fig_width > 20){
      resize_factor <- 20/fig_width
      fig_width %<>% magrittr::multiply_by(resize_factor)
      fig_height %<>% magrittr::multiply_by(resize_factor)
   }

   return(list(width = fig_width, height = fig_height))
   # print
   #suppressMessages(suppressWarnings(p %>% autonomics.support::print2pdf(file, width = fig_width, height = fig_height)))
   #return(invisible(file))
}

old_print_feature_boxes <- function(...){
   old_print_feature_distributions(...)
}

#' @importFrom magrittr %>%
old_print_feature_distributions <- function(p, file, ngroups, nrows, nfvars, nsamples){
   fig_width  <- 3 + 1.5*nrows*((ngroups^0.5)/1.5) # the power ensures saturation at high ngroups
   fig_height <- 3 + nrows*(1.5 + 0.10**nfvars)
   suppressMessages(suppressWarnings(p %>% autonomics.support::print2pdf(file, width = fig_width, height = fig_height)))
   return(invisible(file))
}

#' @importFrom magrittr %>%
old_print_feature_profiles <- function(p, file, ngroups, nrows, nfvars, nsamples){
  fig_width  <- 3 + 1.5*nrows + 3
  fig_height <- 3 + nrows*(1.5 + 0.10*nfvars)
  suppressMessages(suppressWarnings(p %>% autonomics.support::print2pdf(file, width = fig_width, height = fig_height)))
  return(invisible(file))
}

print_feature_bars <- function(p, file, ngroups, nrows, nfvars, nsamples){
   fig_width  <- 3 + 1*nsamples + 3
   fig_height <- 3 + nrows*(1.5 + 0.10*nfvars)
   suppressMessages(suppressWarnings(p %>% autonomics.support::print2pdf(file, width = fig_width, height = fig_height)))
   return(invisible(file))
}


add_annotation <- function(p, title, xlab, ylab, x_text_angle){
  p +
  ggplot2::ggtitle(title) +
  ggplot2::xlab(xlab) +
  ggplot2::ylab(ylab) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = x_text_angle, vjust = 1, hjust = 1))
}

# make_facet_wrap_formula <- function(fvars){
#   facet_def <- sprintf('~ `%s`', fvars[[1]])
#   if (length(fvars) > 1){
#     for (i in seq(2, length(fvars))){
#       facet_def <- sprintf('%s + `%s`', facet_def, fvars[[i]])
#     }
#   }
#   stats::as.formula(facet_def)
# }

make_facet_grid_formula <- function(fvars, facet_var){
   paste0(fvars, collapse = ' + ') %>%
   paste0(' ~ ', facet_var) %>%
   stats::as.formula()
}

# wrap on "feature" to have order correct
# use this labeler to annotate nicely
feature_plot_labeller <- function(plot_df){
   plot_df %>% magrittr::extract2('feature_facet')  %>%
               as.character()                       %>%
               strsplit(' :: ')                     %>%
               do.call(rbind, .)                    %>%
               data.frame(stringsAsFactors = FALSE)
}


#' Plot features
#'
#' @param object           SummarizedExperiment, eSet, or EList
#' @param geom             'point', 'bar', 'violin', 'boxplot'
#' @param x                svar mapped to x
#' @param color_var        svar mapped to color
#' @param color_values     named color vector (names = color_var levels, values = colors)
#' @param shape_var        svar mapped to shape (only relevant for point)
#' @param group_var        svar mapped to group
#' @param txt_var          svar mapped to txt
#' @param facet_var        svar on which to facet plot
#' @param alpha_var        svar (logical) mapped to absence/presence of transparancy
#' @param fvars            fvar used for plot annotation
#' @param title            plot title
#' @param scales           'free'
#' @param x_text_angle     angle of text on x axis
#' @param line             logical: should line be added?
#' @param zero_hline       logical: should y=0 line be added?
#' @param xlab             xlab
#' @param ylab             ylab
#' @param dodge_width      numeric
#' @param verbose          logical
#' @param file             string: file path
#' @param width            numeric
#' @param height           numeric
#' @param theme            ggplot2::theme
#' @param legend.position  position of legend
#' @param ...              only for backward compatibility to deprecated functions
#' @examples
#' require(magrittr)
#' result_dir <- tempdir() %T>% message()
#'
#' # STEM CELL COMPARISON
#' if (require(autonomics.data)){
#'    object <- autonomics.data::stemcomp.proteinratios %>% extract(1:10, )
#'    object %>% autonomics.plot::plot_features(geom = 'violin')
#'    object %>% autonomics.plot::plot_features(geom = 'point')
#'    object %>% autonomics.plot::plot_features(geom = 'boxplot')
#'    object %>% autonomics.plot::plot_features(
#'                  geom = 'boxplot',
#'                  file = paste0(result_dir, '/stemcomp_boxes.pdf'))
#'
#' }
#'
#' # STEM CELL DIFFERENTIATION
#' if (require(autonomics.data)){
#'    object <- autonomics.data::stemdiff.proteinratios %>% extract(1:10, )
#'    object %>% autonomics.plot::plot_features('boxplot')
#' }
#'
#' if (require(billing.differentiation.data)){
#'    object <- billing.differentiation.data::rna.voomcounts %>% extract(1:10, )
#'    object %>% autonomics.plot::plot_features('point',   fvars = 'gene_name')
#'    object %>% autonomics.plot::plot_features('violin',  fvars = 'gene_name')
#'    object %>% autonomics.plot::plot_features('boxplot', fvars = 'gene_name')
#' }
#'
#' # GLUTAMINASE
#' if (require(autonomics.data)){
#'    object <- autonomics.data::glutaminase[1:9, ]
#'    color_values <- c(Control = 'blue', Vehicle = 'green', `Concentration 1` = 'orange',
#'                     `Concentration 2` = 'red')
#'    object$alpha <- object$TREATMENT %in% c(0,1)
#'    object %>% autonomics.plot::plot_features(
#'                 'boxplot',
#'                  x = 'TIME_POINT', color_var = 'GROUP_DESCRIPTION', color_values = color_values,
#'                  alpha_var = 'alpha')
#' }
#'
#' # SUBRAMANIAN.2016
#' if (require(subramanian.2016)){
#'    object <- subramanian.2016::metabolon[1:10, ]
#'    object %>% plot_features('boxplot', color_var = 'condition',
#'                              group_var = 'condition', line = TRUE)
#'
#'    object <- subramanian.2016::rnaseq[1:10, ]
#'    object %>% plot_features('boxplot', color_var = 'condition', group_var = 'condition',
#'                              line = TRUE)
#'
#' }
#'
#' # ATKIN.2014
#' if (require(atkin.2014)){
#'    object <- atkin.2014::soma %>% magrittr::extract(1:10, )
#'    object %>% plot_features('point', x = 'time', color_var = 'condition',
#'                              facet_var = 'subject_id', group_var = 'subject_id', line = TRUE)
#'    object %>% plot_features('point', x = 'time', color_var = 'condition', facet_var = 'subject_id',
#'                              group_var = 'subject_id', line = TRUE, fvar = 'TargetFullName')
#'    object %>% plot_features('point', x = 'time', color_var = 'condition', facet_var = 'subject_id',
#'                              group_var = 'subject_id', line = TRUE, fvar = 'TargetFullName')
#' }
#'
#' @return file path
#' @importFrom magrittr  %>%   %<>%
#' @export
plot_features <- function(
   object,
   geom            = autonomics.plot::default_feature_plots(object)[1],
   x               = autonomics.plot::default_x(object, geom[1]),
   color_var       = autonomics.plot::default_color_var(object),
   color_values    = autonomics.plot::default_color_values(object, color_var),
   shape_var       = autonomics.plot::default_shape_var(object),
   group_var       = autonomics.plot::default_group_var(object),
   txt_var         = autonomics.plot::default_txt_var(object),
   facet_var       = NULL,
   alpha_var       = NULL,
   fvars           = autonomics.plot::default_fvars(object),
   title           = '',
   scales          = 'free_y',
   x_text_angle    = 90, # is better readable than e.g. 60 with many xvalues
   line            = autonomics.plot::default_line(object),
   zero_hline      = autonomics.plot::default_zero_hline(object),
   xlab            = NULL,
   ylab            = NULL,
   dodge_width     = 0,
   verbose         = FALSE,
   file            = NULL,
   width           = NULL,
   height          = NULL,
   theme           = ggplot2::theme_bw(),
   legend.position = 'right'
){
   # Set theme
   old_theme <- ggplot2::theme_get()
   ggplot2::theme_set(theme)

   # Function to plot features
   plot_features_without_printing <- function(object){
      # Defaults
      if (x == 'snames'){        autonomics.import::sdata(object)$sample_id <- autonomics.import::snames(object)
      x <- 'sample_id'}
      if (is.null(fvars)){       autonomics.import::fdata(object)$feature_id <- autonomics.import::fnames(object)
      fvars <- 'feature_id'}

      # Return if no features to plot
      if (nrow(object) == 0){
         message('\t\tno features to plot, no graph created')
         return()
      }

      # Assert
      autonomics.import::assert_is_valid_eset(object)
      assertive.sets::assert_is_subset(geom, autonomics.plot::FEATURE_PLOTS)
      assertive.sets::assert_is_subset(x,         autonomics.import::svars(object))
      assertive.sets::assert_is_subset(color_var, autonomics.import::svars(object))
      if (!is.null(shape_var)){
         assertive.sets::assert_is_subset(shape_var, autonomics.import::svars(object) )
         if (is.numeric(autonomics.import::sdata(object)[[shape_var]])){
            autonomics.import::sdata(object)[[shape_var]] %<>% as.character() %>% autonomics.support::factorify()
         }
      }
      if (!is.null(fvars))     assertive.sets::assert_is_subset(fvars,     autonomics.import::fvars(object))
      if (!is.null(alpha_var)){
         assertive.sets::assert_is_subset(alpha_var, autonomics.import::svars(object))
         assertive.types::assert_is_logical(autonomics.import::sdata(object)[[alpha_var]])
      }

      # Prepare plot df
      plot_df <- autonomics.plot::create_feature_plot_df(object, fvars, verbose = verbose)
      p <- ggplot2::ggplot(plot_df)

      # Add annotation and facet wrap
      if (is.null(facet_var)){ p <- p + ggplot2::facet_wrap(~ feature_facet, scales = scales, labeller = feature_plot_labeller)
      } else {                 p <- p + ggplot2::facet_grid(make_facet_grid_formula(fvars, facet_var), scales = 'free_y', switch = 'y') +
                               ggplot2::theme(strip.text.y = ggplot2::element_text(angle=180))}
      p %<>% add_annotation(title = title, xlab = xlab, ylab = ylab, x_text_angle = x_text_angle)

      # Add zero hline
      if (zero_hline){
         p <- p + ggplot2::geom_hline(yintercept = 0, linetype = 'dashed')
      }

      # Plot either point or distributions
      dodged_position <- ggplot2::position_dodge(width = dodge_width)
      if (       geom == 'point'){   p <- p + ggplot2::geom_point(  ggplot2::aes_string(x = x, y = 'value', color = color_var, alpha = alpha_var, shape = shape_var, group = group_var), na.rm = TRUE, position = dodged_position) # Note: group is required for dodging!
      } else if (geom == 'bar'){     p <- p + ggplot2::geom_bar(    ggplot2::aes_string(x = x, y = 'value', fill  = color_var, alpha = alpha_var), stat = 'identity', na.rm = TRUE)
      } else if (geom == 'violin'){  p <- p + ggplot2::geom_violin( ggplot2::aes_string(x = x, y = 'value', fill  = color_var, alpha = alpha_var),                    na.rm = TRUE)
      p <- p + ggplot2::geom_boxplot(ggplot2::aes_string(x = x, y = 'value', fill  = color_var, alpha = alpha_var), data = plot_df, width = 0.1)
      } else if (geom == 'boxplot'){ p <- p + ggplot2::geom_boxplot(ggplot2::aes_string(x = x, y = 'value', fill  = color_var, alpha = alpha_var), data = plot_df)
      }

      # Take care of alpha transparancy
      if (!is.null(alpha_var)){
         if (all(plot_df[[alpha_var]] == TRUE)){
            p <- p + scale_alpha_discrete(range = c(1,1), guide = FALSE)        # If all subgroups are selected, none should be faded out!
         } else {
            transparent <- if (geom=='boxplot') 0.2 else 0.5
            p <- p + scale_alpha_discrete(range = c(transparent, 1), guide = FALSE)
         }
      }

      # Custom color scale
      if (!is.null(color_values)){
         p <- p + ggplot2::scale_fill_manual(values = color_values)
         p <- p + ggplot2::scale_color_manual(values = color_values)
      }

      # Connect with a line
      if (line){
         # single line with fixed color
         if (group_var == 1){   p <- p + ggplot2::stat_summary(fun.y = 'median', geom = 'line', na.rm = TRUE,  color = 'gray10',
                                                               ggplot2::aes_string(x = x, y = 'value', group = group_var))
         # separate lines with different color
         } else {               p <- p + ggplot2::stat_summary(fun.y = 'median', geom = 'line', na.rm = TRUE,
                                                               ggplot2::aes_string(x = x, y = 'value', group = group_var, color = color_var), position = dodged_position)
         }
      }

      # Add text annotation
      if (!is.null(txt_var)){
         p <- p + ggrepel::geom_text_repel(ggplot2::aes_string(label = txt_var, x=x, y='value', color = color_var),
                                           show.legend = FALSE)
      }

      # Legend position
      p <- p + ggplot2::theme(legend.position = legend.position)
      p

   }


   # Plot and return if no printing required
   if (is.null(file)){
      p <- plot_features_without_printing(object)
      return(p)
   }

   # Print to file
   nperpage <- 9
   eset1 <- object %>% magrittr::extract(seq(1, min(nrow(.), nperpage)), )
   fig_dims <- if (is.null(width) & is.null(height)){
                  compute_feature_plot_dims(plot_features_without_printing(object = eset1))
               } else {
                  if (is.null(width))   width <- height
                  if (is.null(height)) height <- width
                  list(width = width, height = height)
               }
   grDevices::pdf(file, width = fig_dims$width, height = fig_dims$height)
   npage <- ceiling(nrow(object) / nperpage)
   for (ipage in 1:npage){
      first <- (ipage-1)*nperpage + 1
      last <- min(nrow(object), ipage*nperpage)
      eset1 <- object %>% magrittr::extract(first:last, )
      p <- plot_features_without_printing(object = eset1)
      print(p)
   }
   grDevices::dev.off()

   # Return
   ggplot2::theme_set(old_theme)
   return(invisible(file))
}


#' @rdname plot_features
#' @export
plot_feature_bars <- function(...){
   .Deprecated('plot_features')
   plot_features(geom = 'bar', ...)
}


#' @rdname plot_features
#' @export
plot_feature_profiles <- function(...){
   .Deprecated('plot_features')
   plot_features(geom = 'point', ...)
}


#' @rdname plot_features
#' @export
plot_feature_distributions <- function(...){
   .Deprecated('plot_features')
   plot_features(geom = 'violin', ...)
}

#' @rdname plot_features
#' @export
plot_feature_boxes <- function(...){
   .Deprecated('plot_features')
   plot_features(geom = 'boxplot', ...)
}
