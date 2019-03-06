#' Make datatable to plot PCA/PLS/LDA/SMA scores
#' @param projection  output from PCA/PLS/LDA/SMA
#' @param object      SummarizedExperiment
#' @param dims        PCA/LDA dimensions
#' @param color_var   svar mapped to color
#' @param shape_var   svar mapped to shape
#' @param size_var    svar mapped to size
#' @param txt_var     svar mapped to txt
#' @param split_var   svar on which object is splitted prior to lda analysis
#' @param facet_var   svar on which plot is faceted
#' @param group_var   svar mapped to group_var
#' @examples
#' if (require(autonomics.data)){
#'    require(magrittr)
#'    object <- autonomics.data::glutaminase
#'
#'    object %>% pca() %>%
#'               make_projected_samples_df(object,
#'                                         dims = c(1,2),
#'                                         color_var = 'subgroup',
#'                                         shape_var = NULL,
#'                                         size_var  = NULL,
#'                                         txt_var   = NULL,
#'                                         split_var = NULL,
#'                                         facet_var = NULL,
#'                                         group_var = NULL) %>%
#'               print()
#'
#'    object %>% pca(ndim=4) %>%
#'               make_projected_samples_df(object,
#'                                         dims      = 3:4,
#'                                         color_var = 'subgroup',
#'                                         shape_var = NULL,
#'                                         size_var  = NULL,
#'                                         txt_var   = NULL,
#'                                         split_var = NULL,
#'                                         facet_var = NULL,
#'                                         group_var = NULL) %>%
#'               print()
#' }
#' @importFrom data.table  data.table  :=
#' @importFrom magrittr    %<>%
#' @export
make_projected_samples_df <- function(
   projection,
   object,
   dims,
   color_var,
   shape_var,
   size_var,
   txt_var,
   split_var,
   facet_var,
   group_var
){
   # Satisfy CHECK
   color <- shape <- size <- label <- group <- facet <- NULL

   # Add MPA coordinates
   DT <- data.table::data.table(
            sample_id = autonomics.import::snames(object),
            x         = projection$samples[, dims[1]],
            y         = projection$samples[, dims[2]])

   # Add other aesthetics
   DT %<>% magrittr::extract(, color := if(is.null(color_var)) 'default'  else  autonomics.import::sdata(object)[[color_var]] )
   DT %<>% magrittr::extract(, shape := if(is.null(shape_var)) 'default'  else  autonomics.import::sdata(object)[[shape_var]] )
   DT %<>% magrittr::extract(, size  := if(is.null(size_var))  'default'  else  autonomics.import::sdata(object)[[size_var ]] )
   if (!is.null(txt_var))   DT %<>% magrittr::extract(, label := autonomics.import::sdata(object)[[txt_var]])
   if (!is.null(group_var)) DT %<>% magrittr::extract(, group := autonomics.import::sdata(object)[[group_var]])
   if (!is.null(facet_var)) DT %<>% magrittr::extract(, facet := autonomics.import::sdata(object)[[facet_var]])
   if (!is.null(split_var)) DT %<>% magrittr::extract(, facet := sprintf('%s %d%% %d%%', facet, round(projection$var[1]), round(projection$var[2])))
   DT
}




#' Plot PCA/LDA/PLS/SMA sample scores
#' @param object            SummarizedExperiment
#' @param ...               Additional SummarizedExperiment(s)
#' @param listed_objects    Additional SummarizedExperiment(s), wrapped in a \code{link{list}}
#' @param method            'pca', 'lda', 'pls', 'sma'
#' @param implementation    which implementation of pls to use (see \code{\link{project}})
#' @param dims              dimensions
#' @param color_var         svar mapped to color
#' @param color_values      color values vector
#' @param shape_var         svar mapped to shape
#' @param size_var          svar mapped to size
#' @param txt_var           svar mapped to txt
#' @param txt_outliers      logical(1): whether to annotate outliers
#' @param group_var         svar mapped to group
#' @param split_var         svar on which to split data prior to transformation
#' @param facet_var         svar on which to facet sample biplot
#' @param scales            'free' or 'fixed'
#' @param labs              named \code{\link{list}} handed on to \code{\link[ggplot2]{labs}}
#' @param nrow              integer
#' @param mention_method    Whether to use \code{method} in title and/or facet labels
#' @param title             title
#' @param na.impute         TRUE or FALSE
#' @param legend.position   character
#' @examples
#' if (require(autonomics.data)){
#'    require(magrittr)
#'
#'    # GLUTAMINASE
#'       object <- autonomics.data::glutaminase
#'       object %>% plot_pca_samples()
#'       object %>% plot_lda_samples()
#'       object %>% plot_pls_samples()
#'       object %>% plot_projected_samples(
#'                     method = c('pca', 'lda', 'sma', 'pls'),
#'                     facet_var = c('PCA', 'LDA', 'SMA', 'PLS'),
#'                     nrow = 2)
#'
#'    # STEM CELL COMPARISON
#'       object <- autonomics.data::stemcomp.proteinratios
#'       object %>% plot_pca_samples()
#'       object %>% plot_pca_samples(dims = c(3,4))
#'
#'       object %>% plot_projected_samples(
#'                     method = c('pca', 'lda', 'sma', 'pls'),
#'                     facet_var = c('OK PCA', 'Not so nice LDA', 'Nice SMA', 'Very nice PLS'),
#'                     nrow = 2)
#'
#' }
#'
#' @author Aditya Bhagwat, Johannes Graumann
#' @importFrom data.table   data.table   :=
#' @importFrom magrittr %>% %<>%
#' @export
plot_projected_samples <- function(
   object,
   ...,
   listed_objects    = NULL,
   method            = c('pca', 'lda', 'sma', 'pls')[1],
   implementation    = c('mixOmics::plsda', 'mixOmics::splsda', 'ropls::opls')[1],
   dims              = 1:2,
   color_var         = 'subgroup', # default_color_var(object) gives 'block' when present!
   color_values      = default_color_values(object, color_var),
   shape_var         = default_shape_var(object),
   size_var          = NULL,
   txt_var           = default_txt_var(object),
   txt_outliers      = FALSE,
   group_var         = NULL,
   split_var         = NULL,
   facet_var         = NULL,
   scales            = ifelse(is.null(facet_var), 'fixed' ,'free'),
   labs              = list(color = stringi::stri_trans_totitle(color_var), shape = stringi::stri_trans_totitle(shape_var)),
   nrow              = NULL,
   mention_method    = ifelse(is.null(facet_var), TRUE, FALSE),
   title             = NULL,
   na.impute         = FALSE,
   legend.position   = 'right'
){

# Check prerequisites -----------------------------------------------------
   # Objects
   ## Combine all handed in objects
   if(!is.null(listed_objects)) listed_objects %>% assertive.types::assert_is_list()
   obj_list <- list(object) %>% c(list(...), listed_objects)
   obj_list %<>% magrittr::extract(!sapply(obj_list, is.null))

   ## Assert all
   obj_list %>%  lapply(autonomics.import::assert_is_valid_object) %>%
                 lapply(function(x) x %>% autonomics.import::exprs() %>% assertive.properties::assert_is_non_empty())
   obj_list %<>% lapply(function(x) x %>% validify_shape_values(shape_var))
   dims %>% assertive.types::assert_is_numeric()              %>%
            assertive.properties::assert_is_of_length(2)      %>%
            assertive.numbers::assert_all_are_whole_numbers() %>%
            assertive.numbers::assert_all_are_greater_than(0)

   for(var in c(color_var, group_var, shape_var, size_var, split_var, txt_var)){
      obj_list %>% sapply(function(x) x %>% autonomics.import::svars() %>% assertive.sets::assert_is_superset(var))
   }

   if(length(facet_var) == 1){
      obj_list %>% sapply(function(x) x %>% autonomics.import::svars() %>% assertive.sets::assert_is_superset(facet_var))

   } else if(length(facet_var) > 1){
      if(length(obj_list) == 1){ obj_list <- obj_list[[1]] %>% list() %>% rep(times = length(facet_var))
      } else {                   facet_var %>% assertive.properties::assert_are_same_length(obj_list)
      }
   }

   # Must be left AFTER 'facet_var' - as obj_list may change depending on that
   method %<>% match.arg(choices = c('pca', 'lda', 'sma', 'pls'), several.ok = TRUE)
   if(length(method) > 1) method %>% assertive.properties::assert_are_same_length(obj_list)

   # color_values # how to check this?

   scales %<>% match.arg(choices = c('fixed', 'free_x', 'free_y', 'free'), several.ok = FALSE)

   if(!is.null(labs)){
      labs %>% assertive.types::assert_is_list() %>%
               assertive.properties::assert_has_names() # Leaky: 'all_have_names' does not exist
   }

   if(!is.null(nrow)){
      nrow %>% assertive.types::assert_is_a_number() %>%
               assertive.numbers::assert_all_are_whole_numbers() %>%
               assertive.numbers::assert_all_are_greater_than(0)
   }

   mention_method %>% assertive.types::assert_is_a_bool()

   if(!is.null(title)) title %>% assertive.types::assert_is_a_string()

   na.impute %>% assertive.types::assert_is_a_bool()

   if(length(legend.position) == 1){
     legend.position %<>% match.arg(choices    = c("none", "left", "right", "bottom", "top"), several.ok = FALSE)
   } else if(length(legend.position) == 2){
      legend.position %>% assertive.types::assert_is_numeric() %>%
                          assertive.numbers::assert_all_are_in_closed_range(c(0,1))
   } else {
      stop('\'legend.position\' must be of length [1,2].')
   }

   # Transform
   project_list <- seq_along(obj_list) %>%
                   lapply(function(x){ obj_list[[x]] %>%
                                       project(method         = ifelse(length(method) == 1, method, method[x]),
                                               implementation = implementation,
                                               ndim           = max(dims),
                                               na.impute      = na.impute)})

# Processing --------------------------------------------------------------
   # Are we dealing with an on-the-fly generated facetting variable or do we use one from the data?
   ## Augment the facet vars
   if(length(facet_var) > 1){
      facet_var <- seq_along(facet_var) %>%
                   sapply(function(x){ make_sample_scores_title2(
                                          object                = obj_list[[x]],
                                          project_result        = project_list[[x]],
                                          dims                  = dims,
                                          method                = ifelse(length(method) == 1, method, method[x]),
                                          manual_label          = facet_var[x],
                                          mention_method        = mention_method,
                                          mention_feature_count = TRUE,
                                          mention_variance      = TRUE,
                                          separator             = '\n')})
   }

   # Generate data.tables/frames
   plotDF <- seq_along(obj_list) %>%
             lapply(function(x){
                      if(length(facet_var) <= 1){
                         make_projected_samples_df(object     = obj_list[[x]],
                                                   projection = project_list[[x]],
                                                   dims       = dims,
                                                   color_var  = color_var,
                                                   shape_var  = shape_var,
                                                   size_var   = size_var,
                                                   txt_var    = txt_var,
                                                   facet_var  = facet_var,
                                                   split_var  = split_var,
                                                   group_var  = group_var)
                      } else {
                         tmp_df <- make_projected_samples_df(object = obj_list[[x]],
                                                             projection = project_list[[x]],
                                                             dims       = dims,
                                                             color_var  = color_var,
                                                             shape_var  = shape_var,
                                                             size_var   = size_var,
                                                             txt_var    = txt_var,
                                                             facet_var  = NULL,
                                                             split_var  = split_var,
                                                             group_var  = group_var)
                         tmp_df[['facet']] <- facet_var[x]
                         tmp_df  %<>% dplyr::mutate_(facet = ~facet %>% factor(levels = facet_var))
                         tmp_df
                      }}) %>%
             data.table::rbindlist()

   # Initialize plot
   p <- ggplot2::ggplot(plotDF) + ggplot2::theme_bw()

   # Points
   p <- p + ggplot2::geom_point(
            ggplot2::aes_string(
               x     = 'x',
               y     = 'y',
               color = 'color',
               shape = 'shape',
               size  = 'size'))
   if (is.null(shape_var))  p <- p + ggplot2::scale_shape_manual(values = c(default = 15), guide = FALSE)
   if (is.null(size_var))   p <- p + ggplot2::scale_size_manual( values = c(default = 3),  guide = FALSE)
   plot_guide <- if (is.null(color_var)) FALSE  else TRUE
   p <- p + ggplot2::scale_color_manual(values = color_values, guide = plot_guide)
   if (!is.null(color_var)) p <- p + ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(shape = 15, size = 2)))

   # Line
   if (!is.null(group_var)){
      p <- p + ggplot2::geom_path(ggplot2::aes_string(x = 'x', 'y', color = 'color', group = 'group'))
   }

   # txt
   if (!is.null(txt_var)){
      p <- p + ggrepel::geom_text_repel(ggplot2::aes_string(x='x',y='y', label='label', color = 'color'), show.legend = FALSE)
   }
   if (txt_outliers){
      is_an_outlier <- autonomics.support::is_outlier(plotDF$x) | autonomics.support::is_outlier(plotDF$y)
      p <- p + ggrepel::geom_text_repel(data    = plotDF[is_an_outlier, ],
                                        mapping = ggplot2::aes_string(x='x', y='y', color = 'color', label = 'sample_id'))
   }

   # Facet
   if (!is.null(facet_var)){
      p <- p + ggplot2::facet_wrap(
               stats::as.formula('~ facet'),
               scales = scales,
               nrow   = nrow)
      p <- p + ggplot2::theme(strip.text = ggplot2::element_text(hjust = 0.05))
   }

   # Axes
   p <- p + ggplot2::geom_vline(xintercept = 0, linetype = 'dashed') +
            ggplot2::geom_hline(yintercept = 0, linetype = 'dashed')

   # Axis labels
   xlab <- paste0('X', dims[1])
   ylab <- paste0('X', dims[2])
   if(length(facet_var) <= 1){
      if (is.null(split_var)){ var <- project_list[[1]]$var
                               xlab %<>% paste0(' (', var[dims[1]], '%)')
                               ylab %<>% paste0(' (', var[dims[2]], '%)')}
   }
   p <- p + ggplot2::xlab(xlab) + ggplot2::ylab(ylab)

   # Title
   if(length(facet_var) <= 1){
      default_title <- make_sample_scores_title2(
                          object                = obj_list[[1]],
                          method                = method,
                          dims                  = dims,
                          project_result        = project_list[[1]],
                          mention_method        = mention_method,
                          mention_feature_count = TRUE,
                          mention_variance      = FALSE,
                          separator             = '; ')
      title <- if (is.null(title)) default_title   else   sprintf('%s: %s', title, default_title)
      p <- p + ggplot2::ggtitle(title)
   }
   p <- p + ggplot2::theme(legend.position = legend.position)

   # Legends/labs
   if(!is.null(labs)){
      p <- p + do.call(what = ggplot2::labs, args = labs)
   }

   # Left align
   if(length(facet_var) > 1){
      gp <- ggplot2::ggplotGrob(p)

      # Edit the grob
      # Is throwing a warning - fix later - outcomment now
      # gp <- grid::editGrob(grid::grid.force(gp),
      #                      grid::gPath("GRID.stripGrob", "GRID.text"),
      #                      grep = TRUE,
      #                      global = TRUE,
      #                      just = "left",
      #                      x = grid::unit(0, "npc"))
      grid::grid.draw(gp)
      invisible(gp)
   } else {
      p
   }

}


make_sample_scores_title2 <- function(
   object,
   project_result,
   method,
   dims,
   manual_label          = NULL,
   mention_method        = TRUE,
   mention_feature_count = TRUE,
   mention_variance      = FALSE,
   separator             = ' on '
){
   # Manual label
   string1 <- ifelse(is.null(manual_label), NA_character_, manual_label)
   # Method used
   string2 <- ifelse(mention_method, toupper(method), NA_character_)
   # Dimensions used
   if(mention_feature_count){
      selector <- if (method %in% c('pca', 'pls')){
         selector <- matrixStats::rowAnys(!is.na(autonomics.import::exprs(object))) &
            matrixStats::rowAnys(!is.infinite(autonomics.import::exprs(object)))
      } else {
         selector <- matrixStats::rowAlls(!is.na(autonomics.import::exprs(object))) &
            matrixStats::rowAlls(!is.infinite(autonomics.import::exprs(object)))
      }
      string3 <- paste0(sum(selector), ' / ', length(selector), ' features')
   } else {
      string3 <- NA_character_
   }
   # Variances explained
   if(mention_variance){
      var <- project_result[['var']][dims]
      string4 <- 'Variance expl. ' %>%
                 paste0( var     %>%
                         names() %>%
                         stringi::stri_replace_all_regex('.*(\\d+)$', 'X$1') %>%
                         paste0(': ', var, '%', collapse = ', '))
   } else {
      string4 <- NA_character_
   }
   # Assemble
   all_strings <- c(string1, string2, string3, string4) %>%
                  stats::na.omit()
   paste0(all_strings, collapse = separator)
}


#' @rdname plot_projected_samples
#' @export
plot_pca_samples <- function(
   object,
   ...,
   listed_objects    = NULL,
   dims              = 1:2,
   color_var         = 'subgroup', # default_color_var(object) gives 'block' when present!
   color_values      = default_color_values(object, color_var),
   shape_var         = default_shape_var(object),
   size_var          = NULL,
   txt_var           = default_txt_var(object),
   txt_outliers      = FALSE,
   group_var         = NULL,
   split_var         = NULL,
   facet_var         = NULL,
   scales            = ifelse(is.null(facet_var), 'fixed' ,'free'),
   labs              = list(color = stringi::stri_trans_totitle(color_var), shape = stringi::stri_trans_totitle(shape_var)),
   nrow              = NULL,
   mention_method    = ifelse(is.null(facet_var), TRUE, FALSE),
   title             = NULL,
   na.impute         = FALSE,
   legend.position   = 'right'
){
   plot_projected_samples( object,
                           ...,
                           listed_objects    = listed_objects,
                           method            = 'pca',
                           implementation    = character(0),
                           dims              = dims,
                           color_var         = color_var,
                           color_values      = color_values,
                           shape_var         = shape_var,
                           size_var          = size_var,
                           txt_var           = txt_var,
                           txt_outliers      = txt_outliers,
                           group_var         = group_var,
                           split_var         = split_var,
                           facet_var         = facet_var,
                           scales            = scales,
                           labs              = labs,
                           nrow              = nrow,
                           mention_method    = mention_method,
                           title             = title,
                           na.impute         = na.impute,
                           legend.position   = legend.position)
}


#' @rdname plot_projected_samples
#' @export
plot_sma_samples <- function(
   object,
   ...,
   listed_objects    = NULL,
   dims              = 1:2,
   color_var         = 'subgroup', # default_color_var(object) gives 'block' when present!
   color_values      = default_color_values(object, color_var),
   shape_var         = default_shape_var(object),
   size_var          = NULL,
   txt_var           = default_txt_var(object),
   txt_outliers      = FALSE,
   group_var         = NULL,
   split_var         = NULL,
   facet_var         = NULL,
   scales            = ifelse(is.null(facet_var), 'fixed' ,'free'),
   labs              = list(color = stringi::stri_trans_totitle(color_var), shape = stringi::stri_trans_totitle(shape_var)),
   nrow              = NULL,
   mention_method    = ifelse(is.null(facet_var), TRUE, FALSE),
   title             = NULL,
   na.impute         = FALSE,
   legend.position   = 'right'
){
   plot_projected_samples( object,
                           ...,
                           listed_objects    = listed_objects,
                           method            = 'sma',
                           implementation    = character(0),
                           dims              = dims,
                           color_var         = color_var,
                           color_values      = color_values,
                           shape_var         = shape_var,
                           size_var          = size_var,
                           txt_var           = txt_var,
                           txt_outliers      = txt_outliers,
                           group_var         = group_var,
                           split_var         = split_var,
                           facet_var         = facet_var,
                           scales            = scales,
                           labs              = labs,
                           nrow              = nrow,
                           mention_method    = mention_method,
                           title             = title,
                           na.impute         = na.impute,
                           legend.position   = legend.position)
}


#' @rdname plot_projected_samples
#' @export
plot_lda_samples <- function(
   object,
   ...,
   listed_objects    = NULL,
   dims              = 1:2,
   color_var         = 'subgroup', # default_color_var(object) gives 'block' when present!
   color_values      = default_color_values(object, color_var),
   shape_var         = default_shape_var(object),
   size_var          = NULL,
   txt_var           = default_txt_var(object),
   txt_outliers      = FALSE,
   group_var         = NULL,
   split_var         = NULL,
   facet_var         = NULL,
   scales            = ifelse(is.null(facet_var), 'fixed' ,'free'),
   labs              = list(color = stringi::stri_trans_totitle(color_var), shape = stringi::stri_trans_totitle(shape_var)),
   nrow              = NULL,
   mention_method    = ifelse(is.null(facet_var), TRUE, FALSE),
   title             = NULL,
   na.impute         = FALSE,
   legend.position   = 'right'
){
   plot_projected_samples( object,
                           ...,
                           listed_objects    = listed_objects,
                           method            = 'lda',
                           implementation    = character(0),
                           dims              = dims,
                           color_var         = color_var,
                           color_values      = color_values,
                           shape_var         = shape_var,
                           size_var          = size_var,
                           txt_var           = txt_var,
                           txt_outliers      = txt_outliers,
                           group_var         = group_var,
                           split_var         = split_var,
                           facet_var         = facet_var,
                           scales            = scales,
                           labs              = labs,
                           nrow              = nrow,
                           mention_method    = mention_method,
                           title             = title,
                           na.impute         = na.impute,
                           legend.position   = legend.position)
}


#' @rdname plot_projected_samples
#' @export
plot_pls_samples <- function(
   object,
   ...,
   listed_objects    = NULL,
   implementation    = c('mixOmics::plsda', 'mixOmics::splsda', 'ropls::opls')[1],
   dims              = 1:2,
   color_var         = 'subgroup', # default_color_var(object) gives 'block' when present!
   color_values      = default_color_values(object, color_var),
   shape_var         = default_shape_var(object),
   size_var          = NULL,
   txt_var           = default_txt_var(object),
   txt_outliers      = FALSE,
   group_var         = NULL,
   split_var         = NULL,
   facet_var         = NULL,
   scales            = ifelse(is.null(facet_var), 'fixed' ,'free'),
   labs              = list(color = stringi::stri_trans_totitle(color_var), shape = stringi::stri_trans_totitle(shape_var)),
   nrow              = NULL,
   mention_method    = ifelse(is.null(facet_var), TRUE, FALSE),
   title             = NULL,
   na.impute         = FALSE,
   legend.position   = 'right'
){
   plot_projected_samples( object,
                           ...,
                           listed_objects    = listed_objects,
                           method            = 'pls',
                           implementation    = implementation,
                           dims              = dims,
                           color_var         = color_var,
                           color_values      = color_values,
                           shape_var         = shape_var,
                           size_var          = size_var,
                           txt_var           = txt_var,
                           txt_outliers      = txt_outliers,
                           group_var         = group_var,
                           split_var         = split_var,
                           facet_var         = facet_var,
                           scales            = scales,
                           labs              = labs,
                           nrow              = nrow,
                           mention_method    = mention_method,
                           title             = title,
                           na.impute         = na.impute,
                           legend.position   = legend.position)
}

