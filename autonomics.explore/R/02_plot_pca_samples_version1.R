make_sample_scores_title <- function(object, method){
   selector <- if (method %in% c('pca', 'pls')){
                  selector <- matrixStats::rowAnys(!is.na(autonomics.import::exprs(object))) &
                              matrixStats::rowAnys(!is.infinite(autonomics.import::exprs(object)))
               } else {
                  selector <- matrixStats::rowAlls(!is.na(autonomics.import::exprs(object))) &
                              matrixStats::rowAlls(!is.infinite(autonomics.import::exprs(object)))
               }
   sprintf('%s on %d / %d features', toupper(method), sum(selector), length(selector))
}

#' Plot PCA/LDA/PLS sample scores
#' @param object         SummarizedExperiment, eSet, or Elist
#' @param method         'pca', 'lda', 'pls'
#' @param implementation which implementation of the method to use
#' @param dims           dimensions
#' @param color_var      svar mapped to color
#' @param color_values   color values vector
#' @param shape_var      svar mapped to shape
#' @param size_var       svar mapped to size
#' @param txt_var        svar mapped to txt
#' @param group_var      svar mapped to group
#' @param split_var      svar on which to split data prior to transformation
#' @param facet_var      svar on which to facet sample biplot
#' @param scales         'free' or 'fixed'
#' @param nrow           integer
#' @param title          title
#' @param na.impute      TRUE or FALSE
#' @param legend.position character
#' @param ...            passed on to ggplot2::theme(legend.position = .)
#' @examples 
#' require(magrittr)
#' if (require(subramanian.2016)){
#'    object <- subramanian.2016::metabolon
#'    object %>% autonomics.explore::plot_pca_samples1()
#'    object %>% autonomics.explore::plot_sma_samples1()
#'    object %>% autonomics.explore::plot_lda_samples1()
#'    object %>% autonomics.explore::plot_pls_samples1()
#'    object %>% autonomics.explore::plot_pca_samples1(color_var = 'condition')
#'    object %>% autonomics.explore::plot_pca_samples1(color_var = 'condition', size = 'time')
#'    object %>% autonomics.explore::plot_lda_samples1(color_var = 'condition', size = 'time')
#'    object %>% autonomics.explore::plot_pca_samples1(facet_var = 'condition')
#'    object %>% autonomics.explore::plot_pca_samples1(split_var = 'condition')
#'    object %>% autonomics.explore::plot_pca_samples1(
#'                 color_var = 'subgroup',  size_var = 'time', split_var = 'condition')
#'    object %>% autonomics.explore::plot_pca_samples1(
#'                  color_var = 'condition', size_var = 'time', split_var = 'condition')
#' }
#' 
#' # STEM CELL COMPARISON
#' if (require(autonomics.data)){
#'    object <- autonomics.data::stemcomp.proteinratios
#'    object %>% autonomics.explore::plot_pca_samples1()
#'    object %>% autonomics.explore::plot_lda_samples1()
#' }
#' 
#' # STEM CELL DIFFERENTIATION
#' if (require(billing.differentiation.data)){
#'    object <- autonomics.data::stemdiff.proteinratios
#'    object %>% autonomics.explore::plot_pca_samples1()
#'    object %>% autonomics.explore::plot_lda_samples1()
#'    object %>% autonomics.explore::plot_pls_samples1()
#' }
#' if (require(atkin.2014)){
#'    object <- atkin.2014::soma
#'    object %>% autonomics.explore::plot_lda_samples1()
#'    object %>% autonomics.explore::plot_lda_samples1(group_var = 'block')
#'    object %>% autonomics.explore::plot_lda_samples1(group_var = 'block', facet_var = 'condition')
#'    object %>% autonomics.explore::plot_lda_samples1(group_var = 'block', split_var = 'condition')
#' }
#' @importFrom data.table   data.table   :=
#' @export
plot_projected_samples1 <- function(
   object, 
   method         = 'pca',
   implementation = NULL,
   dims           = 1:2, 
   color_var      = 'subgroup', # autonomics.plot::default_color_var(object) gives 'block' when present!
   color_values   = autonomics.plot::default_color_values(object, color_var),
   shape_var      = autonomics.plot::default_shape_var(object), 
   size_var       = NULL, 
   txt_var        = autonomics.plot::default_txt_var(object), 
   group_var      = NULL,
   split_var      = NULL, 
   facet_var      = split_var,
   scales         = if(is.null(split_var)) 'fixed' else 'free',
   nrow           = NULL, 
   title          = NULL, 
   na.impute      = FALSE, 
   legend.position = 'right'
){
   
   # Assert
   autonomics.import::assert_is_valid_object(object)
   assertive.properties::assert_is_non_empty(autonomics.import::exprs(object))
   object %<>% autonomics.plot::validify_shape_values(shape_var)

   # Transform
   object_list <- object %>% autonomics.import::split_by_svar(split_var)
   projection  <- mapply(autonomics.plot::project, 
                         object_list, 
                         MoreArgs = list(method = method, implementation = implementation,  na.impute = na.impute), 
                         SIMPLIFY = FALSE)
   DT <- mapply(autonomics.plot::make_projected_samples_df, 
                object    = object_list, 
                projection   = projection, 
                MoreArgs = list(dims = dims,  color_var = color_var,  shape_var = shape_var,  size_var = size_var,  
                                txt_var = txt_var,  facet_var = facet_var,  split_var = split_var,  group_var = group_var),
                SIMPLIFY = FALSE) %>% 
         data.table::rbindlist()

   # Plot
   p <- ggplot2::ggplot(DT)
      
   # Points
   p <- p + ggplot2::geom_point(ggplot2::aes_string(x = 'x', y = 'y', color = 'color', shape = 'shape', size = 'size'))
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

   # Facet
   if (!is.null(facet_var)) p <- p + ggplot2::facet_wrap(stats::as.formula('~ facet'), scales = scales, nrow = nrow)
   
   # Axes
   p <- p + ggplot2::geom_vline(xintercept = 0, linetype = 'dashed') +
            ggplot2::geom_hline(yintercept = 0, linetype = 'dashed')    
   
   # Axis labels
   xlab <- paste0('X', dims[1])
   ylab <- paste0('X', dims[2])
   if (is.null(split_var)){ var <- projection[[1]]$var[dims]
                            xlab %<>% paste0(' (', var[dims[1]], '%)')
                            ylab %<>% paste0(' (', var[dims[2]], '%)')}
   p <- p + ggplot2::xlab(xlab) + ggplot2::ylab(ylab)

   # Title
   default_title <- make_sample_scores_title(object, method)
   title <- if (is.null(title)) default_title   else   sprintf('%s: %s', title, default_title)
   p <- p + ggplot2::ggtitle(title)
   p <- p + ggplot2::theme(legend.position = legend.position)
   
   # Return
   p
   #suppressWarnings(print(p))
}

#' @rdname plot_projected_samples1
#' @export
plot_pca_samples1 <- function(object, ...){
   plot_projected_samples1(object = object, method = 'pca', ...)
}

#' @rdname plot_projected_samples1
#' @export
plot_sma_samples1 <- function(object, ...){
   plot_projected_samples1(object = object, method = 'sma', ...)
}

#' @rdname plot_projected_samples1
#' @export
plot_lda_samples1 <- function(object, ...){
   plot_projected_samples1(object = object, method = 'lda', ...)
}

#' @rdname plot_projected_samples1
#' @export
plot_pls_samples1 <- function(object, ...){
   plot_projected_samples1(object = object, method = 'pls', ...)
}

