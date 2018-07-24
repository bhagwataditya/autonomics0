get_pca_var <- function(method, dim){
   sprintf('%s%d', method, dim)
}

#' Order on feature loadings
#' @param object    SummarizedExperiment
#' @param method   'pca', 'sma', 'lda', 'pls'
#' @param dim       dimension
#' @param na.impute whether to impute (if applicable)
#' @return re-ordered sumexp
#' @importFrom magrittr %>% 
#' @export
order_on_feature_loadings <- function(object, method, dim, na.impute){
   pca_var <- get_pca_var(method, dim)
   loadings <- autonomics.import::fdata(object) %>% magrittr::extract2(pca_var)
   idx <- order(loadings, na.last = NA)
   object[idx, ]
}

#' Extract top and bottom of sumexp
#' @param object  SummarizedExperiment
#' @param n       number of top features to be selected
#' @return sumex with top n/2 and bottom n/2 rows
#' @importFrom magrittr   %>%
#' @export
extract_top_and_bottom <- function(object, n = 16){

   nfeature <- nrow(object)
   n1 <- ceiling(n/2)
   n2 <- n - n1

   top <- seq(1, n1)
   bottom <- seq(nfeature-n2+1, nfeature) # opposite direction

   object %>% magrittr::extract(c(top, bottom), )
}

#' Plot pca features
#'
#' Plots top pca features.
#' Uses pca factor loadings in fdata when available.
#' When not, add_pca_to_eset is run prior to plotting
#'
#' @param object           SummarizedExperiment, eSet, or Elist
#' @param method           'pca' or 'lda'
#' @param implementation   'character' or NULL
#' @param x                svar mapped to x
#' @param color_var        svar mapped to color
#' @param color_values        color vector (names = color_var levels, values = colors)
#' @param shape_var        svar mapped to shape
#' @param group_var        svar mapped to group
#' @param txt_var          svar mapped to txt
#' @param facet_var        svar used for faceting
#' @param line             whether to connect points with line (logical)
#' @param fvars            fvars used for plot annotation
#' @param dim              principal component dimension
#' @param n                number of top features to plot
#' @param na.impute        TRUE or FALSE
#' @param title            title
#' @param file             file
#' @param feature_plot     value in \code{\link[autonomics.plot]{FEATURE_PLOTS}}
#' @param legend.position  ggplot theme argument
#' @param ...              passed to plot_projection_features
#' @examples
#' require(magrittr)
#' # STEM CELL COMPARISON
#' if (require(autonomics.data)){
#'    object <- autonomics.data::stemcomp.proteinratios
#'    object %>% autonomics.explore::plot_pca_features(n=9)
#'    object %>% autonomics.explore::plot_pca_features(feature_plot = 'bars')
#' }
#' 
#' # STEM CELL DIFFERENTIATION
#' if (require(autonomics.data)){
#'    object <- autonomics.data::stemdiff.proteinratios
#'    object %>% autonomics.explore::plot_pca_features(n=9)
#'    object %>% autonomics.explore::plot_pca_features(dim=2, n=9)
#'    object %>% autonomics.explore::plot_pca_features()
#' }
#' 
#' # GLUTAMINASE
#' if (require(autonomics.data)){
#'    autonomics.data::glutaminase %>% autonomics.explore::plot_pca_features(
#'       color_var    = 'GROUP_DESCRIPTION', 
#'       color_values = c(Control = 'orange', Vehicle = 'red', `Concentration 1` = 'green', 
#'                       `Concentration 2` = 'blue'))
#' }
#' 
#' if (require(subramanian.2016)){
#'    object <- subramanian.2016::metabolon
#'    object %>% autonomics.explore::plot_pca_features(feature_plot = 'boxes')
#'    object %>% autonomics.explore::plot_pls_features(feature_plot = 'boxes')
#'    object %>% autonomics.explore::plot_lda_features(feature_plot = 'boxes')
#' }
#' 
#' if (require(atkin.2014)){
#'    object <- atkin.2014::soma
#'    object %>% autonomics.explore::plot_pca_features()
#' }
#' 
#' @importFrom magrittr  %>%
#' @export
plot_projection_features <- function(
   object,
   method,
   implementation  = NULL,
   feature_plot    = autonomics.plot::default_feature_plots(object)[1],
   x               = autonomics.plot::default_x(object, feature_plot),
   color_var       = autonomics.plot::default_color_var(object),
   color_values    = autonomics.plot::default_color_values(object, color_var),
   shape_var       = autonomics.plot::default_shape_var(object),
   group_var       = autonomics.plot::default_group_var(object),
   txt_var         = autonomics.plot::default_txt_var(object),
   facet_var       = NULL,
   line            = autonomics.plot::default_line(object),
   fvars           = autonomics.plot::default_fvars(object),
   dim             = 1,
   n               = 9,
   na.impute       = FALSE,
   title           = sprintf('X%d', dim),
   file            = NULL,
   legend.position = 'right'
){
   # Check input args
   autonomics.import::assert_is_valid_eset(object)
   if (ncol(object) < 3){
      autonomics.support::cmessage('\tExit %s: only %d samples', method, ncol(object))
      return(invisible(NULL))
   }
   assertive.sets::assert_is_subset(feature_plot, autonomics.plot::FEATURE_PLOTS)

   # Add projection to eset if required
   projection_dim_in_eset <- paste0(method, dim) %in% autonomics.import::fvars(object)
   if (!projection_dim_in_eset){
      # Wipe earlier accounts
      idx <- which(autonomics.import::fvars(object) %>% stringi::stri_detect_regex(paste0('^', method)))
      if (length(idx)>0) autonomics.import::fdata(object) %<>% magrittr::extract(, -idx)
      # Add projection
      object %<>% autonomics.explore::add_projection_to_eset(method, na.impute = na.impute, ndim = dim)
   }
   
   # Order on projection and write to file
   object %<>% autonomics.explore::order_on_feature_loadings(method = method, dim = dim, na.impute = na.impute) %>%
               autonomics.explore::extract_top_and_bottom(n=n)

   # Dispatch to plotfun
   plotfun <- utils::getFromNamespace(paste0('plot_feature_', feature_plot), 'autonomics.plot')

   # Prepare plotargs accordingly
   plotargs <- list(object          = object,
                    x               = x,
                    color_var       = color_var,
                    color_values    = color_values,
                    facet_var       = facet_var,
                    fvars           = fvars,
                    title           = title,
                    file            = file,
                    legend.position = legend.position)
   if (feature_plot!='hbars'){
      plotargs %<>% c(list(shape_var = shape_var,
                           group_var = group_var,
                           txt_var   = txt_var,
                           line      = line))
   }

   # Plot
   do.call(plotfun, plotargs)
}

#' @rdname plot_projection_features
#' @export
plot_pca_features <- function(object, ...){
   plot_projection_features(object, 'pca', ...)
}

#' @rdname plot_projection_features
#' @export
plot_sma_features <- function(object, ...){
   plot_projection_features(object, 'sma', ...)
}

#' @rdname plot_projection_features
#' @export
plot_lda_features <- function(object, ...){
   plot_projection_features(object, 'lda', ...)
}

#' @rdname plot_projection_features
#' @export
plot_pls_features <- function(object, ...){
   plot_projection_features(object, 'pls', ...)
}
