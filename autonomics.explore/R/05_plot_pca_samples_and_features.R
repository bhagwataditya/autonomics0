
#' Create graphics layout for mpa feature + sample plot
#' @param geom feature plot (character)
#' @return graphics layout matrix
#' @importFrom magrittr %>% 
#' @export
layout_sample_projections_and_features <- function(geom){
   if (geom == 'hbar'){
      rbind(c(2,2,2,2,2,2) %>% rep(2) %>% matrix(nrow = 2, byrow = TRUE),
            c(0,1,1,1,0,0) %>% rep(3) %>% matrix(nrow = 3, byrow = TRUE),
            c(3,3,3,3,3,3) %>% rep(2) %>% matrix(nrow = 2, byrow = TRUE)
      )
   } else {
         rbind(c(rep(0,3), rep(1,6), rep(0,4)) %>% rep(4) %>% matrix(nrow = 4, byrow = TRUE),
               c(rep(2,6), rep(0,1), rep(3,6)) %>% rep(7) %>% matrix(nrow = 7, byrow = TRUE))

   }
}

#' Plot pca samples and features
#'
#' Plots pca sample and features results in a combined plot
#'
#' @param object         eset
#' @param method         'pca', 'lda', or 'pls'
#' @param implementation NULL or "package::function"
#' @param dims           pc dimensions to plot
#' @param na.impute      TRUE or FALSE
#' @param geom           value in \code{\link[autonomics.plot]{FEATURE_PLOTS}}
#' @param result_dir     NULL or result directory path
#' @param x              svar mapped to x in feature plots
#' @param color_var      svar mapped to color in sample and feature plots
#' @param color_values   color vector (names = color_var levels, values = colors)
#' @param shape_var      svar mapped to shape in sample and feature plots
#' @param group_var      svar mapped to group in feature plots
#' @param txt_var        svar mapped to txt in feature plots
#' @param line           whether to connect points in feature plot with line (logical)
#' @param n              no of features to be plotted
#' @param ...            passed to plot_projected_samples_and_features
#' @examples
#' require(magrittr)
#' result_dir <- tempdir() %T>% message()
#' \dontrun{
#' 
#'      # STEM CELL COMPARISON
#'      if (require(autonomics.data)){
#'         object <- autonomics.data::stemcomp.proteinratios
#'         object %>% plot_pca_samples_and_features()
#'         object %>% plot_pca_samples_and_features(geom = 'violin')
#'         object %>% plot_pca_samples_and_features(geom = 'violin', na.impute = TRUE)
#'         object %>% plot_pca_samples_and_features(result_dir = result_dir)
#'         object %>% plot_pca_samples_and_features(result_dir = result_dir, geom = 'bar')
#'      }
#'      
#'      # STEM CELL DIFFERENTIATION
#'      if (require(billing.differentiation.data)){
#'         object <- billing.differentiation.data::rna.voomcounts
#'         object %>% plot_pca_samples_and_features()
#'      }
#'      
#'      # GLUTAMINASE
#'      if (require(autonomics.data)){
#'         object <- autonomics.data::glutaminase
#'         object %>% plot_pca_samples_and_features(n=2)
#'      }
#'      
#'      if (require(subramanian.2016)){
#'         object <- subramanian.2016::metabolon
#'         object %>% plot_pca_samples_and_features(n = 4)
#'         object %>% plot_lda_samples_and_features(n = 4)
#'         object %>% plot_pca_samples_and_features(n = 4,
#'                       color_var = 'condition', geom = 'boxplot')
#'      }
#' }
#' @importFrom magrittr   %>%
#' @export
plot_projected_samples_and_features <- function(
   object,
   method,
   implementation = NULL,
   dims         = c(1,2),
   na.impute    = FALSE,
   geom = autonomics.plot::default_feature_plots(object)[1],
   result_dir   = NULL,
   x            = autonomics.plot::default_x(object, geom),
   color_var    = autonomics.plot::default_color_var(object),
   color_values = autonomics.plot::default_color_values(object, color_var),
   shape_var    = NULL,
   group_var    = autonomics.plot::default_group_var(object),
   txt_var      = autonomics.plot::default_txt_var(object),
   line         = autonomics.plot::default_line(object), 
   n            = 9
){
   # Check input args
   autonomics.import::assert_is_valid_eset(object)
   if (ncol(object) < 3){
      autonomics.support::cmessage('\tExit PCA: only %d samples', ncol(object))
      return(invisible(NULL))
   }
   assertive.sets::assert_is_subset(geom, autonomics.plot::FEATURE_PLOTS)

   # Plot
   feature_plot_args <- list(
      object          = object,
      method          = method,
      implementation  = implementation,
      x               = x,
      color_var       = color_var,
      color_values    = color_values,
      shape_var       = shape_var,
      group_var       = group_var,
      line            = line,
      na.impute       = na.impute,
      geom    = geom,
      legend.position = 'none', 
      n               = n
   )
   plotlist <- list(
      samples = autonomics.explore::plot_projected_samples(
                   object, method = method, implementation = implementation, color_var = color_var, color_values = color_values, 
                   shape_var = shape_var, txt_var = txt_var, dims = dims, na.impute = na.impute),
      features1 = autonomics.explore::plot_projection_features %>% do.call(c(feature_plot_args, list(dim = dims[1]))),
      features2 = autonomics.explore::plot_projection_features %>% do.call(c(feature_plot_args, list(dim = dims[2])))
   )

   # Print
   if (!is.null(result_dir)){
      dir.create(sprintf('%s', result_dir, method), showWarnings = FALSE)
      file_name <- sprintf('%s/%s%s%s_%s_%s.pdf',
                           result_dir, method, dims[1], dims[2], ifelse(na.impute, 'all', 'common'), geom)
      grDevices::pdf(file_name, width = 15, height = 12)
   }
   layout <- autonomics.explore::layout_sample_projections_and_features(geom)
   autonomics.plot::multiplot(plotlist = plotlist, layout = layout)
   if (!is.null(result_dir)){
      grDevices::dev.off()
   }

}


#' @rdname plot_projected_samples_and_features
#' @export
plot_pca_samples_and_features <- function(object, ...){
   plot_projected_samples_and_features(object, method = 'pca', ...)
}

#' @rdname plot_projected_samples_and_features
#' @export
plot_sma_samples_and_features <- function(object, ...){
   plot_projected_samples_and_features(object, method = 'sma', ...)
}

#' @rdname plot_projected_samples_and_features
#' @export
plot_lda_samples_and_features <- function(object, ...){
   plot_projected_samples_and_features(object, method = 'lda', ...)
}

#' @rdname plot_projected_samples_and_features
#' @export
plot_pls_samples_and_features <- function(object, ...){
   plot_projected_samples_and_features(object, method = 'pls', ...)
}
