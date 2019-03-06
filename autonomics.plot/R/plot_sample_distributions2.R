#' Plot distribution across features for each sample
#'
#' The distributions are plotted as a combination of a violin plot and a box plot.
#' A separate violin is created for each value of 'x' (rows) and 'facet_var (column).
#' The box is colored per value of color_var
#'
#' @param object             eSet
#' @param ...                additonal eSets, ELists, SummarizedExperiments or matrices AND further parameters for \code{\link{core_ggplot_from_data}}
#' @param x                  sample variable mapped to x axis
#' @param color_var          sample variable mapped to color
#' @param color_values       color values vector (names = subgroups, values = colors)
#' @param displayed_features features to be displayed in the sample distributions (vector of numeric indexes or character \code{feature_id}s)
#' @param jitter_features    whether to use geom_jitter or not for the displayed features
#' @param horizontal         whether to lay the plot out horizontally (measurements on the x, sample names on the y axis)
#' @param xlab               label of y axis (character)
#' @param ylab               label of y axis (character)
#' @param title              title (character)
#' @return ggplot2 object
#' @seealso \code{\link{core_ggplot_from_data}}
#' @author Aditya Bhagwat, Johannes Graumann
#' @importFrom  magrittr  %>%
#' @examples
#' if (require(autonomics.data)){
#'    require(magrittr)
#'    object <- autonomics.data::stemcomp.proteinratios
#'
#'    # Simple plot
#'      object %>% autonomics.plot::plot_sample_distributions2()
#'
#'    # Same thing vertical
#'      object %>% plot_sample_distributions2(horizontal = FALSE)
#'
#'    # More complex using facetting
#'      object %>% plot_sample_distributions2(
#'         2^autonomics.import::exprs(autonomics.data::stemcomp.proteinratios),
#'         facet2_var = c('Logarithmized', 'Raw'))
#'
#'    # Yet more crazy explicitly marking the feature with the maximum value
#'      object %>%
#'      plot_sample_distributions2(2^autonomics.import::exprs(object),
#'                                 facet2_var         = c('Logarithmized', 'Raw'),
#'                                 displayed_features = autonomics.import::exprs(object)        %>%
#'                                                      equals(autonomics.import::exprs(object) %>%
#'                                                      max(na.rm = TRUE))                      %>%
#'                                                      which(arr.ind = TRUE)                   %>%
#'                                                      magrittr::extract(,'row'))
#' }
#' @export
plot_sample_distributions2 <- function(
   object,
   ...,
   x                  =  NULL,
   color_var          =  autonomics.plot::default_color_var(object),
   color_values       =  autonomics.plot::default_color_values(object, color_var),
   displayed_features =  NULL,
   jitter_features    =  FALSE,
   horizontal         =  TRUE,
   xlab               =  '',
   ylab               =  '',
   title              =  sprintf('Per sample distribution of all %d features', nrow(object))
){

   assertive.types::assert_is_a_bool(jitter_features)

   p <- core_ggplot_from_data(object     = object,
                              ...,
                              MARGIN     = 'samples',
                              fill_var   = NULL,
                              marked_ids = displayed_features,
                              horizontal = horizontal,
                              xlab       = xlab,
                              ylab       = ylab,
                              title      = title)

   point_size <- 1

   if(horizontal){
      p <- p + ggstance::geom_violinh(fill = 'gray', na.rm = TRUE) +
               ggstance::geom_boxploth(fill = 'white', width = 0.1, outlier.size = point_size, na.rm = TRUE)
   } else {
      p <- p + ggplot2::geom_violin(fill = 'gray', na.rm = TRUE) +
               ggplot2::geom_boxplot(fill = 'white', width = 0.1, outlier.size = point_size, na.rm = TRUE)
   }

   # Add custom color info
   if (!is.null(color_values)){
      p <- p + ggplot2::scale_color_manual(values = color_values)
   }

   # Add individually plotted features
   if(!is.null(displayed_features)){
      if(jitter_features){
         p <- p + ggplot2::geom_jitter(ggplot2::aes_string(alpha = 'marked_plotvar'),
                                       color = 'black',
                                       size  = point_size)

      } else {
         p <- p + ggplot2::geom_point(ggplot2::aes_string(alpha = 'marked_plotvar'),
                                      color = 'black',
                                      size  = point_size)
      }
      p <- p + ggplot2::scale_alpha_discrete(range = c(0,1), guide = FALSE)
   }

   # Return the object
   return(p)
}
