#' Plot distribution across features for each sample
#'
#' The distributions are plotted as a combination of a violin plot and a box plot.
#' A separate violin is created for each value of 'x' (rows) and 'facet_var (column).
#' The box is colored per value of color_var
#'
#' @param object             SummarizedExperiment
#' @param ...                additonal SummarizedExperiments or matrices AND further parameters for \code{\link{core_ggplot_from_data}}
#' @param x                  string: sample variable mapped to x axis
#' @param color_var          string: sample variable mapped to color
#' @param color_values       string vector: color values vector (names = subgroups, values = colors)
#' @param displayed_features number/string vector: features to be displayed in the sample distributions (vector of numeric indexes or character \code{feature_id}s)
#' @param jitter_features    TRUE/FALSE: whether to use geom_jitter or not for the displayed features
#' @param horizontal         TRUE/FALSE: whether to lay the plot out horizontally (measurements on the x, sample names on the y axis)
#' @param xlab               string: label of y axis
#' @param ylab               string: label of y axis
#' @param title              string: title
#' @return ggplot2 object
#' @seealso \code{\link{core_ggplot_from_data}}
#' @author Aditya Bhagwat, Johannes Graumann
#' @examples
#' if (require(autonomics.data)){
#'    require(magrittr)
#'    object <- autonomics.data::stemcomp.proteinratios
#'
#'    # Simple plot
#'      object %>% plot_sample_distrsibutions()
#'
#'    # Same thing vertical
#'      object %>% plot_sample_distributions(horizontal = FALSE)
#'
#'    # More complex using facetting
#'      object %>%
#'      plot_sample_distributions(
#'         2^autonomics.import::exprs(autonomics.data::stemcomp.proteinratios),
#'         facet2_var = c('Logarithmized', 'Raw'))
#'
#'    # Yet more crazy explicitly marking the feature with the maximum value
#'      object %>%
#'      plot_sample_distributions(2^autonomics.import::exprs(object),
#'                                facet2_var         = c('Logarithmized', 'Raw'),
#'                                displayed_features = autonomics.import::exprs(object)        %>%
#'                                                     equals(autonomics.import::exprs(object) %>%
#'                                                     max(na.rm = TRUE))                      %>%
#'                                                     which(arr.ind = TRUE)                   %>%
#'                                                     magrittr::extract(,'row'))
#' }
#' @importFrom  magrittr  %>%
#' @export
plot_sample_distributions <- function(
   object,
   ...,
   x                  =  NULL,
   color_var          =  default_color_var(object),
   color_values       =  default_color_values(object, color_var),
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






#==============================================
# OLD VERSION
#==============================================



check_args_of_plot_sample_distributions_v1 <- function(object, x, facet_var, color_var, displayed_features, ylab, xlab, title){
   assertive.types::assert_is_any_of(object, c("eSet", 'EList', 'SummarizedExperiment'))
   if (!is.null(x)){
      assertive.types::assert_is_a_string(x)
      assertive.sets::assert_is_subset(x, autonomics.import::svars(object))
   }
   if (!is.null(facet_var)){
      assertive.types::assert_is_a_string(facet_var)
      assertive.sets::assert_is_subset(facet_var, autonomics.import::svars(object))
   }
   if (!is.null(color_var)){
      assertive.types::assert_is_a_string(color_var)
      assertive.sets::assert_is_subset(color_var, autonomics.import::svars(object))
   }
   if (!is.null(displayed_features)){
      autonomics.import::assert_all_are_valid_features(displayed_features, object)
   }
   if (!is.null(xlab)){
      assertive.types::assert_is_a_string(xlab)
   }
   if (!is.null(ylab)){
      assertive.types::assert_is_a_string(ylab)
   }
   if (!is.null(title)){
      assertive.types::assert_is_a_string(title)
   }
}

#' Plot distribution across features for each sample
#'
#' The distributions are plotted as a combination of a violin plot and a box plot.
#' A separate violin is created for each value of 'x' (rows) and 'facet_var (column).
#' The box is colored per value of color_var
#'
#' @param object              SummarizedExperiment
#' @param x                   string: sample variable mapped to x axis
#' @param color_var           string: sample variable mapped to color
#' @param color_values        string vector: names = subgroups, values = colors
#' @param facet_var           string: sample variable for faceting
#' @param displayed_features  features to be displayed in the sample distributions (vector of numeric indexes or character \code{feature_id}s)
#' @param xlab                label of y axis (character)
#' @param ylab                label of y axis (character)
#' @param title               title (character)
#' @param file                file to which to print plot
#' @return ggplot2 object
#' @author Aditya Bhagwat
#' @examples
#' if (require(autonomics.data)){
#'
#'    # ALL
#'       require(magrittr)
#'       object <- autonomics.data::ALL[, 1:30]
#'       object %>% plot_sample_distributions_v1(x = 'cod', color_var = 'BT')
#'
#'       # faceting works, but is not useful here, as there is no block variable
#'       object %>% plot_sample_distributions_v1(x = 'cod', facet_var = 'sex', color_var = 'BT')
#'
#'    # STEM CELL COMPARISON
#'       object <- autonomics.data::stemcomp.proteinratios
#'       object %>% plot_sample_distributions_v1()
#'
#'    # GLUTMINASE
#'       object <- autonomics.data::glutaminase
#'       object %>% plot_sample_distributions_v1()
#' }
#' @importFrom  magrittr  %>%
#' @export
plot_sample_distributions_v1 <- function(
   object,
   x                  =  NULL,
   color_var          =  default_color_var(object),
   color_values       =  default_color_values(object, color_var),
   facet_var          =  NULL,
   displayed_features =  NULL,
   ylab               =  '',
   xlab               =  '',
   title              =  sprintf('Per sample distribution of all %d features', nrow(object)),
   file               =  NULL
){

   # Process arguments
   check_args_of_plot_sample_distributions_v1(object, x, facet_var, color_var, displayed_features, ylab, xlab, title)
   if (!is.null(color_var)){
      if (is.numeric(object[[color_var]])){ # essential to avoid an error in ggplot!
         object[[color_var]] <- factor(as.character(autonomics.import::sdata(object)[[color_var]]))
      }
   }
   if (!is.null(x)){
      if (is.numeric(object[[x]])){    # required to have a separate box per x value
         object[[x]] <- factor(as.character(autonomics.import::sdata(object)[[x]]))
      }
   }
   if (!is.null(facet_var)){
      if (is.numeric(object[[facet_var]])){# essential for faceting
         object[[facet_var]] <- factor(as.character(autonomics.import::sdata(object)[[facet_var]]))
      }
   }

   # Combine sample IDs & exprs and plot distributions
   plotDF <- shape_for_plotting(object, x)

   # Add color info
   if (!is.null(color_var)){
      plotDF <- cbind(plotDF, color = autonomics.import::sdata(object)[[color_var]])
   }
   # Add faceting info
   if (!is.null(facet_var)){
      plotDF <- cbind(facet = autonomics.import::sdata(object)[[facet_var]], plotDF)
   }
   # Draw plot
   p <- ggplot2::ggplot(plotDF) +
        ggplot2::geom_violin(ggplot2::aes_string(x = 'x', y = 'value'), fill = 'gray', na.rm = TRUE)
   p <- p + ggplot2::geom_boxplot(
               ggplot2::aes_string(
                  x = 'x',
                  y = 'value',
                  color = if(is.null(color_var)) NULL else 'color'),
               fill = 'white',
               width = 0.1,
               outlier.size = 1.7,
               na.rm = TRUE)

   # Add custom color info
   if (!is.null(color_values)){
      p <- p + ggplot2::scale_color_manual(values = color_values)
   }

   # Add individually plotted features
   if(!is.null(displayed_features)){
      filtered_object <- object %>%
         magrittr::extract(displayed_features,) %>%
         shape_for_plotting(x)
      if(
         displayed_features %>%
           length() %>%
           magrittr::equals(1))
      {
         p <- p +
            ggplot2::geom_point(
               data = filtered_object,
               ggplot2::aes_string(x = 'x', y = 'value'))
      } else {
         p <- p +
            ggplot2::geom_jitter(
               data = filtered_object,
               ggplot2::aes_string(x = 'x', y = 'value'))
      }
   }

   # Add faceting info
   if (!is.null(facet_var)){
      p <- p + ggplot2::facet_grid(facet ~ .)
   }

   # Add labels & title and return
   p <- p +
      ggplot2::ylab(ylab) +
      ggplot2::xlab(xlab) +
      ggplot2::ggtitle(title) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1))

   if (!is.null(file)){
     p  %>% print_sample_distributions(file)
   }

   return(p)
}

shape_for_plotting <- function(object, x){
   xVec <- autonomics.support::factorify(autonomics.import::snames(object))
   if (!is.null(x)) {
      xVec <- object[[x]]
   }
   prepDF <- data.frame(x = xVec, t(autonomics.import::exprs(object)))
   plotDF <- suppressMessages(
      tidyr::gather(
         prepDF,
         "feature",
         "value",
         setdiff(colnames(prepDF), "x")))
   plotDF %>%
      return()
}
#' Print sample distributions
#' @param ggplot_obj  ggplot object
#' @param file        file to print to
#' @importFrom   magrittr   %>%   %<>%
#' @export
print_sample_distributions <- function(ggplot_obj, file){
   ggplot_obj %>% autonomics.support::print2pdf(file, height = 5, width = 5)
}
