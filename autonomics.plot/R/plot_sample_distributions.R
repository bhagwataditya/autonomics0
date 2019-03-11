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
#'      object %>% plot_sample_distributions()
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

#' Prepare plot datatable
#' @param exprsmat exprs matrix (feature x sample)
#' @param sdata    sample dataframe
#' @return datatable
#' @examples
#' if (require(autonomics.data)){
#'
#'    # STEM CELL COMPARISON
#'    require(magrittr)
#'    prepare_plot_dt(stemcomp.proteinratios)
#'    prepare_plot_dt(autonomics.import::exprs(stemcomp.proteinratios),
#'                    autonomics.import::sdata(stemcomp.proteinratios))
#'
#'    # GLUTAMINASE
#'    prepare_plot_dt(glutaminase)
#'    prepare_plot_dt(autonomics.import::exprs(glutaminase),
#'                    autonomics.import::sdata(glutaminase))
#'
#' }
#' @importFrom magrittr %>%
#' @export
prepare_plot_dt <- function(object, ...){
   UseMethod('prepare_plot_dt', object)
}

#' @rdname prepare_plot_dt
#' @export
prepare_plot_dt.matrix <- function(object, sdata, ...){
   data.table::data.table(feature_id = rownames(object), object)     %>%
   data.table::melt(id.vars = 'feature_id', variable.name = 'sample_id') %>%
   merge(sdata, by = 'sample_id', all.x=TRUE, sort = FALSE)
}

#' @rdname prepare_plot_dt
#' @export
prepare_plot_dt.SummarizedExperiment <- function(object, ...){
   prepare_plot_dt.matrix(
      autonomics.import::exprs(object),
      autonomics.import::sdata(object)
   )
}


#=================================================================================================


#' Plot sample densities
#' @param object  matrix, SummarizedExperiment, or list of SummarizedExperiments
#' @param sdata   dataframe
#' @param color   string (varname). For a list object, an equidimensional stringvector (with varlevels) is also allowed.
#' @param ...     passed to ggplot2::aes()
#' @param facet   string (varname). For a list object, an equidimensional stringvector (with varlevels) is also allowed.
#' @return ggplot2 object
#' @examples
#' if (require(autonomics.data)){
#'
#'    # STEM CELL COMPARISON
#'    require(magrittr)
#'    plot_sample_densities(stemcomp.proteinratios)
#'    plot_sample_densities(object = autonomics.import::exprs(stemcomp.proteinratios),
#'                          sdata  = autonomics.import::sdata(stemcomp.proteinratios))
#'
#'    # GLUTAMINASE
#'    plot_sample_densities(glutaminase, color = 'EXPERIMENT', facet = c('CONCENTRATION', 'TIME_POINT'))
#'    plot_sample_densities(object = autonomics.import::exprs(glutaminase),
#'                          sdata  = autonomics.import::sdata(glutaminase))
#'    plot_sample_densities(object = list(glutaminase %>% magrittr::extract(, 1:12),
#'                                        glutaminase %>% magrittr::extract(, 13:ncol(.))),
#'                          color = c('UT_h10', 'other'),
#'                          facet = c('UT_h10', 'other'))
#' }
#' @export
plot_sample_densities <- function(object, ...){
   UseMethod('plot_sample_densities', object)
}


#' Get sample density aesthetic mappings
#' @param plotdt  plot datable
#' @param x       string
#' @param group   string
#' @param color   string
#' @param ... passed to ggplot2::aes_string()
#' @return
#' @examples
#' if (require(autonomics.data)){
#'    sample_density_aes(prepare_plot_dt(autonomics.import::exprs(glutaminase),
#'                                       autonomics.import::sdata(glutaminase)))
#' }
#' @export
sample_density_aes <- function(
   plotdt,
   x     = 'value',
   group = 'sample_id',
   color = if ('subgroup' %in% names(plotdt)) 'subgroup' else 'sample_id',
   ...
){
   ggplot2::aes_string(x=x, group=group, color=color, ...)
}



#' @rdname plot_sample_densities
#' @export
plot_sample_densities.matrix <- function(
   object,
   sdata,
   mapping = ,
   facet = NULL,
   xlab  = NULL,
   ylab  = NULL,
   title = NULL
){
   assertive.sets::assert_is_subset(color, names(sdata))
   if (!is.null(facet)) assertive.sets::assert_is_subset(facet, names(sdata))

   p <- ggplot2::ggplot(data = prepare_plot_dt(object, sdata),
                        mapping = mapping) +
        ggplot2::theme_bw() +
        ggplot2::geom_line(stat = 'density', na.rm = TRUE) # geom_lined prefered, since geom_density draws a horizontal line

   if (!is.null(facet))                    p <- p + ggplot2::facet_wrap(facet)
   if (length(unique(sdata[[color]]))>25)  p <- p + ggplot2::guides(color = "none")
   p <- p + ggplot2::xlab(xlab) + ggplot2::ylab(ylab) + ggplot2::ggtitle(title)

   p

}


#' @rdname plot_sample_densities
#' @export
plot_sample_densities.SummarizedExperiment <- function(
   object,
   color = if ('subgroup' %in% autonomics.import::svars(object)) 'subgroup' else 'sample_id',
   ...,
   facet = NULL
){
   plot_sample_densities.matrix(object = autonomics.import::exprs(object),
                                sdata  = autonomics.import::sdata(object),
                                color  = color,
                                ...,
                                facet  = facet)
}

#' @rdname plot_sample_densities
#' @export
plot_sample_densities.list <- function(
   object,
   color = 'sample_id',
   ...,
   facet = NULL
){
   # Assert
      objectclass <- Reduce(assertive.base::assert_are_identical, vapply(object, class, character(1)))
      nfeatures   <- Reduce(assertive.sets::assert_are_set_equal, vapply(object, nrow, integer(1)))
      fnames      <- Reduce(assertive.sets::assert_are_set_equal, lapply(object, rownames))


   # Add svars on the fly if required
      for (var in c('color', 'facet', names(list(...)))){
         if (length(get(var)) == length(object)){
            object <- mapply(function(sumexp, value){sumexp[[var]] <- value; sumexp}, sumexp=object, value=get(var), SIMPLIFY=TRUE)
            assign(var, var)
         }
      }

   # Bind cols
      object <- switch(objectclass, SummarizedExperiment = Reduce(SummarizedExperiment::cbind, object),
                                    matrix               = Reduce(cbind, object))
   # Plot
      object %>% plot_sample_densities(color = color, ..., facet = facet)
}

#===================================================================


#' Plot sample violins
#'
#' @param object              SummarizedExperiment or matrix
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
#'       object %>% plot_sample_densities(x = 'cod', color_var = 'BT')
#'
#'       # faceting works, but is not useful here, as there is no block variable
#'       object %>% plot_sample_distributions_v1(x = 'cod', facet_var = 'sex', color_var = 'BT')
#'
#'    # STEM CELL COMPARISON
#'       object <- autonomics.data::stemcomp.proteinratios
#'       object %>% plot_sample_distributions_v1()
#'
#'    # GLUTMINASE
#'       x <- autonomics.data::glutaminase
#'       x %<>% autonomics.import::exprs()
#'       object %>% plot_sample_distributions_v1()
#' }
#' @importFrom  magrittr  %>%
#' @export
plot_sample_violins <- function(object, ...) UseMethod('plot_sample_violins', object)



#' @rdname plot_sample_densities
#' @export
plot_sample_violins.matrix <- function(object, horiz = TRUE, color_var = 'sample_id'){

   plotDF <- prep_column_densities(object)
   p <- ggplot2::ggplot(plotDF) +
        ggplot2::theme_bw() +
        ggplot2::guides(fill = FALSE, color = FALSE) +
        ggplot2::geom_violin( ggplot2::aes_string(x = 'sample_id', y = 'value', fill = color_var),
                              na.rm        = TRUE) +
        ggplot2::geom_boxplot(ggplot2::aes_string(x= 'sample_id', y = 'value', color = color_var),
                              fill         = 'white',
                              width        = 0.1,
                              outlier.size = 1.7,
                              na.rm        = TRUE)
   p <- p + ggplot2::xlab(NULL)
   p <- if (horiz){p + ggplot2::coord_flip()
        } else {   p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1))}

}






#' @rdname plot_sample_densities
#' @export
plot_sample_violins.SummarizedExperiment <- function(
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
   plotDF <- prep_column_densities(object)

   # Draw plot
   p <- ggplot2::ggplot(plotDF) +
        ggplot2::geom_violin( ggplot2::aes_string(x = 'x', y = 'value'), fill = 'gray', na.rm = TRUE) +
        ggplot2::geom_boxplot(ggplot2::aes_string(x = 'x', y = 'value',  fill = 'white',
                                                  color = if(is.null(color_var)) NULL else 'color'),
                                                  width = 0.1,
                                                  outlier.size = 1.7,
                                                  na.rm = TRUE)

   # Add custom color info
   if (!is.null(color_values)) p <- p + ggplot2::scale_color_manual(values = color_values)

   # Add individually plotted features
   if(!is.null(displayed_features)){
      filtered_object <- object %>% magrittr::extract(displayed_features,) %>% shape_for_plotting(x)
      if(length(displayed_features) == 1){
         p <- p + ggplot2::geom_point( data = filtered_object, ggplot2::aes_string(x = 'x', y = 'value'))
      } else {
         p <- p + ggplot2::geom_jitter(data = filtered_object, ggplot2::aes_string(x = 'x', y = 'value'))
      }
   }

   # Add faceting info
   if (!is.null(facet_var)) p <- p + ggplot2::facet_grid(facet ~ .)

   # Add labels & title and return
   p <- p + ggplot2::ylab(ylab) +
            ggplot2::xlab(xlab) +
            ggplot2::ggtitle(title) +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1))

   if (!is.null(file)) p  %>% print_sample_distributions(file)

   # Return
   return(p)
}



#' Print sample distributions
#' @param ggplot_obj  ggplot object
#' @param file        file to print to
#' @importFrom   magrittr   %>%   %<>%
#' @export
print_sample_distributions <- function(ggplot_obj, file){
   ggplot_obj %>% autonomics.support::print2pdf(file, height = 5, width = 5)
}
