#' Plot Venn/Euler diagram
#' @param x             named list of character vectors
#' @param color_values  character vector
#' @param euler         logical(1): whether to use eulerr::euler() rather than VennDiagram::venn.diagram()
#' @param filename      character vector
#' @param width         width  in inches
#' @param height        height in inches
#' @param title         figure title
#' @param scaled        logical: whether or not to scale
#' @param cat.dist      numeric(1): category location (as distance from circle edge)
#' @param main.cex      numeric(1): title font size
#' @param main.pos      numeric(2): title location
#' @param ...           passed to either VennDiagram::venn.diagram(...) or eulerr::plot.euler(...)
#' @examples
#' x <- list(A = c('apple', 'pear'), B = c('pear', 'orange'))
#' filename <- NULL
#' autonomics.plot::plot_venn(x)
#' autonomics.plot::plot_venn(x, euler = TRUE)
#' filename <- tempfile()
#' autonomics.plot::plot_venn(x, filename = filename)
#'
#' @export
plot_venn <- function(
   x,
   color_values = if (assertive.properties::has_names(x)){ autonomics.plot::make_colors(names(x))
                  } else {                                 autonomics.plot::make_gg_colors(1:length(x))},
   euler        = FALSE,
   filename     = NULL,
   width        = 7,
   height       = width,
   title        = NULL,
   scaled       = FALSE,
   cat.dist     = 0.08,
   main.cex     = 1.2,
   main.pos     = c(0.5, 0.97),
   ...
){

   # Restrict to relevant color values
   color_values %<>% magrittr::extract(names(.) %in% intersect(names(.), names(x)))

   # Assert
   assertive.types::assert_is_list(x)
   assertive.types::assert_is_character(color_values)
   assertive.sets::assert_are_set_equal(names(x), names(color_values))

   # Plot Euler or Venn diagram
   p <- if (euler){
           set.seed(1)
           x  %>%
           eulerr::euler() %>% # unlimited number of sets; venn.diagram(scaled = TRUE) is limited
           graphics::plot(quantities = TRUE, fills = color_values, labels = names(x), legend = TRUE, ...)

        } else {
            futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
            VennDiagram::venn.diagram(
                    x,
                    main     = title,
                    filename = NULL,
                    fill     = color_values,
                    col      = color_values,
                    cat.col  = color_values,
                    scaled   = scaled,
                    cat.dist = cat.dist,
                    main.cex = main.cex,
                    main.pos = main.pos,
                    ...
                 )
        }

   # Print
   if (!is.null(filename))   grDevices::pdf(filename, width = width, height = height)
   autonomics.support::cdraw(p)
   if (!is.null(filename))   grDevices::dev.off()
   if (!is.null(filename))   autonomics.support::cmessage(filename)
}

#' @rdname plot_venn
#' @export
plot_venns <- function(...){
   .Deprecated('plot_venn')
   plot_vennn(...)
}

