#' Convenient grid.draw
#' @param x grid object
#' @importFrom magrittr %>%
#' @export
cdraw <- function(x){
  grid::grid.newpage()
  x %>% grid::grid.draw()
}


#' Plot Venn/Euler diagram
#' @param x             named list of character vectors
#' @param color_values  string vector: values = colors, names = x levels
#' @param euler         TRUE/FALSE: whether to use eulerr::euler() rather than VennDiagram::venn.diagram()
#' @param filename      string
#' @param width         number: width  in inches
#' @param height        number: height in inches
#' @param title         string: figure title
#' @param scaled        TRUE/FALSE: whether or not to scale
#' @param cat.dist      number: category location (as distance from circle edge)
#' @param margin        number
#' @param main.cex      number: scaling factor for title size
#' @param main.pos      number: title location
#' @param ...           passed to either VennDiagram::venn.diagram(...) or eulerr::plot.euler(...)
#' @examples
#' x <- list(A = c('apple', 'pear'), B = c('pear', 'orange'))
#' filename <- NULL
#' plot_venn(x)
#' plot_venn(x, euler = TRUE)
#' filename <- tempfile()
#' plot_venn(x, filename = filename)
#' @importFrom magrittr %<>% %>%
#' @export
plot_venn <- function(
   x,
   color_values = if (assertive.properties::has_names(x)){ make_colors(names(x))
                  } else {                                 make_gg_colors(1:length(x))},
   euler        = FALSE,
   filename     = NULL,
   width        = 7,
   height       = width,
   title        = NULL,
   scaled       = FALSE,
   cat.dist     = 0.3,
   margin       = cat.dist,
   main.cex     = 1.5,
   main.pos     = c(0.5, 0.97),
   ...
){

   # Restrict to relevant color values
   color_values %<>% magrittr::extract(names(.) %in% intersect(names(.), names(x)))

   # Assert
   assertive.types::assert_is_list(x)
   assertive.types::assert_is_character(color_values)
   assertive.sets::assert_are_set_equal(names(x), names(color_values))

   names(x) %<>% paste0(' (', vapply(x, length, numeric(1)), ')')

   # Plot Euler or Venn diagram
   p <- if (euler){
           set.seed(1)
           x  %>%
           eulerr::euler() %>% # unlimited number of sets; venn.diagram(scaled = TRUE) is limited
           graphics::plot(quantities = TRUE, fills = color_values, labels = names(x), legend = TRUE, ...)

        } else {
            futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
            if (length(x)>5){
               autonomics.support::cmessage('\t\tLimit to first five sets (restriction of VennDiagram::venn.diagram)')
               x            %<>% magrittr::extract(1:5)
               color_values %<>% magrittr::extract(1:5)
            }
            VennDiagram::venn.diagram(
                    x,
                    main     = title,
                    filename = NULL,
                    fill     = color_values,
                    col      = color_values,
                    cat.col  = color_values,
                    scaled   = scaled,
                    cat.dist = cat.dist,
                    margin   = margin,
                    main.cex = main.cex,
                    main.pos = main.pos,
                    ...
                 )
        }

   # Print
   if (!is.null(filename))   grDevices::pdf(filename, width = width, height = height)
   cdraw(p)
   if (!is.null(filename))   grDevices::dev.off()
   if (!is.null(filename))   autonomics.support::cmessage(filename)
}

