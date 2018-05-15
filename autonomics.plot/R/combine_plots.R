#' Combine multiple ggplots
#'
#' ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
#' - cols:   Number of columns in layout
#' - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#'
#' If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
#' then plot 1 will go in the upper left, 2 will go in the upper right, and
#' 3 will go all the way across the bottom.
#' @param ...            ggplot objects
#' @param plotlist       ggplot objects as a list
#' @param cols           number of columns
#' @param layout         layout specification
#' @param file           file to which to print to
#' @param return_object  whether to return a grob
#' @param width          in inches
#' @param height         in inches
#' @examples
#'    df <- data.frame(x = 1:100, y = 1:100)
#'    p <- ggplot2::ggplot(df) + ggplot2::geom_point(ggplot2::aes(x=x,y=y))
#'    file <- paste0(tempfile(), '.pdf')
#'    autonomics.plot::combine_plots(p, p)
#'    autonomics.plot::combine_plots(p, p, file = file)
#' @author Johannes Graumann
#' @export
combine_plots <- function(
   ...,
   plotlist=NULL,
   cols=1,
   layout=NULL,
   file = NULL,
   return_object = FALSE,
   width = 7,
   height = 7
) {

  # Print to file if requested
   if (!is.null(file)){
      device <- tools::file_ext(file)
      assertive.sets::assert_is_subset(device, c('pdf'))
      get(device)(file)
      combine_plots(..., plotlist = plotlist, cols = cols, layout = layout, file = NULL, return_object = FALSE)
      grDevices::dev.off()
   }

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
     grob_output <- plots[[1]]
  } else {
     grob_output <- gridExtra::arrangeGrob(grobs = plots, layout_matrix = layout)
  }

  if (return_object)
  {
     invisible(grob_output)
  } else {
     # Do not start with grid.newpage(), which adds a blank page to the document
     grid::grid.newpage()
     grid::grid.draw(grob_output)
  }
}

#' @rdname combine_plots
#' @export
multiplot2 <- function(...){
   .Deprecated('combine_plots')
   combine_plots(...)
}
