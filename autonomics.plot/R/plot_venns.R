#' Plot Venn diagrams
#' @param x         named list of character vectors
#' @param filename  character vector
#' @param width     width  in inches
#' @param height    height in inches
#' @param title     figure title
#' @param colors    character vector
#' @param scaled    logical: whether or not to scale
#' @param cat.dist  numeric(1): category location (as distance from circle edge)
#' @param main.cex  numeric(1): title font size
#' @param main.pos  numeric(2): title location
#' @examples
#' x <- list(A = c('apple', 'pear'), B = c('pear', 'orange'))
#' filename <- NULL
#' autonomics.plot::plot_venns(x, filename = NULL)
#' filename <- tempfile()
#' autonomics.plot::plot_venns(x, filename = filename)
#' @export
plot_venns <- function(
   x,
   filename = NULL,
   width    = 7,
   height   = width,
   title    = NULL,
   colors   = autonomics.plot::make_gg_colors(1:length(x)),
   scaled   = FALSE,
   cat.dist = 0.08,
   main.cex = 1.2,
   main.pos = c(0.5, 0.97),
   ...
){
   futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
   p <- VennDiagram::venn.diagram(
           x,
           main     = title,
           filename = NULL,
           fill     = colors,
           col      = colors,
           cat.col  = colors,
           scaled   = scaled,
           cat.dist = cat.dist,
           main.cex = main.cex,
           main.pos = main.pos,
           ...
        )
   if (!is.null(filename))   grDevices::pdf(filename, width = width, height = height)
   autonomics.support::cdraw(p)
   if (!is.null(filename))   grDevices::dev.off()
   if (!is.null(filename))   autonomics.support::cmessage(filename)
}
