#' Plot Venn diagrams
#' @param x named list of character vectors
#' @param filename character vector
#' @param width  width  in inches
#' @param height height in inches
#' @param title figure title
#' @param colors character vector
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
   colors    = autonomics.plot::make_gg_colors(1:length(x))
){
   futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
   colors <- autonomics.plot::make_gg_colors(1:length(x))
   p <- VennDiagram::venn.diagram(
      x,
      main     = title,
      filename = NULL,
      fill     = colors,
      col      = colors,
      cat.col  = colors,
      cat.dist = 0.08,
      main.cex = 1.2,
      main.pos = c(0.5, 0.97))
   if (!is.null(filename))   grDevices::pdf(filename, width = width, height = height)
   autonomics.support::cdraw(p)
   if (!is.null(filename))   grDevices::dev.off()
   if (!is.null(filename))   autonomics.support::cmessage(filename)
}
