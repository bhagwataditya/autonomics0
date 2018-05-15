#' Plot Venn diagrams
#' @param x named list of character vectors
#' @param filename character vector
#' @param width  width  in inches
#' @param height height in inches
#' @param title figure title
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
   width = 7,
   height = 7,
   title = NULL
){
   futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
   p <- VennDiagram::venn.diagram(x, main = title, filename = NULL)
   if (!is.null(filename))   grDevices::pdf(filename, width = width, height = height)
   autonomics.support::cdraw(p)
   if (!is.null(filename))   grDevices::dev.off()
   if (!is.null(filename))   autonomics.support::cmessage(filename)
}
