#' Print graphics object to pdf
#' @param object graphics object
#' @param file file path to write to
#' @param width figure width
#' @param height figure height
#' @param ... additional arguments passed to the pdf() function
#' @export
print2pdf <- function(object, file, width = 7, height = 7, ...){
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  # TODO: replace this with withr::with_pdf once that makes it into the 
  # stable version of that package.
  grDevices::pdf(file, width = width, height = height, ...)
  on.exit(invisible(grDevices::dev.off()))
  suppressWarnings(print(object))
}