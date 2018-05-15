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

#' Print dataframe to txt
#' 
#' Writes a data frame-like object to a tab delimited text file.
#' @param a_df A \code{dta.frame}, or similar.
#' @param file Either a character string naming a file or a connection open for 
#' writing. Passed to \code{\link[utils]{write.table}}.
#' @param quote Passed to \code{\link[utils]{write.table}}.
#' @param sep Passed to \code{\link[utils]{write.table}}.
#' @param row.names Passed to \code{\link[utils]{write.table}}.
#' @param ... Passed to \code{\link[utils]{write.table}}.
#' @seealso \code{\link[data.table]{fwrite}}
#' @examples 
#' require(magrittr)
#' my_file <- tempfile() %T>% message()
#' data.frame(a=1, b=2, c=3) %>% print2txt(file = my_file)
#' @export
print2txt <- function(a_df, file = "", quote = FALSE, sep = '\t', row.names = FALSE, ...){
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(a_df, file = file, quote = quote, sep = sep, row.names = row.names, ...)
  #write.table(a_df, file = file, quote = quote, sep = sep, row.names = row.names, ...)
}
