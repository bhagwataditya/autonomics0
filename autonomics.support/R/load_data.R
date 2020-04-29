#' Load data
#' 
#' Programmatic equivalent of data(.)
#' @param x       character vector
#' @param package package name (character)
#' @return object
#' @examples
#' if (require(autonomics.data)){
#'    load_data('car90', package = 'rpart')
#' }
#' @references http://stackoverflow.com/questions/30951204/load-dataset-from-r-package-using-data-assign-it-directly-to-a-variable
#' @export
load_data <- function(x, package = NULL){
    e <- new.env()
    name <- utils::data(list = x, package = package, envir = e)[1]
    e[[name]]
}