#' Get directory with results for particular limma coef
#' @param result_dir result directory
#' @param contrast_name name of limma contrast
#' @export
get_contrast_subdir <- function(result_dir, contrast_name){
  sprintf('%s/contrasts/%s', result_dir, contrast_name)
}